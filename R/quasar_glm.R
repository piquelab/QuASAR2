#' @import dplyr
#' @import tidyr
#' @import stats
#' @importFrom magrittr %>%
NULL

#' Log-likelihood for Beta-Binomial with Phi reparameterization
#' @param b Coefficients
#' @param phi Dispersion parameter log(1/M)
#' @param X Design matrix
#' @param R Reference reads
#' @param A Alternate reads
#' @param offset Offset
#' @param eps Sequencing error
#' @export
logLikBetaBinomialPhi <- function(b, phi, X, R, A, offset = 0, eps = 0.001) {
  M <- exp(-phi)
  eta <- X %*% b + offset
  pi0 <- plogis(eta)
  p <- pi0 * (1 - eps) + (1 - pi0) * eps
  
  # Avoid p=0 or p=1 exactly if eps=0
  p <- pmin(pmax(p, 1e-10), 1 - 1e-10)
  
  ll <- lgamma(M) + lgamma(R + p * M) + lgamma(A + (1 - p) * M) - 
        lgamma(R + A + M) - lgamma(p * M) - lgamma((1 - p) * M)
  
  return(-sum(ll))
}

#' Gradient of Log-likelihood for Beta-Binomial with Phi reparameterization
#' 
#' @param b Coefficients
#' @param phi Dispersion parameter log(1/M)
#' @param X Design matrix
#' @param R Reference reads
#' @param A Alternate reads
#' @param offset Offset
#' @param eps Sequencing error
#' @export
gLogLikBetaBinomialPhi <- function(b, phi, X, R, A, offset = 0, eps = 0.001) {
  M <- exp(-phi)
  eta <- X %*% b + offset
  pi0 <- plogis(eta)
  p <- pi0 * (1 - eps) + (1 - pi0) * eps
  
  # d_ll / d_p
  dll_dp <- M * (digamma(R + p * M) - digamma(A + (1 - p) * M) - digamma(p * M) + digamma((1 - p) * M))
  
  # d_p / d_eta = (1-2*eps) * pi0 * (1-pi0)
  dp_deta <- (1 - 2 * eps) * pi0 * (1 - pi0)
  
  # d_ll / d_eta
  dll_deta <- dll_dp * dp_deta
  
  grad_b <- -t(X) %*% dll_deta
  return(as.vector(grad_b))
}

#' Fit Beta-Binomial GLM for a single SNP
#' 
#' @param dd Data frame for one SNP
#' @param model Model formula
#' @param phi Dispersion parameter log(1/M)
#' @param eps Sequencing error
#' @export
fitBetaBinomialGLM <- function(dd, model, phi, eps = 0.001, b_init = NULL) {
  X <- model.matrix(model, dd)
  if (is.null(b_init)) b_init <- rep(0, ncol(X))
  
  res <- optim(b_init, 
               fn = logLikBetaBinomialPhi, 
               gr = gLogLikBetaBinomialPhi,
               phi = phi, X = X, R = dd$R, A = dd$A, 
               offset = if("offset" %in% names(dd)) dd$offset else 0,
               eps = eps,
               method = "BFGS", hessian = TRUE)
  
  return(res)
}

#' Cox-Reid Adjusted Profile Log-Likelihood for Phi
#' 
#' @param phi Dispersion parameter log(1/M)
#' @param dd Data frame for one SNP
#' @param model Model formula
#' @param eps Sequencing error
#' @export
adjustedLogLikPhi <- function(phi, dd, model, eps = 0.001) {
  fit <- fitBetaBinomialGLM(dd, model, phi, eps)
  ll <- -fit$value # negative log-likelihood was returned by optim
  
  # Calculate Cox-Reid adjustment
  X <- model.matrix(model, dd)
  eta <- X %*% fit$par + (if("offset" %in% names(dd)) dd$offset else 0)
  pi0 <- plogis(eta)
  p <- pi0 * (1 - eps) + (1 - pi0) * eps
  n <- dd$R + dd$A
  M <- exp(-phi)
  theta <- 1 / (M + 1)
  
  # Expected information weights
  # Var(R) = n * p * (1-p) * (1 + (n-1)*theta)
  # d_mu / d_eta = n * (1-2*eps) * pi0 * (1-pi0)
  weights <- (n^2 * (1 - 2 * eps)^2 * pi0^2 * (1 - pi0)^2) / 
             (n * p * (1 - p) * (1 + (n - 1) * theta))
  
  # Handle potential numerical issues
  weights[is.na(weights) | is.infinite(weights)] <- 0
  
  adj <- 0.5 * determinant(t(X) %*% (as.vector(weights) * X), logarithm = TRUE)$modulus
  
  return(ll - adj)
}

#' Estimate Dispersion for QuASAR2
#' 
#' @param dd Data frame with all SNPs
#' @param model Model formula
#' @param eps Sequencing error
#' @param prior.df Prior degrees of freedom for moderation. If NULL, it will be estimated.
#' @export
estimateDispersionQuasar <- function(dd, model, eps = 0.001, prior.df = NULL) {
  ids <- unique(dd$identifier)
  n_snps <- length(ids)
  
  # Calculate average log-coverage for each SNP
  cat("Calculating average abundances...\n")
  abundance <- dd %>% 
    group_by(identifier) %>% 
    summarize(avg_log_cov = log(mean(R + A)))
  avg_log_cov <- abundance$avg_log_cov
  names(avg_log_cov) <- abundance$identifier
  
  cat("Estimating common dispersion...\n")
  sample_ids <- if(n_snps > 500) sample(ids, 500) else ids
  common_ll <- function(phi) {
    sum(sapply(sample_ids, function(id) {
      adjustedLogLikPhi(phi, dd[dd$identifier == id, ], model, eps)
    }))
  }
  opt_common <- optimize(common_ll, interval = c(-10, 10), maximum = TRUE)
  phi_common <- opt_common$maximum
  cat("Common phi:", phi_common, "(M =", exp(-phi_common), ")\n")
  
  cat("Estimating trended dispersion...\n")
  # Simple local regression (loess) on a sample of tagwise dispersions
  # First get some rough tagwise estimates
  phi_rough <- sapply(sample_ids, function(id) {
    target <- function(phi) adjustedLogLikPhi(phi, dd[dd$identifier == id, ], model, eps)
    opt <- optimize(target, interval = c(-10, 10), maximum = TRUE)
    return(opt$maximum)
  })
  
  # Fit trend
  trend_fit <- loess(phi_rough ~ avg_log_cov[sample_ids], span = 0.7)
  phi_trend <- predict(trend_fit, newdata = data.frame(avg_log_cov = avg_log_cov))
  # Fallback to common if prediction fails
  phi_trend[is.na(phi_trend)] <- phi_common
  names(phi_trend) <- ids
  
  if (is.null(prior.df)) {
    cat("Estimating prior degrees of freedom...\n")
    # Heuristic: base it on the variance of (phi_rough - trend)
    # This is a simplification of the edgeR approach
    resid_var <- var(phi_rough - phi_trend[sample_ids], na.rm = TRUE)
    # We want prior.df such that it represents the "inverse variance" of the prior
    # For a Normal(mu, sigma^2) prior, the penalty is 1/(2*sigma^2)
    # Here we use prior.df as the weight. 
    prior.df <- 1 / max(resid_var, 0.1)
    cat("Estimated prior weight (prior.df):", prior.df, "\n")
  }

  cat("Estimating moderated tagwise dispersions...\n")
  phi_tagwise <- sapply(ids, function(id) {
    snp_dd <- dd[dd$identifier == id, ]
    target <- function(phi) {
      adjustedLogLikPhi(phi, snp_dd, model, eps) - 0.5 * prior.df * (phi - phi_trend[id])^2
    }
    opt <- optimize(target, interval = c(-10, 10), maximum = TRUE)
    return(opt$maximum)
  })
  names(phi_tagwise) <- ids
  
  return(list(common = phi_common, trend = phi_trend, tagwise = phi_tagwise, prior.df = prior.df))
}

#' Refactored fitQuasar2 with GLM and moderated dispersion
#' 
#' @param dd Data frame
#' @param model Model formula
#' @param eps Sequencing error
#' @export
fitQuasar_GLM <- function(dd, model, eps = 0.001, prior.df = 10) {
  # 1. Estimate Dispersions
  disps <- estimateDispersionQuasar(dd, model, eps, prior.df)
  
  # 2. Fit GLMs for each SNP given their moderated dispersion
  cat("Fitting final GLMs...\n")
  ids <- unique(dd$identifier)
  results_list <- lapply(ids, function(id) {
    snp_dd <- dd[dd$identifier == id, ]
    phi <- disps$tagwise[id]
    fit <- fitBetaBinomialGLM(snp_dd, model, phi, eps)
    
    X <- model.matrix(model, snp_dd)
    iH <- solve(fit$hessian)
    se <- sqrt(diag(iH))
    
    tibble(
      identifier = id,
      term = colnames(X),
      estimate = fit$par,
      std.error = se,
      statistic = estimate / std.error,
      p.value = 2 * pnorm(-abs(statistic)),
      phi = phi,
      M = exp(-phi)
    )
  })
  
  results <- bind_rows(results_list)
  
  return(list(
    results = results,
    dispersions = disps,
    eps = eps
  ))
}

#' Likelihood Ratio Test for QuASAR GLM
#' 
#' @param dd Data frame for one SNP
#' @param full_model Full model formula
#' @param reduced_model Reduced model formula
#' @param phi Dispersion
#' @param eps Sequencing error
#' @export
lrtQuasar <- function(dd, full_model, reduced_model, phi, eps = 0.001) {
  fit_full <- fitBetaBinomialGLM(dd, full_model, phi, eps)
  fit_red <- fitBetaBinomialGLM(dd, reduced_model, phi, eps)
  
  lrt_stat <- 2 * (fit_red$value - fit_full$value) # optim returns negative log-lik
  df <- ncol(model.matrix(full_model, dd)) - ncol(model.matrix(reduced_model, dd))
  pval <- pchisq(lrt_stat, df = df, lower.tail = FALSE)
  
  return(list(statistic = lrt_stat, df = df, p.value = pval))
}
