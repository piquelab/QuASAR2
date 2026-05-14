## =============================================================================
##  fitQuasar2CR.R
##
##  A drop-in alternative to fitQuasar2() that:
##    1. Reparametrizes the beta-binomial precision as psi = log(1/M).
##    2. Estimates psi per SNP with the Cox-Reid adjusted profile log-likelihood.
##    3. Moderates psi across SNPs with an edgeR/DESeq2-style trend + EB
##       shrinkage on log-dispersion.
##    4. Exposes a design-matrix interface with first-class LRT and Wald
##       contrast inference.
##
##  Public functions:
##      fitQuasar2CR()           - main fitter
##      estimateDispersionsCR()  - dispersion step alone (callable separately)
##      coef.quasar2CR()         - S3 coef method
##      vcov.quasar2CR()         - S3 vcov method
##      testCoef()               - Wald test for a single named coefficient
##      testContrast()           - Wald test for a general L %*% beta = 0
##      testLRT()                - likelihood ratio test vs. a reduced design
##      makeContrastsQuasar()    - build L from string expressions
##      summary.quasar2CR()      - per-SNP tidy summary
##
##  Internal helpers (prefixed `.bb_*`):
##      .bb_loglik_obs()         - log-likelihood per observation
##      .bb_score_beta()         - score wrt beta_g
##      .bb_info_beta()          - observed information wrt beta_g
##      .bb_fit_beta()           - inner fitter for beta_g at fixed (psi, eps)
##      .bb_apl()                - Cox-Reid APL as a function of psi
##      .bb_fit_psi_gw()         - 1D Brent max of APL (gene-wise psi)
##      .bb_fit_psi_map()        - 1D Brent max of APL + Normal prior (MAP psi)
##      .bb_fit_eps_global()     - 1D pooled max for shared eps
##
##  Depends: dplyr, tidyr, purrr (already QuASAR2 deps). locfit is optional.
## =============================================================================


## ---------- 0. Tiny safety helpers ------------------------------------------

.clip_prob <- function(p, lo = 1e-10, hi = 1 - 1e-10) {
  ## clamp probabilities away from 0 and 1 for numerical safety
  pmin(pmax(p, lo), hi)
}

.logit <- function(p)            log(p) - log1p(-p)
.invlogit <- function(eta)        1 / (1 + exp(-eta))      # equivalent to plogis()


## ---------- 1. Beta-binomial log-likelihood, score, information -------------

#' Beta-binomial log-likelihood per observation.
#'
#' @param R,A integer vectors of reference and alternate counts (length n).
#' @param mu  numeric vector in (0,1) of latent reference probabilities.
#' @param M   scalar beta-binomial precision (sum of Beta shape params), > 0.
#' @param eps scalar sequencing error in [0, 1/2).
#' @return numeric vector of length n; per-observation log-likelihood.
#' @keywords internal
.bb_loglik_obs <- function(R, A, mu, M, eps = 0) {
  mu  <- .clip_prob(mu)
  N   <- R + A
  mut <- (1 - 2 * eps) * mu + eps           # = (1-eps)*mu + eps*(1-mu)
  a   <- M * mut
  b   <- M * (1 - mut)
  lchoose(N, R) +
    lgamma(R + a) + lgamma(A + b) - lgamma(N + M) -
    lgamma(a)     - lgamma(b)     + lgamma(M)
}

#' Score (gradient) of the beta-binomial log-likelihood wrt the logistic
#' regression coefficients beta_g for a single SNP.
#' @keywords internal
.bb_score_beta <- function(beta, R, A, X, offset, M, eps) {
  eta <- as.numeric(X %*% beta) + offset
  mu  <- .invlogit(eta)
  mut <- (1 - 2 * eps) * mu + eps
  a   <- M * mut
  b   <- M * (1 - mut)
  ## d ell / d mut
  d_ll_d_mut <- M * (digamma(R + a) - digamma(A + b) -
                       digamma(a)    + digamma(b))
  ## d mut / d eta = (1 - 2 eps) * mu * (1 - mu)
  d_mut_d_eta <- (1 - 2 * eps) * mu * (1 - mu)
  u <- d_ll_d_mut * d_mut_d_eta             # length n
  drop(crossprod(X, u))                     # length p
}

#' Observed-information matrix (expected-Fisher approximation, dropping the
#' score-curvature cross-term) for beta_g at fixed (M, eps).
#' @keywords internal
.bb_info_beta <- function(beta, R, A, X, offset, M, eps) {
  eta <- as.numeric(X %*% beta) + offset
  mu  <- .invlogit(eta)
  mut <- (1 - 2 * eps) * mu + eps
  a   <- M * mut
  b   <- M * (1 - mut)
  ## - d^2 ell / d mut^2  (>= 0)
  neg_d2_d_mut2 <- M * M * (
    trigamma(a) + trigamma(b) -
      trigamma(R + a) - trigamma(A + b)
  )
  w <- neg_d2_d_mut2 * ((1 - 2 * eps) * mu * (1 - mu))^2     # diagonal weights
  ## X' W X
  crossprod(X * sqrt(pmax(w, 0)))
}


## ---------- 2. Inner fit: beta_g at fixed (psi, eps) ------------------------

#' Maximum-likelihood fit of beta_g for one SNP at fixed psi and eps.
#'
#' @param R,A integer count vectors (length n).
#' @param X   n-by-p design matrix.
#' @param offset numeric offset (length n) or 0.
#' @param psi   scalar log(1/M).
#' @param eps   scalar shared sequencing error.
#' @param beta_start  initial coefficients (length p); defaults to zero.
#' @return list(beta, info, value, convergence).
#' @keywords internal
.bb_fit_beta <- function(R, A, X, offset = 0, psi, eps,
                         beta_start = NULL,
                         control = list(maxit = 200, reltol = 1e-8)) {
  M <- exp(-psi)
  p <- ncol(X)
  if (is.null(beta_start)) beta_start <- rep(0, p)

  fn <- function(b) -sum(.bb_loglik_obs(R, A, .invlogit(as.numeric(X %*% b) + offset),
                                        M, eps))
  gr <- function(b) -.bb_score_beta(b, R, A, X, offset, M, eps)

  opt <- tryCatch(
    optim(par = beta_start, fn = fn, gr = gr, method = "BFGS",
          control = control),
    error = function(e) NULL
  )
  if (is.null(opt) || opt$convergence != 0) {
    ## fall back to Nelder-Mead without gradient
    opt <- optim(par = beta_start, fn = fn, method = "Nelder-Mead",
                 control = list(maxit = 500))
  }
  ## name the coefficient vector so downstream b[coef] / V[coef,coef] work
  beta_hat <- opt$par
  names(beta_hat) <- colnames(X)
  info <- .bb_info_beta(beta_hat, R, A, X, offset, M, eps)
  list(beta = beta_hat, info = info, value = -opt$value,
       convergence = opt$convergence)
}


## ---------- 3. Cox-Reid adjusted profile likelihood for psi -----------------

#' Cox-Reid adjusted profile log-likelihood for psi at one SNP.
#' Profiles beta_g out at the given psi (warm-started from beta_start).
#' @keywords internal
.bb_apl <- function(psi, R, A, X, offset, eps, beta_start = NULL) {
  fit <- .bb_fit_beta(R, A, X, offset, psi, eps, beta_start = beta_start)
  ## log det of observed information; if the matrix is rank-deficient,
  ## drop near-zero eigenvalues to keep the adjustment finite.
  eig <- eigen(fit$info, symmetric = TRUE, only.values = TRUE)$values
  eig <- eig[eig > 1e-12]
  apl <- fit$value - 0.5 * sum(log(eig))
  attr(apl, "fit") <- fit
  apl
}

#' Gene-wise psi: 1-D Brent max of APL_g(psi).
#' @keywords internal
.bb_fit_psi_gw <- function(R, A, X, offset, eps,
                           psi_init = -log(20),
                           interval = c(-15, 15),
                           beta_start = NULL) {
  ## Cache last fitted beta as a warm-start across Brent evaluations.
  cache <- new.env(parent = emptyenv())
  cache$beta <- beta_start
  obj <- function(psi) {
    val <- .bb_apl(psi, R, A, X, offset, eps, beta_start = cache$beta)
    cache$beta <- attr(val, "fit")$beta
    -as.numeric(val)
  }
  opt <- optim(par = psi_init, fn = obj, method = "Brent",
               lower = interval[1], upper = interval[2])
  fit <- .bb_fit_beta(R, A, X, offset, opt$par, eps, beta_start = cache$beta)
  list(psi = opt$par, apl = -opt$value, beta = fit$beta,
       info = fit$info, loglik = fit$value)
}

#' Sampling variance of hat(psi)_gw via second-difference of APL.
#' Used as v_g^samp in the empirical-Bayes variance step.
#' @keywords internal
.bb_psi_var <- function(psi_hat, R, A, X, offset, eps,
                        h = 0.1, beta_start = NULL) {
  f <- function(p) as.numeric(.bb_apl(p, R, A, X, offset, eps,
                                      beta_start = beta_start))
  f0 <- f(psi_hat)
  fp <- f(psi_hat + h)
  fm <- f(psi_hat - h)
  d2 <- (fp - 2 * f0 + fm) / (h^2)
  v  <- if (d2 < -1e-8) -1 / d2 else NA_real_
  v
}

#' MAP estimate of psi with a Normal(psi_trend, sigma2) prior.
#' @keywords internal
.bb_fit_psi_map <- function(R, A, X, offset, eps,
                            psi_trend, sigma2,
                            psi_init = NULL,
                            interval = c(-15, 15),
                            beta_start = NULL) {
  cache <- new.env(parent = emptyenv())
  cache$beta <- beta_start
  if (is.null(psi_init)) psi_init <- psi_trend
  if (!is.finite(sigma2) || sigma2 <= 0) {
    ## degenerate prior - fully shrink to trend
    fit <- .bb_fit_beta(R, A, X, offset, psi_trend, eps,
                        beta_start = beta_start)
    return(list(psi = psi_trend, beta = fit$beta, info = fit$info,
                loglik = fit$value, apl = NA_real_))
  }
  obj <- function(psi) {
    val <- .bb_apl(psi, R, A, X, offset, eps, beta_start = cache$beta)
    cache$beta <- attr(val, "fit")$beta
    -(as.numeric(val) - (psi - psi_trend)^2 / (2 * sigma2))
  }
  opt <- optim(par = psi_init, fn = obj, method = "Brent",
               lower = interval[1], upper = interval[2])
  fit <- .bb_fit_beta(R, A, X, offset, opt$par, eps, beta_start = cache$beta)
  list(psi = opt$par, beta = fit$beta, info = fit$info,
       loglik = fit$value, apl = NA_real_)
}


## ---------- 4. Shared sequencing-error update ------------------------------

#' Update the shared sequencing error eps by pooled max-likelihood, holding all
#' (beta_g, psi_g) fixed.
#' @keywords internal
.bb_fit_eps_global <- function(R, A, mu, M, eps_init = 1e-3,
                               interval = c(1e-6, 0.1)) {
  fn <- function(eps) -sum(.bb_loglik_obs(R, A, mu, M, eps))
  opt <- optim(par = eps_init, fn = fn, method = "Brent",
               lower = interval[1], upper = interval[2])
  list(eps = opt$par, value = -opt$value, convergence = opt$convergence)
}


## ---------- 5. Trend fit + EB prior variance -------------------------------

#' Fit a smooth trend psi_trend(log baseMean).
#' Uses locfit if available, else a smoothing spline.
#' @keywords internal
.bb_fit_trend <- function(log_baseMean, psi_gw) {
  ok <- is.finite(psi_gw) & is.finite(log_baseMean)
  x  <- log_baseMean[ok]
  y  <- psi_gw[ok]
  if (length(x) < 5 || length(unique(x)) < 4) {
    ## not enough x-variation - constant trend at the median
    m <- stats::median(y, na.rm = TRUE)
    return(list(predict = function(z) rep(m, length(z)),
                fit = NULL, type = "constant"))
  }
  if (requireNamespace("locfit", quietly = TRUE)) {
    fit <- tryCatch(
      locfit::locfit(y ~ locfit::lp(x, nn = 0.7)),
      error = function(e) NULL
    )
    if (!is.null(fit))
      return(list(
        predict = function(z) as.numeric(predict(fit, newdata = data.frame(x = z))),
        fit = fit, type = "locfit"))
  }
  ## spline fallback
  fit <- smooth.spline(x, y, df = min(8, length(unique(x)) - 1))
  list(predict = function(z) predict(fit, z)$y,
       fit = fit, type = "smooth.spline")
}

#' Empirical-Bayes prior variance sigma_psi^2.
#' MAD^2 / 0.6745^2  of residuals, minus median sampling variance, truncated at 0.
#' @keywords internal
.bb_prior_var <- function(resid, samp_var) {
  resid    <- resid[is.finite(resid)]
  samp_var <- samp_var[is.finite(samp_var) & samp_var > 0]
  if (length(resid) < 3) return(0)
  mad_var <- (stats::mad(resid, constant = 1.4826))^2
  med_var <- if (length(samp_var)) stats::median(samp_var) else 0
  max(0, mad_var - med_var)
}


## ---------- 6. Public: estimateDispersionsCR -------------------------------

#' Estimate per-SNP dispersion psi = log(1/M) by Cox-Reid APL with
#' empirical-Bayes shrinkage to a smooth trend in log baseMean.
#'
#' This can be called once to get a one-shot dispersion estimate, or driven by
#' fitQuasar2CR() iteratively as part of joint estimation with eps.
#'
#' @param dd       data frame with columns identifier, R, A and any variables
#'                 used in `design`. Optional `offset` column.
#' @param design   one-sided formula. Default: ~ 1.
#' @param eps      current shared sequencing error (scalar).
#' @param interval Brent search interval for psi.
#' @return data frame with one row per identifier and columns
#'         identifier, baseMean, psi_gw, psi_trend, psi_map, var_samp.
#' @export
estimateDispersionsCR <- function(dd, design = ~ 1, eps = 1e-3,
                                  interval = c(-15, 15),
                                  psi_init = -log(20),
                                  verbose = TRUE) {

  stopifnot(all(c("identifier", "R", "A") %in% names(dd)))
  if (is.null(dd$offset)) dd$offset <- 0

  ## Build per-SNP design matrices once (factor levels are taken from the full
  ## data frame, but per-SNP we keep the columns that are non-zero).
  X_full <- stats::model.matrix(design, data = dd)
  dd$.row <- seq_len(nrow(dd))

  groups <- split(dd, dd$identifier)
  if (verbose) message("Cox-Reid gene-wise dispersion for ", length(groups),
                       " identifiers")

  gw <- vector("list", length(groups))
  names(gw) <- names(groups)

  for (k in seq_along(groups)) {
    g  <- groups[[k]]
    X  <- X_full[g$.row, , drop = FALSE]
    ## drop columns that are constant zero within this SNP (rank-deficient)
    keep <- which(apply(X, 2, function(z) length(unique(z)) > 1) |
                    colnames(X) == "(Intercept)")
    Xs <- X[, keep, drop = FALSE]

    fit <- tryCatch(
      .bb_fit_psi_gw(g$R, g$A, Xs, g$offset, eps,
                     psi_init = psi_init, interval = interval),
      error = function(e) NULL
    )
    if (is.null(fit)) {
      gw[[k]] <- data.frame(identifier = unique(g$identifier),
                            baseMean   = mean(g$R + g$A),
                            psi_gw     = NA_real_,
                            var_samp   = NA_real_)
      next
    }
    v <- .bb_psi_var(fit$psi, g$R, g$A, Xs, g$offset, eps,
                     beta_start = fit$beta)
    gw[[k]] <- data.frame(identifier = unique(g$identifier),
                          baseMean   = mean(g$R + g$A),
                          psi_gw     = fit$psi,
                          var_samp   = v)
  }
  gw <- do.call(rbind, gw)

  ## Trend in log baseMean
  trend <- .bb_fit_trend(log(gw$baseMean), gw$psi_gw)
  gw$psi_trend <- trend$predict(log(gw$baseMean))

  ## Prior variance
  resid <- gw$psi_gw - gw$psi_trend
  sigma2 <- .bb_prior_var(resid, gw$var_samp)

  if (verbose) message("EB prior variance sigma^2_psi = ",
                       signif(sigma2, 3),
                       "  (MAD-based, after subtracting median var_samp)")

  attr(gw, "trend")     <- trend
  attr(gw, "sigma2")    <- sigma2
  attr(gw, "eps_used")  <- eps
  gw
}


## ---------- 7. Public: fitQuasar2CR (top-level wrapper) --------------------

#' Beta-binomial logistic regression with Cox-Reid + EB-moderated dispersion.
#'
#' This is a drop-in replacement for fitQuasar2() that uses
#'   1. psi = log(1/M) as the natural log-dispersion,
#'   2. Cox-Reid adjusted profile likelihood to estimate psi per SNP,
#'   3. an edgeR/DESeq2-style smooth trend + EB shrinkage on psi,
#'   4. a model.matrix-based design with proper Wald / LRT inference.
#'
#' @param dd       data frame with identifier, R, A, and any covariates referred
#'                 to in `design`. Optional `offset` column.
#' @param design   one-sided formula, e.g. ~ Batch + Treatment.
#' @param eps_init initial sequencing error (default 1e-3).
#' @param psi_init initial log(1/M) (default -log(20)).
#' @param max_iter outer-loop iteration cap.
#' @param tol_psi,tol_beta,tol_eps convergence tolerances.
#' @param verbose  logical.
#' @return an object of class "quasar2CR" with components
#'   - call, design, formula
#'   - coefficients: named list of beta_g per identifier
#'   - vcov: named list of covariance matrices per identifier
#'   - dispersion: data frame with psi estimates per identifier
#'   - eps, loglik, n_iter, converged
#' @export
fitQuasar2CR <- function(dd, design = ~ 1,
                         eps_init = 1e-3, psi_init = -log(20),
                         max_iter = 10,
                         tol_psi  = 0.01, tol_beta = 0.01, tol_eps = 2e-4,
                         verbose = TRUE) {

  call <- match.call()
  stopifnot(inherits(design, "formula"))
  stopifnot(all(c("identifier", "R", "A") %in% names(dd)))
  if (is.null(dd$offset)) dd$offset <- 0

  ## Stable per-SNP order; row-index lookup for the global model matrix
  dd <- dd[order(dd$identifier), ]
  dd$.row <- seq_len(nrow(dd))
  X_full  <- stats::model.matrix(design, data = dd)
  groups  <- split(dd, dd$identifier)
  if (verbose) message("fitQuasar2CR: ", length(groups), " identifiers, ",
                       ncol(X_full), " design columns.")

  ## Per-SNP design matrices, dropping all-constant non-intercept columns.
  Xs   <- vector("list", length(groups))
  keep <- vector("list", length(groups))
  for (k in seq_along(groups)) {
    g  <- groups[[k]]
    X  <- X_full[g$.row, , drop = FALSE]
    kp <- colnames(X) %in% "(Intercept)" |
          apply(X, 2, function(z) length(unique(z)) > 1)
    Xs[[k]]   <- X[, kp, drop = FALSE]
    keep[[k]] <- kp
  }
  names(Xs) <- names(groups)

  ## ----- initial beta_g at psi_init, eps_init -----
  eps <- eps_init
  psi <- rep(psi_init, length(groups))
  names(psi) <- names(groups)
  betas <- vector("list", length(groups));  names(betas) <- names(groups)
  infos <- vector("list", length(groups));  names(infos) <- names(groups)
  if (verbose) message("Initial coefficient fit (psi=", signif(psi_init, 3),
                       ", eps=", signif(eps_init, 3), ")")
  for (k in seq_along(groups)) {
    g <- groups[[k]]
    fit <- .bb_fit_beta(g$R, g$A, Xs[[k]], g$offset, psi[k], eps)
    betas[[k]] <- fit$beta;  infos[[k]] <- fit$info
  }

  notConverged <- TRUE; it <- 0
  prior_var_path <- numeric(0)
  while (notConverged && it < max_iter) {
    it <- it + 1
    if (verbose) message("\n--- Outer iteration ", it, " ---")

    ## (a) Gene-wise psi via Cox-Reid APL
    psi_old <- psi
    psi_gw  <- numeric(length(groups))
    var_g   <- numeric(length(groups))
    baseMean <- numeric(length(groups))
    for (k in seq_along(groups)) {
      g <- groups[[k]]
      fit <- .bb_fit_psi_gw(g$R, g$A, Xs[[k]], g$offset, eps,
                            psi_init = psi[k], beta_start = betas[[k]])
      psi_gw[k] <- fit$psi
      var_g[k]  <- .bb_psi_var(fit$psi, g$R, g$A, Xs[[k]], g$offset, eps,
                               beta_start = fit$beta)
      baseMean[k] <- mean(g$R + g$A)
      betas[[k]] <- fit$beta;  infos[[k]] <- fit$info
    }
    if (verbose) message("  gene-wise APL: median psi = ",
                         signif(stats::median(psi_gw), 3))

    ## (b) Trend on log baseMean
    trend <- .bb_fit_trend(log(baseMean), psi_gw)
    psi_trend <- trend$predict(log(baseMean))

    ## (c) Prior variance
    resid  <- psi_gw - psi_trend
    sigma2 <- .bb_prior_var(resid, var_g)
    prior_var_path <- c(prior_var_path, sigma2)
    if (verbose) message("  EB sigma^2_psi = ", signif(sigma2, 3))

    ## save for output
    psi_trend_out <- psi_trend
    var_samp_out  <- var_g
    psi_gw_out    <- psi_gw
    sigma2_out    <- sigma2

    ## (d) MAP psi (and updated beta)
    for (k in seq_along(groups)) {
      g <- groups[[k]]
      fit <- .bb_fit_psi_map(g$R, g$A, Xs[[k]], g$offset, eps,
                             psi_trend = psi_trend[k], sigma2 = sigma2,
                             psi_init = psi_gw[k], beta_start = betas[[k]])
      psi[k]     <- fit$psi
      betas[[k]] <- fit$beta
      infos[[k]] <- fit$info
    }

    ## (e) global eps
    eps_old <- eps
    mu_all  <- numeric(nrow(dd))
    M_all   <- numeric(nrow(dd))
    for (k in seq_along(groups)) {
      g <- groups[[k]]
      eta <- as.numeric(Xs[[k]] %*% betas[[k]]) + g$offset
      mu_all[g$.row] <- .invlogit(eta)
      M_all [g$.row] <- exp(-psi[k])
    }
    eps_fit <- .bb_fit_eps_global(dd$R, dd$A, mu_all, M_all, eps_init = eps)
    eps     <- eps_fit$eps
    if (verbose) message("  eps update: ", signif(eps_old, 3), " -> ",
                         signif(eps, 3))

    ## (f) Refit beta at new (psi, eps)
    beta_change <- 0
    for (k in seq_along(groups)) {
      g <- groups[[k]]
      old <- betas[[k]]
      fit <- .bb_fit_beta(g$R, g$A, Xs[[k]], g$offset, psi[k], eps,
                          beta_start = old)
      beta_change <- max(beta_change, max(abs(fit$beta - old)))
      betas[[k]] <- fit$beta
      infos[[k]] <- fit$info
    }

    psi_change <- max(abs(psi - psi_old), na.rm = TRUE)
    eps_change <- abs(eps - eps_old)
    if (verbose) message(sprintf("  max |dpsi|=%.3g  max |dbeta|=%.3g  |deps|=%.3g",
                                 psi_change, beta_change, eps_change))
    if (it > 1 &&
        psi_change  < tol_psi  &&
        beta_change < tol_beta &&
        eps_change  < tol_eps) {
      notConverged <- FALSE
    }
  }

  ## ----- Build output -----
  vcov_list <- lapply(infos, function(I) {
    out <- tryCatch(solve(I), error = function(e) MASS::ginv(I))
    rownames(out) <- colnames(out) <- colnames(I)
    out
  })

  loglik <- numeric(length(groups))
  for (k in seq_along(groups)) {
    g <- groups[[k]]
    eta <- as.numeric(Xs[[k]] %*% betas[[k]]) + g$offset
    loglik[k] <- sum(.bb_loglik_obs(g$R, g$A, .invlogit(eta),
                                    exp(-psi[k]), eps))
  }
  names(loglik) <- names(groups)

  disp_df <- data.frame(
    identifier = names(groups),
    baseMean   = vapply(groups, function(g) mean(g$R + g$A), 0),
    psi        = psi,
    M          = exp(-psi),
    psi_gw     = if (exists("psi_gw_out"))    psi_gw_out    else NA_real_,
    psi_trend  = if (exists("psi_trend_out")) psi_trend_out else NA_real_,
    var_samp   = if (exists("var_samp_out"))  var_samp_out  else NA_real_,
    row.names  = NULL
  )

  fit_obj <- list(
    call         = call,
    formula      = design,
    coefficients = betas,
    vcov         = vcov_list,
    designs      = Xs,
    designs_keep = keep,
    dispersion   = disp_df,
    sigma2_psi   = if (exists("sigma2_out")) sigma2_out else NA_real_,
    eps          = eps,
    loglik       = loglik,
    n_iter       = it,
    converged    = !notConverged,
    prior_var_path = prior_var_path,
    data         = dd
  )
  class(fit_obj) <- "quasar2CR"
  fit_obj
}


## ---------- 8. S3 methods: coef, vcov, print, summary ----------------------

#' @export
coef.quasar2CR <- function(object, id = NULL, ...) {
  if (is.null(id)) return(object$coefficients)
  object$coefficients[[id]]
}

#' @export
vcov.quasar2CR <- function(object, id = NULL, ...) {
  if (is.null(id)) return(object$vcov)
  object$vcov[[id]]
}

#' @export
print.quasar2CR <- function(x, ...) {
  cat("Quasar2 (Cox-Reid + EB) fit\n")
  cat("  Identifiers:    ", length(x$coefficients), "\n")
  cat("  Design:         ", deparse(x$formula), "\n")
  cat("  Coefficients:   ", paste(colnames(x$designs[[1]]), collapse = ", "), "\n")
  cat("  Shared eps:     ", signif(x$eps, 4), "\n")
  cat("  Iterations:     ", x$n_iter, " (converged: ", x$converged, ")\n", sep = "")
  invisible(x)
}

#' Per-SNP tidy summary.
#' Returns a long data frame with one row per (identifier, coefficient).
#' @export
summary.quasar2CR <- function(object, ...) {
  rows <- lapply(seq_along(object$coefficients), function(k) {
    id  <- names(object$coefficients)[k]
    b   <- object$coefficients[[k]]
    V   <- object$vcov[[k]]
    se  <- sqrt(diag(V))
    z   <- b / se
    p   <- 2 * stats::pnorm(-abs(z))
    data.frame(identifier = id,
               term       = names(b),
               estimate   = unname(b),
               std.error  = unname(se),
               statistic  = unname(z),
               p.value    = unname(p),
               row.names  = NULL,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$padj <- stats::p.adjust(out$p.value, method = "BH")
  out
}


## ---------- 9. Wald tests: testCoef, testContrast --------------------------

#' Wald test for a single named coefficient across all SNPs.
#' @param object quasar2CR fit.
#' @param coef   string, the coefficient name (must match colnames of the design
#'               matrix). Use colnames(object$designs[[1]]) to see them.
#' @param ref    reference distribution.  "t" (default) uses Student t; the df
#'               is determined by df_method.  "normal" uses z; only correct
#'               when n_g is large or dispersion is known exactly.
#' @param df_method   how to compute df for the t reference:
#'   "residual" (default; conservative): df = n_g - p_g
#'   "moderated":  df = n_g - p_g + 1/sigma2_psi -- adds the prior df from EB
#'                  shrinkage (analogous to limma's d + d_0).  Closer to
#'                  nominal Type-I error empirically (~0.052 vs ~0.035 for
#'                  "residual" on a coverage x dispersion x n_g grid), at the
#'                  cost of occasionally crossing nominal in pathological
#'                  near-binomial cases.
#'   "manual":      pass a numeric `df_extra` to add to n_g - p_g.
#' @param df_extra    numeric, only used when df_method = "manual".
#' @return data frame with one row per SNP (identifier, estimate, std.error, z, df, pval, padj).
#' @export
testCoef <- function(object, coef,
                     ref        = c("t", "normal"),
                     df_method  = c("residual", "moderated", "manual"),
                     df_extra   = 0) {
  stopifnot(inherits(object, "quasar2CR"))
  ref       <- match.arg(ref)
  df_method <- match.arg(df_method)
  sigma2    <- object$sigma2_psi
  if (is.null(sigma2)) sigma2 <- NA_real_

  rows <- lapply(seq_along(object$coefficients), function(k) {
    id <- names(object$coefficients)[k]
    b  <- object$coefficients[[k]]
    V  <- object$vcov[[k]]
    if (!(coef %in% names(b)))
      return(data.frame(identifier = id, estimate = NA, std.error = NA,
                        statistic = NA, df = NA, p.value = NA))
    est <- unname(b[coef])
    se  <- sqrt(V[coef, coef])
    stat <- est / se
    n_g <- nrow(object$designs[[k]])
    p_g <- length(b)
    df <- switch(df_method,
                 residual  = n_g - p_g,
                 moderated = {
                   d0 <- if (is.finite(sigma2) && sigma2 > 0) 1 / sigma2
                         else 0
                   n_g - p_g + d0
                 },
                 manual    = n_g - p_g + df_extra)
    df  <- max(1, df)
    p   <- if (ref == "t") 2 * stats::pt(-abs(stat), df = df)
           else            2 * stats::pnorm(-abs(stat))
    data.frame(identifier = id, estimate = est, std.error = se,
               statistic = stat, df = df, p.value = p,
               stringsAsFactors = FALSE)
  })
  out <- do.call(rbind, rows)
  out$padj <- stats::p.adjust(out$p.value, method = "BH")
  out
}

#' Wald test for a general contrast L %*% beta = 0.
#' @param object  quasar2CR fit.
#' @param L       a contrast vector (numeric) or matrix (k rows). Names / colnames
#'                must match coefficient names.
#' @param ref     reference distribution.  "F" (default) uses F(k, df2),
#'                "chisq" uses chi-square_k.
#' @param df_method  see ?testCoef for definitions.  Only used when ref="F".
#' @return data frame with statistic, df1, df2, p-value (per SNP).
#' @export
testContrast <- function(object, L,
                         ref       = c("F", "chisq"),
                         df_method = c("residual", "moderated", "manual"),
                         df_extra  = 0) {
  stopifnot(inherits(object, "quasar2CR"))
  ref       <- match.arg(ref)
  df_method <- match.arg(df_method)
  sigma2    <- object$sigma2_psi
  if (is.null(sigma2)) sigma2 <- NA_real_
  if (is.null(dim(L))) L <- matrix(L, nrow = 1, dimnames = list(NULL, names(L)))

  rows <- lapply(seq_along(object$coefficients), function(k) {
    id <- names(object$coefficients)[k]
    b  <- object$coefficients[[k]]
    V  <- object$vcov[[k]]
    cn <- colnames(L);  if (is.null(cn)) cn <- names(b)
    miss <- setdiff(cn, names(b))
    if (length(miss))
      return(data.frame(identifier = id, statistic = NA, df1 = NA, df2 = NA,
                        p.value = NA))
    Lk  <- L[, names(b), drop = FALSE]
    est <- drop(Lk %*% b)
    VV  <- Lk %*% V %*% t(Lk)
    chi2 <- tryCatch(drop(crossprod(est, solve(VV, est))),
                     error = function(e) NA_real_)
    df1 <- nrow(Lk)
    n_g <- nrow(object$designs[[k]])
    p_g <- length(b)
    df2 <- switch(df_method,
                  residual  = n_g - p_g,
                  moderated = {
                    d0 <- if (is.finite(sigma2) && sigma2 > 0) 1 / sigma2
                          else 0
                    n_g - p_g + d0
                  },
                  manual    = n_g - p_g + df_extra)
    df2 <- max(1, df2)
    if (ref == "F") {
      Fstat <- chi2 / df1
      p     <- stats::pf(Fstat, df1, df2, lower.tail = FALSE)
      data.frame(identifier = id, statistic = Fstat,
                 df1 = df1, df2 = df2, p.value = p, stringsAsFactors = FALSE)
    } else {
      p <- stats::pchisq(chi2, df1, lower.tail = FALSE)
      data.frame(identifier = id, statistic = chi2,
                 df1 = df1, df2 = NA, p.value = p, stringsAsFactors = FALSE)
    }
  })
  out <- do.call(rbind, rows)
  out$padj <- stats::p.adjust(out$p.value, method = "BH")
  out
}


## ---------- 10. LRT against a reduced design -------------------------------

#' Likelihood ratio test of `full` vs a reduced design, refitting beta_g per
#' SNP at the **same dispersions and same eps** as `object`.
#'
#' @param object   the full-design fit (class "quasar2CR").
#' @param reduced  a formula. Must be a strict sub-set of object$formula's
#'                 columns after expansion via model.matrix.
#' @param ref      reference distribution.  "F" (default) is the QLF analogue
#'                 F = Lambda / df_diff, "chisq" is the asymptotic chi-square.
#' @param df_method  see ?testCoef.  Only used when ref="F".
#' @return data frame with statistic, df1, df2, p.value (and padj) per SNP.
#' @export
testLRT <- function(object, reduced,
                    ref        = c("F", "chisq"),
                    df_method  = c("residual", "moderated", "manual"),
                    df_extra   = 0) {
  stopifnot(inherits(object, "quasar2CR"))
  ref       <- match.arg(ref)
  df_method <- match.arg(df_method)
  sigma2    <- object$sigma2_psi
  if (is.null(sigma2)) sigma2 <- NA_real_
  dd <- object$data

  X_red_full <- stats::model.matrix(reduced, data = dd)

  rows <- lapply(seq_along(object$coefficients), function(k) {
    id    <- names(object$coefficients)[k]
    g     <- dd[dd$identifier == id, , drop = FALSE]
    X_red <- X_red_full[g$.row, , drop = FALSE]
    kp <- colnames(X_red) %in% "(Intercept)" |
          apply(X_red, 2, function(z) length(unique(z)) > 1)
    X_red <- X_red[, kp, drop = FALSE]

    psi_g <- object$dispersion$psi[k]
    fit_red <- tryCatch(
      .bb_fit_beta(g$R, g$A, X_red, g$offset, psi_g, object$eps),
      error = function(e) NULL
    )
    if (is.null(fit_red))
      return(data.frame(identifier = id, statistic = NA,
                        df1 = NA, df2 = NA, p.value = NA))
    ll_red  <- fit_red$value
    ll_full <- object$loglik[k]
    df_full <- length(object$coefficients[[k]])
    df_red  <- length(fit_red$beta)
    chi2 <- 2 * (ll_full - ll_red)
    chi2 <- pmax(chi2, 0)
    df1  <- df_full - df_red
    n_g  <- nrow(object$designs[[k]])
    df2 <- switch(df_method,
                  residual  = n_g - df_full,
                  moderated = {
                    d0 <- if (is.finite(sigma2) && sigma2 > 0) 1 / sigma2
                          else 0
                    n_g - df_full + d0
                  },
                  manual    = n_g - df_full + df_extra)
    df2 <- max(1, df2)
    if (df1 <= 0)
      return(data.frame(identifier = id, statistic = NA,
                        df1 = df1, df2 = df2, p.value = NA))
    if (ref == "F") {
      Fstat <- chi2 / df1
      p     <- stats::pf(Fstat, df1, df2, lower.tail = FALSE)
      data.frame(identifier = id, statistic = Fstat,
                 df1 = df1, df2 = df2, p.value = p, stringsAsFactors = FALSE)
    } else {
      p <- stats::pchisq(chi2, df1, lower.tail = FALSE)
      data.frame(identifier = id, statistic = chi2,
                 df1 = df1, df2 = NA, p.value = p, stringsAsFactors = FALSE)
    }
  })
  out <- do.call(rbind, rows)
  out$padj <- stats::p.adjust(out$p.value, method = "BH")
  out
}


## ---------- 11. Contrast-builder helper ------------------------------------

#' Build a contrast matrix L from string expressions of coefficient names.
#' Inspired by limma::makeContrasts() but pulls the levels from a quasar2CR fit.
#'
#' Example:
#'   L <- makeContrastsQuasar(c("TreatmentB - TreatmentA",
#'                              "TreatmentC - TreatmentA"), fit)
#'   testContrast(fit, L)
#'
#' @export
makeContrastsQuasar <- function(contrasts, object) {
  stopifnot(inherits(object, "quasar2CR"))
  cn <- colnames(object$designs[[1]])
  L  <- matrix(0, nrow = length(contrasts), ncol = length(cn),
               dimnames = list(contrasts, cn))
  for (i in seq_along(contrasts)) {
    expr <- contrasts[[i]]
    ## tokenize on + and -
    tokens <- regmatches(expr, gregexpr("[+-]?[^+-]+", expr))[[1]]
    tokens <- gsub("\\s", "", tokens)
    for (tk in tokens) {
      sign <- 1
      if (substr(tk, 1, 1) == "-") { sign <- -1; tk <- substring(tk, 2) }
      if (substr(tk, 1, 1) == "+") { tk <- substring(tk, 2) }
      ## handle coefficients like "2*TreatmentA"
      mult <- 1
      if (grepl("\\*", tk)) {
        parts <- strsplit(tk, "\\*")[[1]]
        mult  <- as.numeric(parts[1])
        tk    <- parts[2]
      }
      if (!(tk %in% cn))
        stop("Contrast term '", tk, "' not in coefficient names: ",
             paste(cn, collapse = ", "))
      L[i, tk] <- L[i, tk] + sign * mult
    }
  }
  L
}


## =============================================================================
## End of file.
## =============================================================================
