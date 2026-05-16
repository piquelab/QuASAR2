## ============================================================================
##  fitQuasar2CR.R
##
##  Beta-binomial logistic regression with Cox-Reid + EB-moderated dispersion.
##  See vignettes/quasar2CR.Rmd for the math.
## ============================================================================

#' Fit beta-binomial logistic regression with Cox-Reid + EB-moderated dispersion.
#'
#' A drop-in alternative to \code{\link{fitQuasar2}} that:
#' \enumerate{
#'   \item parameterizes the BB precision as \eqn{\psi = \log(1/M)},
#'   \item estimates \eqn{\psi} per SNP using the Cox-Reid adjusted profile
#'         likelihood (instead of a coverage-binned grid),
#'   \item moderates \eqn{\psi} across SNPs using an edgeR / DESeq2 style
#'         trend + empirical-Bayes shrinkage, and
#'   \item exposes a design-matrix interface so that likelihood ratio tests
#'         and contrasts are first-class operations.
#' }
#'
#' @param dd       data frame.  Must contain columns \code{identifier},
#'                 \code{R}, \code{A}, plus any variables referred to in
#'                 \code{design}.  Optionally a numeric \code{offset} column,
#'                 e.g. \code{qlogis(DNA_proportion)} for MPRA.
#' @param design   one-sided formula, e.g. \code{~ Batch + Treatment}.
#' @param eps_init initial shared sequencing error (default \code{1e-3}).
#' @param psi_init initial log(1/M) (default \code{-log(20)}).
#' @param max_iter outer-loop iteration cap (default 10).
#' @param tol_psi,tol_beta,tol_eps convergence tolerances.
#' @param verbose  logical; print iteration progress.
#'
#' @return An object of class \code{quasar2CR}, a list with:
#'   \describe{
#'     \item{coefficients}{named list: one numeric vector of \eqn{\hat\beta_g} per SNP.}
#'     \item{vcov}{named list: one Wald covariance matrix per SNP.}
#'     \item{dispersion}{data frame with \code{psi}, \code{M}, \code{psi_gw},
#'           \code{psi_trend}, \code{var_samp}, \code{baseMean} per SNP.}
#'     \item{sigma2_psi}{the empirical-Bayes prior variance on \eqn{\psi}.}
#'     \item{eps}{the shared sequencing-error MLE.}
#'     \item{loglik}{vector of per-SNP log-likelihoods at convergence.}
#'     \item{designs}{list of per-SNP design matrices used in fitting.}
#'     \item{formula, data, call, n_iter, converged}{bookkeeping.}
#'   }
#'
#' @seealso \code{\link{testCoef}}, \code{\link{testContrast}},
#'   \code{\link{testLRT}}, \code{\link{estimateDispersionsCR}}.
#'
#' @examples
#' \dontrun{
#' set.seed(1)
#' G <- 50; n_g <- 8
#' trt <- factor(rep(c("A","A","A","A","B","B","B","B"), length.out = n_g))
#' M_g <- exp(runif(G, log(20), log(200)))
#' dat <- do.call(rbind, lapply(seq_len(G), function(g) {
#'   pp <- rbeta(n_g, M_g[g]*0.5, M_g[g]*0.5)
#'   N  <- sample(40:150, n_g, replace = TRUE)
#'   R  <- rbinom(n_g, N, pp); A <- N - R
#'   data.frame(identifier = sprintf("snp_%03d", g),
#'              R = R, A = A, trt = trt)
#' }))
#' fit <- fitQuasar2CR(dat, design = ~ trt)
#' head(testCoef(fit, "trtB"))
#' head(testLRT(fit, reduced = ~ 1))
#' }
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

  dd <- dd[order(dd$identifier), ]
  dd$.row <- seq_len(nrow(dd))
  X_full  <- stats::model.matrix(design, data = dd)
  groups  <- split(dd, dd$identifier)
  if (verbose) message("fitQuasar2CR: ", length(groups), " identifiers, ",
                       ncol(X_full), " design columns.")

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
  psi_gw_out <- psi_trend_out <- var_samp_out <- rep(NA_real_, length(groups))
  sigma2_out <- NA_real_

  while (notConverged && it < max_iter) {
    it <- it + 1
    if (verbose) message("\n--- Outer iteration ", it, " ---")

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

    trend <- .bb_fit_trend(log(baseMean), psi_gw)
    psi_trend <- trend$predict(log(baseMean))

    resid  <- psi_gw - psi_trend
    sigma2 <- .bb_prior_var(resid, var_g)
    prior_var_path <- c(prior_var_path, sigma2)
    if (verbose) message("  EB sigma^2_psi = ", signif(sigma2, 3))

    psi_gw_out    <- psi_gw
    psi_trend_out <- psi_trend
    var_samp_out  <- var_g
    sigma2_out    <- sigma2

    for (k in seq_along(groups)) {
      g <- groups[[k]]
      fit <- .bb_fit_psi_map(g$R, g$A, Xs[[k]], g$offset, eps,
                             psi_trend = psi_trend[k], sigma2 = sigma2,
                             psi_init = psi_gw[k], beta_start = betas[[k]])
      psi[k]     <- fit$psi
      betas[[k]] <- fit$beta
      infos[[k]] <- fit$info
    }

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
    psi_gw     = psi_gw_out,
    psi_trend  = psi_trend_out,
    var_samp   = var_samp_out,
    row.names  = NULL
  )

  fit_obj <- list(
    call           = call,
    formula        = design,
    coefficients   = betas,
    vcov           = vcov_list,
    designs        = Xs,
    designs_keep   = keep,
    dispersion     = disp_df,
    sigma2_psi     = sigma2_out,
    eps            = eps,
    loglik         = loglik,
    n_iter         = it,
    converged      = !notConverged,
    prior_var_path = prior_var_path,
    data           = dd
  )
  class(fit_obj) <- "quasar2CR"
  fit_obj
}


#' Estimate per-SNP dispersion psi = log(1/M) by Cox-Reid APL with
#' empirical-Bayes shrinkage to a smooth trend in log baseMean.
#'
#' This is exposed separately from \code{\link{fitQuasar2CR}} so users can
#' inspect or override the dispersion step in isolation.
#'
#' @param dd       data frame with columns \code{identifier}, \code{R}, \code{A}
#'                 and any variables in \code{design}.
#' @param design   one-sided formula.  Default \code{~ 1}.
#' @param eps      shared sequencing error (scalar; not estimated here).
#' @param interval Brent interval for psi (default \code{c(-15, 15)}).
#' @param psi_init initial psi for Brent.
#' @param verbose  logical.
#' @return data frame with one row per identifier: \code{baseMean},
#'   \code{psi_gw}, \code{psi_trend}, \code{var_samp}.  Trend, prior variance
#'   and eps are attached as attributes.
#' @export
estimateDispersionsCR <- function(dd, design = ~ 1, eps = 1e-3,
                                  interval = c(-15, 15),
                                  psi_init = -log(20),
                                  verbose = TRUE) {

  stopifnot(all(c("identifier", "R", "A") %in% names(dd)))
  if (is.null(dd$offset)) dd$offset <- 0

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

  trend <- .bb_fit_trend(log(gw$baseMean), gw$psi_gw)
  gw$psi_trend <- trend$predict(log(gw$baseMean))

  resid  <- gw$psi_gw - gw$psi_trend
  sigma2 <- .bb_prior_var(resid, gw$var_samp)

  if (verbose) message("EB prior variance sigma^2_psi = ",
                       signif(sigma2, 3))

  attr(gw, "trend")    <- trend
  attr(gw, "sigma2")   <- sigma2
  attr(gw, "eps_used") <- eps
  gw
}
