## ============================================================================
##  bbInternal.R
##
##  Internal helpers for the beta-binomial regression with reparametrized
##  precision psi = log(1/M).  None of these are exported.
##
##  Math reference: vignettes/quasar2CR.Rmd
## ============================================================================

.clip_prob <- function(p, lo = 1e-10, hi = 1 - 1e-10) {
  pmin(pmax(p, lo), hi)
}

.invlogit <- function(eta) 1 / (1 + exp(-eta))


#' Beta-binomial log-likelihood per observation (internal).
#' @keywords internal
#' @noRd
.bb_loglik_obs <- function(R, A, mu, M, eps = 0) {
  mu  <- .clip_prob(mu)
  N   <- R + A
  mut <- (1 - 2 * eps) * mu + eps
  a   <- M * mut
  b   <- M * (1 - mut)
  lchoose(N, R) +
    lgamma(R + a) + lgamma(A + b) - lgamma(N + M) -
    lgamma(a)     - lgamma(b)     + lgamma(M)
}


#' Score of the BB log-likelihood wrt beta_g (internal).
#' @keywords internal
#' @noRd
.bb_score_beta <- function(beta, R, A, X, offset, M, eps) {
  eta <- as.numeric(X %*% beta) + offset
  mu  <- .invlogit(eta)
  mut <- (1 - 2 * eps) * mu + eps
  a   <- M * mut
  b   <- M * (1 - mut)
  d_ll_d_mut <- M * (digamma(R + a) - digamma(A + b) -
                       digamma(a)    + digamma(b))
  d_mut_d_eta <- (1 - 2 * eps) * mu * (1 - mu)
  u <- d_ll_d_mut * d_mut_d_eta
  drop(crossprod(X, u))
}


#' Observed-information matrix wrt beta_g (internal).
#' @keywords internal
#' @noRd
.bb_info_beta <- function(beta, R, A, X, offset, M, eps) {
  eta <- as.numeric(X %*% beta) + offset
  mu  <- .invlogit(eta)
  mut <- (1 - 2 * eps) * mu + eps
  a   <- M * mut
  b   <- M * (1 - mut)
  neg_d2_d_mut2 <- M * M * (
    trigamma(a) + trigamma(b) -
      trigamma(R + a) - trigamma(A + b)
  )
  w <- neg_d2_d_mut2 * ((1 - 2 * eps) * mu * (1 - mu))^2
  crossprod(X * sqrt(pmax(w, 0)))
}


#' Inner per-SNP beta_g optimizer at fixed (psi, eps) (internal).
#' @keywords internal
#' @noRd
.bb_fit_beta <- function(R, A, X, offset = 0, psi, eps,
                         beta_start = NULL,
                         control = list(maxit = 200, reltol = 1e-8)) {
  M <- exp(-psi)
  p <- ncol(X)
  if (is.null(beta_start)) beta_start <- rep(0, p)

  fn <- function(b) -sum(.bb_loglik_obs(R, A,
                          .invlogit(as.numeric(X %*% b) + offset), M, eps))
  gr <- function(b) -.bb_score_beta(b, R, A, X, offset, M, eps)

  opt <- tryCatch(
    stats::optim(par = beta_start, fn = fn, gr = gr, method = "BFGS",
          control = control),
    error = function(e) NULL
  )
  if (is.null(opt) || opt$convergence != 0) {
    opt <- stats::optim(par = beta_start, fn = fn, method = "Nelder-Mead",
                 control = list(maxit = 500))
  }
  beta_hat <- opt$par
  names(beta_hat) <- colnames(X)
  info <- .bb_info_beta(beta_hat, R, A, X, offset, M, eps)
  list(beta = beta_hat, info = info, value = -opt$value,
       convergence = opt$convergence)
}


#' Cox-Reid APL at fixed psi for one SNP (internal).
#' @keywords internal
#' @noRd
.bb_apl <- function(psi, R, A, X, offset, eps, beta_start = NULL) {
  fit <- .bb_fit_beta(R, A, X, offset, psi, eps, beta_start = beta_start)
  eig <- eigen(fit$info, symmetric = TRUE, only.values = TRUE)$values
  eig <- eig[eig > 1e-12]
  apl <- fit$value - 0.5 * sum(log(eig))
  attr(apl, "fit") <- fit
  apl
}


#' Gene-wise psi via 1-D Brent max of APL_g (internal).
#' @keywords internal
#' @noRd
.bb_fit_psi_gw <- function(R, A, X, offset, eps,
                           psi_init = -log(20),
                           interval = c(-15, 15),
                           beta_start = NULL) {
  cache <- new.env(parent = emptyenv())
  cache$beta <- beta_start
  obj <- function(psi) {
    val <- .bb_apl(psi, R, A, X, offset, eps, beta_start = cache$beta)
    cache$beta <- attr(val, "fit")$beta
    -as.numeric(val)
  }
  opt <- stats::optim(par = psi_init, fn = obj, method = "Brent",
               lower = interval[1], upper = interval[2])
  fit <- .bb_fit_beta(R, A, X, offset, opt$par, eps, beta_start = cache$beta)
  list(psi = opt$par, apl = -opt$value, beta = fit$beta,
       info = fit$info, loglik = fit$value)
}


#' Sampling variance of psi_hat via central difference of APL (internal).
#' @keywords internal
#' @noRd
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


#' MAP psi with a Normal(psi_trend, sigma2) prior (internal).
#' @keywords internal
#' @noRd
.bb_fit_psi_map <- function(R, A, X, offset, eps,
                            psi_trend, sigma2,
                            psi_init = NULL,
                            interval = c(-15, 15),
                            beta_start = NULL) {
  cache <- new.env(parent = emptyenv())
  cache$beta <- beta_start
  if (is.null(psi_init)) psi_init <- psi_trend
  if (!is.finite(sigma2) || sigma2 <= 0) {
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
  opt <- stats::optim(par = psi_init, fn = obj, method = "Brent",
               lower = interval[1], upper = interval[2])
  fit <- .bb_fit_beta(R, A, X, offset, opt$par, eps, beta_start = cache$beta)
  list(psi = opt$par, beta = fit$beta, info = fit$info,
       loglik = fit$value, apl = NA_real_)
}


#' Global shared-eps update (internal).
#' @keywords internal
#' @noRd
.bb_fit_eps_global <- function(R, A, mu, M, eps_init = 1e-3,
                               interval = c(1e-6, 0.1)) {
  fn <- function(eps) -sum(.bb_loglik_obs(R, A, mu, M, eps))
  opt <- stats::optim(par = eps_init, fn = fn, method = "Brent",
               lower = interval[1], upper = interval[2])
  list(eps = opt$par, value = -opt$value, convergence = opt$convergence)
}


#' Trend fit: smooth psi_gw vs log baseMean (internal).
#' Uses locfit when available, smooth.spline otherwise.
#' @keywords internal
#' @noRd
.bb_fit_trend <- function(log_baseMean, psi_gw) {
  ok <- is.finite(psi_gw) & is.finite(log_baseMean)
  x  <- log_baseMean[ok]
  y  <- psi_gw[ok]
  if (length(x) < 5 || length(unique(x)) < 4) {
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
        predict = function(z)
          as.numeric(stats::predict(fit, newdata = data.frame(x = z))),
        fit = fit, type = "locfit"))
  }
  fit <- stats::smooth.spline(x, y, df = min(8, length(unique(x)) - 1))
  list(predict = function(z) stats::predict(fit, z)$y,
       fit = fit, type = "smooth.spline")
}


#' Robust prior variance estimator (internal).
#' @keywords internal
#' @noRd
.bb_prior_var <- function(resid, samp_var) {
  resid    <- resid[is.finite(resid)]
  samp_var <- samp_var[is.finite(samp_var) & samp_var > 0]
  if (length(resid) < 3) return(0)
  mad_var <- (stats::mad(resid, constant = 1.4826))^2
  med_var <- if (length(samp_var)) stats::median(samp_var) else 0
  max(0, mad_var - med_var)
}
