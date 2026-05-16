## tests/testthat/test-bb-internals.R
## Validates analytic score and information matrix against numerical
## derivatives at a few SNPs.

context("BB log-likelihood internals")

test_that(".bb_loglik_obs returns numeric vector of expected length", {
  R <- c(10, 20, 5); A <- c(15, 10, 25)
  mu <- c(0.3, 0.5, 0.2); M <- 30; eps <- 1e-3
  ll <- QuASAR2:::.bb_loglik_obs(R, A, mu, M, eps)
  expect_length(ll, 3)
  expect_true(all(is.finite(ll)))
})

test_that(".bb_score_beta matches numerical gradient", {
  skip_if_not_installed("numDeriv")
  set.seed(11)
  n <- 12
  X <- cbind(1, rep(0:1, length.out = n))
  beta <- c(0.1, 0.4)
  N <- rep(80, n)
  R <- rbinom(n, N, plogis(as.numeric(X %*% beta)))
  A <- N - R
  M <- 30; eps <- 1e-3
  fn <- function(b) sum(QuASAR2:::.bb_loglik_obs(
                          R, A, plogis(as.numeric(X %*% b)), M, eps))
  num_grad <- numDeriv::grad(fn, beta)
  ana_grad <- QuASAR2:::.bb_score_beta(beta, R, A, X, 0, M, eps)
  expect_equal(num_grad, ana_grad, tolerance = 1e-5)
})

test_that(".bb_info_beta matches numerical Hessian to 3+ digits", {
  skip_if_not_installed("numDeriv")
  set.seed(22)
  n <- 16
  X <- cbind(1, rep(0:1, length.out = n))
  beta <- c(0.0, 0.5)
  N <- rep(120, n)
  R <- rbinom(n, N, plogis(as.numeric(X %*% beta)))
  A <- N - R
  for (M in c(10, 50, 200)) {
    fn <- function(b) -sum(QuASAR2:::.bb_loglik_obs(
                            R, A, plogis(as.numeric(X %*% b)), M, eps = 1e-3))
    num_H <- numDeriv::hessian(fn, beta)
    ana_I <- QuASAR2:::.bb_info_beta(beta, R, A, X, 0, M, eps = 1e-3)
    expect_equal(det(num_H), det(ana_I), tolerance = 0.03,
                 info = paste("M =", M))
  }
})

test_that(".bb_fit_beta converges from poor starts", {
  set.seed(31)
  n <- 10
  X <- cbind(1, rep(0:1, length.out = n))
  beta_true <- c(0.2, 1.0)
  N <- rep(100, n)
  R <- rbinom(n, N, plogis(as.numeric(X %*% beta_true)))
  A <- N - R
  for (start in list(c(0,0), c(-3,3), c(5,-5))) {
    fit <- QuASAR2:::.bb_fit_beta(R, A, X, 0, psi = -log(50),
                                  eps = 1e-3, beta_start = start)
    expect_true(abs(fit$beta[2] - beta_true[2]) < 0.3,
                info = paste("start =", paste(start, collapse = ",")))
  }
})
