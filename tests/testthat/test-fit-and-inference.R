## tests/testthat/test-fit-and-inference.R

context("fitQuasar2CR end-to-end")

make_bb_data <- function(G = 30, n_g = 8, b_trt = 1.0, M_mid = 50, seed = 1L) {
  set.seed(seed)
  trt <- factor(rep(c("A","A","A","A","B","B","B","B"), length.out = n_g))
  M_g <- exp(runif(G, log(M_mid/3), log(M_mid*3)))
  do.call(rbind, lapply(seq_len(G), function(g) {
    mu <- 1 / (1 + exp(-(b_trt * (trt == "B"))))
    pp <- rbeta(n_g, M_g[g]*mu, M_g[g]*(1-mu))
    N  <- sample(40:150, n_g, replace = TRUE)
    R  <- rbinom(n_g, N, pp); A <- N - R
    data.frame(identifier = sprintf("snp_%03d", g),
               R = R, A = A, trt = trt, stringsAsFactors = FALSE)
  }))
}

test_that("fitQuasar2CR returns a well-formed object", {
  dat <- make_bb_data(G = 20, n_g = 8, b_trt = 1.0)
  fit <- fitQuasar2CR(dat, design = ~ trt, max_iter = 3, verbose = FALSE)
  expect_s3_class(fit, "quasar2CR")
  for (nm in c("coefficients", "vcov", "designs", "dispersion",
               "sigma2_psi", "eps", "loglik"))
    expect_true(nm %in% names(fit), info = nm)
  expect_length(fit$coefficients, 20)
  expect_true(all(c("psi", "M", "psi_gw", "psi_trend", "var_samp",
                    "baseMean") %in% names(fit$dispersion)))
  expect_true(is.finite(fit$sigma2_psi))
})

test_that("coef() and vcov() S3 methods work", {
  dat <- make_bb_data(G = 10, n_g = 8)
  fit <- fitQuasar2CR(dat, design = ~ trt, max_iter = 2, verbose = FALSE)
  one <- names(fit$coefficients)[1]
  expect_named(coef(fit, one), c("(Intercept)", "trtB"))
  V <- vcov(fit, one)
  expect_equal(dim(V), c(2, 2))
  expect_equal(rownames(V), colnames(V))
})

test_that("testCoef gives nominal Type-I on a small null sim", {
  ## Pure null, 250 SNPs.  At alpha=0.05 a binomial(250, 0.05) has SE = 0.014;
  ## empirical rate should be within 0.02-0.10 (~3 SE wide).
  dat <- make_bb_data(G = 250, n_g = 8, b_trt = 0, seed = 99)
  fit <- fitQuasar2CR(dat, design = ~ trt, max_iter = 3, verbose = FALSE)
  p <- testCoef(fit, "trtB")$p.value
  ti <- mean(p <= 0.05, na.rm = TRUE)
  expect_gt(ti, 0.02)
  expect_lt(ti, 0.10)
})

test_that("testLRT against ~ 1 has F-statistic squared agreeing with t-stat", {
  ## For df1 = 1, the F-LRT statistic should be ~equal to t-Wald^2.
  dat <- make_bb_data(G = 15, n_g = 10, b_trt = 1.0, seed = 5)
  fit <- fitQuasar2CR(dat, design = ~ trt, max_iter = 2, verbose = FALSE)
  wald <- testCoef(fit, "trtB")
  lrt  <- testLRT(fit, reduced = ~ 1)
  ## They should be highly correlated (F = chi^2 / 1, ~ t^2 at large n)
  r <- cor(wald$statistic^2, lrt$statistic, use = "complete.obs")
  expect_gt(r, 0.9)
})

test_that("makeContrastsQuasar builds the right L matrix", {
  dat <- make_bb_data(G = 5, n_g = 8, b_trt = 0.5)
  fit <- fitQuasar2CR(dat, design = ~ trt, max_iter = 2, verbose = FALSE)
  L <- makeContrastsQuasar("trtB", fit)
  expect_equal(dim(L), c(1, 2))
  expect_equal(unname(L[1, ]), c(0, 1))
  expect_equal(colnames(L), c("(Intercept)", "trtB"))
})

test_that("df_method='moderated' adds 1/sigma2_psi to df", {
  dat <- make_bb_data(G = 50, n_g = 8, b_trt = 0, seed = 7)
  fit <- fitQuasar2CR(dat, design = ~ trt, max_iter = 3, verbose = FALSE)
  s2  <- fit$sigma2_psi
  if (is.finite(s2) && s2 > 0) {
    r1 <- testCoef(fit, "trtB", df_method = "residual")
    r2 <- testCoef(fit, "trtB", df_method = "moderated")
    expect_equal(r2$df - r1$df, rep(1 / s2, nrow(r1)), tolerance = 1e-8)
  }
})

test_that("ref='normal' gives smaller p-values than ref='t' for the same fit", {
  dat <- make_bb_data(G = 20, n_g = 6, b_trt = 1.0)
  fit <- fitQuasar2CR(dat, design = ~ trt, max_iter = 2, verbose = FALSE)
  pt <- testCoef(fit, "trtB", ref = "t")$p.value
  pn <- testCoef(fit, "trtB", ref = "normal")$p.value
  ok <- is.finite(pt) & is.finite(pn) & pt > 0 & pn > 0
  expect_true(all(pn[ok] <= pt[ok] + 1e-12))
})
