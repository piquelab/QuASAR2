#' Log-likelihood for Beta-Binomial Logistic Regression to Assess ASE
#'
#' This function calculates the log-likelihood for a beta-binomial logistic regression model,
#' which is useful for detecting changes in allele-specific expression (ASE) across different
#' conditions or contexts specified by a design matrix. It adjusts for potential offsets if ASE
#' is not 0.5 under the null hypothesis, which is useful for MPRA applications. Additionally,
#' an L2 regularization term can be applied to control overfitting.
#'
#' @param b Numeric vector of coefficients.
#' @param X Design matrix of predictors.
#' @param M Concentration parameter of the beta distribution ($\alpha + \beta$).
#' @param R Number of reads matching the reference allele.
#' @param A Number of reads matching the alternate allele.
#' @param offset Adjustment for the offset in logistic regression, default is 0.0,
#'        useful when MPRA expected ASE is not 0.5.
#' @param eps Probability of base call error, default is 0.001.
#' @param lambda Regularization parameter for L2 regularization, default is 0.001.
#'
#' @return Numeric value representing the negative log-likelihood.
#' @export
#'
#' @examples
#' b <- c(0.5, -0.3)
#' X <- matrix(c(1, 2, 3, 1, 5, 2), ncol = 2)
#' M <- 10
#' R <- c(3, 6, 2)
#' A <- c(5, 5, 5)
#' logLikBetaBinomialLogistic(b, X, M, R, A)
logLikBetaBinomialLogistic <- function(b, X, M, R, A, offset = 0.0, eps = 0.001, lambda = 0.001) {
  p <- plogis(X %*% b + offset)
  p <- (p * (1 - eps) + (1 - p) * eps)
  aux <- (lgamma(M) + lgamma(R + p * M) + lgamma(A + (1 - p) * M) - lgamma(R + A + M) - lgamma(p * M) - lgamma((1 - p) * M))
  -sum(aux) + lambda * sum(b^2)
}
