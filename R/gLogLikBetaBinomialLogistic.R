#' Gradient of Log-Likelihood for Beta-Binomial Logistic Regression to Assess ASE
#'
#' This function computes the gradient of the log-likelihood for a beta-binomial
#' logistic regression model, useful for detecting allele-specific expression (ASE) changes
#' across different conditions or contexts. This function is often used for optimization purposes
#' in maximum likelihood estimation, helping to find the parameters b values that maximize
#' the log-likelihood function. It also adjusts for potential offsets and includes an L2 
#' regularization term.
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
#' @return A numeric vector representing the gradient of the negative log-likelihood.
#' @export
#'
#' @examples
#' b <- c(0.5, -0.3)
#' X <- matrix(c(1, 2, 3, 1, 5, 2), ncol = 2)
#' M <- 10
#' R <- c(3, 6, 2)
#' A <- c(5, 5, 5)
#' gLogLikBetaBinomialLogistic(b, X, M, R, A)
gLogLikBetaBinomialLogistic <- function(b, X, M, R, A, offset = 0.0, eps = 0.001, lambda = 0.001) {
  p <- plogis(X %*% b + offset)
  p <- (p * (1 - eps) + (1 - p) * eps)
  dp <- (1 - 2 * eps) * (p * (1 - p)) * M
  dlogL <- dp * (digamma(R + p * M) - digamma(A + (1 - p) * M) - digamma(p * M) + digamma((1 - p) * M))
  grad <- -t(X) %*% dlogL + lambda * b
  return(grad)
}


