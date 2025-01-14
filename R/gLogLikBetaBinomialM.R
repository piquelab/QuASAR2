#' Gradient of Log-Likelihood for Beta-Binomial Model M Estimation Using Estimated ASE Proportion
#'
#' This function calculates the gradient of the log-likelihood for a beta-binomial model,
#' which is useful for optimization algorithms aiming to find the best parameter estimates for M.
#' It is designed to work with allele-specific expression (ASE) data where the proportion
#' of allele expression has been estimated. The function helps in adjusting the concentration
#' parameter M based on the gradient descent method or similar optimization approaches.
#'
#' @param M Concentration parameter of the beta distribution ($\alpha + \beta$).
#' @param p Estimated proportion of the reference allele expression, as derived from best fit.
#' @param R Number of reads matching the reference allele.
#' @param A Number of reads matching the alternate allele.
#' @param eps Probability of base call error, default is 0.001.
#'
#' @return Numeric vector representing the gradient of the negative log-likelihood.
#' @export
#'
#' @examples
#' M <- 10
#' p <- 0.5  # Best fit proportion estimate
#' R <- 3
#' A <- 7
#' gLogLikBetaBinomialM(M, p, R, A)
gLogLikBetaBinomialM <- function(M,p,R,A,eps=0.001){
  p <- (p*(1-eps)+(1-p)*eps)
  aux <- (digamma(M) + p * digamma(R + p * M) + (1 - p) * digamma(A + (1 - p) * M) - 
            digamma(R + A + M) - p * digamma(p * M) - (1 - p) * digamma((1 - p) * M ))
  -sum(aux)
}