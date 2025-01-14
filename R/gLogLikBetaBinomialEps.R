#' Gradient of Log-Likelihood for Beta-Binomial Model to Optimize Epsilon Using Estimated ASE Proportion
#'
#' This function calculates the gradient of the log-likelihood for a beta-binomial model, with a
#' focus on optimizing the epsilon parameter. It is useful in optimization algorithms that require 
#' gradient information.
#'
#' @param eps Probability of base call error to be estimated
#' @param M Concentration parameter of the beta distribution ($\alpha + \beta$).
#' @param p Estimated proportion of the reference allele expression
#' @param R Number of reads matching the reference allele.
#' @param A Number of reads matching the alternate allele.
#'
#' @return Gradient of the negative log-likelihood with respect to epsilon.
#' @export
#'
#' @examples
#' M <- 10
#' p <- 0.5  # Estimated proportion
#' R <- 3
#' A <- 7
#' eps <- 0.001
#' gLogLikBetaBinomialEps(eps, M, p, R, A)
gLogLikBetaBinomialEps <- function(eps,M,p,R,A){
  p <- (p*(1-eps)+(1-p)*eps)
  aux <- (1-2*p) * M * (digamma(R+p*M) - digamma(A + (1-p)*M) - digamma(p*M) + digamma((1-p)*M)) 
  -sum(aux)
}
