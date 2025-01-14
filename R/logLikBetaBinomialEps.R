#' Log-Likelihood for Beta-Binomial Model to Optimize Epsilon Using Estimated ASE Proportion
#'
#' This function calculates the log-likelihood for a beta-binomial model, focusing on optimizing
#' the epsilon parameter. Epsilon is probability of base call error when
#' analyzing allele-specific expression (ASE). 
#'
#' @param eps Probability of base call error to be estimated
#' @param M Concentration parameter of the beta distribution ($\alpha + \beta$).
#' @param p Estimated proportion of the reference allele expression
#' @param R Number of reads matching the reference allele.
#' @param A Number of reads matching the alternate allele.
#'
#' @return Numeric value representing the negative log-likelihood.
#' @export
#'
#' @examples
#' M <- 10
#' p <- 0.5  # Estimated proportion
#' R <- 3
#' A <- 7
#' eps <- 0.001
#' logLikBetaBinomialEps(eps, M, p, R, A)
logLikBetaBinomialEps <- function(eps,M,p,R,A){
  p <- (p*(1-eps)+(1-p)*eps)
  aux <- (lgamma(M) + lgamma(R+p*M) + lgamma(A + (1-p)*M) - 
          lgamma(R+A+M) - lgamma(p*M) - lgamma((1-p)*M) )
  -sum(aux)
}

