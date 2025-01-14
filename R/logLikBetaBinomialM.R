#' Log-Likelihood for Beta-Binomial Model To Estimate M Using Estimated ASE Proportion Parameter
#'
#' This function calculates the log-likelihood for a beta-binomial model, suitable for analyzing
#' allele-specific expression (ASE) where the proportion of allele expression is estimated from
#' data. It is designed to estimate the concentration parameter and using best fit and the 
#' the likelihood of observing the given number of reference and alternate allele reads. 
#'
#' @param M Concentration parameter of the beta distribution ($\alpha + \beta$)
#' @param p Estimated proportion of the reference allele expression, as derived from best fit.
#' @param R Number of reads matching the reference allele.
#' @param A Number of reads matching the alternate allele.
#' @param eps Probability of base call error, default is 0.001.
#'
#' @return Numeric value representing the negative log-likelihood.
#' @export
#'
#' @examples
#' M <- 10
#' p <- 0.5  # Best fit proportion estimate
#' R <- 3
#' A <- 7
#' logLikBetaBinomialM(M, p, R, A)
logLikBetaBinomialM <- function(M, p, R, A, eps = 0.001) {
  p <- (p * (1 - eps) + (1 - p) * eps)
  aux <- (lgamma(M) + lgamma(R + p * M) + lgamma(A + (1 - p) * M) - 
            lgamma(R + A + M) - lgamma(p * M) - lgamma((1 - p) * M))
  -sum(aux)
}