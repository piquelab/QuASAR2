% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/logLikBetaBinomialRhoEps.R
\name{logLikBetaBinomialRhoEps}
\alias{logLikBetaBinomialRhoEps}
\title{log-likelihood of beta-binomial model}
\usage{
logLikBetaBinomialRhoEps(rho, eps, D, R, A)
}
\arguments{
\item{rho}{rho value in [0, 1]}

\item{eps}{estimate of error}

\item{D}{dispersion estimate}

\item{R}{reference read count}

\item{A}{alternate read count}
}
\value{
log likelihood
}
\description{
returns the log likelihood of reference & alternate read count data
given rho, dispersion, and error
}
