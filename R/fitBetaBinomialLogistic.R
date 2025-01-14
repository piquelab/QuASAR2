#' Fit Beta-Binomial Logistic Regression Model
#'
#' This function fits a beta-binomial logistic regression model to the data based on a specified
#' model formula. It is ideal for cases with allele-specific expression (ASE) data that exhibits
#' overdispersion relative to a standard binomial model. The function employs gradient descent 
#' optimization techniques to estimate the best fitting model parameters.
#'
#' @param dd Data frame containing the data to be modeled. The data frame should include the
#'           predictor variables specified in the model formula, along with the following columns:
#'           - `M`: Concentration parameter of the beta distribution ($\alpha + \beta$).
#'           - `R`: Number of successes (e.g., reads matching the reference allele).
#'           - `A`: Number of failures (e.g., reads matching the alternate allele).
#'           - `offset`: Offset value if applicable, default is 0. Useful for MPRA (`qlogis(DNA_prop)`)
#' @param model Formula specifying the model to be fitted, like ~ Batch + Treatment. Variables have to
#' present in the `dd` object. 
#' @param b Initial values for the regression coefficients. If NULL, they will be initialized to zero.
#' @param eps Small constant added to probabilities to prevent them from being 0 or 1, default is 0.001.
#'
#' @return A list containing two elements: `res`, a tibble with the fitted coefficients, standard errors,
#'         statistical significance, and other diagnostic measures; and `dd`, the modified data frame
#'         with additional columns for fitted probabilities.
#' @export
#'
#' @examples
#' # Assuming `mydat` is a data frame with required columns and model predictors
#' model_formula <- ~ batch + Tr
#' initial_b <- rep(0, length = ncol(model.matrix(model_formula, mydat)))
#' fitted_model <- fitBetaBinomialLogistic(mydat, model_formula, b = initial_b)
#'
#' @seealso
#' \code{\link{logLikBetaBinomialLogistic}}: Log-likelihood function for beta-binomial logistic regression.
#' \code{\link{gLogLikBetaBinomialLogistic}}: Gradient of the log-likelihood function for beta-binomial logistic regression.
fitBetaBinomialLogistic <- function(dd, model, b = NULL, eps = 0.001) {
  X <- model.matrix(model, dd)
  if (is.null(b)) {
    b <- rep(0, ncol(X))
  }
  stopifnot(ncol(X) == length(b))
  auxLogis <- optim(b,
                    fn = logLikBetaBinomialLogistic,
                    gr = gLogLikBetaBinomialLogistic,
                    X = X,
                    M = dd$M, 
                    R = dd$R,
                    A = dd$A,
                    offset = ifelse("offset" %in% names(dd), dd$offset, 0),
                    eps = eps,
                    method = "L-BFGS-B",
                    hessian = TRUE,
                    lower = -8,
                    upper = 8)
  ##auxLogis$convergence
  ##auxLogis$message
  
  dd$p.fit <- plogis(X %*% auxLogis$par + ifelse("offset" %in% names(dd), dd$offset, 0))
  
  res <- tibble(
    term = colnames(X),
    coeff = auxLogis$par,
    se = 1 / diag(auxLogis$hessian)^0.5,
    convergence = auxLogis$convergence,
    qr.rank = qr(X)$rank,
    df = ncol(X),
    n = nrow(X)
  ) %>% 
    mutate(stat = coeff / se,
           pval = 2 * pnorm(-abs(stat)))
  
  list(res = res, dd = dd)
}