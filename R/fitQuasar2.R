#' fitQuasar2 Fits Beta-Binomial Logistic Regression Model on multiple groups of SNP and conditions
#'
#' This function fits a beta-binomial logistic regression model to the data based on a specified
#' model formula. It is ideal for cases with allele-specific expression (ASE) data that exhibits
#' overdispersion relative to a standard binomial model. The function employs gradient descent 
#' optimization techniques to estimate the best fitting model parameters of the logistic model. 
#' This function is designed to be used on data that has been grouped by SNP:individual sharing
#' the overdispersion parameter M, which is fit on multiple quantile breaks on number of reads. 
#'
#' @param dd Data frame containing the data to be modeled. The data frame should include the
#'           predictor variables specified in the model formula, along with the following columns:
#'           - `identifier` : the betabinomial model will be fit based on grouping variable. 
#'           For MPRA experiments this would be a construct ID, for ASE/cASE on individual SNPs it 
#'           will be the SNP id combined with an individual ID. It is important NOT to mix individuals
#'           unless some phasing is done beforehand flipping R/A accordingly. 
#'           - `R`: Number of successes (e.g., reads matching the reference allele).
#'           - `A`: Number of failures (e.g., reads matching the alternate allele).
#'           - `offset`: Offset value if applicable, default is 0. Useful for MPRA (`qlogis(DNA_prop)`)
#'              This is optional, and it does not need to be set. 
#' @param model Formula specifying the model to be fitted, like ~ Batch + Treatment. Variables have to
#' present in the `dd` object as additional columns and properly formated as factors or continuous 
#' variables.
#' @param M Initial values for the overdispersion parameter (default 20)
#' @param nbreaks Number of breaks to divide the overall expression mean of the `identifier` with `R+A`
#' which is the baseMean. For each break an overdispersion parameter M will be calculated. 
#' @param eps Initial value to account for sequencing errors. Small constant added to probabilities 
#' to prevent them from being 0 or 1, default is 0.001.
#'
#' @return A list containing two elements: `res`, a tibble with the fitted coefficients, standard errors,
#'         statistical significance, and other diagnostic measures; and `dd`, the modified data frame
#'         with additional columns for fitted probabilities.
#' @export
#' @examples
#' \dontrun{
#'   # Assuming `mydat` is a data frame with required columns and model predictors
#'   model_formula <- ~ batch + Tr
#'   initial_b <- rep(0, length = ncol(model.matrix(model_formula, mydat)))
#'   fitted_model <- fitBetaBinomialLogistic(mydat, model_formula, b = initial_b)
#' }
#' @seealso
#' \code{\link{fitBetaBinomialLogistic}}: Log-likelihood function for beta-binomial logistic regression.

fitQuasar2 <- function(dd, model, nbreaks = 20, M = 20, eps = 0.001, 
                       tolM = 0.01, tolCoeff = 0.02, tolEps = 0.0002, maxIt = 100) {
  if (!all(c("identifier", "R", "A") %in% names(dd))) stop("Required columns missing from data frame")
  
  dd <- dd %>% 
    group_by(identifier) %>%
    mutate(MN = mean(R + A))
  
  MN <- dd %>% summarize(MN=mean(R+A)) 
  cov_breaks <- unique(c(0,quantile(MN$MN,(1:nbreaks)/nbreaks)))
  dd <- dd %>% mutate(bin=cut(MN,cov_breaks)) 
  binlab <- levels(dd$bin)
  
  neweps <- eps
  
  dd$M <- M
  Mvec <- rep(M,length(binlab))
  names(Mvec)<- binlab

  cat("Intial model fitting with default M = ",M, "and eps = ",eps,"\n")
  ## Breaks data frame into the groups and fits a first iteration of the BetaBinomialLogistic with default M and eps.  
  aux <- dd %>% group_by(identifier) %>% nest()
  res <- aux %>% mutate(res=map(data,fitBetaBinomialLogistic,model, eps=neweps))
  
  
  res <- res %>% mutate(data=map(res, "dd"),res=map(res,"res"))
  
  res2 <- res %>% select(identifier,res) %>% unnest(cols=c(res))
  res3 <- res %>% select(identifier, data) %>% unnest(cols=c(data))
  
  iteration=1
  notConverged=TRUE
  while(notConverged & iteration < maxIt){
    cat("Iteration: ",iteration,"\n")
    cat("Optimizing M\n")
    ## Consider to separate to a different function. 
    Mres <- sapply(binlab,function(bin){
      selv <- res3$bin==bin;
      auxLogis <- optim(par=Mvec[bin],
                        fn=logLikBetaBinomialM,
                        gr=gLogLikBetaBinomialM,
                        p=res3$p.fit[selv],
                        R=res3$R[selv],
                        A=res3$A[selv],
                        method="L-BFGS-B",#control=c(trace=6),
                        hessian=TRUE,
                        lower=0.1,
                        upper=10000)
      data.frame(M=auxLogis$par,Mse=1/(auxLogis$hessian)^.5,conv=auxLogis$convergence)
    })
    Mres <- t(Mres)
    OldMvec <- Mvec
    Mvec <- unlist(Mres[,1])
    cat("M change max:",max(abs(Mvec-OldMvec)),"\n")
    
    ## Put updated M in the data frame and optimaize again coefficients for the logistic part
    res3$M <- Mvec[res3$bin]
    
    
    cat("Optimizing eps\n")
    oldeps=neweps
    auxLogis <- optim(par=neweps,
                      fn=logLikBetaBinomialEps,
                      gr=gLogLikBetaBinomialEps,
                      p=res3$p.fit,
                      R=res3$R,
                      A=res3$A,
                      M=res3$M,
                      method="L-BFGS-B",#control=c(trace=6),
                      hessian=TRUE,
                      lower=1E-3,
                      upper=0.1)
    auxLogis$convergence
    auxLogis$message
    c(auxLogis$par,1/(auxLogis$hessian)^.5)
    neweps=auxLogis$par
    cat("eps change max:",max(abs(neweps-oldeps)),"\n")
    
    
    cat("Model fitting with eps = ",neweps,"\n")
    aux <- res3 %>% group_by(identifier) %>% nest()
    aux$res <- res$res
    
    res <- aux %>% mutate(res=map2(data,res, ~fitBetaBinomialLogistic(.x,model=model, b=.y$coeff,eps=neweps,control=list(factr=1E10)))) 
    
    res <- res %>% mutate(data=map(res, "dd"),res=map(res,"res"))
    
    res3 <- res %>% select(identifier, data) %>% 
      unnest(cols=c(data))
    
    oldcoeff <- res2$coeff
    res2 <- res %>% select(identifier,res) %>% unnest(cols=c(res))
    
    cat("Max coeff change",max(abs(oldcoeff-res2$coeff)),"\n")
    
    if(max(abs(Mvec - OldMvec)) < tolM & 
       max(abs(oldcoeff - res2$coeff)) < tolCoeff & 
       max(abs(neweps - oldeps)) < tolEps & 
       iteration > 1){
      notConverged = FALSE
      }
    
    iteration=iteration+1;
    }
  
  return(list(
    results = res2,
    Mres = Mres,
    Mvec = Mvec,
    epsilon = neweps,
    iterations = iteration - 1,
    converged = !notConverged
  ))

  }
