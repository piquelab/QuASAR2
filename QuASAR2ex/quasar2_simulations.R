# QuASAR2_simulations.R script
# Purpose: run simulations for different ASE scenarios to see how calibrated QuASAR2 is
# Author: Shreya Nirmalan
# Date Last Edited: 1-5-2025
suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(qqman)
  library(QuASAR2)
})

sim_quasar2_df <- function(
  n_snps = 20000, # Number of SNPs tested
  n_ctrl = 5, # Number of control samples
  n_trt  = 5, # Number of treatment samples
  N = 100,        # Read coverage for all SNPs × samples
  M = 80,         # beta-binomial concentration (overdispersion)
  # fractions of SNP types
  frac_ASE_only  = 0.0, # Fraction of SNPs with ASE
  frac_cASE_only = 0.0, # Fraction of SNPs with cASE (treatment v control)
  frac_both      = 0.0, # Fraction of SNPs with ASE in both conditions and cASE
  # proportions
  delta_ASE  = 0.10,   # ASE SNPs: ref proportion = 0.5 ± delta_ASE for instances of ASE
  delta_cASE = 0.10,   # cASE SNPs: ref proportion = 0.5 ± delta_cASE in treatment only samples for instances of cASE
  seed = 1 # for reproducibility
) {
  set.seed(seed)

  ## sample covariates
  n_samp <- n_ctrl + n_trt
  samples <- data.frame(
    SampleID   = paste0("S", seq_len(n_samp)),
    Treatment = factor(
      c(rep("control", n_ctrl), rep("treatment", n_trt)),
      levels = c("control", "treatment")
    )
  )

  ## SNP truth assignment
  snps <- paste0("rsSIM", seq_len(n_snps))
  class <- rep("null", n_snps)

  n_both <- round(frac_both * n_snps)
  if (n_both > 0) class[sample.int(n_snps, n_both)] <- "both"

  remaining <- which(class == "null")
  n_ase <- round(frac_ASE_only * n_snps)
  if (n_ase > 0) class[sample(remaining, min(n_ase, length(remaining)))] <- "ASE_only"

  remaining <- which(class == "null")
  n_case <- round(frac_cASE_only * n_snps)
  if (n_case > 0) class[sample(remaining, min(n_case, length(remaining)))] <- "cASE_only"

  truth <- data.frame(
    identifier = snps,
    class = class,
    is_ASE  = class %in% c("ASE_only", "both"),
    is_cASE = class %in% c("cASE_only", "both")
  )

  ## randomly assigns if ref or alt allele would have the higher counts
  dir <- sample(c(-1, 1), n_snps, replace = TRUE)
  truth$direction <- dir

  ## make dataframe input forquasar2 with sample and snp information
  
  dd <- merge(
    truth,
    samples,
    all = TRUE
  )


  snp_idx <- match(dd$identifier, snps)
  is_treated <- dd$Treatment == "treatment"

  ## set mean ref proportion p

  p <- rep(0.5, nrow(dd))   # null baseline

  # ASE only: imbalance in ALL samples
  ase_rows <- dd$is_ASE
  p[ase_rows] <- 0.5 + dd$direction[ase_rows] * delta_ASE

  # cASE only: imbalance ONLY in treated
  case_rows <- dd$is_cASE & is_treated
  p[case_rows] <- 0.5 + dd$direction[case_rows] * delta_cASE

  # both: baseline ASE + extra treated shift
  both_all <- dd$class == "both"
  p[both_all] <- 0.5 + dd$direction[both_all] * delta_ASE
  both_trt <- dd$class == "both" & is_treated
  p[both_trt] <- p[both_trt] + dd$direction[both_trt] * delta_cASE

  # ensure proportion never becomes 0 or 1 to ensure rbeta works correctly
  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  dd$prop_true <- p


  ## fixed coverage + sampling
  dd$N <- N[match(dd$identifier, snps)]              # one N per SNP, broadcast correctly
  theta <- rbeta(nrow(dd), p * M, (1 - p) * M)
  dd$R  <- rbinom(nrow(dd), size = dd$N, prob = theta)
  dd$A  <- dd$N - dd$R
  
  
  #mean coverage per SNP
  dd$mean_RA <- ave(dd$R + dd$A, dd$identifier, FUN = mean)

  dd[order(dd$identifier, dd$SampleID), ]
}

Nvec=round(runif(20000,60,1000))

##Nvec=100

dd <- sim_quasar2_df(n_snps = 20000, n_ctrl = 5, n_trt = 5, M = 200, N = Nvec, seed=1,frac_ASE_only=0.05,frac_cASE_only=0.05)

## ---------- 2. Fit ----------
fit <- fitQuasar2CR(dd, design = ~ Treatment,
                    max_iter = 6, verbose = TRUE)

print(fit)

## ---------- 3. Recovery of trt effect ----------
resCR.Tr <- testCoef(fit, coef = "Treatmenttreatment",df_method = "moderated")

resCR.Int <- testCoef(fit, coef = "(Intercept)",df_method = "moderated")


sum(resCR.Int$padj<0.1)

sum(resCR.Tr$padj<0.1)

qq(resCR.Tr$p.value)

######
resObj <- fitQuasar2(dd,~Treatment)


head(resObj)

aseInt <- resObj$results %>% dplyr::filter(term=="(Intercept)") ##%>% left_join(dd)

aseTr <- resObj$results %>% dplyr::filter(term=="Treatmenttreatment") ##%>% left_join(dd)

aseInt$padj <- p.adjust(aseInt$pval)

aseTr$padj <- p.adjust(aseTr$pval)

sum(aseTr$padj<0.1)

sum(aseInt$padj<0.1)

qq1 <- qq(aseInt$pval)

qq2 <- qq(aseTr$pval)



## Using exact M just in fitBetaBinomialLogistic. 

aux <- dd %>% mutate(M=200) %>% group_by(identifier) %>% nest()

res <- aux %>% mutate(res=map(data,fitBetaBinomialLogistic,~Treatment, eps=0.001))

res <- res %>% mutate(data=map(res, "dd"),res=map(res,"res"))

res2 <- res %>% select(identifier,res) %>% unnest(cols=c(res))

aseInt3 <- res2 %>% dplyr::filter(term=="(Intercept)") ##%>% left_join(dd)
aseTr3 <- res2 %>% dplyr::filter(term=="Treatmenttreatment") ##%>% left_join(dd)



qq(aseInt3$pval)

qq(aseTr3$pval)

plot(-log10(aseInt$pval),-log10(aseInt3$pval))
abline(0,1)


plot(-log10(resCR.Int$p.value),-log10(aseInt3$pval))
abline(0,1)

plot(-log10(aseTr$pval),-log10(aseTr3$pval))
abline(0,1)


plot(-log10(resCR.Tr$p.value),-log10(aseTr3$pval))
abline(0,1)

disp <- fit$dispersion

plot(disp$baseMean,disp$M)

hist(disp$M[disp$baseMean>400],breaks=100)


  mean(dd$R == 0 | dd$A == 0)
# [1] 0
sum(dd$R == 0 | dd$A == 0)
# [1] 0


debugonce(fitQuasar2)
resObj <- fitQuasar2(dd, ~Treatment)