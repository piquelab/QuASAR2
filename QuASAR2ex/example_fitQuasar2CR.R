## =============================================================================
##  example_fitQuasar2CR.R
##
##  Worked example & sanity check for fitQuasar2CR().
##  Run with:    Rscript example_fitQuasar2CR.R
##
##  Demonstrates:
##    - simulating beta-binomial ASE-like data with batch + treatment
##    - fitting with fitQuasar2CR(dat, ~ batch + trt)
##    - recovering coefficients (testCoef)
##    - building & testing arbitrary contrasts (makeContrastsQuasar, testContrast)
##    - LRT vs a reduced model (testLRT)
## =============================================================================

## Sanity check: simulate beta-binomial ASE-like data and fit it
suppressPackageStartupMessages({
  source("fitQuasar2CR.R")
})

set.seed(20260512)

## ---------- 1. Simulate ----------
G   <- 10000              # SNPs
n_g <- 8               # replicates per SNP (4 in treatment A, 4 in B)
batch <- factor(rep(c("X","Y"), each = 4))
trt   <- factor(rep(c("A","A","B","B"), times = 2))
M_true_global <- 30    # base precision

dat <- do.call(rbind, lapply(seq_len(G), function(g) {
  ## per-SNP true coefficients
  b0 <- 0                                    # intercept (logit p = 0 -> p = 0.5)
  b_batch <- if (g %% 2 == 0) 0.15 else -0.1
  b_trt   <- if (g <= 100) 1.2 else if (g <= 200) -0.8 else rnorm(1, 0, 0.005)

  eta <- b0 + b_batch * (batch == "Y") + b_trt * (trt == "B")
  mu  <- 1 / (1 + exp(-eta))

  ## per-SNP variation in M
  M_g  <- M_true_global * exp(rnorm(1, 0, 0.4))

  ## draw beta-binomial: p ~ Beta(M*mu, M*(1-mu)); R ~ Binom(N, p)
  N    <- sample(20:200, n_g, replace = TRUE)
  p_g  <- rbeta(n_g, M_g * mu, M_g * (1 - mu))
  R    <- rbinom(n_g, N, p_g)
  A    <- N - R
  data.frame(
    identifier = sprintf("snp_%03d", g),
    R = R, A = A, batch = batch, trt = trt,
    true_b_trt   = b_trt,
    true_b_batch = b_batch,
    true_M       = M_g,
    stringsAsFactors = FALSE
  )
}))
cat("Simulated dataset:", nrow(dat), "rows /", G, "SNPs\n")

## ---------- 2. Fit ----------
fit <- fitQuasar2CR(dat, design = ~ batch + trt,
                    max_iter = 6, verbose = TRUE)

print(fit)


resObj <- fitQuasar2(dat,~ batch + trt)

aseInt <- resObj$results %>% dplyr::filter(term=="(Intercept)") ##%>% left_join(dd)
aseTr <- resObj$results %>% dplyr::filter(term=="trtB") ##%>% left_join(dd)

Mvec <- resObj$Mvec

dd <- dat %>% 
  group_by(identifier) %>%
  mutate(MN = mean(R + A))

nbreaks=20
MN <- dd %>% summarize(MN=mean(R+A)) 
cov_breaks <- unique(c(0,quantile(MN$MN,(1:nbreaks)/nbreaks)))
dd <- dd %>% mutate(bin=cut(MN,cov_breaks)) 
binlab <- levels(dd$bin)

dd$estM <- Mvec[dd$bin]

plot(log(dd$true_M),log(dd$estM))
abline(0,1)

plot((dd$true_M),(dd$estM))
abline(0,1)

cat("\nCorrelation log M_hat vs log M_true: ",
    round(cor(log(dd$true_M), log(dd$estM)), 3), "\n")




## ---------- 3. Recovery of trt effect ----------
res <- testCoef(fit, coef = "trtB")
true_trt <- unique(dat[, c("identifier","true_b_trt")])
res2 <- merge(res, true_trt, by = "identifier")
cat("\nCorrelation between estimated and true trt effect: ",
    round(cor(res2$estimate, res2$true_b_trt), 3), "\n")
cat("RMSE:                                              ",
    round(sqrt(mean((res2$estimate - res2$true_b_trt)^2)), 3), "\n")

plot( res2$true_b_trt, res2$estimate)
abline(0,1)

plot(res2$true_b_trt,aseTr$coeff)
abline(0,1)


plot(res2$estimate,aseTr$coeff)
abline(0,1)


res3 <- merge(res,aseTr) %>% merge(true_trt) %>% mutate(myP=abs(true_b_trt)>0.5)

res3$padj2 <- p.adjust(res3$pval)

library(tidyverse)



res3 %>% arrange(myP) %>%
ggplot(aes(x=-log10(p.value),y=-log10(pval),color=myP)) + 
  geom_point() +
  geom_abline(slope=1,intercept=0,color="red") + 
  theme_bw()

table(res3$myP,res3$padj<0.1)

table(res3$myP,res3$padj2<0.1)

table(res3$myP,res3$padj<0.01)


## ---------- 4. Recovery of dispersion ----------
disp <- fit$dispersion
true_disp <- unique(dat[, c("identifier", "true_M")])
disp2 <- merge(disp, true_disp, by = "identifier")
cat("\nCorrelation log M_hat vs log M_true: ",
    round(cor(log(disp2$M), log(disp2$true_M)), 3), "\n")

plot(log(disp2$true_M),log(disp2$M))
abline(0,1)



## ---------- 5. LRT against ~ batch alone ----------
lrt <- testLRT(fit, reduced = ~ batch)
cat("\nLRT (drop trt): how many SNPs with padj < 0.05?",
    sum(lrt$padj < 0.05, na.rm = TRUE), "\n")
cat("                            of the 200 truly DE:",
    sum(lrt$padj[as.numeric(sub("snp_", "", lrt$identifier)) <= 200] < 0.05,
        na.rm = TRUE), "\n")


sum(lrt$padj < 0.05, na.rm = TRUE)

## ---------- 6. Contrast: trtB - 0 (same as testCoef, just to exercise machinery) ----------
L <- makeContrastsQuasar("trtB", fit)
print(L)
ctrast <- testContrast(fit, L)
cat("\nWald chi2 from contrast vs testCoef^2 (first 5):\n")
print(head(data.frame(coefSq = res2$statistic^2, contrast = ctrast$chi2), 5))

cat("\nALL TESTS RAN OK.\n")
