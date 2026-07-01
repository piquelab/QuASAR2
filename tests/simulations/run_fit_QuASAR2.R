#!/usr/bin/env Rscript
# run_fit_QuASAR2.R
# Fits fitQuasar2 (binned shared overdispersion) on a saved simulated dataset.
#
# Usage:
#   Rscript run_fit_QuASAR2.R <rds_path> <out_dir>


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("Usage: Rscript run_fit_QuASAR2.R <rds_path> <out_dir>")

rds_path <- args[1]
out_dir  <- args[2]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

devtools::load_all("/rs/rs_grp_scaipgenetic/QuASAR2")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(tibble)
})

# ---- load data ---------------------------------------------------------------
obj        <- readRDS(rds_path)
dd         <- obj$dd
truth      <- obj$truth        # pre-computed in 01_generate_data.R
sim_params <- obj$sim_params
stem       <- tools::file_path_sans_ext(basename(rds_path))

# ---- fit ---------------------------------------------------------------------
message("Fitting fitQuasar2 ...")
t0          <- proc.time()
res_q2      <- fitQuasar2(dd, ~ Treatment)
elapsed_sec <- (proc.time() - t0)[["elapsed"]]
message(sprintf("fitQuasar2 finished in %.1f seconds", elapsed_sec))

# ---- tidy results for each coefficient ---------------------------------------
tidy_coef <- function(term_str, test_label, truth_col) {
  raw <- res_q2$results %>%
    filter(term == term_str) %>%
    mutate(padj = p.adjust(pval, "BH"))

  raw <- as.data.frame(raw)
  raw <- merge(
    raw[, c("identifier", "pval", "padj")],
    as.data.frame(truth)[, c("identifier", "is_ASE", "is_cASE")],
    by = "identifier", all.x = TRUE, sort = FALSE
  )

  tibble(
    method    = "QuASAR2",
    test      = test_label,
    truth_col = truth_col,
    pvalue    = as.numeric(raw$pval),
    padj      = as.numeric(raw$padj),
    truth_pos = as.logical(raw[[truth_col]])
  )
}

res <- bind_rows(
  tidy_coef("Treatmenttreatment", "cASE (trt coef)",  "is_cASE"),
  tidy_coef("(Intercept)",        "ASE (intercept)",  "is_ASE")
)

# ---- M estimates (bin-level, mapped back to per-SNP) -------------------------
# fitQuasar2 estimates one M per coverage bin via fit_q2$Mvec.
# We reconstruct the same binning used internally (nbreaks = 20 quantiles)
# and look up each SNP's bin M.
nbreaks_q2  <- 20
cov_per_snp <- dd %>%
  group_by(identifier) %>%
  summarise(mean_RA = mean(R + A), .groups = "drop")

cov_breaks <- unique(c(
  0,
  quantile(cov_per_snp$mean_RA,
           probs = seq_len(nbreaks_q2) / nbreaks_q2,
           na.rm = TRUE)
))

m_est <- cov_per_snp %>%
  mutate(
    bin   = cut(mean_RA, breaks = cov_breaks, include.lowest = TRUE),
    M_est = as.numeric(res_q2$Mvec[as.character(bin)])
  ) %>%
  transmute(
    method         = "QuASAR2",
    identifier,
    M_est,
    estimate_level = "coverage_bin",
    mean_RA
  ) %>%
  mutate(M_true = sim_params$M)

# ---- save as named list ------------------------------------------------------
out_file <- file.path(out_dir, paste0(stem, "_Q2.rds"))
saveRDS(
  list(
    results     = res,
    M_estimates = m_est,
    sim_params  = sim_params,
    elapsed_sec = elapsed_sec
  ),
  out_file
)
message("Saved: ", out_file)