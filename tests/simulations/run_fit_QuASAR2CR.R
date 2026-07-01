#!/usr/bin/env Rscript
# /run_fit_QuASAR2CR.R
# Fits fitQuasar2CR on a saved simulated dataset and saves raw results.
#
# Usage:
#   Rscript /run_fit_QuASAR2CR.R <rds_path> <out_dir>
#
# Output:
#   <out_dir>/<stem>_CR.rds   — tibble with columns:
#     identifier, term, estimate, std.error, statistic, p.value, padj,
#     is_ASE, is_cASE, cov_bin   + sim_params as an attribute


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("Usage: Rscript run_fit_QuASAR2CR.R <rds_path> <out_dir>")

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

# ---- load data ----------------------------------------------------------------
obj        <- readRDS(rds_path)
dd         <- obj$dd
truth      <- obj$truth        # pre-computed in 01_generate_data.R
sim_params <- obj$sim_params
stem       <- tools::file_path_sans_ext(basename(rds_path))

# ---- fit ----------------------------------------------------------------------
message("Fitting fitQuasar2CR ...")
t0          <- proc.time()
fit         <- fitQuasar2CR(dd, design = ~ Treatment, max_iter = 6, verbose = FALSE)
elapsed_sec <- (proc.time() - t0)[["elapsed"]]
message(sprintf("fitQuasar2CR finished in %.1f seconds", elapsed_sec))

# ---- extract results for both coefficients ------------------------------------
extract_coef <- function(coef_name, test_label, truth_col) {
  raw <- testCoef(fit, coef = coef_name, df_method = "moderated")

  # Defensive: testCoef() may return matrix/list-columns (e.g. covariance terms).
  # Join only the identifier + flat truth columns we actually need, then pull
  # everything out as plain atomic vectors before building the tibble.
  raw <- as.data.frame(raw)  # strip any exotic tibble subclassing
  raw <- merge(
    raw[, c("identifier", "p.value", "padj")],
    as.data.frame(truth)[, c("identifier", "is_ASE", "is_cASE")],
    by = "identifier", all.x = TRUE, sort = FALSE
  )

  tibble(
    method    = "QuASAR2CR",
    test      = test_label,
    truth_col = truth_col,
    pvalue    = as.numeric(raw$p.value),
    padj      = as.numeric(raw$padj),
    truth_pos = as.logical(raw[[truth_col]])
  )
}

res <- bind_rows(
  extract_coef("Treatmenttreatment", "cASE (trt coef)", "is_cASE"),
  extract_coef("(Intercept)",        "ASE (intercept)", "is_ASE")
)

# ---- M estimates (per-SNP from dispersion slot) ------------------------------
m_est <- fit$dispersion %>%
  transmute(
    method         = "QuASAR2CR",
    identifier,
    M_est          = M,
    estimate_level = "SNP"
  ) %>%
  left_join(truth %>% select(identifier, mean_RA), by = "identifier") %>%
  mutate(M_true = sim_params$M)

# ---- save as named list ------------------------------------------------------
out_file <- file.path(out_dir, paste0(stem, "_CR.rds"))
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