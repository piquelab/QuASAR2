#!/usr/bin/env Rscript
# run_fit_QuASAR2_GLM.R
# Fits fitQuasar_GLM on a saved simulated dataset and saves raw results.
#
# Usage:
#   Rscript run_fit_QuASAR2_GLM.R <rds_path> <out_dir>


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("Usage: Rscript run_fit_QuASAR2_GLM.R <rds_path> <out_dir>")

rds_path <- args[1]
out_dir  <- args[2]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

devtools::load_all("/rs/rs_grp_scaipgenetic/QuASAR2")
source("/rs/rs_grp_scaipgenetic/QuASAR2/R/quasar_glm.R")

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
message("Fitting fitQuasar_GLM ...")
t0          <- proc.time()
res_glm     <- fitQuasar_GLM(dd, ~ Treatment)
elapsed_sec <- (proc.time() - t0)[["elapsed"]]
message(sprintf("fitQuasar_GLM finished in %.1f seconds", elapsed_sec))

# ---- tidy results for each coefficient ---------------------------------------
tidy_coef <- function(term_str, test_label, truth_col) {
  raw <- res_glm$results %>%
    filter(term == term_str) %>%
    mutate(padj = p.adjust(p.value, "BH"))

  raw <- as.data.frame(raw)
  raw <- merge(
    raw[, c("identifier", "p.value", "padj")],
    as.data.frame(truth)[, c("identifier", "is_ASE", "is_cASE")],
    by = "identifier", all.x = TRUE, sort = FALSE
  )

  tibble(
    method    = "QuASAR_GLM",
    test      = test_label,
    truth_col = truth_col,
    pvalue    = as.numeric(raw$p.value),
    padj      = as.numeric(raw$padj),
    truth_pos = as.logical(raw[[truth_col]])
  )
}

res <- bind_rows(
  tidy_coef("Treatmenttreatment", "cASE (trt coef)",  "is_cASE"),
  tidy_coef("(Intercept)",        "ASE (intercept)",  "is_ASE")
)

# ---- M estimates (per-SNP from results slot) ---------------------------------
m_est <- res_glm$results %>%
  distinct(identifier, M) %>%
  transmute(
    method         = "QuASAR_GLM",
    identifier,
    M_est          = M,
    estimate_level = "SNP"
  ) %>%
  left_join(truth %>% select(identifier, mean_RA), by = "identifier") %>%
  mutate(M_true = sim_params$M)

# ---- save as named list ------------------------------------------------------
out_file <- file.path(out_dir, paste0(stem, "_GLM.rds"))
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