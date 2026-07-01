#!/usr/bin/env Rscript
# run_fit_LM.R
# Fits a standard per-SNP linear model log((R+1)/(A+1)) ~ Treatment.
#
# Usage:
#   Rscript run_fit_LM.R <rds_path> <out_dir>


args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("Usage: Rscript run_fit_LM.R <rds_path> <out_dir>")

rds_path <- args[1]
out_dir  <- args[2]
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

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

# ---- per-SNP LM helper -------------------------------------------------------
mylm <- function(data) {
  obj     <- lm(beta1 ~ Treatment, data = data)
  mysum   <- summary(obj)
  mycoeff <- as.data.frame(mysum$coefficients)
  colnames(mycoeff) <- c("estimate", "std.error", "statistic", "p.value")
  mycoeff %>%
    rownames_to_column("term") %>%
    mutate(deg.f = obj$df.residual, rank = obj$qr$rank)
}

# ---- fit ---------------------------------------------------------------------
message("Fitting LM ...")
t0    <- proc.time()
dd_lm <- dd %>%
  mutate(
    beta1 = log((R + 1) / (A + 1))
  )

all_coefs <- dd_lm %>%
  group_by(identifier) %>%
  nest() %>%
  mutate(res = map(data, mylm)) %>%
  select(identifier, data = res) %>%
  unnest(cols = c(data)) %>%
  as.data.frame()
elapsed_sec <- (proc.time() - t0)[["elapsed"]]
message(sprintf("LM finished in %.1f seconds", elapsed_sec))

# ---- tidy results per coefficient --------------------------------------------
tidy_coef <- function(term_str, test_label, truth_col) {
  raw <- all_coefs %>%
    filter(term == term_str) %>%
    mutate(padj = p.adjust(p.value, "BH"))

  raw <- as.data.frame(raw)
  raw <- merge(
    raw[, c("identifier", "p.value", "padj")],
    as.data.frame(truth)[, c("identifier", "is_ASE", "is_cASE")],
    by = "identifier", all.x = TRUE, sort = FALSE
  )

  tibble(
    method    = "LM",
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

# ---- M estimates: LM does not estimate overdispersion -----------------------
m_est <- NULL

# ---- save as named list ------------------------------------------------------
out_file <- file.path(out_dir, paste0(stem, "_LM.rds"))
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