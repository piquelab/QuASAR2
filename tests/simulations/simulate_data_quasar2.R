#!/usr/bin/env Rscript
# 01_generate_data.R
# Generates simulated ASE/cASE datasets and saves each as an RDS file.
#
# Usage (single condition, one seed at a time — called from bash):
#   Rscript 01_generate_data.R <seed> <out_dir> [N_lo] [N_hi] [M] [delta]
#
# Default condition (if optional args omitted):
#   N_lo=60, N_hi=300, M=100, delta=0.10
#
# Output:
#   <out_dir>/sim_data_seed<seed>_N<N_lo>-<N_hi>_M<M>_d<delta>.rds

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("Usage: Rscript 01_generate_data.R <seed> <out_dir> [N_lo N_hi M delta]")

seed    <- as.integer(args[1])
out_dir <- args[2]
N_lo    <- if (length(args) >= 4) as.numeric(args[3]) else 60
N_hi    <- if (length(args) >= 4) as.numeric(args[4]) else 300
M       <- if (length(args) >= 5) as.numeric(args[5]) else 100
delta   <- if (length(args) >= 6) as.numeric(args[6]) else 0.10

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# ============================================================
# Simulation function (self-contained copy so this script has
# no dependency on the QuASAR2 package)
# ============================================================
sim_quasar2_df <- function(
  n_snps         = 5000,
  n_ctrl         = 5,
  n_trt          = 5,
  N_range        = c(60, 300),
  M              = 100,
  frac_ASE_only  = 0.05,
  frac_cASE_only = 0.05,
  frac_both      = 0.0,
  delta_ASE_range  = c(0.10, 0.10),
  delta_cASE_range = c(0.10, 0.10),
  seed           = 1
) {
  set.seed(seed)

  n_samp <- n_ctrl + n_trt

  samples <- data.frame(
    SampleID  = paste0("S", seq_len(n_samp)),
    Treatment = factor(
      c(rep("control", n_ctrl), rep("treatment", n_trt)),
      levels = c("control", "treatment")
    )
  )

  snps  <- paste0("rsSIM", seq_len(n_snps))
  class <- rep("null", n_snps)

  n_both <- round(frac_both * n_snps)
  if (n_both > 0) class[sample.int(n_snps, n_both)] <- "both"

  remaining <- which(class == "null")
  n_ase     <- round(frac_ASE_only * n_snps)
  if (n_ase  > 0) class[sample(remaining, min(n_ase, length(remaining)))] <- "ASE_only"

  remaining <- which(class == "null")
  n_case    <- round(frac_cASE_only * n_snps)
  if (n_case > 0) class[sample(remaining, min(n_case, length(remaining)))] <- "cASE_only"

  dir_snp        <- sample(c(-1L, 1L), n_snps, replace = TRUE)
  delta_ASE_snp  <- runif(n_snps, min = delta_ASE_range[1],  max = delta_ASE_range[2])
  delta_cASE_snp <- runif(n_snps, min = delta_cASE_range[1], max = delta_cASE_range[2])

  truth <- data.frame(
    identifier  = snps,
    class       = class,
    is_ASE      = class %in% c("ASE_only", "both"),
    is_cASE     = class %in% c("cASE_only", "both"),
    direction   = dir_snp,
    delta_ASE   = delta_ASE_snp,
    delta_cASE  = delta_cASE_snp
  )

  dd <- merge(truth, samples, all = TRUE)
  dd <- dd[order(dd$identifier, dd$SampleID), ]

  is_treated <- dd$Treatment == "treatment"

  p <- rep(0.5, nrow(dd))

  ase_rows  <- dd$is_ASE
  p[ase_rows] <- 0.5 + dd$direction[ase_rows] * dd$delta_ASE[ase_rows]

  case_rows <- dd$is_cASE & is_treated
  p[case_rows] <- 0.5 + dd$direction[case_rows] * dd$delta_cASE[case_rows]

  both_all <- dd$class == "both"
  p[both_all] <- 0.5 + dd$direction[both_all] * dd$delta_ASE[both_all]

  both_trt <- dd$class == "both" & is_treated
  p[both_trt] <- p[both_trt] + dd$direction[both_trt] * dd$delta_cASE[both_trt]

  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  dd$prop_true <- p

  dd$N     <- round(runif(nrow(dd), min = N_range[1], max = N_range[2]))
  theta    <- rbeta(nrow(dd), p * M, (1 - p) * M)
  dd$R     <- rbinom(nrow(dd), size = dd$N, prob = theta)
  dd$A     <- dd$N - dd$R
  dd$mean_RA <- ave(dd$R + dd$A, dd$identifier, FUN = mean)

  dd
}

# ============================================================
# Generate and save
# ============================================================
message(sprintf(
  "Generating: seed=%d  N=[%g,%g]  M=%g  delta=%g",
  seed, N_lo, N_hi, M, delta
))

dd <- sim_quasar2_df(
  n_snps           = 10000,
  n_ctrl           = 5,
  n_trt            = 5,
  N_range          = c(N_lo, N_hi),
  M                = M,
  frac_ASE_only    = 0.05,
  frac_cASE_only   = 0.05,
  delta_ASE_range  = c(delta, delta),
  delta_cASE_range = c(delta, delta),
  seed             = seed
)

suppressPackageStartupMessages(library(dplyr))

# ---- truth table (one row per SNP, computed once here for all methods) -------
# cov_bin boundaries are fixed here so they are identical across all method scripts.
truth <- dd %>%
  group_by(identifier) %>%
  summarise(
    is_ASE  = first(is_ASE),
    is_cASE = first(is_cASE),
    mean_RA = mean(R + A),
    .groups = "drop"
  )

# ---- condition metadata ------------------------------------------------------
sim_params <- list(
  seed  = seed,
  N_lo  = N_lo,
  N_hi  = N_hi,
  M     = M,
  delta = delta
)

# ---- save as named list so fit scripts load dd, truth, and params together ---
out_file <- file.path(
  out_dir,
  sprintf("sim_data_seed%02d_N%g-%g_M%g_d%g.rds", seed, N_lo, N_hi, M, delta)
)

saveRDS(
  list(dd = dd, truth = truth, sim_params = sim_params),
  out_file
)
message("Saved: ", out_file)