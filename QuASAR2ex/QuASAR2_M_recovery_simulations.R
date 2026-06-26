# QuASAR2_M_recovery_simulations.R
# Purpose: Simulation analysis comparing estimated overdispersion M across
#          fitQuasar2CR, fitQuasar2, and fitQuasar_GLM.
# Author:  Shreya Nirmalan
# Date:    6/24/2026

# Load QuASAR2 locally
devtools::load_all("/rs/rs_grp_scaipgenetic/QuASAR2")

# fitQuasar_GLM is not yet exported from the QuASAR2 package, load it directly.
source("/rs/rs_grp_scaipgenetic/QuASAR2/R/quasar_glm.R")

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(tibble)
  library(ggplot2)
})


# ============================================================
# 1. SIMULATION DATA GENERATOR
# ============================================================

sim_quasar2_df <- function(
  n_snps         = 2000,
  n_ctrl         = 5,
  n_trt          = 5,
  N_range        = c(60, 300),
  M              = 100,
  frac_ASE_only  = 0.05,
  frac_cASE_only = 0.05,
  frac_both      = 0.00,
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
  n_ase <- round(frac_ASE_only * n_snps)
  if (n_ase > 0) class[sample(remaining, min(n_ase, length(remaining)))] <- "ASE_only"

  remaining <- which(class == "null")
  n_case <- round(frac_cASE_only * n_snps)
  if (n_case > 0) class[sample(remaining, min(n_case, length(remaining)))] <- "cASE_only"

  truth <- data.frame(
    identifier = snps,
    class = class,
    is_ASE = class %in% c("ASE_only", "both"),
    is_cASE = class %in% c("cASE_only", "both"),
    direction = sample(c(-1L, 1L), n_snps, replace = TRUE),
    delta_ASE = runif(n_snps, delta_ASE_range[1], delta_ASE_range[2]),
    delta_cASE = runif(n_snps, delta_cASE_range[1], delta_cASE_range[2])
  )

  dd <- merge(truth, samples, all = TRUE)
  dd <- dd[order(dd$identifier, dd$SampleID), ]

  is_treated <- dd$Treatment == "treatment"

  p <- rep(0.5, nrow(dd))

  ase_rows <- dd$is_ASE
  p[ase_rows] <- 0.5 + dd$direction[ase_rows] * dd$delta_ASE[ase_rows]

  case_rows <- dd$is_cASE & is_treated
  p[case_rows] <- 0.5 + dd$direction[case_rows] * dd$delta_cASE[case_rows]

  both_all <- dd$class == "both"
  p[both_all] <- 0.5 + dd$direction[both_all] * dd$delta_ASE[both_all]

  both_trt <- dd$class == "both" & is_treated
  p[both_trt] <- p[both_trt] + dd$direction[both_trt] * dd$delta_cASE[both_trt]

  p <- pmin(pmax(p, 1e-6), 1 - 1e-6)
  dd$prop_true <- p

  dd$N <- round(runif(nrow(dd), min = N_range[1], max = N_range[2]))

  theta <- rbeta(nrow(dd), p * M, (1 - p) * M)
  dd$R <- rbinom(nrow(dd), size = dd$N, prob = theta)
  dd$A <- dd$N - dd$R

  dd$mean_RA <- ave(dd$R + dd$A, dd$identifier, FUN = mean)

  dd
}


# ============================================================
# 2. SINGLE SIMULATION RUN → M ESTIMATES
# ============================================================

run_one_M_sim <- function(sim_params, verbose = FALSE) {

  dd <- do.call(sim_quasar2_df, sim_params)
  truth_cov <- dd %>%
  group_by(identifier) %>%
  summarise(
    mean_RA = mean(R + A),
    .groups = "drop"
  )

  M_true <- sim_params$M

  results <- list()

  # ---- fitQuasar2CR ----
  tryCatch({
    fit_cr <- fitQuasar2CR(
      dd,
      design = ~ Treatment,
      max_iter = 6,
      verbose = verbose
    )

    results$CR <- fit_cr$dispersion %>%
      transmute(
        method = "QuASAR2CR",
        identifier,
        M_est = M,
        estimate_level = "SNP"
      )
  }, error = function(e) message("fitQuasar2CR error: ", conditionMessage(e)))


  # ---- fitQuasar_GLM ----
  tryCatch({
    fit_glm <- fitQuasar_GLM(dd, ~ Treatment)

    results$GLM <- fit_glm$results %>%
      distinct(identifier, M) %>%
      transmute(
        method = "QuASAR_GLM",
        identifier,
        M_est = M,
        estimate_level = "SNP"
      )
  }, error = function(e) message("fitQuasar_GLM error: ", conditionMessage(e)))


  # ---- fitQuasar2 original ----
tryCatch({
  nbreaks_q2 <- 20

  fit_q2 <- fitQuasar2(
    dd,
    ~ Treatment,
    nbreaks = nbreaks_q2
  )

  q2_bins <- dd %>%
    group_by(identifier) %>%
    summarise(
      mean_RA = mean(R + A),
      .groups = "drop"
    )

  cov_breaks <- unique(c(
    0,
    quantile(
      q2_bins$mean_RA,
      probs = (1:nbreaks_q2) / nbreaks_q2,
      na.rm = TRUE
    )
  ))

  q2_bins <- q2_bins %>%
    mutate(
      bin = cut(mean_RA, breaks = cov_breaks),
      M_est = as.numeric(fit_q2$Mvec[as.character(bin)])
    ) %>%
    transmute(
      method = "QuASAR2",
      identifier,
      M_est,
      estimate_level = "coverage_bin"
    )

  results$Q2 <- q2_bins

}, error = function(e) message("fitQuasar2 error: ", conditionMessage(e)))

  if (length(results) == 0) return(NULL)

  bind_rows(results) %>%
  left_join(truth_cov, by = "identifier") %>%
  mutate(M_true = M_true)
}


# ============================================================
# 3. M RECOVERY GRID
# ============================================================

build_M_recovery_grid <- function(
  N_scenarios = list(
    mid = c(60, 300)
  ),
  M_scenarios = c(
    overdispersed = 20,
    moderate = 100,
    tight = 500
  ),
  n_snps = 2000,
  n_ctrl = 5,
  n_trt = 5,
  frac_ASE_only = 0.05,
  frac_cASE_only = 0.05,
  delta = 0.10,
  seed_base = 42,
  verbose = FALSE
) {

  grid <- expand.grid(
    N_name = names(N_scenarios),
    M_name = names(M_scenarios),
    stringsAsFactors = FALSE
  )

  message(sprintf("Grid: %d cells × %d SNPs each.", nrow(grid), n_snps))

  all_rows <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {

    row <- grid[i, ]

    Nrng <- N_scenarios[[row$N_name]]
    Mval <- M_scenarios[row$M_name]

    message("\nCell ", i, "/", nrow(grid),
            ": N = ", row$N_name,
            ", M = ", row$M_name,
            " (", Mval, ")")

    sim_p <- list(
      n_snps = n_snps,
      n_ctrl = n_ctrl,
      n_trt = n_trt,
      N_range = Nrng,
      M = Mval,
      frac_ASE_only = frac_ASE_only,
      frac_cASE_only = frac_cASE_only,
      delta_ASE_range = c(delta, delta),
      delta_cASE_range = c(delta, delta),
      seed = seed_base + i
    )

    res <- run_one_M_sim(sim_p, verbose = verbose)

    if (!is.null(res)) {
      all_rows[[i]] <- res %>%
        mutate(
          N_scenario = row$N_name,
          M_scenario = row$M_name,
          N_lo = Nrng[1],
          N_hi = Nrng[2],
          cell_id = i
        )
    }
  }

  bind_rows(all_rows)
}


# ============================================================
# 4. PLOTTING: M RECOVERY BOXPLOTS
# ============================================================

plot_M_recovery <- function(
  M_df,
  out_pdf = "QuASAR2_M_recovery_boxplots.pdf"
) {

  M_df <- M_df %>%
    filter(!is.na(M_est), is.finite(M_est), M_est > 0) %>%
    mutate(
      M_scenario = factor(
        M_scenario,
        levels = c("overdispersed", "moderate", "tight"),
        labels = c("M overdispersed", "M moderate", "M tight")
      ),
      method = factor(
        method,
        levels = c("QuASAR2", "QuASAR2CR", "QuASAR_GLM")
      )
    )

  p <- ggplot(M_df, aes(x = method, y = M_est)) +
    geom_boxplot(outlier.alpha = 0.15) +
    geom_hline(
      aes(yintercept = M_true),
      linetype = "dashed",
      color = "grey40"
    ) +
    facet_wrap(~ M_scenario, nrow = 1) +
    scale_y_log10() +
    labs(
      x = "Method",
      y = "Estimated M",
      title = "Estimated overdispersion by true M scenario",
      subtitle = "Dashed line = true simulated M"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "grey92"),
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1)
    )

  pdf(out_pdf, width = 11, height = 4.5)
  print(p)
  dev.off()

  message("Saved: ", out_pdf)

  p
}

plot_M_vs_coverage <- function(
  M_df,
  out_pdf = "QuASAR2_M_vs_coverage.pdf"
) {

  M_df <- M_df %>%
    filter(
      !is.na(M_est),
      is.finite(M_est),
      M_est > 0,
      !is.na(mean_RA),
      mean_RA > 0
    ) %>%
    mutate(
      log10_mean_RA = log10(mean_RA),
      method = factor(
        method,
        levels = c("QuASAR2", "QuASAR2CR", "QuASAR_GLM")
      )
    )

  p <- ggplot(M_df, aes(x = log10_mean_RA, y = M_est)) +
    geom_point(alpha = 0.25, size = 0.7) +
    geom_hline(
      aes(yintercept = M_true),
      linetype = "dashed",
      color = "grey40"
    ) +
    facet_wrap(~ method, nrow = 1) +
    scale_y_log10() +
    labs(
      x = "log10(mean coverage per SNP)",
      y = "Estimated M",
      title = "Estimated M as a function of coverage",
      subtitle = "Dashed line = true simulated M"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "grey92"),
      panel.grid.minor = element_blank()
    )

  pdf(out_pdf, width = 11, height = 4)
  print(p)
  dev.off()

  message("Saved: ", out_pdf)

  p
}
# ============================================================
# 5. SUMMARY TABLE
# ============================================================

summarize_M_recovery <- function(M_df) {

  M_df %>%
    filter(!is.na(M_est), is.finite(M_est), M_est > 0) %>%
    group_by(M_scenario, M_true, method, estimate_level) %>%
    summarise(
      n_estimates = n(),
      median_M_est = median(M_est),
      mean_M_est = mean(M_est),
      q25 = quantile(M_est, 0.25),
      q75 = quantile(M_est, 0.75),
      median_ratio = median_M_est / first(M_true),
      .groups = "drop"
    )
}


# ============================================================
# 6. RUN
# ============================================================

if (TRUE) {

  M_df <- build_M_recovery_grid(
    N_scenarios = list(
      mid = c(60, 300)
    ),
    M_scenarios = c(
      overdispersed = 20,
      moderate = 100,
      tight = 500
    ),
    n_snps = 10000,
    n_ctrl = 5,
    n_trt = 5,
    frac_ASE_only = 0.05,
    frac_cASE_only = 0.05,
    delta = 0.10,
    seed_base = 42,
    verbose = FALSE
  )

  write_csv(M_df, paste0("QuASAR2_M_recovery_estimates_", Sys.Date(), ".csv"))

  M_summary <- summarize_M_recovery(M_df)
  write_csv(M_summary, paste0("QuASAR2_M_recovery_summary_", Sys.Date(), ".csv"))

  plot_M_recovery(
    M_df,
    out_pdf = "QuASAR2_M_recovery_boxplots.pdf"
  )

  print(M_summary)
}


# M versus coverage plot
M_df <- build_M_recovery_grid(
  N_scenarios = list(
    wide = c(20, 1000)
  ),
  M_scenarios = c(
    moderate = 100
  ),
  n_snps = 10000,
  n_ctrl = 5,
  n_trt = 5,
  frac_ASE_only = 0.05,
  frac_cASE_only = 0.05,
  delta = 0.10,
  seed_base = 42,
  verbose = FALSE
)
plot_M_vs_coverage(
  M_df,
  out_pdf = paste0("QuASAR2_M_vs_coverage", Sys.Date(), ".pdf")
)