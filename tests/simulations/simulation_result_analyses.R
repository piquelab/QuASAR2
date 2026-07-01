#!/usr/bin/env Rscript
# 06_aggregate_and_plot.R
# Reads all per-method RDS results, computes power/FPR/empirical-FDR averaged
# across seeds for each (N, M, delta) condition, and saves plots to PDF.
#
# Usage:
#   Rscript 06_aggregate_and_plot.R <results_root_dir> <out_pdf>
#
# <results_root_dir> should contain subdirectories CR/, GLM/, Q2/, LM/ (or
# you can just point it at a flat directory if you saved everything there).
# The script finds ALL .rds files recursively and binds them.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2)
  stop("Usage: Rscript 06_aggregate_and_plot.R <results_root_dir> <out_pdf>")

res_dir <- args[1]
out_pdf <- args[2]
fdr_alpha <- 0.1   # must match what was used at fit time

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(patchwork)
})

# ============================================================
# 1.  LOAD ALL RESULTS
# ============================================================
rds_files <- list.files(res_dir, pattern = "\\.rds$",
                         recursive = TRUE, full.names = TRUE)

# exclude the raw sim data files (they contain a data.frame without a 'method'
# column — skip gracefully)
message(sprintf("Found %d RDS files under %s", length(rds_files), res_dir))

all_res <- lapply(rds_files, function(f) {
  obj <- readRDS(f)

  # Raw sim data files are lists with dd/truth/sim_params — skip them
  if (!is.list(obj) || !"results" %in% names(obj)) return(NULL)

  res    <- obj$results
  params <- obj$sim_params

  if (!is.null(params)) {
    res$seed        <- params$seed
    res$N_lo        <- params$N_lo
    res$N_hi        <- params$N_hi
    res$M_val       <- params$M
    res$delta_val   <- params$delta
    res$elapsed_sec <- obj$elapsed_sec   # NA if missing (old files)
  }

  list(results = res, M_estimates = obj$M_estimates)
})

# Split into results and M_estimates streams
all_res        <- Filter(Negate(is.null), all_res)
res_long       <- bind_rows(lapply(all_res, `[[`, "results"))
M_estimates_all <- bind_rows(Filter(Negate(is.null),
                                    lapply(all_res, `[[`, "M_estimates")))

if (nrow(res_long) == 0)
  stop("No valid method result files found. Check your results directory.")

message(sprintf(
  "Loaded %d rows from %d files. Methods: %s",
  nrow(res_long),
  length(all_res),
  paste(unique(res_long$method), collapse = ", ")
))

# ============================================================
# 2.  COMPUTE POWER / FPR PER (method, test, seed, condition)
# ============================================================
compute_metrics <- function(padj, truth_pos, alpha = fdr_alpha) {
  called <- !is.na(padj) & padj < alpha
  tp     <- sum( called &  truth_pos, na.rm = TRUE)
  fp     <- sum( called & !truth_pos, na.rm = TRUE)
  fn     <- sum(!called &  truth_pos, na.rm = TRUE)
  tn     <- sum(!called & !truth_pos, na.rm = TRUE)
  data.frame(
    power         = if ((tp + fn) > 0) tp / (tp + fn) else NA_real_,
    fpr           = if ((fp + tn) > 0) fp / (fp + tn) else NA_real_,
    empirical_fdr = if ((tp + fp) > 0) fp / (tp + fp) else NA_real_,
    n_called      = sum(called),
    tp = tp, fp = fp
  )
}

per_seed <- res_long %>%
  group_by(method, test, seed, N_lo, N_hi, M_val, delta_val) %>%
  summarise(
    compute_metrics(padj, truth_pos),
    .groups = "drop"
  )

# ============================================================
# 3.  AVERAGE ACROSS SEEDS  (mean ± sd)
# ============================================================
summary_df <- per_seed %>%
  group_by(method, test, N_lo, N_hi, M_val, delta_val) %>%
  summarise(
    n_seeds       = n(),
    power_mean    = mean(power,         na.rm = TRUE),
    power_sd      = sd(power,           na.rm = TRUE),
    fpr_mean      = mean(fpr,           na.rm = TRUE),
    fpr_sd        = sd(fpr,             na.rm = TRUE),
    efdr_mean     = mean(empirical_fdr, na.rm = TRUE),
    efdr_sd       = sd(empirical_fdr,   na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    # human-readable coverage label
    N_label = paste0("N [", N_lo, ",", N_hi, "]")
  )

saveRDS(summary_df, sub("\\.pdf$", "_summary.rds", out_pdf))
message("Summary saved: ", sub("\\.pdf$", "_summary.rds", out_pdf))

# ============================================================
# 3b.  TIMING SUMMARY
# ============================================================
# elapsed_sec is per-seed, so average across seeds per (method, condition).
# Each result file has one elapsed_sec value covering all SNPs — one row per
# (seed x method), so deduplicate before summarising.

timing_df <- res_long %>%
  select(method, seed, N_lo, N_hi, M_val, delta_val, elapsed_sec) %>%
  distinct() %>%
  group_by(method, N_lo, N_hi, M_val, delta_val) %>%
  summarise(
    n_seeds          = n(),
    elapsed_mean_sec = mean(elapsed_sec, na.rm = TRUE),
    elapsed_sd_sec   = sd(elapsed_sec,   na.rm = TRUE),
    elapsed_min_sec  = min(elapsed_sec,  na.rm = TRUE),
    elapsed_max_sec  = max(elapsed_sec,  na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(N_label = paste0("N [", N_lo, ",", N_hi, "]"))

timing_csv <- sub("\\.pdf$", "_timing.csv", out_pdf)
write.csv(timing_df, timing_csv, row.names = FALSE)
message("Timing saved: ", timing_csv)

print(timing_df)

# ============================================================
# 3c.  M RECOVERY SUMMARY
# ============================================================
if (nrow(M_estimates_all) > 0) {

  M_summary <- M_estimates_all %>%
    filter(!is.na(M_est), is.finite(M_est), M_est > 0) %>%
    group_by(method, M_true, estimate_level) %>%
    summarise(
      n_estimates  = n(),
      median_M_est = median(M_est),
      mean_M_est   = mean(M_est),
      q25          = quantile(M_est, 0.25),
      q75          = quantile(M_est, 0.75),
      median_ratio = median(M_est) / first(M_true),
      .groups = "drop"
    )

  M_csv <- sub("\\.pdf$", "_M_recovery.csv", out_pdf)
  write.csv(M_summary, M_csv, row.names = FALSE)
  message("M recovery summary saved: ", M_csv)
  print(M_summary)

} else {
  message("No M estimates found — skipping M recovery summary.")
  M_summary <- NULL
}
method_colors <- c(
  "QuASAR2CR"  = "#E41A1C",
  "QuASAR_GLM" = "#377EB8",
  "QuASAR2"    = "#4DAF4A",
  "LM"         = "#984EA3"
)

base_theme <- theme_bw(base_size = 10) +
  theme(
    strip.background = element_rect(fill = "grey92"),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

# Generic panel: mean line + ±1 SD ribbon
make_panel <- function(df, x_var, x_lab,
                       y_mean, y_sd, y_lab,
                       facet_var, facet_lab,
                       h_line = NULL, reverse_x = FALSE,
                       pct_x = FALSE) {
  df$xv   <- df[[x_var]]
  df$ymn  <- df[[y_mean]]
  df$ysd  <- df[[y_sd]]
  df$fv   <- df[[facet_var]]

  p <- ggplot(df, aes(x = xv, y = ymn, color = method,
                      fill = method, group = method)) +
    geom_ribbon(aes(ymin = pmax(ymn - ysd, 0),
                    ymax = pmin(ymn + ysd, 1)),
                alpha = 0.12, color = NA) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    facet_wrap(~ fv, nrow = 1,
               labeller = labeller(fv = function(x) paste(facet_lab, x))) +
    { if (!is.null(h_line))
        geom_hline(yintercept = h_line, linetype = "dashed", color = "grey50") } +
    { if (reverse_x) scale_x_reverse() } +
    { if (pct_x)
        scale_x_continuous(labels = scales::percent_format(accuracy = 1)) } +
    scale_color_manual(values = method_colors) +
    scale_fill_manual(values  = method_colors) +
    scale_y_continuous(
      labels = scales::percent_format(accuracy = 1),
      limits = c(0, 1)
    ) +
    labs(x = x_lab, y = y_lab, color = "Method", fill = "Method",
         caption = paste0("Mean ± 1 SD across seeds  |  FDR threshold = ", fdr_alpha)) +
    base_theme
  p
}

pdf(out_pdf, width = 11, height = 13)

for (test_label in unique(summary_df$test)) {
  for (N_lab in unique(summary_df$N_label)) {

    sub <- summary_df %>%
      filter(test == test_label, N_label == N_lab)
    if (nrow(sub) == 0) next

    page_title <- sprintf("%s  |  Coverage: %s", test_label, N_lab)


    # ---- Plot B: x = M (overdispersion), facets = delta ----------------
    pB_pow <- make_panel(
      sub, "M_val", "M (overdispersion)  ←  noisier",
      "power_mean", "power_sd", "Power (TPR)",
      "delta_val", "delta =", reverse_x = TRUE
    ) + ggtitle(paste("B — Power by overdispersion |", page_title))


    pB_fdr <- make_panel(
      sub, "M_val", "M (overdispersion)  ←  noisier",
      "efdr_mean", "efdr_sd", "Empirical FDR",
      "delta_val", "delta =", h_line = fdr_alpha, reverse_x = TRUE
    )

    print(
      (pB_pow / pB_fdr) +
        plot_layout(guides = "collect") &
        theme(legend.position = "bottom", legend.title = element_blank())
    )
  }
}

# ---- M recovery boxplot (appended as final page) ----------------------------
if (nrow(M_estimates_all) > 0) {

  method_colors_M <- c(
    "QuASAR2CR"  = "#E41A1C",
    "QuASAR_GLM" = "#377EB8",
    "QuASAR2"    = "#4DAF4A"
  )

  M_plot_df <- M_estimates_all %>%
    filter(!is.na(M_est), is.finite(M_est), M_est > 0) %>%
    mutate(
      method   = factor(method, levels = names(method_colors_M)),
      M_label  = paste0("True M = ", M_true)
    )

  p_M <- ggplot(M_plot_df, aes(x = method, y = M_est, fill = method)) +
    geom_boxplot(outlier.alpha = 0.15, outlier.size = 0.8) +
    geom_hline(aes(yintercept = M_true),
               linetype = "dashed", color = "grey40") +
    facet_wrap(~ M_label, nrow = 1) +
    scale_y_log10() +
    scale_fill_manual(values = method_colors_M) +
    labs(
      x        = "Method",
      y        = "Estimated M  (log scale)",
      title    = "M recovery — estimated vs true overdispersion",
      subtitle = "Dashed line = true simulated M  |  LM does not estimate M"
    ) +
    theme_bw(base_size = 11) +
    theme(
      strip.background = element_rect(fill = "grey92"),
      panel.grid.minor = element_blank(),
      axis.text.x      = element_text(angle = 35, hjust = 1),
      legend.position  = "none"
    )

  print(p_M)
}

dev.off()
message("Plots saved: ", out_pdf)