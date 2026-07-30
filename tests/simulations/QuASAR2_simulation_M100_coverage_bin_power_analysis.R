
# Purpose:
#   For simulations with true M = 100 and coverage range 20–1000:
#     1. Load results from QuASAR2CR, QuASAR_GLM, QuASAR2, and LM.
#     2. Match each result to its simulated SNP-level mean coverage.
#     3. Bin SNPs into approximately log2-spaced coverage bins.
#     4. Compute power and empirical FDR within each bin.
#     5. Average performance across simulation seeds.
#     6. Save summary tables and plots.
#
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

res_dir <- "/rs/rs_grp_scaipgenetic/QuASAR2/tests/simulations/results"
out_pdf     <- "/rs/rs_grp_scaipgenetic/QuASAR2/tests/simulations/results/M_100_coverage_bin_power_FDR_plots.pdf"
alpha <- 0.10

# ============================================================
# 1. LOAD TRUTH TABLES FOR M = 100, N = 20–1000
# ============================================================

truth_all <- bind_rows(lapply(
  list.files(file.path(res_dir, "data"), "\\.rds$", full.names = TRUE),
  function(f) {
    x <- readRDS(f)
    p <- x$sim_params

    if (is.null(p) || p$M != 100 || p$N_lo != 20 || p$N_hi != 1000)
      return(NULL)

    x$truth %>%
      select(identifier, mean_RA) %>%
      mutate(
        seed = p$seed,
        delta_val = p$delta
      )
  }
))

# log2-style coverage bins
breaks <- c(20, 40, 80, 160, 320, 640, 1001)
labels <- c("20–40", "40–80", "80–160",
            "160–320", "320–640", "640–1000")
midpoints <- c(30, 60, 120, 240, 480, 800)

truth_all <- truth_all %>%
  mutate(
    coverage_bin = cut(
      mean_RA,
      breaks = breaks,
      labels = labels,
      include.lowest = TRUE,
      right = FALSE
    ),
    coverage_midpoint = midpoints[as.integer(coverage_bin)]
  )


# ============================================================
# 2. LOAD ALL FOUR METHOD RESULTS
# ============================================================

result_files <- unlist(lapply(
  file.path(res_dir, c("CR", "GLM", "Q2", "LM")),
  list.files,
  pattern = "\\.rds$",
  full.names = TRUE
))

res_long <- bind_rows(lapply(result_files, function(f) {
  x <- readRDS(f)
  p <- x$sim_params

  if (is.null(p) || p$M != 100 || p$N_lo != 20 || p$N_hi != 1000)
    return(NULL)

  x$results %>%
    select(identifier, method, test, padj, truth_pos) %>%
    mutate(
      seed = p$seed,
      delta_val = p$delta
    )
}))

# Join SNP coverage onto method results
dat <- res_long %>%
  left_join(
    truth_all,
    by = c("identifier", "seed", "delta_val")
  ) %>%
  filter(!is.na(coverage_bin))


# ============================================================
# 3. COMPUTE POWER, FPR, EMPIRICAL FDR PER SEED/BIN
# ============================================================

metrics <- dat %>%
  mutate(called = !is.na(padj) & padj < alpha) %>%
  group_by(method, test, seed, delta_val,
           coverage_bin, coverage_midpoint) %>%
  summarise(
    tp = sum(called & truth_pos, na.rm = TRUE),
    fp = sum(called & !truth_pos, na.rm = TRUE),
    fn = sum(!called & truth_pos, na.rm = TRUE),
    tn = sum(!called & !truth_pos, na.rm = TRUE),

    power = ifelse(tp + fn > 0, tp / (tp + fn), NA_real_),
    efdr  = ifelse(tp + fp > 0, fp / (tp + fp), NA_real_),

    .groups = "drop"
  )

summary_df <- metrics %>%
  group_by(method, test, delta_val,
           coverage_bin, coverage_midpoint) %>%
  summarise(
    power_mean = mean(power, na.rm = TRUE),
    power_sd   = sd(power, na.rm = TRUE),
    efdr_mean  = mean(efdr, na.rm = TRUE),
    efdr_sd    = sd(efdr, na.rm = TRUE),
    n_seeds    = n(),
    .groups = "drop"
  )

write.csv(
  summary_df,
  sub("\\.pdf$", "_summary.csv", out_pdf),
  row.names = FALSE
)


# ============================================================
# 4. PLOT
# ============================================================

method_colors <- c(
  "QuASAR2CR"  = "#E41A1C",
  "QuASAR_GLM" = "#377EB8",
  "QuASAR2"    = "#4DAF4A",
  "LM"         = "#984EA3"
)

make_plot <- function(df, mean_col, sd_col, ylab, hline = NULL) {
  df$y <- df[[mean_col]]
  df$s <- df[[sd_col]]

  p <- ggplot(
    df,
    aes(
      x = coverage_midpoint,
      y = y,
      color = method,
      fill = method,
      group = method
    )
  ) +
    geom_ribbon(
      aes(
        ymin = pmax(y - s, 0),
        ymax = pmin(y + s, 1)
      ),
      alpha = 0.12,
      color = NA
    ) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2) +
    facet_wrap(~ delta_val, nrow = 1,
               labeller = label_both) +
    scale_x_log10(
      breaks = midpoints,
      labels = labels
    ) +
    scale_y_continuous(
      limits = c(0, 1),
      labels = scales::percent_format()
    ) +
    scale_color_manual(values = method_colors) +
    scale_fill_manual(values = method_colors) +
    labs(
      x = "Mean coverage per SNP",
      y = ylab,
      color = "Method",
      fill = "Method"
    ) +
    theme_bw() +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      axis.text.x = element_text(angle = 35, hjust = 1)
    )

  if (!is.null(hline))
    p <- p + geom_hline(
      yintercept = hline,
      linetype = "dashed",
      color = "grey50"
    )

  p
}

pdf(out_pdf, width = 12, height = 9)

for (tt in unique(summary_df$test)) {

  sub <- filter(summary_df, test == tt)

  p_power <- make_plot(
    sub,
    "power_mean",
    "power_sd",
    "Power"
  ) + ggtitle(paste(tt, "— power by coverage"))

  p_fdr <- make_plot(
    sub,
    "efdr_mean",
    "efdr_sd",
    "Empirical FDR",
    hline = alpha
  )

  print(
    (p_power / p_fdr) +
      plot_layout(guides = "collect") &
      theme(
        legend.position = "bottom",
        legend.title = element_blank()
      )
  )
}

dev.off()

message("Saved: ", out_pdf)