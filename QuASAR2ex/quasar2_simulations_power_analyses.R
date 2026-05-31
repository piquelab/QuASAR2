# QuASAR2_simulations_v2.R
# Purpose: Simulation and power/FPR analysis comparing fitQuasar2CR, fitQuasar2,
#          fitQuasar_GLM, and a standard linear model (LM) for ASE/cASE detection.
# Author:  Shreya Nirmalan (updated)
# Date:    5/28/2026

# Load QuASAR2 locally
devtools::load_all("/rs/rs_grp_scaipgenetic/QuASAR2")

# fitQuasar_GLM is not yet exported from the QuASAR2 package, load it directly.
# TODO: remove this line once fitQuasar_GLM is integrated into QuASAR2.
source("/rs/rs_grp_scaipgenetic/QuASAR2/R/quasar_glm.R")   # adjust path as needed, e.g. source("R/fitQuasar_GLM.R")


suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
  library(readr)
  library(tibble)      
  library(qqman)
  #library(QuASAR2)    #uncomment when package is updated
  library(ggplot2)
  library(patchwork)   
})


# ============================================================
# 1.  SIMULATION DATA GENERATOR  (updated)
# ============================================================
# Key changes 
# -----------------
# * N_range        : per-replicate coverage drawn from Uniform(N_range[1], N_range[2])
#                    so each SNP × sample cell gets its own read depth.
#                    (Pass N_range = c(x, x) to keep coverage fixed, as before.)
# * delta_ASE_range / delta_cASE_range
#                  : effect sizes drawn per-SNP from Uniform(lo, hi).
#                    Pass c(x, x) for a fixed delta, as before.

sim_quasar2_df <- function(
  n_snps         = 20000,
  n_ctrl         = 5,
  n_trt          = 5,
  N_range        = c(60, 1000),   # per-replicate coverage range  [NEW]
  M              = 80,
  frac_ASE_only  = 0.0,
  frac_cASE_only = 0.0,
  frac_both      = 0.0,
  delta_ASE_range  = c(0.10, 0.10),  # [lo, hi] for per-SNP ASE effect  [NEW]
  delta_cASE_range = c(0.10, 0.10),  # [lo, hi] for per-SNP cASE effect [NEW]
  seed           = 1
) {
  set.seed(seed)

  n_samp <- n_ctrl + n_trt

  ## ---- sample metadata ----
  samples <- data.frame(
    SampleID  = paste0("S", seq_len(n_samp)),
    Treatment = factor(
      c(rep("control", n_ctrl), rep("treatment", n_trt)),
      levels = c("control", "treatment")
    )
  )

  ## ---- SNP truth assignment ----
  snps  <- paste0("rsSIM", seq_len(n_snps))
  class <- rep("null", n_snps)

  n_both <- round(frac_both * n_snps)
  if (n_both > 0) class[sample.int(n_snps, n_both)] <- "both"

  remaining <- which(class == "null")
  n_ase     <- round(frac_ASE_only * n_snps)
  if (n_ase  > 0) class[sample(remaining, min(n_ase,  length(remaining)))] <- "ASE_only"

  remaining <- which(class == "null")
  n_case    <- round(frac_cASE_only * n_snps)
  if (n_case > 0) class[sample(remaining, min(n_case, length(remaining)))] <- "cASE_only"

  ## ---- per-SNP direction and effect sizes ----
  dir <- sample(c(-1L, 1L), n_snps, replace = TRUE)

  delta_ASE_snp  <- runif(n_snps,
                          min = delta_ASE_range[1],
                          max = delta_ASE_range[2])
  delta_cASE_snp <- runif(n_snps,
                          min = delta_cASE_range[1],
                          max = delta_cASE_range[2])

  truth <- data.frame(
    identifier    = snps,
    class         = class,
    is_ASE        = class %in% c("ASE_only", "both"),
    is_cASE       = class %in% c("cASE_only", "both"),
    direction     = dir,
    delta_ASE     = delta_ASE_snp,
    delta_cASE    = delta_cASE_snp
  )

  ## ---- cross SNP × sample ----
  dd <- merge(truth, samples, all = TRUE)
  # merge() sorts by first key; re-sort for clarity
  dd <- dd[order(dd$identifier, dd$SampleID), ]

  snp_idx    <- match(dd$identifier, snps)
  is_treated <- dd$Treatment == "treatment"

  ## ---- mean ref proportion p ----
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

  ## ---- per-replicate coverage (NEW) ----
  dd$N <- round(runif(nrow(dd), min = N_range[1], max = N_range[2]))

  ## ---- draw allele counts ----
  theta <- rbeta(nrow(dd), p * M, (1 - p) * M)
  dd$R  <- rbinom(nrow(dd), size = dd$N, prob = theta)
  dd$A  <- dd$N - dd$R

  ## ---- mean coverage per SNP (for binning / diagnostics) ----
  dd$mean_RA <- ave(dd$R + dd$A, dd$identifier, FUN = mean)

  dd
}


# ============================================================
# 2.  STANDARD LINEAR MODEL 
# ============================================================
# Per-SNP log(ref+1 / alt+1) ~ Treatment.
# Returns tidy results with p.value, padj for intercept and treatment term.

mylm <- function(data) {
  # beta1 = log-ratio of ref to alt counts (pseudocount +1 on each)
  obj    <- lm(beta1 ~ Treatment, data = data)
  mysum  <- summary(obj)
  mycoeff <- as.data.frame(mysum$coefficients)
  colnames(mycoeff) <- c("estimate", "std.error", "statistic", "p.value")
  mycoeff %>%
    rownames_to_column("term") %>%
    mutate(
      deg.f = obj$df.residual,
      rank  = obj$qr$rank
    )
}

fit_lm_ase <- function(dd, fdr_method = "BH") {
  # pre-compute beta1 = log((R+1)/(A+1)) before nesting
  dd_lm <- dd %>%
    mutate(
      ref.reads1 = R + 1,
      alt.reads1 = A + 1,
      beta1      = log(ref.reads1 / alt.reads1)
    )

  aux   <- dd_lm %>% group_by(identifier) %>% nest()
  resc  <- aux %>% mutate(res = map(data, mylm))
  resc2 <- resc %>%
    select(identifier, data = res) %>%
    unnest(cols = c(data)) %>%
    as.data.frame()

  # split into intercept and treatment results, add FDR adjustment
  lm_trt <- resc2 %>%
    filter(term == "Treatmenttreatment") %>%
    mutate(padj = p.adjust(p.value, method = fdr_method))

  lm_int <- resc2 %>%
    filter(term == "(Intercept)") %>%
    mutate(padj = p.adjust(p.value, method = fdr_method))

  list(trt = lm_trt, int = lm_int)
}


# ============================================================
# 3.  POWER / FPR EXTRACTION HELPER
# ============================================================
# Given a results table joined with truth columns is_cASE / is_ASE,
# returns power (TPR) and FPR at a given FDR threshold.

compute_power_fpr <- function(padj, truth_positive, alpha = 0.1) {
  called     <- !is.na(padj) & padj < alpha
  tp         <- sum( called &  truth_positive, na.rm = TRUE)
  fp         <- sum( called & !truth_positive, na.rm = TRUE)
  fn         <- sum(!called &  truth_positive, na.rm = TRUE)
  tn         <- sum(!called & !truth_positive, na.rm = TRUE)

  power <- if ((tp + fn) > 0) tp / (tp + fn) else NA_real_
  fpr   <- if ((fp + tn) > 0) fp / (fp + tn) else NA_real_

  list(power = power, fpr = fpr, n_called = sum(called))
}


# ============================================================
# 4.  SINGLE SIMULATION RUN → METRICS  (one row per method)
# ============================================================

run_one_sim <- function(sim_params, fdr_alpha = 0.1, verbose = FALSE) {

  dd <- do.call(sim_quasar2_df, sim_params)

  ## truth table (one row per SNP)
  truth <- dd %>%
    group_by(identifier) %>%
    summarise(
      is_ASE  = first(is_ASE),
      is_cASE = first(is_cASE),
      mean_N  = mean(N),
      mean_RA = mean(R + A),
      .groups = "drop"
    )

  ## coverage bin labels for the truth table (used in stratified power curves)
  truth$cov_bin <- cut(
    truth$mean_RA,
    breaks = quantile(truth$mean_RA, probs = c(0, 0.33, 0.67, 1)),
    labels = c("low", "mid", "high"),
    include.lowest = TRUE
  )

  ## ---- fit models ----
  results <- list()

  # 4a. fitQuasar2CR
  tryCatch({
    fit_cr <- fitQuasar2CR(dd, design = ~ Treatment, max_iter = 6, verbose = verbose)
    res_cr_trt <- testCoef(fit_cr, coef = "Treatmenttreatment", df_method = "moderated") %>%
      left_join(truth, by = "identifier")
    res_cr_int <- testCoef(fit_cr, coef = "(Intercept)",        df_method = "moderated") %>%
      left_join(truth, by = "identifier")

    results$CR_trt <- tibble(
      method    = "QuASAR2CR",
      test      = "cASE (trt coef)",
      truth_col = "is_cASE",
      pvalue    = res_cr_trt$p.value,
      padj      = res_cr_trt$padj,
      truth_pos = res_cr_trt$is_cASE,
      cov_bin   = res_cr_trt$cov_bin
    )
    results$CR_int <- tibble(
      method    = "QuASAR2CR",
      test      = "ASE (intercept)",
      truth_col = "is_ASE",
      pvalue    = res_cr_int$p.value,
      padj      = res_cr_int$padj,
      truth_pos = res_cr_int$is_ASE,
      cov_bin   = res_cr_int$cov_bin
    )
  }, error = function(e) message("fitQuasar2CR error: ", conditionMessage(e)))

  # 4b. fitQuasar_GLM
  tryCatch({
    res_glm  <- fitQuasar_GLM(dd, ~ Treatment)
    trt_glm  <- res_glm$results %>%
      filter(term == "Treatmenttreatment") %>%
      mutate(padj = p.adjust(p.value, "BH")) %>%
      left_join(truth, by = "identifier")
    int_glm  <- res_glm$results %>%
      filter(term == "(Intercept)") %>%
      mutate(padj = p.adjust(p.value, "BH")) %>%
      left_join(truth, by = "identifier")

    results$GLM_trt <- tibble(
      method    = "QuASAR_GLM",
      test      = "cASE (trt coef)",
      truth_col = "is_cASE",
      pvalue    = trt_glm$p.value,
      padj      = trt_glm$padj,
      truth_pos = trt_glm$is_cASE,
      cov_bin   = trt_glm$cov_bin
    )
    results$GLM_int <- tibble(
      method    = "QuASAR_GLM",
      test      = "ASE (intercept)",
      truth_col = "is_ASE",
      pvalue    = int_glm$p.value,
      padj      = int_glm$padj,
      truth_pos = int_glm$is_ASE,
      cov_bin   = int_glm$cov_bin
    )
  }, error = function(e) message("fitQuasar_GLM error: ", conditionMessage(e)))

  # 4c. fitQuasar2 (original, no Cox-Reid)
  tryCatch({
    res_q2  <- fitQuasar2(dd, ~ Treatment)
    trt_q2  <- res_q2$results %>%
      filter(term == "Treatmenttreatment") %>%
      mutate(padj = p.adjust(pval, "BH")) %>%
      left_join(truth, by = "identifier")
    int_q2  <- res_q2$results %>%
      filter(term == "(Intercept)") %>%
      mutate(padj = p.adjust(pval, "BH")) %>%
      left_join(truth, by = "identifier")

    results$Q2_trt <- tibble(
      method    = "QuASAR2",
      test      = "cASE (trt coef)",
      truth_col = "is_cASE",
      pvalue    = trt_q2$p.value,
      padj      = trt_q2$padj,
      truth_pos = trt_q2$is_cASE,
      cov_bin   = trt_q2$cov_bin
    )
    results$Q2_int <- tibble(
      method    = "QuASAR2",
      test      = "ASE (intercept)",
      truth_col = "is_ASE",
      pvalue    = int_q2$p.value,
      padj      = int_q2$padj,
      truth_pos = int_q2$is_ASE,
      cov_bin   = int_q2$cov_bin
    )
  }, error = function(e) message("fitQuasar2 error: ", conditionMessage(e)))

  # 4d. Standard linear model  (log((R+1)/(A+1)) ~ Treatment, per-SNP)
  tryCatch({
    res_lm <- fit_lm_ase(dd)
    lm_trt <- res_lm$trt %>% left_join(truth, by = "identifier")
    lm_int <- res_lm$int %>% left_join(truth, by = "identifier")

    results$LM_trt <- tibble(
      method    = "LM",
      test      = "cASE (trt coef)",
      truth_col = "is_cASE",
      pvalue    = lm_trt$p.value,
      padj      = lm_trt$padj,
      truth_pos = lm_trt$is_cASE,
      cov_bin   = lm_trt$cov_bin
    )
    results$LM_int <- tibble(
      method    = "LM",
      test      = "ASE (intercept)",
      truth_col = "is_ASE",
      pvalue    = lm_int$p.value,
      padj      = lm_int$padj,
      truth_pos = lm_int$is_ASE,
      cov_bin   = lm_int$cov_bin
    )
  }, error = function(e) message("LM error: ", conditionMessage(e)))

  if (length(results) == 0) return(NULL)

  bind_rows(results)
}


# ============================================================
# 5.  POWER ANALYSIS GRID
# ============================================================
# Sweeps over:
#   - N_range  : coverage scenarios (low / medium / high)
#   - M        : overdispersion scenarios (tight / moderate / overdispersed)
#   - delta    : effect sizes (subtle / moderate / large)
# For each cell: simulate once, compute power & FPR per method × test × coverage bin.

build_power_grid <- function(
  N_scenarios    = list(low  = c(20,  80),
                        mid  = c(60,  300),
                        high = c(200, 1000)),
  M_scenarios    = c(tight = 500, moderate = 100, overdispersed = 20),
  delta_scenarios= c(subtle = 0.05, moderate = 0.10, large = 0.20),
  n_snps         = 20000,   # fewer SNPs to keep run time manageable
  n_ctrl         = 5,
  n_trt          = 5,
  frac_ASE_only  = 0.05,
  frac_cASE_only = 0.05,
  fdr_alpha      = 0.1,
  seed_base      = 42,
  verbose        = FALSE
) {
  grid <- expand.grid(
    N_name     = names(N_scenarios),
    M_name     = names(M_scenarios),
    delta_name = names(delta_scenarios),
    stringsAsFactors = FALSE
  )

  pb <- txtProgressBar(min = 0, max = nrow(grid), style = 3)

  all_rows <- vector("list", nrow(grid))

  for (i in seq_len(nrow(grid))) {
    setTxtProgressBar(pb, i)
    row  <- grid[i, ]
    Nrng <- N_scenarios[[row$N_name]]
    Mval <- M_scenarios[row$M_name]
    dval <- delta_scenarios[row$delta_name]

    sim_p <- list(
      n_snps         = n_snps,
      n_ctrl         = n_ctrl,
      n_trt          = n_trt,
      N_range        = Nrng,
      M              = Mval,
      frac_ASE_only  = frac_ASE_only,
      frac_cASE_only = frac_cASE_only,
      delta_ASE_range  = c(dval, dval),
      delta_cASE_range = c(dval, dval),
      seed           = seed_base + i
    )

    res_long <- run_one_sim(sim_p, fdr_alpha = fdr_alpha, verbose = verbose)
    if (is.null(res_long)) next

    # Compute power & FPR per method × test × coverage bin
    summary_rows <- res_long %>%
      group_by(method, test, cov_bin) %>%
      summarise(
        power    = compute_power_fpr(padj, truth_pos, alpha = fdr_alpha)$power,
        fpr      = compute_power_fpr(padj, truth_pos, alpha = fdr_alpha)$fpr,
        n_called = compute_power_fpr(padj, truth_pos, alpha = fdr_alpha)$n_called,
        .groups  = "drop"
      ) %>%
      mutate(
        N_scenario    = row$N_name,
        M_scenario    = row$M_name,
        delta_scenario= row$delta_name,
        N_lo          = Nrng[1],
        N_hi          = Nrng[2],
        M_val         = Mval,
        delta_val     = dval
      )

    all_rows[[i]] <- summary_rows
  }
  close(pb)

  bind_rows(all_rows)
}


# ============================================================
# 6.  PLOTTING: POWER & FPR CURVES  →  PDF
# ============================================================

plot_power_fpr <- function(
  power_df,
  out_pdf  = "QuASAR2_power_analysis.pdf",
  fdr_alpha = 0.1
) {
  method_colors <- c(
    "QuASAR2CR"  = "#E41A1C",
    "QuASAR_GLM" = "#377EB8",
    "QuASAR2"    = "#4DAF4A",
    "LM"         = "#984EA3"
  )

  # helper: facet over delta (x-axis) × M (facet col) per N scenario
  make_curve_plot <- function(df, y_var, y_lab, title_suffix, h_line = NULL) {
    df$y_val <- df[[y_var]]
    ggplot(df, aes(x = delta_val, y = y_val,
                   color = method, group = method)) +
      geom_line(linewidth = 0.9) +
      geom_point(size = 2) +
      facet_grid(cov_bin ~ M_scenario,
                 labeller = labeller(
                   cov_bin    = function(x) paste("coverage:", x),
                   M_scenario = function(x) paste("M:", x)
                 )) +
      { if (!is.null(h_line))
          geom_hline(yintercept = h_line, linetype = "dashed", color = "grey50") } +
      scale_color_manual(values = method_colors) +
      scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                         limits = c(0, 1)) +
      labs(
        title  = paste(y_lab, "—", title_suffix),
        x      = "Effect size (delta)",
        y      = y_lab,
        color  = "Method",
        caption= paste0("FDR threshold = ", fdr_alpha)
      ) +
      theme_bw(base_size = 10) +
      theme(
        strip.background = element_rect(fill = "grey92"),
        legend.position  = "bottom"
      )
  }

  pdf(out_pdf, width = 11, height = 8.5)

  for (test_label in unique(power_df$test)) {
    for (N_lab in unique(power_df$N_scenario)) {
      sub <- power_df %>%
        filter(test == test_label, N_scenario == N_lab)

      if (nrow(sub) == 0) next

      title_str <- sprintf("%s | Coverage: %s", test_label, N_lab)

      p_power <- make_curve_plot(sub, "power", "Power (TPR)", title_str)
      p_fpr   <- make_curve_plot(sub, "fpr",   "False Positive Rate", title_str,
                                  h_line = fdr_alpha)

      print(p_power / p_fpr)   # patchwork: stack vertically
    }
  }

  dev.off()
  message("Saved: ", out_pdf)
}


# ============================================================
# 7.  QQ-PLOT COMPARISON  →  PDF
# ============================================================
# Overlay QQ plots for all four methods on the same axes (null SNPs only).

plot_qq_comparison <- function(
  sim_params,
  out_pdf = "QuASAR2_qq_comparison.pdf",
  verbose = FALSE
) {
  res_long <- run_one_sim(sim_params, verbose = verbose)
    qq_data <- res_long %>%
        group_by(method, test) %>%
        arrange(pvalue) %>%
        mutate(
          expected = -log10(ppoints(n())),
          observed = -log10(pmax(pvalue, 1e-300))
        ) %>%
        ungroup()
    
      p_qq <- ggplot(
        qq_data,
        aes(x = expected, y = observed, color = method)
      ) +
        geom_abline(
          intercept = 0,
          slope = 1,
          linetype = "dashed",
          color = "grey50"
        ) +
        geom_line(alpha = 0.8) +
        facet_wrap(~ test, ncol = 2) +
        scale_color_manual(values = c(
          "QuASAR2CR"  = "#E41A1C",
          "QuASAR_GLM" = "#377EB8",
          "QuASAR2"    = "#4DAF4A",
          "LM"         = "#984EA3"
        )) +
        labs(
          title = "QQ plot",
          x = "Expected -log10(p)",
          y = "Observed -log10(p)",
          color = "Method"
        ) +
        theme_bw(base_size = 11)
    
      pdf(out_pdf, width = 10, height = 6)
      print(p_qq)
      dev.off()
  message("Saved: ", out_pdf)
}


# ============================================================
# 8.  EXAMPLE / QUICK-RUN  (edit as needed)
# ============================================================

if (TRUE) {   # set to TRUE to execute

  ## --- 8a. Single simulation (sanity check) ---
  dd <- sim_quasar2_df(
    n_snps          = 10000,
    n_ctrl          = 5,
    n_trt           = 5,
    N_range         = c(60, 1000),        # per-replicate, varies by sample
    M               = 200,
    frac_ASE_only   = 0.05,
    frac_cASE_only  = 0.05,
    delta_ASE_range  = c(0.05, 0.20),     # per-SNP effect varies between 5–20%
    delta_cASE_range = c(0.05, 0.20),
    seed            = 1
  )

  # coverage distribution check
  hist(dd$N, breaks = 60, main = "Per-replicate coverage", xlab = "N")

  # fit QuASAR2CR
  fit <- fitQuasar2CR(dd, design = ~ Treatment, max_iter = 6, verbose = TRUE)
  print(fit)

  resCR.Tr  <- testCoef(fit, coef = "Treatmenttreatment", df_method = "moderated")
  resCR.Int <- testCoef(fit, coef = "(Intercept)",        df_method = "moderated")

  qq(resCR.Tr$p.value[!dd %>% group_by(identifier) %>%
                         summarise(is_cASE = first(is_cASE)) %>% pull(is_cASE)],
     main = "QQ - cASE (null SNPs only)")

  ## --- 8b. QQ-plot comparison across methods ---
  plot_qq_comparison(
    sim_params = list(
      n_snps          = 10000,
      n_ctrl          = 5,
      n_trt           = 5,
      N_range         = c(60, 300),
      M               = 100,
      frac_ASE_only   = 0.05,
      frac_cASE_only  = 0.05,
      delta_ASE_range  = c(0.10, 0.10),
      delta_cASE_range = c(0.10, 0.10),
      seed            = 99
    ),
    out_pdf = "QuASAR2_qq_comparison.pdf"
  )

  ## --- 8c. Full power analysis grid ---
  power_df <- build_power_grid(
    N_scenarios     = list(low  = c(20,  80),
                           mid  = c(60,  300),
                           high = c(200, 1000)),
    M_scenarios     = c(tight = 500, moderate = 100, overdispersed = 20),
    delta_scenarios = c(subtle = 0.05, moderate = 0.10, large = 0.20),
    n_snps          = 10000,
    n_ctrl          = 5,
    n_trt           = 5,
    frac_ASE_only   = 0.05,
    frac_cASE_only  = 0.05,
    fdr_alpha       = 0.1,
    seed_base       = 42
  )

  saveRDS(power_df, "QuASAR2_power_grid.rds")  # cache results

  plot_power_fpr(power_df, out_pdf = "QuASAR2_power_analysis.pdf")
}