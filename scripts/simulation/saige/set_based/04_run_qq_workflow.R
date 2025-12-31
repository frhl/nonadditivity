#!/usr/bin/env Rscript
# author: frederik lassen
# note: this script creates QQ plots comparing genetic models at variant and gene levels
#       faceted by heritability and replicate

library(bravastring)
library(data.table)
library(latex2exp)
library(ggplot2)
library(ggrastr)

# Define paths
data_dir <- "/Users/flassen/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/derived/saige_simulation_results"
output_dir <- "/Users/flassen/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/derived/saige_simulation_results/plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Function to create confidence band ribbon
create_confidence_ribbon <- function(n_tests, ribbon_p = 0.95) {
  n <- round(n_tests * 2)
  seq.p.value <- (1:n) / n

  data.table(
    p.value.expt = -log10(seq.p.value),
    clower = -log10(qbeta(p = (1 + ribbon_p) / 2, shape1 = 1:n, shape2 = n:1)),
    cupper = -log10(qbeta(p = (1 - ribbon_p) / 2, shape1 = 1:n, shape2 = n:1))
  )
}

# Extract mode from filename
extract_mode <- function(filename) {
  if (grepl("\\.additive\\.", filename)) {
    return("Additive")
  } else if (grepl("\\.recessive\\.", filename)) {
    return("Recessive")
  } else if (grepl("\\.dominance\\.", filename)) {
    return("Dominance")
  } else {
    return(NA_character_)
  }
}

# Extract heritability from filename
extract_h2 <- function(filename) {
  match <- regmatches(filename, regexpr("p_[0-9.]+_continuous", filename))
  if (length(match) > 0) {
    h2 <- gsub("p_", "", gsub("_continuous.*", "", match))
    return(paste0("h²=", h2))
  }
  return(NA_character_)
}

# Extract replicate from filename
extract_replicate <- function(filename) {
  match <- regmatches(filename, regexpr("continuous_[0-9]+", filename))
  if (length(match) > 0) {
    rep <- gsub("continuous_", "", match)
    return(paste0("Rep ", rep))
  }
  return(NA_character_)
}

# ============================================
# Process variant-level results by mode
# ============================================
cat("============================================\n")
cat("Processing variant-level SAIGE results\n")
cat("============================================\n\n")

variant_files <- list.files(data_dir, pattern = "\\.(additive|recessive|dominance)\\..*\\.txt\\.gz$", full.names = TRUE)
variant_files <- variant_files[!grepl("chr", variant_files)]  # Exclude gene-based files

if (length(variant_files) == 0) {
  cat("WARNING: No variant-level result files found!\n\n")
  dt_variant_all <- NULL
} else {
  cat("Found", length(variant_files), "variant-level result files\n\n")

  # Read all files and tag with mode, heritability, and replicate
  dt_variant_list <- lapply(variant_files, function(vf) {
    cat("Reading:", basename(vf), "\n")
    dt <- fread(vf)
    dt[, mode := extract_mode(vf)]
    dt[, heritability := extract_h2(vf)]
    dt[, replicate := extract_replicate(vf)]
    dt[, file := basename(vf)]
    dt[, test_level := "Variant"]
    return(dt)
  })

  dt_variant_all <- rbindlist(dt_variant_list, fill = TRUE)

  # Remove any rows with missing metadata
  dt_variant_all <- dt_variant_all[!is.na(mode) & !is.na(heritability) & !is.na(replicate)]

  # Calculate expected p-values
  dt_variant_all[, p.value.expt := bravastring::get_expected_p(p.value), by = .(mode, heritability, replicate)]

  cat("\nSummary by mode (variant-level):\n")
  variant_summary <- dt_variant_all[, .(
    n_tests = .N,
    n_files = length(unique(file)),
    lambda_gc = bravastring::calc_inflation(p.value)
  ), by = mode]
  print(variant_summary)

  cat("\nSummary by heritability:\n")
  h2_summary <- dt_variant_all[, .(
    n_tests = .N,
    lambda_gc = bravastring::calc_inflation(p.value)
  ), by = .(mode, heritability)]
  print(h2_summary)

  # ============================================
  # Create variant-level QQ plots by mode
  # ============================================
  cat("\n============================================\n")
  cat("Creating variant-level QQ plots\n")
  cat("============================================\n\n")

  # 1. Overall comparison by mode
  p_variant_modes <- ggplot(dt_variant_all, aes(x = -log10(p.value.expt), y = -log10(p.value), color = mode)) +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
    geom_point_rast(size = 0.5, alpha = 0.4) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "bottom",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10)
    ) +
    scale_color_manual(
      name = "Genetic Model",
      values = c("Additive" = "#E41A1C", "Recessive" = "#377EB8", "Dominance" = "#4DAF4A")
    ) +
    ggtitle(
      "Variant-Level Tests: Genetic Model Comparison",
      TeX(sprintf(
        "$\\lambda_{GC}$: Add=%.3f, Rec=%.3f, Dom=%.3f",
        variant_summary[mode == "Additive", lambda_gc],
        variant_summary[mode == "Recessive", lambda_gc],
        variant_summary[mode == "Dominance", lambda_gc]
      ))
    ) +
    ylab(TeX("$-\\log_{10}(P_{observed})$")) +
    xlab(TeX("$-\\log_{10}(P_{expected})$"))

  ggsave(file.path(output_dir, "qq_plot_variant_modes_overall.png"), p_variant_modes, width = 10, height = 8, dpi = 300)
  cat("Saved: qq_plot_variant_modes_overall.png\n")

  # 2. Faceted grid: heritability (rows) x replicate (columns)
  cat("Creating plots faceted by heritability (rows) and replicate (columns)...\n")

  p_h2_rep_grid <- ggplot(dt_variant_all, aes(x = -log10(p.value.expt), y = -log10(p.value), color = mode)) +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
    geom_point_rast(size = 0.3, alpha = 0.5) +
    facet_grid(heritability ~ replicate) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 8),
      legend.position = "bottom",
      legend.title = element_text(size = 11, face = "bold"),
      legend.text = element_text(size = 10),
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(size = 9, face = "bold")
    ) +
    scale_color_manual(
      name = "Genetic Model",
      values = c("Additive" = "#E41A1C", "Recessive" = "#377EB8", "Dominance" = "#4DAF4A")
    ) +
    ggtitle("Variant-Level Tests: Heritability (rows) × Replicate (columns)") +
    ylab(TeX("$-\\log_{10}(P_{observed})$")) +
    xlab(TeX("$-\\log_{10}(P_{expected})$"))

  ggsave(file.path(output_dir, "qq_plot_variant_h2_by_replicate.png"), p_h2_rep_grid, width = 16, height = 18, dpi = 300)
  cat("Saved: qq_plot_variant_h2_by_replicate.png\n")

  # 4. Faceted by mode (separate panels for each genetic model)
  cat("Creating plots faceted by mode...\n")

  max_n <- max(dt_variant_all[, .N, by = mode]$N)
  ribbon_dt <- create_confidence_ribbon(max_n, ribbon_p = 0.95)

  p_mode_facet <- ggplot(dt_variant_all, aes(x = -log10(p.value.expt), y = -log10(p.value))) +
    geom_ribbon(
      data = ribbon_dt,
      aes(x = p.value.expt, ymin = clower, ymax = cupper, y = NULL),
      fill = "grey80", alpha = 0.5, inherit.aes = FALSE
    ) +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
    geom_point_rast(aes(color = mode), size = 0.5, alpha = 0.5) +
    facet_wrap(~ mode, nrow = 1) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "none",
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(size = 11, face = "bold")
    ) +
    scale_color_manual(
      values = c("Additive" = "#E41A1C", "Recessive" = "#377EB8", "Dominance" = "#4DAF4A")
    ) +
    ggtitle("Variant-Level Tests by Genetic Model") +
    ylab(TeX("$-\\log_{10}(P_{observed})$")) +
    xlab(TeX("$-\\log_{10}(P_{expected})$"))

  ggsave(file.path(output_dir, "qq_plot_variant_by_mode.png"), p_mode_facet, width = 15, height = 5, dpi = 300)
  cat("Saved: qq_plot_variant_by_mode.png\n")

  # 5. Individual plots for each mode with confidence bands
  cat("\nCreating individual plots for each mode...\n")
  for (gm in c("Additive", "Recessive", "Dominance")) {
    dt_mode <- dt_variant_all[mode == gm]
    if (nrow(dt_mode) == 0) next

    ribbon_mode <- create_confidence_ribbon(nrow(dt_mode), ribbon_p = 0.95)
    lambda_mode <- bravastring::calc_inflation(dt_mode$p.value)
    n_tests_mode <- nrow(dt_mode)

    p_mode <- ggplot(dt_mode, aes(x = -log10(p.value.expt), y = -log10(p.value))) +
      geom_ribbon(
        data = ribbon_mode,
        aes(x = p.value.expt, ymin = clower, ymax = cupper, y = NULL),
        fill = "grey80", alpha = 0.5, inherit.aes = FALSE
      ) +
      geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
      geom_point_rast(size = 0.8, alpha = 0.6) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 14, face = "bold"),
        plot.subtitle = element_text(size = 11),
        axis.title = element_text(size = 12),
        axis.text = element_text(size = 10)
      ) +
      ggtitle(
        sprintf("Variant-Level: %s Model", gm),
        TeX(sprintf("$\\lambda_{GC} = %.4f$, $n_{tests} = %s$", lambda_mode, format(n_tests_mode, big.mark = ",")))
      ) +
      ylab(TeX("$-\\log_{10}(P_{observed})$")) +
      xlab(TeX("$-\\log_{10}(P_{expected})$"))

    ggsave(file.path(output_dir, sprintf("qq_plot_variant_%s.png", tolower(gm))), p_mode, width = 8, height = 8, dpi = 300)
    cat("  Saved: qq_plot_variant_", tolower(gm), ".png\n", sep = "")
  }
}

# ============================================
# Process gene-level results
# ============================================
cat("\n============================================\n")
cat("Processing gene-level SAIGE results\n")
cat("============================================\n\n")

gene_files <- list.files(data_dir, pattern = "^UKB\\.auto\\.chr.*\\.txt\\.gz$", full.names = TRUE)
gene_files <- gene_files[!grepl("singleAssoc", gene_files)]  # Exclude single-variant within gene results

if (length(gene_files) == 0) {
  cat("WARNING: No gene-level result files found!\n\n")
  dt_gene_all <- NULL
} else {
  cat("Found", length(gene_files), "gene-level result files\n\n")

  # Read all files and tag with heritability and replicate
  dt_gene_list <- lapply(gene_files, function(gf) {
    cat("Reading:", basename(gf), "\n")
    dt <- fread(gf)
    dt[, heritability := extract_h2(gf)]
    dt[, replicate := extract_replicate(gf)]
    dt[, file := basename(gf)]
    dt[, test_level := "Gene"]
    return(dt)
  })

  dt_gene_all <- rbindlist(dt_gene_list, fill = TRUE)

  # Use burden test p-value
  if ("Pvalue_Burden" %in% names(dt_gene_all)) {
    dt_gene_all[, p.value := Pvalue_Burden]
  } else if ("Pvalue" %in% names(dt_gene_all)) {
    dt_gene_all[, p.value := Pvalue]
  } else {
    stop("Cannot find p-value column in gene-based results!")
  }

  # Remove NA p-values
  dt_gene_all <- dt_gene_all[!is.na(p.value)]
  dt_gene_all <- dt_gene_all[!is.na(heritability) & !is.na(replicate)]

  # Calculate expected p-values
  dt_gene_all[, p.value.expt := bravastring::get_expected_p(p.value), by = .(heritability, replicate)]

  cat("\nSummary (gene-level):\n")
  gene_summary <- dt_gene_all[, .(
    n_tests = .N,
    n_files = length(unique(file)),
    n_genes = length(unique(Region)),
    lambda_gc = bravastring::calc_inflation(p.value)
  )]
  print(gene_summary)

  cat("\nSummary by heritability (gene-level):\n")
  gene_h2_summary <- dt_gene_all[, .(
    n_tests = .N,
    lambda_gc = bravastring::calc_inflation(p.value)
  ), by = heritability]
  print(gene_h2_summary)

  # ============================================
  # Create gene-level QQ plots
  # ============================================
  cat("\n============================================\n")
  cat("Creating gene-level QQ plots\n")
  cat("============================================\n\n")

  # 1. Overall gene-level plot
  ribbon_gene <- create_confidence_ribbon(nrow(dt_gene_all), ribbon_p = 0.95)

  p_gene <- ggplot(dt_gene_all, aes(x = -log10(p.value.expt), y = -log10(p.value))) +
    geom_ribbon(
      data = ribbon_gene,
      aes(x = p.value.expt, ymin = clower, ymax = cupper, y = NULL),
      fill = "grey80", alpha = 0.5, inherit.aes = FALSE
    ) +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
    geom_point_rast(size = 1.2, alpha = 0.6, color = "#984EA3") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      plot.subtitle = element_text(size = 11),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10)
    ) +
    ggtitle(
      "Gene-Level Set-Based Tests",
      TeX(sprintf(
        "$\\lambda_{GC} = %.4f$, $n_{tests} = %s$, $n_{genes} = %s$",
        gene_summary$lambda_gc,
        format(gene_summary$n_tests, big.mark = ","),
        format(gene_summary$n_genes, big.mark = ",")
      ))
    ) +
    ylab(TeX("$-\\log_{10}(P_{observed})$")) +
    xlab(TeX("$-\\log_{10}(P_{expected})$"))

  ggsave(file.path(output_dir, "qq_plot_gene_overall.png"), p_gene, width = 8, height = 8, dpi = 300)
  cat("Saved: qq_plot_gene_overall.png\n")

  # 2. Faceted grid: heritability (rows) x replicate (columns)
  cat("Creating gene-level plots faceted by heritability (rows) and replicate (columns)...\n")

  p_gene_h2_rep <- ggplot(dt_gene_all, aes(x = -log10(p.value.expt), y = -log10(p.value))) +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
    geom_point_rast(size = 0.8, alpha = 0.6, color = "#984EA3") +
    facet_grid(heritability ~ replicate) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 8),
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(size = 9, face = "bold")
    ) +
    ggtitle("Gene-Level Tests: Heritability (rows) × Replicate (columns)") +
    ylab(TeX("$-\\log_{10}(P_{observed})$")) +
    xlab(TeX("$-\\log_{10}(P_{expected})$"))

  ggsave(file.path(output_dir, "qq_plot_gene_h2_by_replicate.png"), p_gene_h2_rep, width = 16, height = 18, dpi = 300)
  cat("Saved: qq_plot_gene_h2_by_replicate.png\n")
}

# ============================================
# Create combined variant vs gene comparison
# ============================================
if (!is.null(dt_variant_all) && !is.null(dt_gene_all)) {
  cat("\n============================================\n")
  cat("Creating variant vs gene level comparison\n")
  cat("============================================\n\n")

  # Prepare data - take additive mode for variant-level comparison
  dt_variant_add <- dt_variant_all[mode == "Additive", .(p.value, p.value.expt)]
  dt_variant_add[, comparison := "Variant-Level (Additive)"]

  dt_gene_comp <- dt_gene_all[, .(p.value, p.value.expt)]
  dt_gene_comp[, comparison := "Gene-Level (Set-Based)"]

  dt_comparison <- rbind(dt_variant_add, dt_gene_comp)

  p_comparison <- ggplot(dt_comparison, aes(x = -log10(p.value.expt), y = -log10(p.value), color = comparison)) +
    geom_abline(intercept = 0, slope = 1, color = "black", linetype = "dashed") +
    geom_point_rast(size = 0.5, alpha = 0.4) +
    facet_wrap(~ comparison, scales = "free") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, face = "bold"),
      axis.title = element_text(size = 12),
      axis.text = element_text(size = 10),
      legend.position = "none",
      strip.background = element_rect(fill = "grey90"),
      strip.text = element_text(size = 11, face = "bold")
    ) +
    scale_color_manual(values = c("Variant-Level (Additive)" = "#E41A1C", "Gene-Level (Set-Based)" = "#984EA3")) +
    ggtitle("SAIGE Association Tests: Variant-Level vs Gene-Level") +
    ylab(TeX("$-\\log_{10}(P_{observed})$")) +
    xlab(TeX("$-\\log_{10}(P_{expected})$"))

  ggsave(file.path(output_dir, "qq_plot_variant_vs_gene.png"), p_comparison, width = 12, height = 6, dpi = 300)
  cat("Saved: qq_plot_variant_vs_gene.png\n")
}

cat("\n============================================\n")
cat("All QQ plots generated successfully!\n")
cat("============================================\n")
cat("Output directory:", output_dir, "\n")
cat("\nGenerated plots:\n")
cat("  VARIANT-LEVEL:\n")
cat("    - qq_plot_variant_modes_overall.png (all modes overlay)\n")
cat("    - qq_plot_variant_by_mode.png (faceted by mode)\n")
cat("    - qq_plot_variant_h2_by_replicate.png (grid: h² rows × rep columns)\n")
cat("    - qq_plot_variant_additive.png (additive only)\n")
cat("    - qq_plot_variant_recessive.png (recessive only)\n")
cat("    - qq_plot_variant_dominance.png (dominance only)\n")
cat("  GENE-LEVEL:\n")
cat("    - qq_plot_gene_overall.png (all combined)\n")
cat("    - qq_plot_gene_h2_by_replicate.png (grid: h² rows × rep columns)\n")
cat("  COMPARISON:\n")
cat("    - qq_plot_variant_vs_gene.png (variant vs gene)\n")
