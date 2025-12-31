#!/usr/bin/env Rscript

options(scipen = 999)

# Required packages
packages <- c('data.table', 'ggplot2', 'ggrastr', 'latex2exp', 'optparse', 'stringr')

for (p in packages) {
  if (!require(p, character.only = TRUE, quietly = TRUE)) {
    cat(paste0("Installing package: ", p, "\n"))
    if (p == "bravastring") {
      if (!require("devtools", quietly = TRUE)) {
        install.packages("devtools", repos = "http://cran.us.r-project.org")
      }
      devtools::install_github("frhl/bravastring")
    } else {
      install.packages(p, repos = "http://cran.us.r-project.org")
    }
    library(p, character.only = TRUE)
  }
}

# Load bravastring separately
if (!require("bravastring", quietly = TRUE)) {
  if (!require("devtools", quietly = TRUE)) {
    install.packages("devtools", repos = "http://cran.us.r-project.org")
    library(devtools)
  }
  devtools::install_github("frhl/bravastring")
  library(bravastring)
}

cat("=============================================================\n")
cat("Comparing SAIGE and REGENIE results\n")
cat("=============================================================\n\n")

opt <- list(
  ancestry = "eur",
  annotation = "pLoF_damaging_missense",
  saige_dir = "~/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/derived/simulation/saige/",
  regenie_dir = "~/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/derived/simulation/regenie/",
  output_dir = "~/Desktop/null_simulation_plots/"
)

# List SAIGE files
saige_pattern <- paste0(".*", opt$ancestry, ".*", opt$annotation, "\\.combined\\.tsv\\.gz$")
saige_files <- list.files(opt$saige_dir, pattern = saige_pattern, full.names = TRUE)

cat(paste0("Found ", length(saige_files), " SAIGE files:\n"))
for (f in saige_files) {
  cat(paste0("  - ", basename(f), "\n"))
}
cat("\n")

if (length(saige_files) == 0) {
  stop("No SAIGE files found matching pattern: ", saige_pattern, call.=FALSE)
}

# List REGENIE files
regenie_pattern <- paste0(".*", opt$ancestry, ".*", opt$annotation, "\\.combined\\.regenie\\.tsv\\.gz$")
regenie_files <- list.files(opt$regenie_dir, pattern = regenie_pattern, full.names = TRUE)
regenie_files <- regenie_files[!grepl(regenie_files, pattern="0\\.1")]

cat(paste0("Found ", length(regenie_files), " REGENIE files:\n"))
for (f in regenie_files) {
  cat(paste0("  - ", basename(f), "\n"))
}
cat("\n")

if (length(regenie_files) == 0) {
  stop("No REGENIE files found matching pattern: ", regenie_pattern, call.=FALSE)
}

# Load SAIGE files
cat("Loading SAIGE files...\n")
saige_list <- lapply(saige_files, function(f){
  cat(paste0("  Loading: ", basename(f), "\n"))
  d <- fread(f)
  d[, heritability := stringr::str_extract(basename(f), "(?<=p_)\\d*\\.?\\d+")]
  d[, id := stringr::str_extract(basename(f), "(?<=continuous_)\\d+")]
  d[, annotation := stringr::str_extract(f, "(?<=continuous_\\d\\.).+?(?=\\.combined)")]

  # Restrict to variants with at least one biallelic variant
  d <- d[!is.na(d$DOM.P),]

  cat(paste0("    Heritability: ", unique(d$heritability),
             ", ID: ", unique(d$id),
             ", Variants: ", nrow(d), "\n"))
  return(d)
})

saige_data <- rbindlist(saige_list)
cat(paste0("\nTotal SAIGE variants: ", nrow(saige_data), "\n"))

# Load REGENIE files
cat("\nLoading REGENIE files...\n")
regenie_list <- lapply(regenie_files, function(f){
  cat(paste0("  Loading: ", basename(f), "\n"))
  d <- fread(f)
  d[, heritability := stringr::str_extract(basename(f), "(?<=p_)\\d*\\.?\\d+")]
  d[, id := stringr::str_extract(basename(f), "(?<=continuous_)\\d+")]
  d[, annotation := stringr::str_extract(f, "(?<=continuous_\\d\\.).+?(?=\\.combined)")]

  # Restrict to variants with at least one biallelic variant
  d <- d[!is.na(d$DOM.P),]

  cat(paste0("    Heritability: ", unique(d$heritability),
             ", ID: ", unique(d$id),
             ", Variants: ", nrow(d), "\n"))
  return(d)
})

regenie_data <- rbindlist(regenie_list)
cat(paste0("\nTotal REGENIE variants: ", nrow(regenie_data), "\n"))

# Create unique variant identifiers for merging
# Using CHR/CHROM, POS/GENPOS, and variant ID
cat("\nCreating variant identifiers for merging...\n")
saige_data[, variant_id := paste(CHR, POS, MarkerID, sep=":")]
regenie_data[, variant_id := paste(CHROM, GENPOS, ID, sep=":")]

# Filter for N_HOMALT >= 5 before merging
cat("Filtering for N_HOMALT >= 10...\n")
cat(paste0("SAIGE variants before filtering: ", nrow(saige_data), "\n"))
saige_data <- saige_data[N_HOMALT >= 10,]
cat(paste0("SAIGE variants after filtering: ", nrow(saige_data), "\n"))

cat(paste0("REGENIE variants before filtering: ", nrow(regenie_data), "\n"))
regenie_data <- regenie_data[N_HOMALT >= 10,]
cat(paste0("REGENIE variants after filtering: ", nrow(regenie_data), "\n"))

# Merge SAIGE and REGENIE data
cat("\nMerging SAIGE and REGENIE data...\n")
merged <- merge(
  saige_data,
  regenie_data,
  by = c("variant_id", "heritability", "id"),
  suffixes = c(".saige", ".regenie")
)

cat(paste0("Merged variants: ", nrow(merged), "\n"))
cat(paste0("Unique heritabilities: ", paste(unique(merged$heritability), collapse=", "), "\n"))
cat(paste0("Unique IDs: ", paste(unique(merged$id), collapse=", "), "\n\n"))

if (nrow(merged) == 0) {
  stop("No overlapping variants found between SAIGE and REGENIE!", call.=FALSE)
}

# Calculate 95% confidence intervals
merged[, ADD.CI_lower.saige := ADD.BETA.saige - 1.96 * ADD.SE.saige]
merged[, ADD.CI_upper.saige := ADD.BETA.saige + 1.96 * ADD.SE.saige]
merged[, ADD.CI_lower.regenie := ADD.BETA.regenie - 1.96 * ADD.SE.regenie]
merged[, ADD.CI_upper.regenie := ADD.BETA.regenie + 1.96 * ADD.SE.regenie]

merged[, REC.CI_lower.saige := REC.BETA.saige - 1.96 * REC.SE.saige]
merged[, REC.CI_upper.saige := REC.BETA.saige + 1.96 * REC.SE.saige]
merged[, REC.CI_lower.regenie := REC.BETA.regenie - 1.96 * REC.SE.regenie]
merged[, REC.CI_upper.regenie := REC.BETA.regenie + 1.96 * REC.SE.regenie]

merged[, DOM.CI_lower.saige := DOM.BETA.saige - 1.96 * DOM.SE.saige]
merged[, DOM.CI_upper.saige := DOM.BETA.saige + 1.96 * DOM.SE.saige]
merged[, DOM.CI_lower.regenie := DOM.BETA.regenie - 1.96 * DOM.SE.regenie]
merged[, DOM.CI_upper.regenie := DOM.BETA.regenie + 1.96 * DOM.SE.regenie]

# Calculate plot dimensions
n_heritabilities <- length(unique(merged$heritability))
n_ids <- length(unique(merged$id))
n_facets <- n_heritabilities * n_ids
width <- min(24, max(16, ceiling(sqrt(n_facets)) * 4))
height <- min(24, max(12, ceiling(n_facets / ceiling(sqrt(n_facets))) * 3))

width = 12
height = 10

# ========== P-VALUE COMPARISON PLOTS ==========
cat("\n=== Creating P-value comparison plots ===\n")

# Additive model P-value comparison
cat("Creating Additive P-value comparison plot...\n")
p_add_pval <- ggplot(merged, aes(x = -log10(ADD.P.saige), y = -log10(ADD.P.regenie))) +
  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = "dashed", size = 1) +
  geom_point_rast(alpha = 0.5, size = 1.5) +
  facet_wrap(~ heritability + id, labeller = label_both, scales = "free") +
  ggtitle("Additive Model: P-value Comparison") +
  xlab(TeX("SAIGE: $-\\log_{10}(P)$")) +
  ylab(TeX("REGENIE: $-\\log_{10}(P)$")) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey90", color = "grey90"),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(margin = ggplot2::margin(t = 15)),
    axis.title.y = element_text(margin = ggplot2::margin(r = 15)),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
  )

output_add_pval <- paste0(opt$output_dir, "saige_regenie_comparison_additive_pvalue.pdf")
cat(paste0("Saving to: ", output_add_pval, "\n"))
ggsave(filename = output_add_pval, plot = p_add_pval, width = width, height = height, units = "in", dpi = 300)

# Recessive model P-value comparison
cat("Creating Recessive P-value comparison plot...\n")
p_rec_pval <- ggplot(merged, aes(x = -log10(REC.P.saige), y = -log10(REC.P.regenie))) +
  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = "dashed", size = 1) +
  geom_point_rast(alpha = 0.5, size = 1.5) +
  facet_wrap(~ heritability + id, labeller = label_both, scales = "free") +
  ggtitle("Recessive Model: P-value Comparison") +
  xlab(TeX("SAIGE: $-\\log_{10}(P)$")) +
  ylab(TeX("REGENIE: $-\\log_{10}(P)$")) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey90", color = "grey90"),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(margin = ggplot2::margin(t = 15)),
    axis.title.y = element_text(margin = ggplot2::margin(r = 15)),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
  )

output_rec_pval <- paste0(opt$output_dir, "saige_regenie_comparison_recessive_pvalue.pdf")
cat(paste0("Saving to: ", output_rec_pval, "\n"))
ggsave(filename = output_rec_pval, plot = p_rec_pval, width = width, height = height, units = "in", dpi = 300)

# Dominance model P-value comparison
cat("Creating Dominance P-value comparison plot...\n")
p_dom_pval <- ggplot(merged, aes(x = -log10(DOM.P.saige), y = -log10(DOM.P.regenie))) +
  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = "dashed", size = 1) +
  geom_point_rast(alpha = 0.5, size = 1.5) +
  facet_wrap(~ heritability + id, labeller = label_both, scales = "free") +
  ggtitle("Nonadditive Model: P-value Comparison") +
  xlab(TeX("SAIGE: $-\\log_{10}(P)$")) +
  ylab(TeX("REGENIE: $-\\log_{10}(P)$")) +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey90", color = "grey90"),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(margin = ggplot2::margin(t = 15)),
    axis.title.y = element_text(margin = ggplot2::margin(r = 15)),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
  )

output_dom_pval <- paste0(opt$output_dir, "saige_regenie_comparison_dominance_pvalue.pdf")
cat(paste0("Saving to: ", output_dom_pval, "\n"))
ggsave(filename = output_dom_pval, plot = p_dom_pval, width = width, height = height, units = "in", dpi = 300)

# ========== BETA COMPARISON PLOTS ==========
cat("\n=== Creating BETA comparison plots ===\n")

# Additive model BETA comparison
cat("Creating Additive BETA comparison plot...\n")
p_add_beta <- ggplot(merged, aes(x = ADD.BETA.saige, y = ADD.BETA.regenie)) +
  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = "dashed", size = 1) +
  geom_point_rast(alpha = 0.5, size = 1.5) +
  geom_errorbar(aes(ymin = ADD.CI_lower.regenie, ymax = ADD.CI_upper.regenie),
                alpha = 0.2, width = 0) +
  geom_errorbarh(aes(xmin = ADD.CI_lower.saige, xmax = ADD.CI_upper.saige),
                 alpha = 0.2, height = 0) +
  facet_wrap(~ heritability + id, labeller = label_both, scales = "free") +
  ggtitle("Additive Model: Effect Size (BETA) Comparison with 95% CI") +
  xlab("SAIGE: BETA") +
  ylab("REGENIE: BETA") +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey90", color = "grey90"),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(margin = ggplot2::margin(t = 15)),
    axis.title.y = element_text(margin = ggplot2::margin(r = 15)),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
  )

output_add_beta <- paste0(opt$output_dir, "saige_regenie_comparison_additive_beta.pdf")
cat(paste0("Saving to: ", output_add_beta, "\n"))
ggsave(filename = output_add_beta, plot = p_add_beta, width = width, height = height, units = "in", dpi = 300)

# Recessive model BETA comparison
cat("Creating Recessive BETA comparison plot...\n")
p_rec_beta <- ggplot(merged, aes(x = REC.BETA.saige, y = REC.BETA.regenie)) +
  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = "dashed", size = 1) +
  geom_point_rast(alpha = 0.5, size = 1.5) +
  geom_errorbar(aes(ymin = REC.CI_lower.regenie, ymax = REC.CI_upper.regenie),
                alpha = 0.2, width = 0) +
  geom_errorbarh(aes(xmin = REC.CI_lower.saige, xmax = REC.CI_upper.saige),
                 alpha = 0.2, height = 0) +
  facet_wrap(~ heritability + id, labeller = label_both, scales = "free") +
  ggtitle("Recessive Model: Effect Size (BETA) Comparison with 95% CI") +
  xlab("SAIGE: BETA") +
  ylab("REGENIE: BETA") +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey90", color = "grey90"),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(margin = ggplot2::margin(t = 15)),
    axis.title.y = element_text(margin = ggplot2::margin(r = 15)),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
  )

output_rec_beta <- paste0(opt$output_dir, "saige_regenie_comparison_recessive_beta.pdf")
cat(paste0("Saving to: ", output_rec_beta, "\n"))
ggsave(filename = output_rec_beta, plot = p_rec_beta, width = width, height = height, units = "in", dpi = 300)

# Dominance model BETA comparison
cat("Creating Dominance BETA comparison plot...\n")
p_dom_beta <- ggplot(merged, aes(x = DOM.BETA.saige, y = DOM.BETA.regenie)) +
  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = "dashed", size = 1) +
  geom_point_rast(alpha = 0.5, size = 1.5) +
  geom_errorbar(aes(ymin = DOM.CI_lower.regenie, ymax = DOM.CI_upper.regenie),
                alpha = 0.2, width = 0) +
  geom_errorbarh(aes(xmin = DOM.CI_lower.saige, xmax = DOM.CI_upper.saige),
                 alpha = 0.2, height = 0) +
  facet_wrap(~ heritability + id, labeller = label_both, scales = "free") +
  ggtitle("Nonadditive Model: Effect Size (BETA) Comparison with 95% CI") +
  xlab("SAIGE: BETA") +
  ylab("REGENIE: BETA") +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey90", color = "grey90"),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(margin = ggplot2::margin(t = 15)),
    axis.title.y = element_text(margin = ggplot2::margin(r = 15)),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
  )

output_dom_beta <- paste0(opt$output_dir, "saige_regenie_comparison_dominance_beta.pdf")
cat(paste0("Saving to: ", output_dom_beta, "\n"))
ggsave(filename = output_dom_beta, plot = p_dom_beta, width = width, height = height, units = "in", dpi = 300)

# ========== STANDARD ERROR COMPARISON PLOTS ==========
cat("\n=== Creating Standard Error comparison plots ===\n")

# Additive model SE comparison
cat("Creating Additive SE comparison plot...\n")
p_add_se <- ggplot(merged, aes(x = ADD.SE.saige, y = ADD.SE.regenie)) +
  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = "dashed", size = 1) +
  geom_point_rast(alpha = 0.5, size = 1.5) +
  facet_wrap(~ heritability + id, labeller = label_both, scales = "free") +
  ggtitle("Additive Model: Standard Error Comparison") +
  xlab("SAIGE: SE") +
  ylab("REGENIE: SE") +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey90", color = "grey90"),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(margin = ggplot2::margin(t = 15)),
    axis.title.y = element_text(margin = ggplot2::margin(r = 15)),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
  )

output_add_se <- paste0(opt$output_dir, "saige_regenie_comparison_additive_se.pdf")
cat(paste0("Saving to: ", output_add_se, "\n"))
ggsave(filename = output_add_se, plot = p_add_se, width = width, height = height, units = "in", dpi = 300)

# Recessive model SE comparison
cat("Creating Recessive SE comparison plot...\n")
p_rec_se <- ggplot(merged, aes(x = REC.SE.saige, y = REC.SE.regenie)) +
  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = "dashed", size = 1) +
  geom_point_rast(alpha = 0.5, size = 1.5) +
  facet_wrap(~ heritability + id, labeller = label_both, scales = "free") +
  ggtitle("Recessive Model: Standard Error Comparison") +
  xlab("SAIGE: SE") +
  ylab("REGENIE: SE") +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey90", color = "grey90"),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(margin = ggplot2::margin(t = 15)),
    axis.title.y = element_text(margin = ggplot2::margin(r = 15)),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
  )

output_rec_se <- paste0(opt$output_dir, "saige_regenie_comparison_recessive_se.pdf")
cat(paste0("Saving to: ", output_rec_se, "\n"))
ggsave(filename = output_rec_se, plot = p_rec_se, width = width, height = height, units = "in", dpi = 300)

# Dominance model SE comparison
cat("Creating Dominance SE comparison plot...\n")
p_dom_se <- ggplot(merged, aes(x = DOM.SE.saige, y = DOM.SE.regenie)) +
  geom_abline(intercept = 0, slope = 1, color = 'red', linetype = "dashed", size = 1) +
  geom_point_rast(alpha = 0.5, size = 1.5) +
  facet_wrap(~ heritability + id, labeller = label_both, scales = "free") +
  ggtitle("Nonadditive Model: Standard Error Comparison") +
  xlab("SAIGE: SE") +
  ylab("REGENIE: SE") +
  theme_classic() +
  theme(
    strip.background = element_rect(fill = "grey90", color = "grey90"),
    strip.text = element_text(size = 12, face = "bold"),
    axis.text = element_text(size = 10),
    axis.title = element_text(size = 16, face = "bold"),
    axis.title.x = element_text(margin = ggplot2::margin(t = 15)),
    axis.title.y = element_text(margin = ggplot2::margin(r = 15)),
    plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
  )

output_dom_se <- paste0(opt$output_dir, "saige_regenie_comparison_dominance_se.pdf")
cat(paste0("Saving to: ", output_dom_se, "\n"))
ggsave(filename = output_dom_se, plot = p_dom_se, width = width, height = height, units = "in", dpi = 300)

cat("\n=============================================================\n")
cat("Done! Created the following comparison plots:\n")
cat("P-value comparisons:\n")
cat(paste0("  - ", output_add_pval, "\n"))
cat(paste0("  - ", output_rec_pval, "\n"))
cat(paste0("  - ", output_dom_pval, "\n"))
cat("\nBETA comparisons (with 95% CI):\n")
cat(paste0("  - ", output_add_beta, "\n"))
cat(paste0("  - ", output_rec_beta, "\n"))
cat(paste0("  - ", output_dom_beta, "\n"))
cat("\nStandard Error comparisons:\n")
cat(paste0("  - ", output_add_se, "\n"))
cat(paste0("  - ", output_rec_se, "\n"))
cat(paste0("  - ", output_dom_se, "\n"))
cat("=============================================================\n")
