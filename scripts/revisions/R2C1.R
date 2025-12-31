#!/usr/bin/env Rscript
# Manhattan plot for Olink variant analysis
# Description: Aggregate Manhattan plot for selected proteins from olink_variant analysis

# Load required libraries
library(data.table)
library(ggplot2)
library(dplyr)
library(ggrastr)
library(latex2exp)
library(ggrepel)

# Parameters
data_dir <- "/Users/flassen/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/data/olink_variant_results"
output_dir <- "/Users/flassen/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/reviewer_comments"

# Specific proteins to include
proteins <- c(
  "MMP10",
  "VIT",
  "GIMAP7",
  "RBKS",
  "CCL26",
  "MASP1",
  "LGALS3",
  "MMP8",
  "CD63",
  "ASAH2"
)

# Model types to include (additive, dominance, and recessive)
model_types <- c("additive", "dominance", "recessive")

# P-value threshold for genome-wide significance
gwas_threshold <- 5e-8
suggestive_threshold <- 1e-5

# Create output directory if it doesn't exist
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Function to calculate expected p-values for QQ plot
get_expected_p <- function(pvalues) {
  n <- length(pvalues)
  rank <- rank(pvalues, ties.method = "first")
  expected <- (rank - 0.5) / n
  return(expected)
}

# Function to read and process a single result file
read_variant_file <- function(file_path, protein, chromosome, model_type) {
  cat(sprintf("Reading: %s\n", basename(file_path)))

  # Read the file
  dt <- fread(file_path, header = TRUE)

  # Add metadata columns
  dt[, `:=`(
    protein = protein,
    chr_num = as.integer(gsub("chr", "", CHR))
  )]

  # Select relevant columns including AC_Allele2
  dt <- dt[, .(chr_num, CHR, POS, MarkerID, p.value, AC_Allele2, BETA, SE, protein)]

  return(dt)
}

# Find all files for the specified proteins
cat("Searching for variant files...\n")
all_files <- list.files(data_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)

# Filter files for our proteins and model types
protein_pattern <- paste(proteins, collapse = "|")
model_pattern <- paste(paste0("variant_", model_types), collapse = "|")

filtered_files <- all_files[
  grepl(protein_pattern, all_files) &
  grepl(model_pattern, all_files)
]

cat(sprintf("Found %d files matching criteria\n", length(filtered_files)))

if (length(filtered_files) == 0) {
  stop("No files found for the specified proteins. Please check that files have been downloaded.")
}

# Organize files by model type
results_by_model <- list(
  additive = list(),
  dominance = list(),
  recessive = list()
)

for (file in filtered_files) {
  # Extract metadata from filename
  # Format: UKB.wes.chrX.PROTEIN.eur.af05.variant_MODEL.txt.gz
  basename_file <- basename(file)
  parts <- strsplit(basename_file, "\\.")[[1]]

  chrom <- parts[3]  # e.g., chr1
  protein <- parts[4]  # e.g., MMP10
  model_full <- parts[7]  # e.g., variant_additive, variant_dominance, or variant_recessive
  model_type <- gsub("variant_", "", model_full)

  # Only process if protein is in our list
  if (protein %in% proteins && model_type %in% model_types) {
    tryCatch({
      dt <- read_variant_file(file, protein, chrom, model_type)
      results_by_model[[model_type]][[length(results_by_model[[model_type]]) + 1]] <- dt
    }, error = function(e) {
      cat(sprintf("Error reading %s: %s\n", basename(file), e$message))
    })
  }
}

# Combine results within each model type
cat("Combining results by model type...\n")
additive_combined <- if (length(results_by_model$additive) > 0) {
  rbindlist(results_by_model$additive, use.names = TRUE, fill = TRUE)
} else NULL

dominance_combined <- if (length(results_by_model$dominance) > 0) {
  rbindlist(results_by_model$dominance, use.names = TRUE, fill = TRUE)
} else NULL

recessive_combined <- if (length(results_by_model$recessive) > 0) {
  rbindlist(results_by_model$recessive, use.names = TRUE, fill = TRUE)
} else NULL

# Rename columns with model-specific suffixes
if (!is.null(additive_combined)) {
  setnames(additive_combined,
           c("p.value", "AC_Allele2", "BETA", "SE"),
           c("p.value.add", "ac.add", "beta.add", "se.add"))
}

if (!is.null(dominance_combined)) {
  setnames(dominance_combined,
           c("p.value", "AC_Allele2", "BETA", "SE"),
           c("p.value.dom", "ac.dom", "beta.dom", "se.dom"))
}

if (!is.null(recessive_combined)) {
  setnames(recessive_combined,
           c("p.value", "AC_Allele2", "BETA", "SE"),
           c("p.value.rec", "ac.rec", "beta.rec", "se.rec"))
}

# Merge all three model types
cat("Merging additive, dominance, and recessive results...\n")
combined_results <- NULL

# Start with additive if available
if (!is.null(additive_combined)) {
  combined_results <- additive_combined
}

# Merge dominance
if (!is.null(dominance_combined)) {
  if (is.null(combined_results)) {
    combined_results <- dominance_combined
  } else {
    combined_results <- merge(
      combined_results,
      dominance_combined,
      by = c("chr_num", "CHR", "POS", "MarkerID", "protein"),
      all = TRUE
    )
  }
}

# Merge recessive
if (!is.null(recessive_combined)) {
  if (is.null(combined_results)) {
    combined_results <- recessive_combined
  } else {
    combined_results <- merge(
      combined_results,
      recessive_combined,
      by = c("chr_num", "CHR", "POS", "MarkerID", "protein"),
      all = TRUE
    )
  }
}

#### READY TO FILTER
combined_results 

# Filter: Keep only variants where ac.rec >= 10
if ("ac.rec" %in% names(combined_results)) {
  cat(sprintf("Variants before ac.rec filter: %d\n", nrow(combined_results)))
  combined_results <- combined_results[!is.na(ac.rec) & ac.rec >= 10]
  cat(sprintf("Variants after ac.rec >= 10 filter: %d\n", nrow(combined_results)))
}

cat(sprintf("Total variants: %d\n", nrow(combined_results)))
cat(sprintf("Unique proteins: %d\n", length(unique(combined_results$protein))))
cat(sprintf("Chromosomes represented: %s\n", paste(sort(unique(combined_results$chr_num)), collapse = ", ")))

# Remove variants where ALL p-values are missing
combined_results <- combined_results[
  !is.na(p.value.add) | !is.na(p.value.dom) | !is.na(p.value.rec)
]

# For plotting, we need to reshape to long format
cat("Reshaping data for plotting...\n")
plot_data_list <- list()

if ("p.value.add" %in% names(combined_results)) {
  add_data <- combined_results[!is.na(p.value.add) & p.value.add > 0 & p.value.add <= 1,
                                .(chr_num, CHR, POS, MarkerID, protein, p.value = p.value.add, model_type = "additive")]
  plot_data_list[[length(plot_data_list) + 1]] <- add_data
}

if ("p.value.dom" %in% names(combined_results)) {
  dom_data <- combined_results[!is.na(p.value.dom) & p.value.dom > 0 & p.value.dom <= 1,
                                .(chr_num, CHR, POS, MarkerID, protein, p.value = p.value.dom, model_type = "dominance")]
  plot_data_list[[length(plot_data_list) + 1]] <- dom_data
}

if ("p.value.rec" %in% names(combined_results)) {
  rec_data <- combined_results[!is.na(p.value.rec) & p.value.rec > 0 & p.value.rec <= 1,
                                .(chr_num, CHR, POS, MarkerID, protein, p.value = p.value.rec, model_type = "recessive")]
  plot_data_list[[length(plot_data_list) + 1]] <- rec_data
}

plot_data_long <- rbindlist(plot_data_list, use.names = TRUE, fill = TRUE)

# Calculate -log10(p-value)
plot_data_long[, neglog10p := -log10(p.value)]

# Create cumulative position for Manhattan plot
# Calculate chromosome lengths and cumulative positions
chr_info <- plot_data_long[, .(
  chr_len = max(POS),
  n_variants = .N
), by = chr_num][order(chr_num)]

chr_info[, cumul_start := cumsum(c(0, chr_len[-.N]))]
chr_info[, chr_center := cumul_start + chr_len / 2]

# Merge cumulative positions back to plot data
plot_data_long <- merge(
  plot_data_long,
  chr_info[, .(chr_num, cumul_start)],
  by = "chr_num"
)

plot_data_long[, pos_cumul := POS + cumul_start]

# Print summary statistics
cat("\n=== Summary Statistics ===\n")

# Find significant variants across all models
sig_variants_add <- if ("p.value.add" %in% names(combined_results)) {
  combined_results[!is.na(p.value.add) & p.value.add < gwas_threshold]
} else data.table()

sig_variants_dom <- if ("p.value.dom" %in% names(combined_results)) {
  combined_results[!is.na(p.value.dom) & p.value.dom < gwas_threshold]
} else data.table()

sig_variants_rec <- if ("p.value.rec" %in% names(combined_results)) {
  combined_results[!is.na(p.value.rec) & p.value.rec < gwas_threshold]
} else data.table()

cat(sprintf("Genome-wide significant variants (p < %.0e):\n", gwas_threshold))
cat(sprintf("  Additive: %d\n", nrow(sig_variants_add)))
cat(sprintf("  Dominance: %d\n", nrow(sig_variants_dom)))
cat(sprintf("  Recessive: %d\n", nrow(sig_variants_rec)))

# Print top variants from plot_data_long
if (nrow(plot_data_long) > 0) {
  cat("\nTop significant variants (across all models):\n")
  top_variants <- plot_data_long[order(p.value)][1:min(10, nrow(plot_data_long))]
  print(top_variants[, .(protein, CHR, POS, MarkerID, p.value, model_type)])
}

# Create Manhattan plot
cat("\nGenerating Manhattan plot...\n")

# Create color palette for chromosomes
chr_colors <- rep(c("#2E4057", "#66A182"), length.out = length(unique(plot_data_long$chr_num)))

# Manhattan plot
dominance_dt <- plot_data_long[plot_data_long$model_type=="dominance",]
dominance_dt[, expt.p := get_expected_p(p.value)]
dominance_dt[, FDR := p.adjust(p.value, method="fdr")]
bonf_p <- 0.05/nrow(dominance_dt)
#dominance_dt[, label := ifelse(p.value<bonf_p, MarkerID, NA )]
dominance_dt[, label := ifelse(FDR<0.001, MarkerID, NA )]
dominance_dt[p.value<bonf_p,]

dominance_dt$label[!is.na(dominance_dt$label)]

ggplot(dominance_dt, aes(x = pos_cumul, y = neglog10p, label=label)) +
  geom_point(aes(color = as.factor(chr_num)), size = 1.2) +
  geom_label_repel() +
  scale_color_manual(values = chr_colors, guide = "none") +
  scale_x_continuous(
    label = chr_info$chr_num,
    breaks = chr_info$chr_center,
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(expand = c(0.01, 0)) +
  geom_hline(yintercept = -log10(bonf_p), color = "red", linetype = "dashed", linewidth = 0.5) +
  labs(
    #title = "Manhattan Plot: Olink Variant Analysis",
    x = "Chromosome",
    y = expression(-log[10](italic(p)))
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40")
  )

# Save the plot
output_file <- file.path(output_dir, "manhattan_plot_olink_variant.pdf")
ggsave(
  output_file,
  width = 6,
  height = 3,
  dpi = 300
)

ggplot(dominance_dt, aes(x = -log10(expt.p), y=-log10(p.value), label=label)) + 
  geom_point() +
  geom_label_repel(box.padding = 0.8, point.padding = 1.5, 
                   size=3, max.overlaps = Inf, alpha=0.6) +
  geom_abline(linetype="dashed") +
  xlab(TeX("$-log_{10}$(Expt. $\\textit{P}$-value)")) +
  ylab(TeX("$-log_{10}$(Obs. $\\textit{P}$-value)")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size=10),
    axis.text.y = element_text(angle = 0, hjust = 0.5, size=10),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40")
  )

# Save the plot
output_file <- file.path(output_dir, "qq_plot_olink_variant.pdf")
ggsave(
  output_file,
  width = 4,
  height = 2.5,
  dpi = 300
)

#### recessive

# Manhattan plot
recessive_dt <- plot_data_long[plot_data_long$model_type=="recessive",]
recessive_dt[, expt.p := get_expected_p(p.value)]
recessive_dt[, FDR := p.adjust(p.value, method="fdr")]
bonf_p_add <- 0.05/nrow(recessive_dt)
#recessive_dt[, label := ifelse(p.value<bonf_p_add, MarkerID, NA )]
recessive_dt[, label := ifelse(FDR<0.001, MarkerID, NA )]

recessive_dt$label[!is.na(recessive_dt$label)]

recessive_dt$p.value[recessive_dt$p.value < 1e-20] <- 1e-20

ggplot(recessive_dt, aes(x = pos_cumul, y = neglog10p, label=label)) +
  geom_point(aes(color = as.factor(chr_num)), size = 1.2) +
  #geom_label_repel() +
  scale_color_manual(values = chr_colors, guide = "none") +
  scale_x_continuous(
    label = chr_info$chr_num,
    breaks = chr_info$chr_center,
    expand = c(0.01, 0.01)
  ) +
  scale_y_continuous(expand = c(0.01, 0)) +
  geom_hline(yintercept = -log10(bonf_p_add), color = "red", linetype = "dashed", linewidth = 0.5) +
  labs(
    #title = "Manhattan Plot: Olink Variant Analysis",
    x = "Chromosome",
    y = expression(-log[10](italic(p)))
  ) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.x = element_blank(),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40")
  )

# Save the plot
output_file <- file.path(output_dir, "manhattan_plot_olink_variant_add.pdf")
ggsave(
  output_file,
  width = 6,
  height = 3,
  dpi = 300
)

ggplot(recessive_dt, aes(x = -log10(expt.p), y=-log10(p.value), label=label)) + 
  geom_point() +
  geom_label_repel(box.padding = 0.8, point.padding = 1.5, 
                   size=3, max.overlaps = Inf, alpha=0.6) +
  geom_abline(linetype="dashed") +
  xlab(TeX("$-log_{10}$(Expt. $\\textit{P}$-value)")) +
  ylab(TeX("$-log_{10}$(Obs. $\\textit{P}$-value)")) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5, size=10),
    axis.text.y = element_text(angle = 0, hjust = 0.5, size=10),
    plot.title = element_text(face = "bold", size = 14),
    plot.subtitle = element_text(size = 10, color = "gray40")
  )

# Save the plot
output_file <- file.path(output_dir, "qq_plot_olink_variant_add.pdf")
ggsave(
  output_file,
  width = 4,
  height = 2.5,
  dpi = 300
)




