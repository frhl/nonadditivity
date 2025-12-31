#!/usr/bin/env Rscript
# Simple Power Analysis - Step-by-step QQ Plots for Reviewer
# Author: Frederik Lassen

library(data.table)
library(ggplot2)
library(dplyr)

# ============================================================================
# STEP 1: CONFIGURATION
# ============================================================================

results_dir <- "/Users/flassen/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/data/saige_power_analysis"
output_dir <- "/Users/flassen/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/data/power_simulation_plots"

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Choose parameters for plots
H2 <- 0.01      # heritability (0.01 or 0.20)
N_CAUSAL <- 10  # number of causal variants (10 or 20) - NOTE: M=20 has very sparse data
REPLICATE <- 4  # which replicate to use (set to NULL to combine all replicates)
                # NOTE: For h2=0.01, M=10: use replicate 4 (most complete)

# ============================================================================
# STEP 2: LOAD ALL DATA
# ============================================================================

cat("Loading data files...\n")

# Get all result files
files <- list.files(results_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)
cat("Found", length(files), "files\n")

# Load all files and combine
all_data_list <- lapply(files, function(f) {

  # Read file
  df <- fread(f, sep = "\t")

  # Parse filename to get metadata
  # Example: UKB.auto.arch_additive_h2_001_M_10_rep_3.additive.eur.txt.gz
  filename <- basename(f)
  parts <- strsplit(filename, "\\.")[[1]]

  # Extract parameters from the 3rd part
  param_string <- parts[3]  # e.g., "arch_additive_h2_001_M_10_rep_3"
  param_parts <- strsplit(param_string, "_")[[1]]

  # Get architecture
  arch_idx <- which(param_parts == "arch")
  architecture <- param_parts[arch_idx + 1]

  # Get heritability
  h2_idx <- which(param_parts == "h2")
  heritability <- as.numeric(param_parts[h2_idx + 1]) / 100

  # Get number of causal variants
  M_idx <- which(param_parts == "M")
  n_causal <- as.numeric(param_parts[M_idx + 1])

  # Get replicate
  rep_idx <- which(param_parts == "rep")
  replicate <- as.numeric(param_parts[rep_idx + 1])

  # Get model type (additive, dominance, recessive)
  model <- parts[4]

  # Add metadata columns
  df$architecture <- architecture
  df$heritability <- heritability
  df$n_causal <- n_causal
  df$replicate <- replicate
  df$model <- model

  return(df)
})

# Combine all data
all_data <- rbindlist(all_data_list, fill = TRUE)

# Convert p.value to numeric (in case it was read as character)
all_data$p.value <- as.numeric(all_data$p.value)

cat("Loaded", nrow(all_data), "total results\n")
cat("\nArchitectures:", paste(unique(all_data$architecture), collapse = ", "), "\n")
cat("Heritabilities:", paste(unique(all_data$heritability), collapse = ", "), "\n")
cat("Models:", paste(unique(all_data$model), collapse = ", "), "\n")

# ============================================================================
# STEP 3: HELPER FUNCTION FOR QQ PLOTS
# ============================================================================

make_qq_plot <- function(data, title) {

  # Remove NA p-values
  data <- data %>% filter(!is.na(p.value) & p.value > 0)

  # Calculate observed -log10(p)
  data$obs_log10p <- -log10(data$p.value)

  # For each model, calculate expected -log10(p)
  # Sort by ascending p-value (smallest first) for proper QQ plot
  data <- data %>%
    group_by(model) %>%
    arrange(p.value) %>%
    mutate(
      rank = row_number(),
      n_total = n(),
      exp_log10p = -log10(rank / (n_total + 1))
    ) %>%
    ungroup()

  # Get axis limits
  max_val <- max(c(data$obs_log10p, data$exp_log10p), na.rm = TRUE)

  # Create plot
  p <- ggplot(data, aes(x = exp_log10p, y = obs_log10p, color = model)) +
    geom_point(size = 1.5, alpha = 0.6) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
    labs(
      title = title,
      x = expression(Expected~~-log[10](italic(P))),
      y = expression(Observed~~-log[10](italic(P)))
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      legend.position = "bottom",
      panel.grid.minor = element_blank()
    ) +
    #coord_cartesian(xlim = c(0, max_val * 1.05), ylim = c(0, max_val * 1.05)) +
    scale_color_brewer(palette = "Set1", name = "Model")

  return(p)
}

# ============================================================================
# STEP 4: FIGURE R7 - ADDITIVE ARCHITECTURE
# ============================================================================

cat("\n=== Creating Figure R7: Additive Architecture ===\n")

# Filter data for additive architecture
additive_data <- all_data %>%
  filter(
    architecture == "additive",
    heritability == H2,
    n_causal == N_CAUSAL,
    model %in% c("additive", "dominance", "recessive")
  )

# Filter by replicate if specified
if (!is.null(REPLICATE)) {
  additive_data <- additive_data %>% filter(replicate == REPLICATE)
  cat("Using replicate", REPLICATE, "\n")
}

cat("Using", nrow(additive_data), "results\n")

# Create plot
plot_title <- sprintf("Additive Architecture (h²=%.1f%%, %d causal variants)", H2*100, N_CAUSAL)
p_additive <- make_qq_plot(additive_data, plot_title)

# Save plot
output_file <- file.path(output_dir, sprintf("figure_R7_additive_h2_%03d_M_%02d.pdf", H2*100, N_CAUSAL))
#ggsave(output_file, p_additive, width = 8, height = 6)
cat("Saved:", output_file, "\n")

# Display plot
print(p_additive)

# ============================================================================
# STEP 5: FIGURE R8 - RECESSIVE ARCHITECTURE
# ============================================================================

cat("\n=== Creating Figure R8: Recessive Architecture ===\n")

# Filter data for recessive architecture
recessive_data <- all_data %>%
  filter(
    architecture == "recessive",
    heritability == H2,
    n_causal == N_CAUSAL,
    model %in% c("additive", "dominance", "recessive")
  )

# Filter by replicate if specified
if (!is.null(REPLICATE)) {
  recessive_data <- recessive_data %>% filter(replicate == REPLICATE)
  cat("Using replicate", REPLICATE, "\n")
}

cat("Using", nrow(recessive_data), "results\n")

# Create plot
plot_title <- sprintf("Recessive Architecture (h²=%.1f%%, %d causal variants)", H2*100, N_CAUSAL)
p_recessive <- make_qq_plot(recessive_data, plot_title)

# Save plot
output_file <- file.path(output_dir, sprintf("figure_R8_recessive_h2_%03d_M_%02d.pdf", H2*100, N_CAUSAL))
ggsave(output_file, p_recessive, width = 8, height = 6)
cat("Saved:", output_file, "\n")

# Display plot
print(p_recessive)

# ============================================================================
# STEP 6: FIGURE R9 - PARTIALLY RECESSIVE ARCHITECTURE
# ============================================================================

cat("\n=== Creating Figure R9: Partially Recessive Architecture ===\n")

# Filter data for partially recessive architecture
partially_data <- all_data %>%
  filter(
    architecture == "partially",
    heritability == H2,
    n_causal == N_CAUSAL,
    model %in% c("additive", "dominance", "recessive")
  )

# Filter by replicate if specified
if (!is.null(REPLICATE)) {
  partially_data <- partially_data %>% filter(replicate == REPLICATE)
  cat("Using replicate", REPLICATE, "\n")
}

cat("Using", nrow(partially_data), "results\n")

# Create plot
plot_title <- sprintf("Partially Recessive Architecture (h²=%.1f%%, %d causal variants)", H2*100, N_CAUSAL)
p_partially <- make_qq_plot(partially_data, plot_title)

# Save plot
output_file <- file.path(output_dir, sprintf("figure_R9_partially_recessive_h2_%03d_M_%02d.pdf", H2*100, N_CAUSAL))
#ggsave(output_file, p_partially, width = 8, height = 6)
cat("Saved:", output_file, "\n")

# Display plot
print(p_partially)

# ============================================================================
# STEP 7: FIGURE R10 - OVERDOMINANT ARCHITECTURE
# ============================================================================

cat("\n=== Creating Figure R10: Overdominant Architecture ===\n")

# Filter data for overdominant architecture
overdominant_data <- all_data %>%
  filter(
    architecture == "overdominant",
    heritability == H2,
    n_causal == N_CAUSAL,
    model %in% c("additive", "dominance", "recessive")
  )

# Filter by replicate if specified
if (!is.null(REPLICATE)) {
  overdominant_data <- overdominant_data %>% filter(replicate == REPLICATE)
  cat("Using replicate", REPLICATE, "\n")
}

cat("Using", nrow(overdominant_data), "results\n")

# Create plot
plot_title <- sprintf("Overdominant Architecture (h²=%.1f%%, %d causal variants)", H2*100, N_CAUSAL)
p_overdominant <- make_qq_plot(overdominant_data, plot_title)

# Save plot
output_file <- file.path(output_dir, sprintf("figure_R10_overdominant_h2_%03d_M_%02d.pdf", H2*100, N_CAUSAL))
#ggsave(output_file, p_overdominant, width = 8, height = 6)
cat("Saved:", output_file, "\n")

# Display plot
print(p_overdominant)

# ============================================================================
# DONE
# ============================================================================

cat("\n=== All plots created! ===\n")
cat("Output directory:", output_dir, "\n")
