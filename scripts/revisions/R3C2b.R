#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(dplyr)

# Configuration
DATE <- "251225"
results_dir <- paste0("/Users/flassen/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/data/", DATE, "/power_simulation/")
output_dir <- "/Users/flassen/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/data/power_simulation_plots"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

H2 <- 0.02  # Using highest h2 value for testing
N_CAUSAL <- 100
REPLICATE <- 1  # Set to NULL to combine all replicates

# Load and parse data files
cat("Loading data files...\n")
files <- list.files(results_dir, pattern = "\\.txt\\.gz$", full.names = TRUE)
cat("Found", length(files), "files\n")

all_data_list <- lapply(files, function(f) {
  df <- fread(f, sep = "\t")
  filename <- basename(f)
  
  if (!grepl("arch_", filename)) {
    cat(sprintf("WARNING: Could not find 'arch_' in filename: %s\n", filename))
    return(NULL)
  }
  
  # Extract architecture and model using regex
  param_string <- sub("^.*\\.(arch_.*?)\\.(additive|dominance|recessive)\\..*$", "\\1", filename)
  model <- sub("^.*\\.(additive|dominance|recessive)\\..*$", "\\1", filename)
  param_parts <- strsplit(param_string, "_")[[1]]
  
  # Parse architecture (everything between "arch" and "h2")
  arch_idx <- which(param_parts == "arch")
  h2_idx <- which(param_parts == "h2")
  
  if (length(arch_idx) == 0 || length(h2_idx) == 0) {
    cat(sprintf("WARNING: Could not parse file: %s\n", filename))
    return(NULL)
  }
  
  architecture <- paste(param_parts[(arch_idx + 1):(h2_idx - 1)], collapse = "_")
  heritability <- as.numeric(param_parts[h2_idx + 1]) / 1000000  # 6 decimal places (e.g., 0000050 -> 0.00005)
  n_causal <- as.numeric(param_parts[which(param_parts == "M") + 1])
  replicate <- as.numeric(param_parts[which(param_parts == "rep") + 1])
  
  df$architecture <- architecture
  df$heritability <- heritability
  df$n_causal <- n_causal
  df$replicate <- replicate
  df$model <- model
  
  return(df)
})

all_data <- rbindlist(all_data_list[!sapply(all_data_list, is.null)], fill = TRUE)
all_data$p.value <- as.numeric(all_data$p.value)

cat("Loaded", nrow(all_data), "total results\n")
cat("Architectures:", paste(unique(all_data$architecture), collapse = ", "), "\n")
cat("Heritabilities:", paste(unique(all_data$heritability), collapse = ", "), "\n")

# QQ plot function
make_qq_plot <- function(data, title) {
  data <- data %>%
    filter(!is.na(p.value) & p.value > 0) %>%
    mutate(obs_log10p = -log10(p.value)) %>%
    group_by(model) %>%
    arrange(p.value) %>%
    mutate(
      rank = row_number(),
      exp_log10p = -log10(rank / (n() + 1))
    ) %>%
    ungroup()
  
  ggplot(data, aes(x = exp_log10p, y = obs_log10p, color = model)) +
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
    scale_color_brewer(palette = "Set1", name = "Model")
}

# Helper to filter and plot
plot_architecture <- function(arch, fig_name, title_prefix) {
  data <- all_data %>%
    filter(
      architecture == arch,
      heritability == H2,
      n_causal == N_CAUSAL,
      model %in% c("additive", "dominance", "recessive")
    )
  
  if (!is.null(REPLICATE)) {
    data <- data %>% filter(replicate == REPLICATE)
  }
  
  if (nrow(data) == 0) {
    cat(sprintf("WARNING: No data for %s\n", arch))
    return(NULL)
  }
  
  cat(sprintf("\n=== %s (%d results) ===\n", title_prefix, nrow(data)))
  
  # Build title
  title <- sprintf("%s (h²=%.4f%%, %d causal variants)", title_prefix, H2*100, N_CAUSAL)
  if (grepl("partially_recessive", arch)) {
    het <- as.numeric(sub("^partially_recessive_", "", arch))
    title <- sprintf("Partially Recessive (het=%.1f, h²=%.4f%%, %d causal)", het, H2*100, N_CAUSAL)
  }
  
  p <- make_qq_plot(data, title)
  
  # Save plot
  output_file <- file.path(output_dir, sprintf("%s_h2_%06d_M_%02d.pdf", fig_name, H2*1000000, N_CAUSAL))
  ggsave(output_file, p, width = 8, height = 6)
  cat("Saved:", output_file, "\n")
  
  print(p)
  return(p)
}

# Generate all figures
plot_architecture("additive", "figure_R7_additive", "Additive Architecture")
plot_architecture("recessive", "figure_R8_recessive", "Recessive Architecture")
plot_architecture("partially_recessive_0.1", "figure_R9_partially_recessive_0.1", "Partially Recessive")
plot_architecture("partially_recessive_0.2", "figure_R9b_partially_recessive_0.2", "Partially Recessive")
plot_architecture("partially_recessive_0.3", "figure_R9c_partially_recessive_0.3", "Partially Recessive")
# plot_architecture("overdominant", "figure_R10_overdominant", "Overdominant Architecture")  # Not simulated

cat("\n=== Done ===\n")
cat("Output directory:", output_dir, "\n")
