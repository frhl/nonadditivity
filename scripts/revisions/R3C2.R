#!/usr/bin/env Rscript

library(data.table)
library(ggplot2)
library(dplyr)
library(tidyr)
library(scales)
library(ggrastr)
library(latex2exp)
library(patchwork)

# Configuration
#DATE <- "251226" # eur
DATE <- "251227" # afr
results_dir <- paste0("/Users/flassen/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/data/", DATE, "/power_simulation/")
#output_dir <- "~/Downloads/simulation/" # eur
output_dir <- "~/Downloads/simulation/afr/" # afr
output_dir <- path.expand(output_dir)  # Expand ~ to full path
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Power curve settings
N_CAUSAL <- 100
REPLICATE <- 1  # Set to NULL to combine all replicates
ALPHA <- c(0.05, 0.01, 0.001, 5e-8)  # Significance thresholds to test

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

# QQ plot function (faceted by architecture)
make_qq_plots_faceted <- function(data, h2_val, filename_suffix) {

  plot_data <- data %>%
    filter(
      heritability == h2_val,
      n_causal == N_CAUSAL,
      model %in% c("additive", "dominance", "recessive"),
      !is.na(p.value) & p.value > 0
    ) %>%
    mutate(
      arch_label = case_when(
        architecture == "additive" ~ "Additive [0,1,2]",
        architecture == "recessive" ~ "Recessive [0,0,2]",
        architecture == "partially_recessive_0.05" ~ "Partially Recessive [0,0.05,2]",
        architecture == "partially_recessive_0.1" ~ "Partially Recessive [0,0.1,2]",
        architecture == "partially_recessive_0.2" ~ "Partially Recessive [0,0.2,2]",
        architecture == "partially_recessive_0.3" ~ "Partially Recessive [0,0.3,2]",
        TRUE ~ architecture
      ),
      # Order the facets: Recessive -> Partially Recessive (low to high) -> Additive
      arch_label = factor(arch_label, levels = c(
        "Recessive [0,0,2]",
        "Partially Recessive [0,0.05,2]",
        "Partially Recessive [0,0.1,2]",
        "Partially Recessive [0,0.2,2]",
        "Partially Recessive [0,0.3,2]",
        "Additive [0,1,2]"
      )),
      obs_log10p = -log10(p.value)
    ) %>%
    group_by(arch_label, model) %>%
    arrange(p.value) %>%
    mutate(
      rank = row_number(),
      exp_log10p = -log10(rank / (n() + 1))
    ) %>%
    ungroup()

  if (!is.null(REPLICATE)) {
    plot_data <- plot_data %>% filter(replicate == REPLICATE)
  }

  if (nrow(plot_data) == 0) {
    cat("WARNING: No data for QQ plots\n")
    return(NULL)
  }

  p <- ggplot(plot_data, aes(x = exp_log10p, y = obs_log10p, color = model)) +
    geom_point_rast(size = 1.2, shape = 19, raster.dpi = 300) +
    geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray30", linewidth = 0.5) +
    facet_wrap(~ arch_label, ncol = 3, scales = "free") +
    scale_color_manual(
      name = "Encoding",
      values = c("additive" = "#E31A1C", "dominance" = "#1F78B4", "recessive" = "#33A02C"),
      labels = c("additive" = "Additive", "dominance" = "Nonadditive", "recessive" = "Recessive")
    ) +
    labs(
      title = "QQ Plots for Simulated Traits",
      subtitle = TeX(sprintf("$h^2 = %.2f\\%%$, %d causal variants per trait", h2_val * 100, N_CAUSAL)),
      x = TeX("Expected $-\\log_{10}(P)$"),
      y = TeX("Observed $-\\log_{10}(P)$")
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 3)),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40", margin = margin(b = 12)),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10),
      legend.key.width = unit(1.5, "cm"),
      panel.grid.major = element_line(color = "gray92", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 10.5, margin = margin(3, 3, 3, 3)),
      axis.title.x = element_text(face = "bold", size = 11, margin = margin(t = 8)),
      axis.title.y = element_text(face = "bold", size = 11, margin = margin(r = 8)),
      axis.text.x = element_text(size = 9.5, color = "gray20"),
      axis.text.y = element_text(size = 9.5, color = "gray20"),
      axis.ticks = element_line(color = "gray70", linewidth = 0.3),
      axis.ticks.length = unit(2, "pt"),
      plot.margin = margin(8, 8, 8, 8)
    )

  # Save plot
  output_file <- file.path(output_dir, sprintf("qq_plots_%s_h2_%06d_M_%02d.pdf", filename_suffix, h2_val*1000000, N_CAUSAL))
  ggsave(output_file, p, width = 11, height = 6.5, dpi=300)
  cat(sprintf("Saved QQ plot: %s\n", output_file))

  return(p)
}

# QQ plot grid function (multiple heritabilities in a grid layout)
make_qq_grid <- function(data, h2_values, architecture_filter, filename_suffix) {

  plot_list <- lapply(h2_values, function(h2_val) {

    plot_data <- data %>%
      filter(
        heritability == h2_val,
        architecture == architecture_filter,
        n_causal == N_CAUSAL,
        model %in% c("additive", "dominance", "recessive"),
        !is.na(p.value) & p.value > 0
      ) %>%
      mutate(obs_log10p = -log10(p.value)) %>%
      group_by(model) %>%
      arrange(p.value) %>%
      mutate(
        rank = row_number(),
        exp_log10p = -log10(rank / (n() + 1))
      ) %>%
      ungroup()

    if (!is.null(REPLICATE)) {
      plot_data <- plot_data %>% filter(replicate == REPLICATE)
    }

    if (nrow(plot_data) == 0) {
      return(NULL)
    }

    # Create individual QQ plot
    p <- ggplot(plot_data, aes(x = exp_log10p, y = obs_log10p, color = model)) +
      geom_point_rast(size = 0.8, shape = 19, raster.dpi = 300) +
      geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "gray30", linewidth = 0.3) +
      scale_color_manual(
        name = "Encoding",
        values = c("additive" = "#E31A1C", "dominance" = "#1F78B4", "recessive" = "#33A02C"),
        labels = c("additive" = "Additive", "dominance" = "Nonadditive", "recessive" = "Recessive")
      ) +
      labs(
        title = TeX(sprintf("$h^2 = %.2f\\%%$", h2_val * 100)),
        x = TeX("Expected $-\\log_{10}(P)$"),
        y = TeX("Observed $-\\log_{10}(P)$")
      ) +
      theme_bw(base_size = 8) +
      theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 9, margin = margin(b = 2)),
        legend.position = "none",  # Remove individual legends
        panel.grid.major = element_line(color = "gray92", linewidth = 0.25),
        panel.grid.minor = element_blank(),
        panel.border = element_rect(color = "gray70", linewidth = 0.4),
        axis.title.x = element_text(size = 7, margin = margin(t = 3)),
        axis.title.y = element_text(size = 7, margin = margin(r = 3)),
        axis.text.x = element_text(size = 6.5, color = "gray20"),
        axis.text.y = element_text(size = 6.5, color = "gray20"),
        axis.ticks = element_line(color = "gray70", linewidth = 0.2),
        axis.ticks.length = unit(1.5, "pt"),
        plot.margin = margin(3, 3, 3, 3)
      )

    return(p)
  })

  # Remove NULL plots
  plot_list <- plot_list[!sapply(plot_list, is.null)]

  if (length(plot_list) == 0) {
    cat("WARNING: No plots generated for QQ grid\n")
    return(NULL)
  }

  # Determine grid layout based on number of plots
  n_plots <- length(plot_list)
  ncol <- min(4, n_plots)  # Max 4 columns
  nrow <- ceiling(n_plots / ncol)

  # Create architecture label for title
  arch_label <- case_when(
    architecture_filter == "additive" ~ "Additive [0,1,2]",
    architecture_filter == "recessive" ~ "Recessive [0,0,2]",
    grepl("partially_recessive", architecture_filter) ~ {
      het <- sub("partially_recessive_", "", architecture_filter)
      sprintf("Partially Recessive [0,%s,2]", het)
    },
    TRUE ~ architecture_filter
  )

  # Combine plots with patchwork
  combined_plot <- wrap_plots(plot_list, ncol = ncol) +
    plot_layout(guides = "collect") +
    plot_annotation(
      title = sprintf("QQ Plots: %s Architecture", arch_label),
      subtitle = sprintf("%d causal variants per trait", N_CAUSAL),
      theme = theme(
        plot.title = element_text(hjust = 0.5, face = "bold", size = 12, margin = margin(b = 2)),
        plot.subtitle = element_text(hjust = 0.5, size = 9, color = "gray40", margin = margin(b = 8))
      )
    ) &
    theme(legend.position = "bottom")

  # Save plot (A4 width, adjusted height for square plots)
  # With 4 columns and typical margins, each plot width ≈ 2 inches
  # Setting height to make plots roughly square
  plot_height <- (nrow * 2) + 1.5  # 2 inches per row + space for title/legend
  output_file <- file.path(output_dir, sprintf("qq_grid_%s_M_%02d.pdf", filename_suffix, N_CAUSAL))
  ggsave(output_file, combined_plot, width = 8.27, height = plot_height, dpi = 300)
  cat(sprintf("Saved QQ grid: %s\n", output_file))

  return(combined_plot)
}

# =============================================================================
# QQ PLOTS (faceted by architecture)
# =============================================================================

cat("\n=============================================================================\n")
cat("GENERATING QQ PLOTS\n")
cat("=============================================================================\n")

# Generate QQ plots for highest heritability value
qq_h2_val <- 0.02  # You can change this to any h2 value in your data
qq1 <- make_qq_plots_faceted(all_data, h2_val = qq_h2_val, filename_suffix = "all_architectures")

# Optionally generate for a lower heritability too
# qq2 <- make_qq_plots_faceted(all_data, h2_val = 0.005, filename_suffix = "all_architectures")

print(qq1)

# =============================================================================
# QQ PLOTS GRID (multiple h2 values in grid layout - fits on A4 page)
# =============================================================================

cat("\n=============================================================================\n")
cat("GENERATING QQ PLOT GRIDS\n")
cat("=============================================================================\n")

# Define h2 values to include in grid
#h2_values_for_grid <- c(0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.0075, 0.01, 0.015, 0.02)
h2_values_for_grid <- c(0.0005, 0.001, 0.002, 0.003, 0.005, 0.0075, 0.01, 0.02)



# Generate grid for recessive architecture
qq_grid_recessive <- make_qq_grid(
  data = all_data,
  h2_values = h2_values_for_grid,
  architecture_filter = "recessive",
  filename_suffix = "recessive"
)

# Generate grid for additive architecture
qq_grid_additive <- make_qq_grid(
  data = all_data,
  h2_values = h2_values_for_grid,
  architecture_filter = "additive",
  filename_suffix = "additive"
)

# Optionally generate for partially recessive
 qq_grid_partial <- make_qq_grid(
   data = all_data,
   h2_values = h2_values_for_grid,
   architecture_filter = "partially_recessive_0.1",
   filename_suffix = "partially_recessive_0.1"
 )
 
 # Optionally generate for partially recessive
 qq_grid_partial <- make_qq_grid(
   data = all_data,
   h2_values = h2_values_for_grid,
   architecture_filter = "partially_recessive_0.2",
   filename_suffix = "partially_recessive_0.2"
 )

print(qq_grid_recessive)

# =============================================================================
# POWER CURVE ANALYSIS
# =============================================================================

cat("\n=============================================================================\n")
cat("POWER CURVE ANALYSIS\n")
cat("=============================================================================\n\n")

# Load metadata to identify causal variants
metadata_file <- file.path(results_dir, "metadata.tsv.gz")
if (!file.exists(metadata_file)) {
  stop(sprintf("Metadata file not found: %s\nPlease download it from DNAnexus.", metadata_file))
}

cat("Loading causal variant metadata...\n")
metadata <- fread(metadata_file)
cat(sprintf("  Loaded %d causal variant entries\n", nrow(metadata)))

# Filter to data for power analysis (only causal variants)
# Join with metadata to get only causal variants
power_data_raw <- all_data %>%
  filter(
    n_causal == N_CAUSAL,
    model %in% c("additive", "dominance", "recessive")
  )

if (!is.null(REPLICATE)) {
  power_data_raw <- power_data_raw %>% filter(replicate == REPLICATE)
}

cat(sprintf("Total tests (all variants): %d\n", nrow(power_data_raw)))

# Create a lookup of phenotype -> causal variants
causal_lookup <- metadata %>%
  select(phenotype_name, variant_id) %>%
  distinct()

# For each result, check if it's a causal variant for that phenotype
# Extract phenotype from filename parsing
# NOTE: h2 is formatted with 7 digits (6 decimal places * 1000000, then remove decimal)
power_data_raw$phenotype_name <- paste0(
  "arch_", power_data_raw$architecture,
  "_h2_", sprintf("%07d", round(power_data_raw$heritability * 1000000)),
  "_M_", power_data_raw$n_causal,
  "_rep_", power_data_raw$replicate
)

# Mark causal variants
power_data <- power_data_raw %>%
  left_join(
    causal_lookup %>% mutate(is_causal = TRUE),
    by = c("phenotype_name", "MarkerID" = "variant_id")
  ) %>%
  mutate(is_causal = ifelse(is.na(is_causal), FALSE, is_causal))

n_causal_tests <- sum(power_data$is_causal)
n_null_tests <- sum(!power_data$is_causal)

cat(sprintf("\nCausal variant tests: %d (%.1f%%)\n", n_causal_tests, 100*n_causal_tests/nrow(power_data)))
cat(sprintf("Null variant tests: %d (%.1f%%)\n", n_null_tests, 100*n_null_tests/nrow(power_data)))
cat(sprintf("Architectures: %s\n", paste(unique(power_data$architecture), collapse = ", ")))
cat(sprintf("Heritabilities: %s\n", paste(sort(unique(power_data$heritability)), collapse = ", ")))
cat(sprintf("Models tested: %s\n\n", paste(unique(power_data$model), collapse = ", ")))

# Calculate power at different significance thresholds
calculate_power <- function(pvals, alpha) {
  sum(pvals < alpha, na.rm = TRUE) / length(pvals[!is.na(pvals)])
}

# Compute power for CAUSAL variants only
power_summary_causal <- power_data %>%
  filter(is_causal == TRUE) %>%
  group_by(architecture, heritability, model) %>%
  summarise(
    n_tests = n(),
    power_0.05 = calculate_power(p.value, 0.05),
    power_0.01 = calculate_power(p.value, 0.01),
    power_0.001 = calculate_power(p.value, 0.001),
    power_5e8 = calculate_power(p.value, 5e-8),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = starts_with("power_"),
    names_to = "alpha_level",
    values_to = "power",
    names_prefix = "power_"
  ) %>%
  mutate(
    alpha = case_when(
      alpha_level == "0.05" ~ 0.05,
      alpha_level == "0.01" ~ 0.01,
      alpha_level == "0.001" ~ 0.001,
      alpha_level == "5e8" ~ 5e-8
    ),
    alpha_label = case_when(
      alpha_level == "0.05" ~ "$\\alpha = 0.05$",
      alpha_level == "0.01" ~ "$\\alpha = 0.01$",
      alpha_level == "0.001" ~ "$\\alpha = 0.001$",
      alpha_level == "5e8" ~ "$\\alpha = 5 \\times 10^{-8}$"
    )
  )

# Also compute Type I error rate for NULL variants
type1_error <- power_data %>%
  filter(is_causal == FALSE) %>%
  group_by(architecture, heritability, model) %>%
  summarise(
    n_tests = n(),
    type1_0.05 = calculate_power(p.value, 0.05),
    type1_0.01 = calculate_power(p.value, 0.01),
    type1_0.001 = calculate_power(p.value, 0.001),
    type1_5e8 = calculate_power(p.value, 5e-8),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = starts_with("type1_"),
    names_to = "alpha_level",
    values_to = "type1_error",
    names_prefix = "type1_"
  ) %>%
  mutate(
    alpha = case_when(
      alpha_level == "0.05" ~ 0.05,
      alpha_level == "0.01" ~ 0.01,
      alpha_level == "0.001" ~ 0.001,
      alpha_level == "5e8" ~ 5e-8
    )
  )

# Print power summary table (causal variants)
cat("\n=== POWER (Causal Variants, α = 0.05) ===\n")
power_0.05_table <- power_summary_causal %>%
  filter(alpha == 0.05) %>%
  select(architecture, heritability, model, power) %>%
  pivot_wider(names_from = model, values_from = power)
print(power_0.05_table, n = Inf)

# Print Type I error summary (null variants)
cat("\n=== TYPE I ERROR (Null Variants, α = 0.05) ===\n")
type1_0.05_table <- type1_error %>%
  filter(alpha == 0.05) %>%
  select(architecture, heritability, model, type1_error) %>%
  pivot_wider(names_from = model, values_from = type1_error)
print(type1_0.05_table, n = Inf)

# Create power curves faceted by architecture
make_power_curves <- function(data, alpha_val, alpha_label, filename_suffix) {

  plot_data <- data %>%
    filter(alpha == alpha_val) %>%
    mutate(
      arch_label = case_when(
        architecture == "additive" ~ "Additive [0,1,2]",
        architecture == "recessive" ~ "Recessive [0,0,2]",
        architecture == "partially_recessive_0.05" ~ "Partially Recessive [0,0.05,2]",
        architecture == "partially_recessive_0.1" ~ "Partially Recessive [0,0.1,2]",
        architecture == "partially_recessive_0.2" ~ "Partially Recessive [0,0.2,2]",
        architecture == "partially_recessive_0.3" ~ "Partially Recessive [0,0.3,2]",
        TRUE ~ architecture
      ),
      # Order the facets: Recessive -> Partially Recessive (low to high) -> Additive
      arch_label = factor(arch_label, levels = c(
        "Recessive [0,0,2]",
        "Partially Recessive [0,0.05,2]",
        "Partially Recessive [0,0.1,2]",
        "Partially Recessive [0,0.2,2]",
        "Partially Recessive [0,0.3,2]",
        "Additive [0,1,2]"
      ))
    )

  p <- ggplot(plot_data, aes(x = heritability * 100, y = power, color = model, group = model)) +
    geom_line(linewidth = 1, alpha = 0.85) +
    geom_point(size = 2.5, shape = 19) +
    facet_wrap(~ arch_label, ncol = 3, scales = "free_x") +
    scale_x_continuous(
      name = "Heritability (%)",
      labels = function(x) {
        ifelse(x >= 1, sprintf("%.1f", x), sprintf("%.2f", x))
      },
      expand = expansion(mult = c(0.05, 0.05))
    ) +
    scale_y_continuous(
      name = "Power",
      limits = c(0, 1),
      breaks = c(0, 0.25, 0.5, 0.75, 1),
      labels = c("0%", "25%", "50%", "75%", "100%"),
      expand = expansion(mult = c(0.02, 0.02))
    ) +
    scale_color_manual(
      name = "Encoding",
      values = c("additive" = "#E31A1C", "dominance" = "#1F78B4", "recessive" = "#33A02C"),
      labels = c("additive" = "Additive", "dominance" = "Nonadditive", "recessive" = "Recessive")
    ) +
    geom_hline(yintercept = 0.8, linetype = "dashed", color = "gray40", linewidth = 0.5, alpha = 0.6) +
    labs(
      title = TeX(sprintf("Power to Detect Causal Variants (%s)", alpha_label)),
      subtitle = sprintf("%d causal variants per trait", N_CAUSAL)
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold", size = 14, margin = margin(b = 3)),
      plot.subtitle = element_text(hjust = 0.5, size = 10, color = "gray40", margin = margin(b = 12)),
      legend.position = "bottom",
      legend.title = element_text(face = "bold", size = 11),
      legend.text = element_text(size = 10),
      legend.key.width = unit(1.5, "cm"),
      panel.grid.major = element_line(color = "gray92", linewidth = 0.35),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = 10.5, margin = margin(3, 3, 3, 3)),
      axis.title.x = element_text(face = "bold", size = 11, margin = margin(t = 8)),
      axis.title.y = element_text(face = "bold", size = 11, margin = margin(r = 8)),
      axis.text.x = element_text(size = 9.5, color = "gray20"),
      axis.text.y = element_text(size = 9.5, color = "gray20"),
      axis.ticks = element_line(color = "gray70", linewidth = 0.3),
      axis.ticks.length = unit(2, "pt"),
      plot.margin = margin(8, 8, 8, 8)
    )

  # Save plot
  output_file <- file.path(output_dir, sprintf("power_curves_%s_M_%02d.pdf", filename_suffix, N_CAUSAL))
  ggsave(output_file, p, width = 11, height = 6.5, dpi=300)
  cat(sprintf("Saved: %s\n", output_file))

  return(p)
}

# Generate power curves for different alpha levels (causal variants only)
cat("\n=== Generating Power Curves (Causal Variants Only) ===\n")
p1 <- make_power_curves(power_summary_causal, 0.05, "$\\alpha = 0.05$", "causal_alpha_0.05")
p2 <- make_power_curves(power_summary_causal, 0.01, "$\\alpha = 0.01$", "causal_alpha_0.01")
p3 <- make_power_curves(power_summary_causal, 0.001, "$\\alpha = 0.001$", "causal_alpha_0.001")
p4 <- make_power_curves(power_summary_causal, 5e-8, "$\\alpha = 5 \\times 10^{-8}$", "causal_alpha_5e8")

# Print final plot
print(p3)

# =============================================================================
# SUMMARY: Best performing model for each architecture
# =============================================================================

cat("\n=== Best Performing Model (Highest Power) ===\n")
best_model <- power_summary_causal %>%
  filter(alpha == 0.05) %>%
  group_by(architecture, heritability) %>%
  slice_max(power, n = 1, with_ties = FALSE) %>%
  select(architecture, heritability, best_model = model, power) %>%
  arrange(architecture, heritability)

print(best_model, n = Inf)

# =============================================================================
# HEATMAP (OPTIONAL - commented out)
# =============================================================================
#
# # Create a heatmap showing power by architecture x model
# cat("\n=== Generating Power Heatmap ===\n")
# heatmap_data <- power_summary_causal %>%
#   filter(alpha == 0.05) %>%
#   mutate(
#     arch_label = case_when(
#       architecture == "additive" ~ "Additive",
#       architecture == "recessive" ~ "Recessive",
#       grepl("partially_recessive", architecture) ~ {
#         het <- sub("partially_recessive_", "", architecture)
#         sprintf("Partially Recessive (h = %s)", het)
#       },
#       TRUE ~ architecture
#     ),
#     model_label = case_when(
#       model == "additive" ~ "Additive",
#       model == "dominance" ~ "Dominance",
#       model == "recessive" ~ "Recessive"
#     )
#   )
#
# p_heatmap <- ggplot(heatmap_data, aes(x = heritability * 100, y = model_label, fill = power)) +
#   geom_tile(color = "white", size = 0.5) +
#   geom_text(aes(label = sprintf("%.2f", power)), size = 3, color = "black") +
#   facet_wrap(~ arch_label, ncol = 2) +
#   scale_fill_gradient2(
#     low = "white", mid = "yellow", high = "darkred",
#     midpoint = 0.5,
#     limits = c(0, 1),
#     name = "Power"
#   ) +
#   scale_x_continuous(
#     name = "Heritability (%)",
#     labels = function(x) sprintf("%.2f", x)
#   ) +
#   scale_y_discrete(name = "Test Model") +
#   labs(
#     title = sprintf("Power Heatmap (α = 0.05, M = %d)", N_CAUSAL),
#     caption = sprintf("Replicate: %s", ifelse(is.null(REPLICATE), "All", REPLICATE))
#   ) +
#   theme_bw(base_size = 11) +
#   theme(
#     plot.title = element_text(hjust = 0.5, face = "bold", size = 14),
#     plot.caption = element_text(hjust = 0.5, size = 9, color = "gray40"),
#     legend.position = "right",
#     panel.grid = element_blank(),
#     strip.background = element_rect(fill = "gray90"),
#     strip.text = element_text(face = "bold", size = 10),
#     axis.text.x = element_text(angle = 45, hjust = 1)
#   )
#
# heatmap_file <- file.path(output_dir, sprintf("power_heatmap_alpha_0.05_M_%02d.pdf", N_CAUSAL))
# ggsave(heatmap_file, p_heatmap, width = 11, height = 8)
# cat(sprintf("Saved: %s\n", heatmap_file))
#
# print(p_heatmap)

cat("\n=== Done ===\n")
cat("Output directory:", output_dir, "\n")
