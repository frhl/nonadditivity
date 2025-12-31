#!/usr/bin/env Rscript
# Interactive analysis of SAIGE null simulation results
# Author: Frederik Lassen

options(scipen = 999)


# Required packages
packages <- c('data.table', 'ggplot2', 'ggrastr', 'latex2exp', 'dplyr', 'stringr')

cat("Installing/loading required packages...\n")
for (p in packages) {
    if (!require(p, character.only = TRUE, quietly = TRUE)) {
        cat(paste0("Installing package: ", p, "\n"))
        install.packages(p, repos = "http://cran.us.r-project.org")
        library(p, character.only = TRUE)
    }
}

# Load bravastring (optional - provides calc_inflation and get_expected_p functions)
if (!require("bravastring", quietly = TRUE)) {
    cat("bravastring not found. Installing from GitHub...\n")
    if (!require("devtools", quietly = TRUE)) {
        install.packages("devtools", repos = "http://cran.us.r-project.org")
        library(devtools)
    }
    devtools::install_github("frhl/bravastring")
    library(bravastring)
}


# =============================================================================
# SECTION 2: Define Helper Functions
# =============================================================================

# Function to calculate genomic inflation factor (lambda)
calc_inflation <- function(p_values) {
    # Remove NAs and values that are 0 or 1
    p_values <- p_values[!is.na(p_values) & p_values > 0 & p_values < 1]

    if (length(p_values) == 0) return(NA)

    # Calculate lambda using median chi-square method
    chisq <- qchisq(1 - p_values, 1)
    lambda <- median(chisq, na.rm = TRUE) / qchisq(0.5, 1)
    return(lambda)
}

# Function to calculate expected p-values
get_expected_p <- function(p_values) {
    # Remove NAs
    p_values <- p_values[!is.na(p_values)]
    n <- length(p_values)

    if (n == 0) return(rep(NA, length(p_values)))

    # Rank the p-values
    ranks <- rank(p_values, ties.method = "first")

    # Calculate expected p-values
    expected_p <- ranks / (n + 1)

    return(expected_p)
}

# Function to calculate Type I error rate and 95% CI
calc_type1_error <- function(p_values, alpha) {
    # Remove NAs
    p_values <- p_values[!is.na(p_values)]
    n <- length(p_values)

    if (n == 0) {
        return(list(rate = NA, ci_lower = NA, ci_upper = NA, n = 0))
    }

    # Calculate empirical rate
    empirical_rate <- mean(p_values < alpha)

    # Calculate 95% Confidence Interval (Standard Error of proportion)
    se <- sqrt((empirical_rate * (1 - empirical_rate)) / n)
    ci_lower <- empirical_rate - 1.96 * se
    ci_upper <- empirical_rate + 1.96 * se

    # Clamp CI to [0, 1]
    ci_lower <- max(0, ci_lower)
    ci_upper <- min(1, ci_upper)

    return(list(
        rate = empirical_rate,
        ci_lower = ci_lower,
        ci_upper = ci_upper,
        n = n
    ))
}

cat("Helper functions defined.\n\n")

# =============================================================================
# SECTION 3: Set Parameters - MODIFY THESE AS NEEDED
# =============================================================================

# Local data directory
data_dir <- "/Users/flassen/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/data/null_simulation_saige"

# Ancestries to analyze (modify as needed)
ancestries <- c("eur", "sas", "afr", "eas")

# Annotation pattern to filter files
annotation <- "pLoF_damaging_missense"

# Significance thresholds for Type I error assessment
alphas <- c(0.05, 0.01, 0.001, 0.0001, 0.00001)

# Output directory for plots (will be created if it doesn't exist)
output_dir <- file.path("~/Desktop", "null_simulation_plots")
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

cat("Parameters set:\n")
cat(paste0("  Data directory: ", data_dir, "\n"))
cat(paste0("  Ancestries: ", paste(ancestries, collapse = ", "), "\n"))
cat(paste0("  Annotation: ", annotation, "\n"))
cat(paste0("  Output directory: ", output_dir, "\n\n"))

# =============================================================================
# SECTION 4: Load Data for Selected Ancestry
# =============================================================================

# Select which ancestry to analyze (change this to analyze different ancestries)
selected_ancestry <- "eur"  # Change to "sas", "afr", "eas" as needed

cat(paste0("Loading data for ancestry: ", selected_ancestry, "\n"))

# Find all files matching the ancestry and annotation
pattern <- paste0("saige_step2\\.", selected_ancestry, "\\..*", annotation, "\\.combined\\.tsv\\.gz$")
files <- list.files(data_dir, pattern = pattern, full.names = TRUE)

cat(paste0("Found ", length(files), " files:\n"))
for (f in files) {
    cat(paste0("  - ", basename(f), "\n"))
}

if (length(files) == 0) {
    stop(paste0("No files found matching pattern: ", pattern))
}

# Load all files and combine
cat("\nLoading files...\n")
data_list <- lapply(files, function(f) {
    cat(paste0("  Loading: ", basename(f), "\n"))
    d <- fread(f)

    # Extract metadata from filename
    d[, heritability := stringr::str_extract(basename(f), "(?<=p_)\\d*\\.?\\d+")]
    d[, replicate := stringr::str_extract(basename(f), "(?<=continuous_)\\d+")]
    d[, annotation := stringr::str_extract(basename(f), "(?<=continuous_\\d\\.).+?(?=\\.combined)")]
    d[, ancestry := selected_ancestry]

    # Filter to variants with at least one homozygous alternate (i.e., present in nonadditive model)
    d <- d[!is.na(DOM.P)]

    cat(paste0("    H2: ", unique(d$heritability),
               ", Rep: ", unique(d$replicate),
               ", Variants: ", nrow(d), "\n"))
    return(d)
})

# Combine all data
dt <- rbindlist(data_list, fill = TRUE)
cat(paste0("\nTotal variants loaded: ", nrow(dt), "\n"))

# Apply N_HOMALT >= 5 filter (standard for recessive/nonadditive tests)
cat("\nFiltering for N_HOMALT >= 5...\n")
cat(paste0("  Before: ", nrow(dt), " variants\n"))
dt <- dt[N_HOMALT >= 5]
cat(paste0("  After: ", nrow(dt), " variants\n\n"))

cat(paste0("Unique heritabilities: ", paste(unique(dt$heritability), collapse = ", "), "\n"))
cat(paste0("Unique replicates: ", paste(unique(dt$replicate), collapse = ", "), "\n\n"))

# =============================================================================
# SECTION 5: Calculate Expected P-values and Create Bins
# =============================================================================

cat("Calculating expected p-values for all models...\n")
dt[, ADD.PEXPT := get_expected_p(ADD.P), by = .(heritability, replicate)]
dt[, REC.PEXPT := get_expected_p(REC.P), by = .(heritability, replicate)]
dt[, DOM.PEXPT := get_expected_p(DOM.P), by = .(heritability, replicate)]

cat("Creating allele count bins...\n")
# AC_Allele2 bins for additive model
dt[, AC2_BIN := cut(
    AC_Allele2,
    breaks = c(0, 5, 10, 25, Inf),
    labels = c("<5", "[5, 10)", "[10, 25)", ">25"),
    include.lowest = TRUE
)]

# N_HOMALT bins for recessive and nonadditive models
dt[, N_HOMALT_BIN := cut(
    N_HOMALT,
    breaks = c(5, 10, 25, Inf),
    labels = c("[5, 10)", "[10, 25)", ">25"),
    include.lowest = TRUE
)]

cat("Data preparation complete!\n\n")

# Quick summary
cat("=== Data Summary ===\n")
cat(paste0("Ancestry: ", selected_ancestry, "\n"))
cat(paste0("Total variants: ", nrow(dt), "\n"))
cat(paste0("Heritabilities: ", paste(sort(as.numeric(unique(dt$heritability))), collapse = ", "), "\n"))
cat(paste0("Replicates: ", paste(sort(as.numeric(unique(dt$replicate))), collapse = ", "), "\n\n"))

# =============================================================================
# SECTION 6: Calculate Type I Error Rates
# =============================================================================

cat("=============================================================================\n")
cat("CALCULATING TYPE I ERROR RATES\n")
cat("=============================================================================\n\n")

type1_results <- list()

for (alpha in alphas) {
    cat(paste0("Calculating Type I error at alpha = ", alpha, "...\n"))

    # Additive model
    add_res <- dt[, calc_type1_error(ADD.P, alpha), by = .(heritability, replicate)]
    add_res[, model := "Additive"]
    add_res[, alpha := alpha]

    # Recessive model
    rec_res <- dt[, calc_type1_error(REC.P, alpha), by = .(heritability, replicate)]
    rec_res[, model := "Recessive"]
    rec_res[, alpha := alpha]

    # nonadditive model
    dom_res <- dt[, calc_type1_error(DOM.P, alpha), by = .(heritability, replicate)]
    dom_res[, model := "nonadditive"]
    dom_res[, alpha := alpha]

    # Combine
    type1_results[[length(type1_results) + 1]] <- rbind(add_res, rec_res, dom_res)
}

type1_summary <- rbindlist(type1_results)
type1_summary[, ancestry := selected_ancestry]

# Calculate average across replicates
type1_avg <- type1_summary[, .(
    mean_rate = mean(rate, na.rm = TRUE),
    mean_ci_lower = mean(ci_lower, na.rm = TRUE),
    mean_ci_upper = mean(ci_upper, na.rm = TRUE),
    n_replicates = .N
), by = .(ancestry, model, heritability, alpha)]

cat("\n=== Type I Error Summary (averaged across replicates) ===\n")
print(type1_avg[order(alpha, model, heritability)])

# Check calibration at alpha = 0.05
cat("\n=== Calibration Check at alpha = 0.05 ===\n")
calib_05 <- type1_avg[alpha == 0.05]
calib_05[, deviation := abs(mean_rate - alpha)]
calib_05[, well_calibrated := mean_ci_lower <= alpha & mean_ci_upper >= alpha]
print(calib_05[order(model, heritability)])

# =============================================================================
# SECTION 7: Calculate Genomic Inflation Factors (Lambda)
# =============================================================================

cat("\n=============================================================================\n")
cat("CALCULATING GENOMIC INFLATION FACTORS (LAMBDA)\n")
cat("=============================================================================\n\n")

lambda_add <- dt[, .(lambda = calc_inflation(ADD.P)), by = .(heritability, replicate)]
lambda_add[, model := "Additive"]

lambda_rec <- dt[, .(lambda = calc_inflation(REC.P)), by = .(heritability, replicate)]
lambda_rec[, model := "Recessive"]

lambda_dom <- dt[, .(lambda = calc_inflation(DOM.P)), by = .(heritability, replicate)]
lambda_dom[, model := "nonadditive"]

lambda_summary <- rbind(lambda_add, lambda_rec, lambda_dom)
lambda_summary[, ancestry := selected_ancestry]

# Average across replicates
lambda_avg <- lambda_summary[, .(
    mean_lambda = mean(lambda, na.rm = TRUE),
    sd_lambda = sd(lambda, na.rm = TRUE)
), by = .(ancestry, model, heritability)]

cat("=== Lambda Summary ===\n")
print(lambda_avg[order(model, heritability)])

# =============================================================================
# SECTION 8: Create QQ Plots
# =============================================================================

cat("\n=============================================================================\n")
cat("CREATING QQ PLOTS\n")
cat("=============================================================================\n\n")

# Define colors
colors <- list(
    red = "#E35278",
    orange = "#E2A98C",
    green = "#9EC0A6",
    blue = "#009894"
)

# Calculate plot dimensions
n_replicates <- length(unique(dt$replicate))
n_facets <- n_heritabilities * n_replicates
(width <- min(24, max(16, ceiling(sqrt(n_facets)) * 4)))
(height <- min(24, max(12, ceiling(n_facets / ceiling(sqrt(n_facets))) * 3)))

# width
width <- 12
height <- 10


# --- ADDITIVE MODEL ---
cat("Creating ADDITIVE model QQ plot...\n")

# Add lambda labels
lambda_add[, lambda_label := paste0("λ = ", round(lambda, 3))]

# Calculate 95% confidence ribbon
ribbon_p <- 0.95
dt_ribbon_add <- dt[, {
    n <- round(.N * 2)
    seq.p.value <- (1:n) / n
    list(
        ADD.PEXPT = -log10(seq.p.value),
        clower = -log10(qbeta(p = (1 + ribbon_p) / 2, shape1 = 1:n, shape2 = n:1)),
        cupper = -log10(qbeta(p = (1 - ribbon_p) / 2, shape1 = 1:n, shape2 = n:1))
    )
}, by = .(heritability, replicate)]

p_add <- ggplot(dt, aes(x = -log10(ADD.PEXPT), y = -log10(ADD.P))) +
    geom_ribbon(
        data = dt_ribbon_add,
        aes(x = ADD.PEXPT, ymin = clower, ymax = cupper, y = NULL),
        fill = "grey80", alpha = 0.5, inherit.aes = FALSE
    ) +
    geom_abline(intercept = 0, slope = 1, color = 'black', linetype = "dashed") +
    geom_point_rast(size = 2, alpha = 0.9) +
    scale_color_manual(
        name = "AC_Allele2",
        values = c("<5" = colors$red, "[5, 10)" = colors$orange,
                  "[10, 25)" = colors$green, ">25" = colors$blue)
    ) +
    geom_text(
        data = lambda_add,
        aes(x = -Inf, y = Inf, label = lambda_label),
        hjust = -0.2, vjust = 1.5, inherit.aes = FALSE, size = 4
    ) +
    facet_wrap(~ heritability + replicate, labeller = label_both) +
    ggtitle(paste0("Additive Model - ", toupper(selected_ancestry))) +
    ylab(TeX("$-\\log_{10}(P_{observed})$")) +
    xlab(TeX("$-\\log_{10}(P_{expected})$")) +
    theme_classic() +
    theme(
        legend.position = "right",
        strip.background = element_rect(fill = "grey90", color = "grey90"),
        strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(margin = ggplot2::margin(t = 15)),
        axis.title.y = element_text(margin = ggplot2::margin(r = 15)),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
    )

# Save plot
output_add <- file.path(output_dir, paste0("qq_plot_", selected_ancestry, "_additive.pdf"))
ggsave(filename = output_add, plot = p_add, width = width, height = height, units = "in", dpi = 300)
cat(paste0("  Saved to: ", output_add, "\n"))

# Display plot (for RStudio)
print(p_add)

# --- RECESSIVE MODEL ---
cat("\nCreating RECESSIVE model QQ plot...\n")

lambda_rec[, lambda_label := paste0("λ = ", round(lambda, 3))]

dt_ribbon_rec <- dt[, {
    n <- round(.N * 2)
    seq.p.value <- (1:n) / n
    list(
        REC.PEXPT = -log10(seq.p.value),
        clower = -log10(qbeta(p = (1 + ribbon_p) / 2, shape1 = 1:n, shape2 = n:1)),
        cupper = -log10(qbeta(p = (1 - ribbon_p) / 2, shape1 = 1:n, shape2 = n:1))
    )
}, by = .(heritability, replicate)]

p_rec <- ggplot(dt, aes(x = -log10(REC.PEXPT), y = -log10(REC.P))) +
    geom_ribbon(
        data = dt_ribbon_rec,
        aes(x = REC.PEXPT, ymin = clower, ymax = cupper, y = NULL),
        fill = "grey80", alpha = 0.5, inherit.aes = FALSE
    ) +
    geom_abline(intercept = 0, slope = 1, color = 'black', linetype = "dashed") +
    geom_point_rast(size = 2, alpha = 0.9) +
    scale_color_manual(
        name = "N_HOMALT",
        values = c("[5, 10)" = colors$orange,
                  "[10, 25)" = colors$green, ">25" = colors$blue)
    ) +
    geom_text(
        data = lambda_rec,
        aes(x = -Inf, y = Inf, label = lambda_label),
        hjust = -0.2, vjust = 1.5, inherit.aes = FALSE, size = 4
    ) +
    facet_wrap(~ heritability + replicate, labeller = label_both) +
    ggtitle(paste0("Recessive Model - ", toupper(selected_ancestry))) +
    ylab(TeX("$-\\log_{10}(P_{observed})$")) +
    xlab(TeX("$-\\log_{10}(P_{expected})$")) +
    theme_classic() +
    theme(
        legend.position = "right",
        strip.background = element_rect(fill = "grey90", color = "grey90"),
        strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(margin = ggplot2::margin(t = 15)),
        axis.title.y = element_text(margin = ggplot2::margin(r = 15)),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
    )

output_rec <- file.path(output_dir, paste0("qq_plot_", selected_ancestry, "_recessive.pdf"))
ggsave(filename = output_rec, plot = p_rec, width = width, height = height, units = "in", dpi = 300)
cat(paste0("  Saved to: ", output_rec, "\n"))

print(p_rec)

# --- nonadditive MODEL ---
cat("\nCreating nonadditive model QQ plot...\n")

lambda_dom[, lambda_label := paste0("λ = ", round(lambda, 3))]

dt_ribbon_dom <- dt[, {
    n <- round(.N * 2)
    seq.p.value <- (1:n) / n
    list(
        DOM.PEXPT = -log10(seq.p.value),
        clower = -log10(qbeta(p = (1 + ribbon_p) / 2, shape1 = 1:n, shape2 = n:1)),
        cupper = -log10(qbeta(p = (1 - ribbon_p) / 2, shape1 = 1:n, shape2 = n:1))
    )
}, by = .(heritability, replicate)]

p_dom <- ggplot(dt, aes(x = -log10(DOM.PEXPT), y = -log10(DOM.P))) +
    geom_ribbon(
        data = dt_ribbon_dom,
        aes(x = DOM.PEXPT, ymin = clower, ymax = cupper, y = NULL),
        fill = "grey80", alpha = 0.5, inherit.aes = FALSE
    ) +
    geom_abline(intercept = 0, slope = 1, color = 'black', linetype = "dashed") +
    geom_point_rast(size = 2, alpha = 0.9) +
    scale_color_manual(
        name = "N_HOMALT",
        values = c("[5, 10)" = colors$orange,
                  "[10, 25)" = colors$green, ">25" = colors$blue)
    ) +
    geom_text(
        data = lambda_dom,
        aes(x = -Inf, y = Inf, label = lambda_label),
        hjust = -0.2, vjust = 1.5, inherit.aes = FALSE, size = 4
    ) +
    facet_wrap(~ heritability + replicate, labeller = label_both) +
    ggtitle(paste0("Nonadditive ", toupper(selected_ancestry))) +
    ylab(TeX("$-\\log_{10}(P_{observed})$")) +
    xlab(TeX("$-\\log_{10}(P_{expected})$")) +
    theme_classic() +
    theme(
        legend.position = "right",
        strip.background = element_rect(fill = "grey90", color = "grey90"),
        strip.text = element_text(size = 12, face = "bold"),
        axis.text = element_text(size = 10),
        axis.title = element_text(size = 16, face = "bold"),
        axis.title.x = element_text(margin = ggplot2::margin(t = 15)),
        axis.title.y = element_text(margin = ggplot2::margin(r = 15)),
        plot.title = element_text(hjust = 0.5, size = 18, face = "bold")
    )

output_dom <- file.path(output_dir, paste0("qq_plot_", selected_ancestry, "_nonadditive.pdf"))
ggsave(filename = output_dom, plot = p_dom, width = width, height = height, units = "in", dpi = 300)
cat(paste0("  Saved to: ", output_dom, "\n"))

print(p_dom)

# =============================================================================
# SECTION 9: Save Summary Tables
# =============================================================================

cat("\n=============================================================================\n")
cat("SAVING SUMMARY TABLES\n")
cat("=============================================================================\n\n")

# Save type 1 error summary
type1_file <- file.path(output_dir, paste0("type1_error_", selected_ancestry, ".tsv"))
fwrite(type1_avg, type1_file, sep = "\t")
cat(paste0("Type I error summary saved to: ", type1_file, "\n"))

# Save lambda summary
lambda_file <- file.path(output_dir, paste0("lambda_", selected_ancestry, ".tsv"))
fwrite(lambda_avg, lambda_file, sep = "\t")
cat(paste0("Lambda summary saved to: ", lambda_file, "\n"))

cat("\n=============================================================================\n")
cat("ANALYSIS COMPLETE!\n")
cat("=============================================================================\n\n")

cat("To analyze a different ancestry, change 'selected_ancestry' in SECTION 4\n")
cat("and re-run sections 4-9.\n\n")

cat("Output files location: ", output_dir, "\n\n")
