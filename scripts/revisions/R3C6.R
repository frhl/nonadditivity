#!/usr/bin/env Rscript
# Interactive QQ plot creation script
# You can step through this script line by line in RStudio or R console

options(scipen = 999)

# ============================================================
# SETUP: Install and load packages
# ============================================================

packages <- c('data.table', 'ggplot2', 'ggrastr', 'latex2exp', 'stringr')

for (p in packages) {
    if (!require(p, character.only = TRUE, quietly = TRUE)) {
        cat(paste0("Installing package: ", p, "\n"))
        install.packages(p, repos = "http://cran.us.r-project.org")
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

# ============================================================
# PARAMETERS: Modify these as needed
# ============================================================

# Input directory with downloaded data
combined_dir <- "/Users/flassen/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/data/rint_comparison"

# Output directory for plots
out_dir <- "/Users/flassen/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/data/rint_comparison/qq_plots"
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Filtering parameters
ancestry <- "eur"
annotation <- "pLoF_damaging_missense"  # or "synonymous"
min_biallelic <- 5  # Minimum N_HOMALT

# Output filename base
output_base <- file.path(out_dir, paste0("qq_rint_comparison_", ancestry, "_", annotation, "_minbi", min_biallelic))

# ============================================================
# STEP 1: Find and list input files
# ============================================================

cat("=============================================================\n")
cat("Creating QQ plots for RINT vs no-RINT comparison\n")
cat("=============================================================\n\n")

# List combined files for RINT
rint_pattern <- paste0(".*", ancestry, ".*", annotation, "\\.rint\\.combined\\.tsv\\.gz$")
rint_files <- list.files(combined_dir, pattern = rint_pattern, full.names = TRUE)

cat(paste0("Found ", length(rint_files), " RINT combined files:\n"))
for (f in rint_files) {
    cat(paste0("  - ", basename(f), "\n"))
}
cat("\n")

# List combined files for standard
standard_pattern <- paste0(".*", ancestry, ".*", annotation, "\\.standard\\.combined\\.tsv\\.gz$")
standard_files <- list.files(combined_dir, pattern = standard_pattern, full.names = TRUE)

cat(paste0("Found ", length(standard_files), " standard combined files:\n"))
for (f in standard_files) {
    cat(paste0("  - ", basename(f), "\n"))
}
cat("\n")

if (length(rint_files) == 0 || length(standard_files) == 0) {
    stop("No combined files found for comparison.", call.=FALSE)
}

# ============================================================
# STEP 2: Load data
# ============================================================

# Load RINT files
cat("Loading RINT files...\n")
rint_list <- lapply(rint_files, function(f) {
    cat(paste0("  Loading: ", basename(f), "\n"))
    d <- fread(f)
    d[, phenotype := sub("\\..*", "", basename(f))]
    d[, version := "RINT"]
    return(d)
})
d_rint <- rbindlist(rint_list, fill = TRUE)

# Load standard files
cat("Loading standard files...\n")
standard_list <- lapply(standard_files, function(f) {
    cat(paste0("  Loading: ", basename(f), "\n"))
    d <- fread(f)
    d[, phenotype := sub("\\..*", "", basename(f))]
    d[, version := "standard"]
    return(d)
})
d_standard <- rbindlist(standard_list, fill = TRUE)

# Combine both versions
d <- rbindlist(list(d_rint, d_standard), fill = TRUE)

cat(paste0("\nTotal variants: ", nrow(d), "\n"))
cat(paste0("Unique phenotypes: ", paste(unique(d$phenotype), collapse=", "), "\n\n"))

# ============================================================
# STEP 3: Convert p-values to numeric and apply filters
# ============================================================

cat("Converting p-value columns to numeric...\n")
d[, ADD.P := as.numeric(ADD.P)]
d[, REC.P := as.numeric(REC.P)]
d[, DOM.P := as.numeric(DOM.P)]

cat("Applying filters...\n")
cat(paste0("  Initial variants: ", nrow(d), "\n"))

# Filter by minimum biallelic carriers (N_HOMALT >= min_biallelic)
cat(paste0("  Filtering variants with N_HOMALT >= ", min_biallelic, "\n"))
d <- d[N_HOMALT >= min_biallelic]
cat(paste0("  After filtering: ", nrow(d), " variants\n"))

# ============================================================
# STEP 4: Calculate expected p-values
# ============================================================

cat("\nCalculating expected p-values...\n")
d[, ADD.PEXPT := get_expected_p(ADD.P), by = .(phenotype, version)]
d[, REC.PEXPT := get_expected_p(REC.P), by = .(phenotype, version)]
d[, DOM.PEXPT := get_expected_p(DOM.P), by = .(phenotype, version)]

# ============================================================
# STEP 5: Ensure stratification bins exist
# ============================================================

if (!"AC2_interval" %in% colnames(d)) {
    cat("Creating AC2 interval bins...\n")
    d[, AC2_interval := cut(
        AC_Allele2,
        breaks = c(0, 5, 10, 25, Inf),
        labels = c("<5", "[5, 10)", "[10, 25)", ">25"),
        include.lowest = TRUE
    )]
}

if (!"N_HOMALT_BIN" %in% colnames(d)) {
    cat("Creating N_HOMALT bins...\n")
    d[, N_HOMALT_BIN := cut(
        N_HOMALT,
        breaks = c(-1, 0, 1, 5, 10, 25, Inf),
        labels = c("0", "1", "2-4", "5-9", "10-24", ">=25"),
        include.lowest = TRUE
    )]
}

# ============================================================
# STEP 6: Define plot settings
# ============================================================

# Color scheme
colors <- list(
    red = c("#E35278"),
    orange = c("#E2A98C"),
    green = c("#9EC0A6"),
    blue = c("#009894")
)

# Calculate plot dimensions
n_phenotypes <- length(unique(d$phenotype))
n_versions <- 2  # RINT and standard
n_facets <- n_phenotypes * n_versions
width <- min(24, max(16, ceiling(sqrt(n_facets)) * 4))
height <- min(24, max(12, ceiling(n_facets / ceiling(sqrt(n_facets))) * 3))

# Ribbon confidence level
ribbon_p <- 0.95

# ============================================================
# STEP 7: Create ADDITIVE model QQ plot
# ============================================================

cat("\n=== Creating ADDITIVE model QQ plot ===\n")

# Calculate lambda (genomic inflation factor)
lambda_add <- d[, .(lambda = calc_inflation(ADD.P)), by = .(phenotype, version)]
lambda_add[, lambda_label := paste0("λ = ", sprintf("%.3f", lambda))]
print(lambda_add)

# Calculate confidence ribbon
dt_ribbon_add <- d[, {
    n <- round(.N * 2)
    seq.p.value <- (1:n) / n
    list(
        ADD.PEXPT = -log10(seq.p.value),
        clower = -log10(qbeta(p = (1 + ribbon_p) / 2, shape1 = 1:n, shape2 = n:1)),
        cupper = -log10(qbeta(p = (1 - ribbon_p) / 2, shape1 = 1:n, shape2 = n:1))
    )
}, by = .(phenotype, version)]

# Create plot
p_add <- ggplot(d, aes(x = -log10(ADD.PEXPT), y = -log10(ADD.P), color = AC2_interval)) +
    geom_ribbon(
        data = dt_ribbon_add,
        aes(x = ADD.PEXPT, ymin = clower, ymax = cupper, y = NULL),
        fill = "grey80", alpha = 0.5, inherit.aes = FALSE
    ) +
    geom_abline(intercept = 0, slope = 1, color = 'black', linetype = "dashed") +
    geom_point_rast(size = 2, alpha = 0.7) +
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
    facet_wrap(~ phenotype + version, labeller = label_both) +
    ggtitle("Additive Model") +
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
output_add <- paste0(output_base, "_additive.pdf")
cat(paste0("Saving additive plot to: ", output_add, "\n"))
ggsave(filename = output_add, plot = p_add, width = width, height = height, units = "in", dpi = 300)

# ============================================================
# STEP 8: Create RECESSIVE model QQ plot
# ============================================================

cat("\n=== Creating RECESSIVE model QQ plot ===\n")

# Calculate lambda
lambda_rec <- d[, .(lambda = calc_inflation(REC.P)), by = .(phenotype, version)]
lambda_rec[, lambda_label := paste0("λ = ", sprintf("%.3f", lambda))]
print(lambda_rec)

# Calculate confidence ribbon
dt_ribbon_rec <- d[, {
    n <- round(.N * 2)
    seq.p.value <- (1:n) / n
    list(
        REC.PEXPT = -log10(seq.p.value),
        clower = -log10(qbeta(p = (1 + ribbon_p) / 2, shape1 = 1:n, shape2 = n:1)),
        cupper = -log10(qbeta(p = (1 - ribbon_p) / 2, shape1 = 1:n, shape2 = n:1))
    )
}, by = .(phenotype, version)]

# Create plot
p_rec <- ggplot(d, aes(x = -log10(REC.PEXPT), y = -log10(REC.P), color = N_HOMALT_BIN)) +
    geom_ribbon(
        data = dt_ribbon_rec,
        aes(x = REC.PEXPT, ymin = clower, ymax = cupper, y = NULL),
        fill = "grey80", alpha = 0.5, inherit.aes = FALSE
    ) +
    geom_abline(intercept = 0, slope = 1, color = 'black', linetype = "dashed") +
    geom_point_rast(size = 2, alpha = 0.7) +
    scale_color_manual(
        name = "N_HOMALT",
        values = c("0" = "grey50", "1" = colors$red, "2-4" = colors$orange,
                  "5-9" = colors$green, "10-24" = colors$blue, ">=25" = "darkblue")
    ) +
    geom_text(
        data = lambda_rec,
        aes(x = -Inf, y = Inf, label = lambda_label),
        hjust = -0.2, vjust = 1.5, inherit.aes = FALSE, size = 4
    ) +
    facet_wrap(~ phenotype + version, labeller = label_both) +
    ggtitle("Recessive Model") +
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
output_rec <- paste0(output_base, "_recessive.pdf")
cat(paste0("Saving recessive plot to: ", output_rec, "\n"))
ggsave(filename = output_rec, plot = p_rec, width = width, height = height, units = "in", dpi = 300)

# ============================================================
# STEP 9: Create DOMINANCE model QQ plot
# ============================================================

cat("\n=== Creating DOMINANCE model QQ plot ===\n")

# Calculate lambda
lambda_dom <- d[, .(lambda = calc_inflation(DOM.P)), by = .(phenotype, version)]
lambda_dom[, lambda_label := paste0("λ = ", sprintf("%.3f", lambda))]
print(lambda_dom)

# Calculate confidence ribbon
dt_ribbon_dom <- d[, {
    n <- round(.N * 2)
    seq.p.value <- (1:n) / n
    list(
        DOM.PEXPT = -log10(seq.p.value),
        clower = -log10(qbeta(p = (1 + ribbon_p) / 2, shape1 = 1:n, shape2 = n:1)),
        cupper = -log10(qbeta(p = (1 - ribbon_p) / 2, shape1 = 1:n, shape2 = n:1))
    )
}, by = .(phenotype, version)]

# Create plot
p_dom <- ggplot(d, aes(x = -log10(DOM.PEXPT), y = -log10(DOM.P), color = N_HOMALT_BIN)) +
    geom_ribbon(
        data = dt_ribbon_dom,
        aes(x = DOM.PEXPT, ymin = clower, ymax = cupper, y = NULL),
        fill = "grey80", alpha = 0.5, inherit.aes = FALSE
    ) +
    geom_abline(intercept = 0, slope = 1, color = 'black', linetype = "dashed") +
    geom_point_rast(size = 2, alpha = 0.7) +
    scale_color_manual(
        name = "N_HOMALT",
        values = c("0" = "grey50", "1" = colors$red, "2-4" = colors$orange,
                  "5-9" = colors$green, "10-24" = colors$blue, ">=25" = "darkblue")
    ) +
    geom_text(
        data = lambda_dom,
        aes(x = -Inf, y = Inf, label = lambda_label),
        hjust = -0.2, vjust = 1.5, inherit.aes = FALSE, size = 4
    ) +
    facet_wrap(~ phenotype + version, labeller = label_both) +
    ggtitle("Dominance Model") +
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
output_dom <- paste0(output_base, "_dominance.pdf")
cat(paste0("Saving dominance plot to: ", output_dom, "\n"))
ggsave(filename = output_dom, plot = p_dom, width = width, height = height, units = "in", dpi = 300)

# ============================================================
# DONE!
# ============================================================

cat("\n=============================================================\n")
cat("Done!\n")
cat("=============================================================\n")
cat("Output files:\n")
cat(paste0("  - ", output_add, "\n"))
cat(paste0("  - ", output_rec, "\n"))
cat(paste0("  - ", output_dom, "\n"))
