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

main <- function(opt)
{
    # Check required inputs
    if (is.null(opt$combined_dir) || is.null(opt$ancestry) || is.null(opt$annotation) || is.null(opt$output)) {
        print_help(opt_parser)
        stop("All parameters (combined_dir, ancestry, annotation, output) must be supplied.", call.=FALSE)
    }

    cat("=============================================================\n")
    cat("Creating QQ plots from combined REGENIE results\n")
    cat("=============================================================\n\n")

    # List all combined files matching ancestry and annotation
    pattern <- paste0(".*", opt$ancestry, ".*", opt$annotation, "\\.combined\\.regenie\\.tsv\\.gz$")
    files <- list.files(opt$combined_dir, pattern = pattern, full.names = TRUE)

    cat(paste0("Found ", length(files), " combined files:\n"))
    for (f in files) {
        cat(paste0("  - ", basename(f), "\n"))
    }
    cat("\n")

    if (length(files) == 0) {
        stop("No combined files found matching pattern: ", pattern, call.=FALSE)
    }

    # Load and process all files
    cat("Loading and processing files...\n")
    l <- lapply(files, function(f){
        cat(paste0("  Loading: ", basename(f), "\n"))
        d <- fread(f)
        d[, heritability := stringr::str_extract(basename(f), "(?<=p_)\\d*\\.?\\d+")]
        d[, id := stringr::str_extract(basename(f), "(?<=continuous_)\\d+")]
        d[, annotation := stringr::str_extract(f, "(?<=continuous_\\d\\.).+?(?=\\.combined)")]

        # Restrict to variants with at least one biallelic variant
        # (i.e., variants present in dominance model)
        d <- d[!is.na(d$DOM.P),]

        cat(paste0("    Heritability: ", unique(d$heritability),
                   ", ID: ", unique(d$id),
                   ", Variants: ", nrow(d), "\n"))
        return(d)
    })

    # Combine all data
    d <- rbindlist(l)
    cat(paste0("\nTotal variants across all files: ", nrow(d), "\n"))

    # Filter for N_HOMALT >= 5
    cat(paste0("Filtering for N_HOMALT >= 5...\n"))
    d <- d[N_HOMALT >= 5,]
    cat(paste0("Variants after filtering: ", nrow(d), "\n"))

    cat(paste0("Unique heritabilities: ", paste(unique(d$heritability), collapse=", "), "\n"))
    cat(paste0("Unique IDs: ", paste(unique(d$id), collapse=", "), "\n\n"))

    # Calculate expected p-values for all models
    cat("Calculating expected p-values...\n")
    d[, ADD.PEXPT := get_expected_p(ADD.P), by = .(heritability, id)]
    d[, REC.PEXPT := get_expected_p(REC.P), by = .(heritability, id)]
    d[, DOM.PEXPT := get_expected_p(DOM.P), by = .(heritability, id)]

    # Create stratification bins
    cat("Creating stratification bins...\n")
    # AC2 bins for additive model
    d[, AC2_BIN := cut(
        AC2,
        breaks = c(0, 5, 10, 25, Inf),
        labels = c("<5", "[5, 10)", "[10, 25)", ">25"),
        include.lowest = TRUE
    )]
    # N_HOMALT bins for recessive and dominance models
    d[, N_HOMALT_BIN := cut(
        N_HOMALT,
        breaks = c(0, 5, 10, 25, Inf),
        labels = c("<5", "[5, 10)", "[10, 25)", ">25"),
        include.lowest = TRUE
    )]

    # Filter out variants in the "<5" N_HOMALT_BIN category
    cat(paste0("Filtering out N_HOMALT_BIN '<5' variants...\n"))
    cat(paste0("Variants before filtering: ", nrow(d), "\n"))
    d <- d[N_HOMALT_BIN != "<5",]
    cat(paste0("Variants after filtering: ", nrow(d), "\n"))

    # Define color scale
    colors <- list(red=c("#E35278"), orange=c("#E2A98C"), green=c("#9EC0A6"), blue=c("#009894"))

    # Create base output filename (without extension)
    output_base <- sub("\\.pdf$", "", opt$output)

    # Calculate plot dimensions
    n_heritabilities <- length(unique(d$heritability))
    n_ids <- length(unique(d$id))
    n_facets <- n_heritabilities * n_ids
    width <- min(24, max(16, ceiling(sqrt(n_facets)) * 4))
    height <- min(24, max(12, ceiling(n_facets / ceiling(sqrt(n_facets))) * 3))

    # ========== ADDITIVE MODEL PLOT ==========
    cat("\n=== Creating ADDITIVE model QQ plot ===\n")
    lambda_add <- d[, .(lambda = calc_inflation(ADD.P)), by = .(heritability, id)]
    lambda_add[, lambda_label := paste0("λ = ", round(lambda, 3))]
    print(lambda_add)

    ribbon_p <- 0.95
    dt_ribbon_add <- d[, {
        n <- round(.N * 2)
        seq.p.value <- (1:n) / n
        list(
            ADD.PEXPT = -log10(seq.p.value),
            clower = -log10(qbeta(p = (1 + ribbon_p) / 2, shape1 = 1:n, shape2 = n:1)),
            cupper = -log10(qbeta(p = (1 - ribbon_p) / 2, shape1 = 1:n, shape2 = n:1))
        )
    }, by = .(heritability, id)]

    p_add <- ggplot(d, aes(x = -log10(ADD.PEXPT), y = -log10(ADD.P), color = AC2_BIN)) +
        geom_ribbon(
            data = dt_ribbon_add,
            aes(x = ADD.PEXPT, ymin = clower, ymax = cupper, y = NULL),
            fill = "grey80", alpha = 0.5, inherit.aes = FALSE
        ) +
        geom_abline(intercept = 0, slope = 1, color = 'black', linetype = "dashed") +
        geom_point_rast(size = 2, alpha = 0.7) +
        scale_color_manual(
            name = "AC2",
            values = c("<5" = colors$red, "[5, 10)" = colors$orange,
                      "[10, 25)" = colors$green, ">25" = colors$blue)
        ) +
        geom_text(
            data = lambda_add,
            aes(x = -Inf, y = Inf, label = lambda_label),
            hjust = -0.2, vjust = 1.5, inherit.aes = FALSE, size = 4
        ) +
        facet_wrap(~ heritability + id, labeller = label_both) +
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

    output_add <- paste0(output_base, "_additive.pdf")
    cat(paste0("Saving additive plot to: ", output_add, "\n"))
    ggsave(filename = output_add, plot = p_add, width = width, height = height, units = "in", dpi = 300)

    # ========== RECESSIVE MODEL PLOT ==========
    cat("\n=== Creating RECESSIVE model QQ plot ===\n")
    lambda_rec <- d[, .(lambda = calc_inflation(REC.P)), by = .(heritability, id)]
    lambda_rec[, lambda_label := paste0("λ = ", round(lambda, 3))]
    print(lambda_rec)

    dt_ribbon_rec <- d[, {
        n <- round(.N * 2)
        seq.p.value <- (1:n) / n
        list(
            REC.PEXPT = -log10(seq.p.value),
            clower = -log10(qbeta(p = (1 + ribbon_p) / 2, shape1 = 1:n, shape2 = n:1)),
            cupper = -log10(qbeta(p = (1 - ribbon_p) / 2, shape1 = 1:n, shape2 = n:1))
        )
    }, by = .(heritability, id)]

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
            values = c("<5" = colors$red, "[5, 10)" = colors$orange,
                      "[10, 25)" = colors$green, ">25" = colors$blue)
        ) +
        geom_text(
            data = lambda_rec,
            aes(x = -Inf, y = Inf, label = lambda_label),
            hjust = -0.2, vjust = 1.5, inherit.aes = FALSE, size = 4
        ) +
        facet_wrap(~ heritability + id, labeller = label_both) +
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

    output_rec <- paste0(output_base, "_recessive.pdf")
    cat(paste0("Saving recessive plot to: ", output_rec, "\n"))
    ggsave(filename = output_rec, plot = p_rec, width = width, height = height, units = "in", dpi = 300)

    # ========== DOMINANCE MODEL PLOT ==========
    cat("\n=== Creating DOMINANCE model QQ plot ===\n")
    lambda_dom <- d[, .(lambda = calc_inflation(DOM.P)), by = .(heritability, id)]
    lambda_dom[, lambda_label := paste0("λ = ", round(lambda, 3))]
    print(lambda_dom)

    dt_ribbon_dom <- d[, {
        n <- round(.N * 2)
        seq.p.value <- (1:n) / n
        list(
            DOM.PEXPT = -log10(seq.p.value),
            clower = -log10(qbeta(p = (1 + ribbon_p) / 2, shape1 = 1:n, shape2 = n:1)),
            cupper = -log10(qbeta(p = (1 - ribbon_p) / 2, shape1 = 1:n, shape2 = n:1))
        )
    }, by = .(heritability, id)]

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
            values = c("<5" = colors$red, "[5, 10)" = colors$orange,
                      "[10, 25)" = colors$green, ">25" = colors$blue)
        ) +
        geom_text(
            data = lambda_dom,
            aes(x = -Inf, y = Inf, label = lambda_label),
            hjust = -0.2, vjust = 1.5, inherit.aes = FALSE, size = 4
        ) +
        facet_wrap(~ heritability + id, labeller = label_both) +
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

    output_dom <- paste0(output_base, "_dominance.pdf")
    cat(paste0("Saving dominance plot to: ", output_dom, "\n"))
    ggsave(filename = output_dom, plot = p_dom, width = width, height = height, units = "in", dpi = 300)

    cat("\n=============================================================\n")
    cat("Done!\n")
    cat("=============================================================\n")
}

# Add arguments
option_list = list(
    make_option(c("-c", "--combined_dir"), type="character", default=NULL,
        help="directory containing combined REGENIE files", metavar="character"),
    make_option(c("-a", "--ancestry"), type="character", default=NULL,
        help="ancestry to filter files", metavar="character"),
    make_option(c("-n", "--annotation"), type="character", default=NULL,
        help="annotation to filter files", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
        help="output PDF file for plots", metavar="character")
)

opt_parser <- OptionParser(add_help_option=TRUE, option_list=option_list)
opt <- parse_args(opt_parser)

main(opt)
