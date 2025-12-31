#!/usr/bin/env Rscript

options(scipen = 999)

# Required packages
packages <- c('data.table', 'ggplot2', 'ggrastr', 'latex2exp', 'optparse', 'stringr', 'gridExtra')

for (p in packages) {
    if (!require(p, character.only = TRUE, quietly = TRUE)) {
        cat(paste0("Installing package: ", p, "\n"))
        install.packages(p, repos = "http://cran.us.r-project.org")
        library(p, character.only = TRUE)
    }
}

main <- function(opt)
{
    # Check required inputs
    if (is.null(opt$saige_dir) || is.null(opt$regenie_dir) ||
        is.null(opt$ancestry) || is.null(opt$annotation) || is.null(opt$output)) {
        print_help(opt_parser)
        stop("All parameters must be supplied.", call.=FALSE)
    }

    cat("=============================================================\n")
    cat("Comparing SAIGE and REGENIE results\n")
    cat("=============================================================\n\n")

    # List all combined files for SAIGE
    saige_pattern <- paste0(".*", opt$ancestry, ".*", opt$annotation, "\\.combined\\.tsv\\.gz$")
    saige_files <- list.files(opt$saige_dir, pattern = saige_pattern, full.names = TRUE)

    cat(paste0("Found ", length(saige_files), " SAIGE files:\n"))
    for (f in saige_files) {
        cat(paste0("  - ", basename(f), "\n"))
    }
    cat("\n")

    # List all combined files for REGENIE
    regenie_pattern <- paste0(".*", opt$ancestry, ".*", opt$annotation, "\\.combined\\.regenie\\.tsv\\.gz$")
    regenie_files <- list.files(opt$regenie_dir, pattern = regenie_pattern, full.names = TRUE)

    cat(paste0("Found ", length(regenie_files), " REGENIE files:\n"))
    for (f in regenie_files) {
        cat(paste0("  - ", basename(f), "\n"))
    }
    cat("\n")

    if (length(saige_files) == 0 || length(regenie_files) == 0) {
        stop("No files found for comparison.", call.=FALSE)
    }

    # Load SAIGE data
    cat("Loading SAIGE data...\n")
    saige_list <- lapply(saige_files, function(f){
        cat(paste0("  Loading: ", basename(f), "\n"))
        d <- fread(f)
        d[, heritability := stringr::str_extract(basename(f), "(?<=p_)\\d*\\.?\\d+")]
        d[, id := stringr::str_extract(basename(f), "(?<=continuous_)\\d+")]
        d[, source := "SAIGE"]
        return(d)
    })
    saige_data <- rbindlist(saige_list, fill = TRUE)

    # Load REGENIE data
    cat("Loading REGENIE data...\n")
    regenie_list <- lapply(regenie_files, function(f){
        cat(paste0("  Loading: ", basename(f), "\n"))
        d <- fread(f)
        d[, heritability := stringr::str_extract(basename(f), "(?<=p_)\\d*\\.?\\d+")]
        d[, id := stringr::str_extract(basename(f), "(?<=continuous_)\\d+")]
        d[, source := "REGENIE"]
        return(d)
    })
    regenie_data <- rbindlist(regenie_list, fill = TRUE)

    cat(paste0("\nSAIGE variants: ", nrow(saige_data), "\n"))
    cat(paste0("REGENIE variants: ", nrow(regenie_data), "\n\n"))

    # Standardize column names for merging
    # SAIGE uses MarkerID, REGENIE uses ID
    if ("MarkerID" %in% colnames(saige_data)) {
        setnames(saige_data, "MarkerID", "ID")
    }

    # Create list to store plots
    plot_list <- list()

    # Loop through each heritability and replicate
    heritabilities <- sort(unique(c(saige_data$heritability, regenie_data$heritability)))
    ids <- sort(unique(c(saige_data$id, regenie_data$id)))

    cat("Creating comparison plots...\n")
    cat(paste0("Heritabilities: ", paste(heritabilities, collapse=", "), "\n"))
    cat(paste0("IDs: ", paste(ids, collapse=", "), "\n\n"))

    # For each genetic model (additive, recessive, dominance)
    models <- c("ADD", "REC", "DOM")

    for (model in models) {
        cat(paste0("Processing model: ", model, "\n"))

        for (h2 in heritabilities) {
            for (rep_id in ids) {

                # Filter data for this heritability and replicate
                saige_subset <- saige_data[heritability == h2 & id == rep_id]
                regenie_subset <- regenie_data[heritability == h2 & id == rep_id]

                if (nrow(saige_subset) == 0 || nrow(regenie_subset) == 0) {
                    cat(paste0("  Skipping h2=", h2, ", id=", rep_id, " (no data)\n"))
                    next
                }

                # Select relevant columns for this model
                saige_cols <- c("ID", "heritability", "id", paste0(model, ".P"))
                regenie_cols <- c("ID", "heritability", "id", paste0(model, ".P"))

                saige_model <- saige_subset[, ..saige_cols]
                regenie_model <- regenie_subset[, ..regenie_cols]

                # Rename P-value columns
                setnames(saige_model, paste0(model, ".P"), "SAIGE_P")
                setnames(regenie_model, paste0(model, ".P"), "REGENIE_P")

                # Merge by gene ID
                merged <- merge(saige_model, regenie_model, by = c("ID", "heritability", "id"))

                # Remove missing values
                merged <- merged[!is.na(SAIGE_P) & !is.na(REGENIE_P)]
                merged <- merged[is.finite(SAIGE_P) & is.finite(REGENIE_P)]
                merged <- merged[SAIGE_P > 0 & REGENIE_P > 0]

                if (nrow(merged) == 0) {
                    cat(paste0("  Skipping h2=", h2, ", id=", rep_id, ", model=", model, " (no overlapping data)\n"))
                    next
                }

                cat(paste0("  Creating plot for h2=", h2, ", id=", rep_id, ", model=", model, " (n=", nrow(merged), ")\n"))

                # Calculate correlation
                cor_pearson <- cor(merged$SAIGE_P, merged$REGENIE_P, method = "pearson")
                cor_spearman <- cor(merged$SAIGE_P, merged$REGENIE_P, method = "spearman")

                # Transform to -log10 scale
                merged[, SAIGE_logP := -log10(SAIGE_P)]
                merged[, REGENIE_logP := -log10(REGENIE_P)]

                # Create scatter plot
                p <- ggplot(merged, aes(x = SAIGE_logP, y = REGENIE_logP)) +
                    geom_point_rast(dpi = 300, size = 1.5) +
                    geom_abline(intercept = 0, slope = 1, color = "red", linetype = "dashed") +
                    annotate("text", x = -Inf, y = Inf,
                             label = sprintf("Pearson r = %.3f\nSpearman ρ = %.3f\nn = %d",
                                           cor_pearson, cor_spearman, nrow(merged)),
                             hjust = -0.1, vjust = 1.2, size = 3.5) +
                    coord_equal() +
                    labs(
                        title = sprintf("%s Model | h²=%.2f | Rep %s", model, as.numeric(h2), rep_id),
                        x = TeX("$-\\log_{10}(P_{SAIGE})$"),
                        y = TeX("$-\\log_{10}(P_{REGENIE})$")
                    ) +
                    theme_classic() +
                    theme(
                        plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
                        axis.text = element_text(size = 9),
                        axis.title = element_text(size = 11, face = "bold"),
                        plot.margin = margin(t = 5, r = 5, b = 10, l = 10, unit = "pt")
                    )

                plot_list[[paste0(model, "_h2_", h2, "_id_", rep_id)]] <- p
            }
        }
    }

    # Save all plots to a single PDF
    cat(paste0("\nSaving plots to: ", opt$output, "\n"))

    if (length(plot_list) == 0) {
        stop("No plots were generated. Check that data is available for comparison.", call.=FALSE)
    }

    # Arrange plots in grid - 3 columns, multiple pages if needed
    n_plots <- length(plot_list)
    n_per_page <- 9  # 3x3 grid per page
    n_cols <- 3

    # A4 landscape dimensions (in inches)
    width <- 11.69
    height <- 8.27

    cat(paste0("Creating PDF with ", n_plots, " plots\n"))
    cat(paste0("PDF dimensions: ", width, " x ", height, " inches per page\n"))
    cat(paste0("Number of pages: ", ceiling(n_plots / n_per_page), "\n"))

    pdf(opt$output, width = width, height = height)

    # Arrange plots in grid, filling empty spaces with blank plots
    for (i in seq(1, length(plot_list), by = n_per_page)) {
        end_idx <- min(i + n_per_page - 1, length(plot_list))
        plots_on_page <- plot_list[i:end_idx]

        # If this is the last page and it's not full, add empty plots
        n_plots_on_page <- length(plots_on_page)
        if (n_plots_on_page < n_per_page) {
            # Create empty/blank plots to fill the grid
            for (j in (n_plots_on_page + 1):n_per_page) {
                empty_plot <- ggplot() + theme_void()
                plots_on_page[[j]] <- empty_plot
            }
        }

        grid.arrange(grobs = plots_on_page, ncol = n_cols)
    }

    dev.off()

    cat("\n=============================================================\n")
    cat("Done!\n")
    cat("=============================================================\n")
}

# Add arguments
option_list = list(
    make_option(c("-s", "--saige_dir"), type="character", default=NULL,
        help="directory containing SAIGE combined files", metavar="character"),
    make_option(c("-r", "--regenie_dir"), type="character", default=NULL,
        help="directory containing REGENIE combined files", metavar="character"),
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
