#!/usr/bin/env Rscript

options(scipen = 999)

# Required packages
packages <- c('data.table', 'dplyr', 'optparse')

for (p in packages) {
    if (!require(p, character.only = TRUE, quietly = TRUE)) {
        cat(paste0("Installing package: ", p, "\n"))
        install.packages(p, repos = "http://cran.us.r-project.org")
        library(p, character.only = TRUE)
    }
}

main <- function(opt)
{
    if (is.null(opt$input)) {
        print_help(opt_parser)
        stop("Input file must be supplied.", call.=FALSE)
    }

    cat(paste0("Reading input file: ", opt$input, "\n"))
    dt <- fread(opt$input)

    # REGENIE requires FID and IID columns
    # Assuming the ID column is 'eid' (from the phenotype generation)
    id_col <- NULL
    for (col in c("eid", "IID", "participant_eid")) {
        if (col %in% colnames(dt)) {
            id_col <- col
            break
        }
    }

    if (is.null(id_col)) {
        stop("Could not find ID column (eid, IID, or participant_eid)")
    }

    cat(paste0("Using ID column: ", id_col, "\n"))

    # Create output with FID IID format
    output_dt <- data.table(
        FID = dt[[id_col]],
        IID = dt[[id_col]]
    )

    # Add all other columns (phenotypes and covariates)
    other_cols <- setdiff(colnames(dt), id_col)
    output_dt <- cbind(output_dt, dt[, ..other_cols])

    cat(paste0("Output dimensions: ", nrow(output_dt), " rows, ", ncol(output_dt), " columns\n"))
    cat(paste0("Columns: ", paste(colnames(output_dt), collapse=", "), "\n"))

    # Verify FID and IID are the first two columns
    if (colnames(output_dt)[1] != "FID" || colnames(output_dt)[2] != "IID") {
        stop("First two columns must be FID and IID")
    }

    # Write output with space separator (REGENIE requirement)
    cat(paste0("Writing to: ", opt$output, "\n"))
    fwrite(output_dt, file=opt$output, sep=" ", na="NA", quote=FALSE, col.names=TRUE)

    # Verify output
    cat("First few lines of output:\n")
    system(paste0("head -3 ", opt$output))

    cat("Done!\n")
}

# Add arguments
option_list = list(
    make_option(c("-i", "--input"), type="character", default=NULL,
        help="input phenotype file", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
        help="output file for REGENIE", metavar="character")
)

opt_parser <- OptionParser(add_help_option=TRUE, option_list=option_list)
opt <- parse_args(opt_parser)

main(opt)
