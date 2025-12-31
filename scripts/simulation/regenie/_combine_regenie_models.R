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
    # Check required inputs
    if (is.null(opt$additive) || is.null(opt$recessive) || is.null(opt$dominance)) {
        print_help(opt_parser)
        stop("All three model files (additive, recessive, dominance) must be supplied.", call.=FALSE)
    }

    if (is.null(opt$output)) {
        print_help(opt_parser)
        stop("Output file must be supplied.", call.=FALSE)
    }

    cat("=============================================================\n")
    cat("Combining REGENIE results from three genetic models\n")
    cat("=============================================================\n\n")

    # Load the files
    cat(paste0("Loading additive model: ", opt$additive, "\n"))
    d_add <- fread(opt$additive)
    cat(paste0("  Variants: ", nrow(d_add), "\n"))

    cat(paste0("Loading recessive model: ", opt$recessive, "\n"))
    d_rec <- fread(opt$recessive)
    cat(paste0("  Variants: ", nrow(d_rec), "\n"))

    cat(paste0("Loading dominance model: ", opt$dominance, "\n"))
    d_dom <- fread(opt$dominance)
    cat(paste0("  Variants: ", nrow(d_dom), "\n\n"))

    # Select and rename columns for each model
    cat("Selecting and renaming columns...\n")

    # Check which columns are available in each file
    cat("Additive columns:", paste(colnames(d_add), collapse=", "), "\n")
    cat("Recessive columns:", paste(colnames(d_rec), collapse=", "), "\n")
    cat("Dominance columns:", paste(colnames(d_dom), collapse=", "), "\n\n")

    # Build additive selection dynamically based on available columns
    add_base_cols <- c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N")
    if ("INFO" %in% colnames(d_add)) {
        add_base_cols <- c(add_base_cols[1:6], "INFO", "N")
    }

    add_cols <- d_add[, c(add_base_cols, "BETA", "SE", "LOG10P", "CHISQ"), with=FALSE]
    setnames(add_cols, c("BETA", "SE", "LOG10P", "CHISQ"), c("ADD.BETA", "ADD.SE", "ADD.P", "ADD.CHISQ"))

    # For recessive: calculate N_HOMALT from A1FREQ in recessive encoding
    # In recessive BGEN encoding [0,1,2]->[0,0,2]: only homozygous alternates have dosage=2
    # A2FREQ in recessive model = frequency of homozygous alternate individuals
    d_rec[, N_HOMALT := round((1 - A1FREQ) * N)]
    rec_cols <- d_rec[, .(ID, N_HOMALT, REC.BETA = BETA, REC.SE = SE, REC.P = LOG10P, REC.CHISQ = CHISQ)]

    dom_cols <- d_dom[, .(ID, DOM.BETA = BETA, DOM.SE = SE, DOM.P = LOG10P, DOM.CHISQ = CHISQ)]

    # Merge with inner joins (only keep variants present in all three files)
    # This matches the behavior of SAIGE combine script
    cat("Merging additive and recessive models...\n")
    combined <- merge(add_cols, rec_cols, by = "ID", all = FALSE)

    cat("Merging with dominance model...\n")
    combined <- merge(combined, dom_cols, by = "ID", all = FALSE)

    cat(paste0("Combined variants: ", nrow(combined), "\n\n"))

    # Calculate allele counts and frequencies
    cat("Calculating allele counts and frequencies...\n")
    combined[, A2FREQ := 1 - A1FREQ]
    combined[, AC1 := round(A1FREQ * N * 2)]
    combined[, AC2 := round(A2FREQ * N * 2)]  # AC2 = ALTERNATE/VARIANT allele count

    # N_HOMALT was already calculated from recessive file before merging
    # Calculate N_HET (heterozygotes): AC_rare = 2*N_HOMALT + N_HET, so N_HET = AC_rare - 2*N_HOMALT
    # Use the minor allele count (whichever is smaller: AC1 or AC2)
    cat("Calculating heterozygote counts...\n")
    combined[, AC_rare := pmin(AC1, AC2)]
    combined[, N_HET := AC_rare - (2 * N_HOMALT)]

    # Create AC2 interval bins
    cat("Creating AC2 interval bins...\n")
    combined[, AC2_interval := cut(
      AC2,
      breaks = c(0, 5, 10, 25, Inf),
      labels = c("<5", "[5, 10)", "[10, 25)", ">25"),
      include.lowest = TRUE
    )]

    # Create N_HOMALT bins
    cat("Creating N_HOMALT bins...\n")
    combined[, N_HOMALT_BIN := cut(
      N_HOMALT,
      breaks = c(-1, 0, 1, 5, 10, 25, Inf),
      labels = c("0", "1", "2-4", "5-9", "10-24", ">=25"),
      include.lowest = TRUE
    )]

    # ensure that we are not on log10 scale anymore
    combined[ , REC.P := 10^(-REC.P)]
    combined[ , DOM.P := 10^(-DOM.P)]
    combined[ , ADD.P := 10^(-ADD.P)]


    # Reorder columns (dynamically handle INFO column)
    cat("Reordering columns...\n")
    col_order <- c("CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1",
                   "A1FREQ", "A2FREQ", "AC1", "AC2", "AC2_interval",
                   "N_HET", "N_HOMALT", "N_HOMALT_BIN")

    if ("INFO" %in% colnames(combined)) {
        col_order <- c(col_order, "INFO")
    }

    col_order <- c(col_order, "N",
                   "ADD.BETA", "ADD.SE", "ADD.P", "ADD.CHISQ",
                   "REC.BETA", "REC.SE", "REC.P", "REC.CHISQ",
                   "DOM.BETA", "DOM.SE", "DOM.P", "DOM.CHISQ")

    setcolorder(combined, col_order)

    # Summary statistics
    cat("\n=============================================================\n")
    cat("Summary Statistics\n")
    cat("=============================================================\n")
    cat(paste0("Total variants: ", nrow(combined), "\n"))
    cat(paste0("Columns: ", ncol(combined), "\n\n"))

    cat("AC2 interval distribution:\n")
    print(table(combined$AC2_interval, useNA = "ifany"))
    cat("\n")

    cat("N_HOMALT bin distribution:\n")
    print(table(combined$N_HOMALT_BIN, useNA = "ifany"))
    cat("\n")

    cat("Variants with missing data:\n")
    cat(paste0("  Missing additive results: ", sum(is.na(combined$ADD.BETA)), "\n"))
    cat(paste0("  Missing recessive results: ", sum(is.na(combined$REC.BETA)), "\n"))
    cat(paste0("  Missing dominance results: ", sum(is.na(combined$DOM.BETA)), "\n\n"))

    # Write output
    cat(paste0("Writing output to: ", opt$output, "\n"))
    fwrite(combined, file=opt$output, sep="\t", na="NA", quote=FALSE, col.names=TRUE)

    # Verify output
    cat("\nFirst few lines of output:\n")
    cat(paste0(capture.output(head(combined, 3)), collapse="\n"), "\n")

    cat("\n=============================================================\n")
    cat("Done!\n")
    cat("=============================================================\n")
}

# Add arguments
option_list = list(
    make_option(c("-a", "--additive"), type="character", default=NULL,
        help="additive model REGENIE file", metavar="character"),
    make_option(c("-r", "--recessive"), type="character", default=NULL,
        help="recessive model REGENIE file", metavar="character"),
    make_option(c("-d", "--dominance"), type="character", default=NULL,
        help="dominance model REGENIE file", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
        help="output file for combined results", metavar="character")
)

opt_parser <- OptionParser(add_help_option=TRUE, option_list=option_list)
opt <- parse_args(opt_parser)

main(opt)
