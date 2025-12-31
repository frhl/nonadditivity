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
    cat("Combining SAIGE results from three genetic models\n")
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

    # Check column names
    cat("Additive columns:", paste(colnames(d_add), collapse=", "), "\n")
    cat("Recessive columns:", paste(colnames(d_rec), collapse=", "), "\n")
    cat("Dominance columns:", paste(colnames(d_dom), collapse=", "), "\n\n")

    # SAIGE columns: CHR, POS, MarkerID, Allele1, Allele2, AC_Allele2, AF_Allele2, MissingRate, BETA, SE, Tstat, var, p.value, N

    # Select and rename columns for additive model
    cat("Selecting and renaming columns...\n")
    add_cols <- d_add[, .(
        CHR, POS, MarkerID, Allele1, Allele2, AC_Allele2, AF_Allele2, MissingRate, N,
        ADD.BETA = BETA,
        ADD.SE = SE,
        ADD.Tstat = Tstat,
        ADD.var = var,
        ADD.P = p.value
    )]

    # For recessive model, we need to calculate N_HOMALT
    # In SAIGE recessive mode, AC_Allele2 represents the count of homozygous alternate genotypes
    # Because recessive encoding treats only homozygous alts as "cases"
    # Need to divide by 2 because AC_Allele2=4 whill indicate 4 varants, i.e. 2 biallelic variants.
    d_rec[, N_HOMALT := round(AC_Allele2) / 2]
    rec_cols <- d_rec[, .(
        MarkerID,
        N_HOMALT,
        REC.BETA = BETA,
        REC.SE = SE,
        REC.Tstat = Tstat,
        REC.var = var,
        REC.P = p.value
    )]

    # Dominance model
    dom_cols <- d_dom[, .(
        MarkerID,
        DOM.BETA = BETA,
        DOM.SE = SE,
        DOM.Tstat = Tstat,
        DOM.var = var,
        DOM.P = p.value
    )]

    # Merge all three models
    cat("Merging additive and recessive models...\n")
    combined <- merge(add_cols, rec_cols, by = "MarkerID", all = FALSE)

    cat("Merging with dominance model...\n")
    combined <- merge(combined, dom_cols, by = "MarkerID", all = FALSE)

    cat(paste0("Combined variants: ", nrow(combined), "\n\n"))

    # Calculate allele counts and frequencies
    cat("Calculating allele counts and frequencies...\n")
    combined[, AF_Allele1 := 1 - AF_Allele2]
    combined[, AC_Allele1 := round(AF_Allele1 * N * 2)]

    # Calculate N_HET (heterozygotes)
    # AC_Allele2 = 2*N_HOMALT + N_HET (for the alternate allele)
    # So N_HET = AC_Allele2 - 2*N_HOMALT
    cat("Calculating heterozygote counts...\n")
    combined[, N_HET := AC_Allele2 - (2 * N_HOMALT)]

    # Create AC2 interval bins (AC_Allele2 is the alternate allele count)
    cat("Creating AC2 interval bins...\n")
    combined[, AC2_interval := cut(
      AC_Allele2,
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

    # Reorder columns for consistency with REGENIE format
    cat("Reordering columns...\n")
    col_order <- c("CHR", "POS", "MarkerID", "Allele1", "Allele2",
                   "AF_Allele1", "AF_Allele2", "AC_Allele1", "AC_Allele2", "AC2_interval",
                   "N_HET", "N_HOMALT", "N_HOMALT_BIN", "MissingRate", "N",
                   "ADD.BETA", "ADD.SE", "ADD.Tstat", "ADD.var", "ADD.P",
                   "REC.BETA", "REC.SE", "REC.Tstat", "REC.var", "REC.P",
                   "DOM.BETA", "DOM.SE", "DOM.Tstat", "DOM.var", "DOM.P")

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
        help="additive model SAIGE file", metavar="character"),
    make_option(c("-r", "--recessive"), type="character", default=NULL,
        help="recessive model SAIGE file", metavar="character"),
    make_option(c("-d", "--dominance"), type="character", default=NULL,
        help="dominance model SAIGE file", metavar="character"),
    make_option(c("-o", "--output"), type="character", default=NULL,
        help="output file for combined results", metavar="character")
)

opt_parser <- OptionParser(add_help_option=TRUE, option_list=option_list)
opt <- parse_args(opt_parser)

main(opt)
