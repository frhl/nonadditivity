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

# Beta weight function from SAIGE
beta_weight <- function(AAF, weights.beta=c(1,25)){
    n <- length(AAF)
    weights <- rep(0, n)
    IDX_0 <- which(AAF == 0)
    if(length(IDX_0) == n){
      stop("No polymorphic SNPs")
    } else if(length(IDX_0) == 0){
      weights <- dbeta(AAF, weights.beta[1], weights.beta[2])
    } else {
      weights[-IDX_0] <- dbeta(AAF[-IDX_0], weights.beta[1], weights.beta[2])
    }
    return(weights)
}

main <- function(opt)
{
    cat("=============================================================\n")
    cat("Creating REGENIE annotation and set list files from VEP data\n")
    cat("=============================================================\n\n")

    # Load VEP canonical annotation file
    cat(paste0("Loading VEP file: ", opt$vep_file, "\n"))
    vep <- fread(opt$vep_file)
    cat(paste0("  Total variants in VEP: ", nrow(vep), "\n"))
    cat(paste0("  Columns: ", paste(colnames(vep), collapse=", "), "\n\n"))

    # Load markers file (variants to keep)
    cat(paste0("Loading markers file: ", opt$markers_file, "\n"))
    markers <- fread(opt$markers_file, header=FALSE)$V1
    cat(paste0("  Markers to keep: ", length(markers), "\n\n"))

    # Load allele count file for AAF calculation
    cat(paste0("Loading AC file: ", opt$ac_file, "\n"))
    ac <- fread(opt$ac_file)
    cat(paste0("  AC variants: ", nrow(ac), "\n"))
    cat(paste0("  Columns: ", paste(colnames(ac), collapse=", "), "\n\n"))

    # Calculate AAF and identify singletons from AC data
    cat("Calculating AAF and identifying singletons...\n")
    ac[, AC_A1 := 2*`C(HOM A1)` + `C(HET)`]
    ac[, AC_A2 := 2*`C(HOM A2)` + `C(HET)`]
    ac[, AN := AC_A1 + AC_A2]

    # AAF is for the alternate allele (A1 in the frqx file)
    ac[, AAF := AC_A1 / AN]

    # Singleton: AC_A1 == 1
    ac[, IS_SINGLETON := ifelse(AC_A1 == 1, 1, 0)]

    # Calculate beta weights
    ac$WEIGHT <- beta_weight(ac$AAF)

    cat(paste0("  Singletons: ", sum(ac$IS_SINGLETON), "\n"))
    cat(paste0("  AAF range: ", min(ac$AAF, na.rm=TRUE), " - ", max(ac$AAF, na.rm=TRUE), "\n\n"))

    # Filter VEP data to only markers in dominance encoding
    cat("Filtering VEP data to markers...\n")
    vep_filtered <- vep[vep$SNP_ID %in% markers, ]
    cat(paste0("  Variants after filtering: ", nrow(vep_filtered), "\n\n"))

    if (nrow(vep_filtered) == 0) {
        stop("No variants remaining after filtering!")
    }

    # Merge VEP with AC data to get weights
    cat("Merging VEP with AC data...\n")
    vep_filtered <- merge(vep_filtered, ac[, .(SNP, AAF, WEIGHT, IS_SINGLETON)],
                          by.x="SNP_ID", by.y="SNP", all.x=TRUE)
    cat(paste0("  Variants with AC data: ", sum(!is.na(vep_filtered$AAF)), "\n"))
    cat(paste0("  Variants missing AC data: ", sum(is.na(vep_filtered$AAF)), "\n\n"))

    # Create REGENIE annotation file with weights
    # Extended format: VARIANT_ID GENE_NAME ANNOTATION WEIGHT
    cat("Creating REGENIE annotation file with weights...\n")
    anno <- vep_filtered[, .(SNP_ID, GENE, ANNOTATION, WEIGHT)]

    # Remove any rows with missing data
    anno <- anno[!is.na(SNP_ID) & !is.na(GENE) & !is.na(ANNOTATION) & !is.na(WEIGHT), ]

    cat(paste0("  Annotation entries: ", nrow(anno), "\n"))
    cat(paste0("  Unique genes: ", length(unique(anno$GENE)), "\n"))
    cat(paste0("  Unique variants: ", length(unique(anno$SNP_ID)), "\n"))
    cat(paste0("  Weight range: ", min(anno$WEIGHT, na.rm=TRUE), " - ", max(anno$WEIGHT, na.rm=TRUE), "\n\n"))

    # Show annotation distribution
    cat("Annotation distribution:\n")
    print(table(anno$ANNOTATION))
    cat("\n")

    # Write annotation file (space-separated, no header)
    cat(paste0("Writing annotation file: ", opt$out_anno, "\n"))
    fwrite(anno, file=opt$out_anno, sep=" ", col.names=FALSE, quote=FALSE)

    # Create REGENIE set list file
    # Format: GENE_NAME CHR POS VARIANT1,VARIANT2,...
    cat("\nCreating REGENIE set list file...\n")

    genes <- unique(anno$GENE)
    cat(paste0("  Processing ", length(genes), " genes...\n"))

    setlist <- lapply(genes, function(gene) {
        # Get all variants for this gene
        gene_variants <- anno$SNP_ID[anno$GENE == gene]

        # Extract chromosome and position from first variant
        # Format: chr21:10538674:C:A
        first_var <- gene_variants[1]
        var_parts <- strsplit(first_var, ":")[[1]]
        chr <- var_parts[1]
        pos <- var_parts[2]

        # Create comma-separated variant list
        variant_list <- paste(gene_variants, collapse=",")

        # Return gene row: GENE CHR POS VARIANT_LIST
        return(data.table(
            GENE = gene,
            CHR = chr,
            POS = pos,
            VARIANTS = variant_list
        ))
    })

    setlist_dt <- rbindlist(setlist)

    cat(paste0("  Set list entries: ", nrow(setlist_dt), "\n"))
    cat(paste0("  Chromosomes: ", paste(unique(setlist_dt$CHR), collapse=", "), "\n\n"))

    # Write set list file (space-separated, no header)
    cat(paste0("Writing set list file: ", opt$out_setlist, "\n"))
    fwrite(setlist_dt, file=opt$out_setlist, sep=" ", col.names=FALSE, quote=FALSE)

    # Create AAF file with singleton indicators
    # Format: VARIANT_ID AAF IS_SINGLETON
    cat("\nCreating AAF file with singleton indicators...\n")
    aaf_dt <- vep_filtered[, .(SNP_ID, AAF, IS_SINGLETON)]

    # Remove duplicates and missing data
    aaf_dt <- unique(aaf_dt)
    aaf_dt <- aaf_dt[!is.na(SNP_ID) & !is.na(AAF) & !is.na(IS_SINGLETON), ]

    cat(paste0("  AAF entries: ", nrow(aaf_dt), "\n"))
    cat(paste0("  Singletons: ", sum(aaf_dt$IS_SINGLETON), "\n"))
    cat(paste0("  AAF range: ", min(aaf_dt$AAF, na.rm=TRUE), " - ", max(aaf_dt$AAF, na.rm=TRUE), "\n\n"))

    # Write AAF file (space-separated, no header)
    cat(paste0("Writing AAF file: ", opt$out_aaf, "\n"))
    fwrite(aaf_dt, file=opt$out_aaf, sep=" ", col.names=FALSE, quote=FALSE)

    # Show first few lines of each file
    cat("\n=============================================================\n")
    cat("Sample output\n")
    cat("=============================================================\n")
    cat("\nAnnotation file (first 10 lines):\n")
    print(head(anno, 10))

    cat("\nSet list file (first 5 entries):\n")
    print(head(setlist_dt, 5))

    cat("\nAAF file (first 10 lines):\n")
    print(head(aaf_dt, 10))

    cat("\n=============================================================\n")
    cat("Done!\n")
    cat("=============================================================\n")
}

# Add arguments
option_list = list(
    make_option(c("--vep_file"), type="character", default=NULL,
        help="VEP canonical annotation file", metavar="character"),
    make_option(c("--markers_file"), type="character", default=NULL,
        help="Markers file (variants to keep)", metavar="character"),
    make_option(c("--ac_file"), type="character", default=NULL,
        help="Allele count file (frqx format)", metavar="character"),
    make_option(c("--out_anno"), type="character", default=NULL,
        help="Output annotation file for REGENIE", metavar="character"),
    make_option(c("--out_setlist"), type="character", default=NULL,
        help="Output set list file for REGENIE", metavar="character"),
    make_option(c("--out_aaf"), type="character", default=NULL,
        help="Output AAF file for REGENIE", metavar="character"),
    make_option(c("--chr"), type="character", default=NULL,
        help="Chromosome (for logging)", metavar="character")
)

opt_parser <- OptionParser(add_help_option=TRUE, option_list=option_list)
opt <- parse_args(opt_parser)

main(opt)
