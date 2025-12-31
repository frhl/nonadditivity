#!/usr/bin/env Rscript
# author: frederik lassen
# note: Create SAIGE group/weight files with filtering for minimum homozygous carriers
# This ensures we have sufficient homozygous variants for nonadditive encoding

library(argparse)
library(data.table)

null_omit <- function(lst) {
    lst[!vapply(lst, is.null, logical(1))]
}

# Beta weighting function from SAIGE
beta_weight <- function(MAF, weights.beta=c(1,25)){
    n <- length(MAF)
    weights <- rep(0, n)
    IDX_0 <- which(MAF == 0)
    if(length(IDX_0) == n){
      stop("No polymorphic SNPs")
    } else if(length(IDX_0) == 0){
      weights <- dbeta(MAF, weights.beta[1], weights.beta[2])
    } else {
      weights[-IDX_0] <- dbeta(MAF[-IDX_0], weights.beta[1], weights.beta[2])
    }
    return(weights)
}

main <- function(args){

    cat("=============================================================\n")
    cat("Creating SAIGE group/weight file with homozygous carrier filter\n")
    cat("=============================================================\n\n")

    print(args)

    min_mac <- as.numeric(args$min_mac)
    min_hom_carriers <- as.numeric(args$min_hom_carriers)
    max_af <- as.numeric(args$max_af)

    cat(paste0("Minimum MAC: ", min_mac, "\n"))
    cat(paste0("Minimum homozygous carriers: ", min_hom_carriers, "\n"))
    cat(paste0("Maximum AF: ", max_af, "\n\n"))

    # Load annotation file
    cat(paste0("Loading annotations: ", args$input_path, "\n"))
    d <- fread(args$input_path)
    cat(paste0("  Variants: ", nrow(d), "\n"))

    M <- d[,c("GENE", "SNP_ID", "ANNOTATION")]
    colnames(M) <- c('gene','variant','anno')
    M$weight <- NA

    # Load popmax exclusion list
    cat(paste0("\nLoading popmax exclusion list: ", args$popmax_exclude, "\n"))
    popmax_exclude <- fread(args$popmax_exclude, header=FALSE)$V1
    cat(paste0("  Variants to exclude: ", length(popmax_exclude), "\n"))

    # Load dominance-encoded markers
    cat(paste0("\nLoading dominance markers: ", args$markers_path, "\n"))
    markers <- fread(args$markers_path, header=FALSE)$V1
    cat(paste0("  Markers: ", length(markers), "\n"))

    # Subset to markers in dominance file
    cat(paste0("\nFiltering to dominance markers...\n"))
    cat(paste0("  Before: ", nrow(M), " variants\n"))
    M <- M[M$variant %in% markers,]
    stopifnot(nrow(M) > 0)
    cat(paste0("  After: ", nrow(M), " variants\n"))

    # Load allele counts and calculate weights
    cat(paste0("\nLoading allele counts: ", args$ac_path, "\n"))
    dt_AC <- fread(args$ac_path)
    cat(paste0("  Variants in AC file: ", nrow(dt_AC), "\n"))

    # Calculate allele counts
    cat("\nCalculating allele counts and frequencies...\n")
    dt_AC[, AC_A1 := 2*`C(HOM A1)` + `C(HET)`]
    dt_AC[, AC_A2 := 2*`C(HOM A2)` + `C(HET)`]
    dt_AC[, MAC := pmin(AC_A1, AC_A2)]
    dt_AC[, MAF := MAC/(AC_A1 + AC_A2)]

    # Identify homozygous carriers
    # For each variant, determine which allele is minor and count homozygotes
    dt_AC[, N_HOM_MINOR := ifelse(AC_A1 < AC_A2, `C(HOM A1)`, `C(HOM A2)`)]

    cat("\nHomozygous carrier distribution:\n")
    cat("  Total variants:", nrow(dt_AC), "\n")
    cat("  0 hom carriers:", sum(dt_AC$N_HOM_MINOR == 0), "\n")
    cat("  1-4 hom carriers:", sum(dt_AC$N_HOM_MINOR >= 1 & dt_AC$N_HOM_MINOR < 5), "\n")
    cat("  5-9 hom carriers:", sum(dt_AC$N_HOM_MINOR >= 5 & dt_AC$N_HOM_MINOR < 10), "\n")
    cat("  10+ hom carriers:", sum(dt_AC$N_HOM_MINOR >= 10), "\n")

    # Apply filters
    cat(paste0("\nApplying filters...\n"))
    cat(paste0("  Before filtering: ", nrow(dt_AC), " variants\n"))

    # Remove monomorphic variants
    dt_AC <- dt_AC[dt_AC$MAC > 0,]
    cat(paste0("  After removing monomorphic (MAC=0): ", nrow(dt_AC), " variants\n"))

    # Apply MAF filter
    dt_AC <- dt_AC[dt_AC$MAF < max_af,]
    cat(paste0("  After MAF < ", max_af, " filter: ", nrow(dt_AC), " variants\n"))

    # Exclude popmax variants
    dt_AC <- dt_AC[!dt_AC$SNP %in% popmax_exclude,]
    cat(paste0("  After excluding popmax variants: ", nrow(dt_AC), " variants\n"))

    # Apply MAC filter
    dt_AC <- dt_AC[dt_AC$MAC >= min_mac,]
    cat(paste0("  After MAC >= ", min_mac, " filter: ", nrow(dt_AC), " variants\n"))

    # Apply homozygous carrier filter
    dt_AC <- dt_AC[dt_AC$N_HOM_MINOR >= min_hom_carriers,]
    cat(paste0("  After N_HOM >= ", min_hom_carriers, " filter: ", nrow(dt_AC), " variants\n"))

    # Calculate beta weights
    cat("\nCalculating beta weights...\n")
    dt_AC$weight <- beta_weight(dt_AC$MAF)

    # Filter annotation table to variants passing QC
    cat(paste0("\nFiltering annotations to QC-passing variants...\n"))
    cat(paste0("  Before: ", nrow(M), " variants\n"))
    M <- M[M$variant %in% unique(dt_AC$SNP), ]
    cat(paste0("  After: ", nrow(M), " variants\n"))

    # Create weight mapping
    snp_weights <- dt_AC$weight
    names(snp_weights) <- dt_AC$SNP

    # Generate SAIGE group file format
    cat("\nGenerating SAIGE group file...\n")
    genes <- unique(M$gene)
    cat(paste0("  Genes: ", length(genes), "\n"))

    out <- lapply(genes, function(g){
      variants <- M$variant[M$gene %in% g]
      annotations <- M$anno[M$gene %in% g]
      nas <- (is.na(variants) | is.na(annotations))
      variants <- variants[!nas]
      annotations <- annotations[!nas]
      accepted <- annotations %in% c('pLoF','damaging_missense', "other_missense","synonymous", "non_coding")
      if (sum(accepted) > 0){
          variants <- variants[accepted]
          annotations <- annotations[accepted]
          row1 <- paste(c(g,'var', variants), collapse = " ")
          row2 <- paste(c(g, 'anno', annotations), collapse = " ")
          row3 <- paste(c(g, 'weight', abs(snp_weights[variants])), collapse = " ")
          return(paste0(c(row1, row2, row3), collapse = '\n'))
      }
    })
    out <- null_omit(out)

    # Write output
    cat(paste0("\nWriting to: ", args$output_path, "\n"))
    writeLines(paste(out, collapse = '\n'), args$output_path)

    # Summary statistics
    cat("\n=============================================================\n")
    cat("Summary\n")
    cat("=============================================================\n")
    cat(paste0("Genes with variants: ", length(out), "\n"))
    cat(paste0("Total variants in output: ", nrow(M), "\n"))
    cat(paste0("All variants have >= ", min_hom_carriers, " homozygous carriers\n"))
    cat("=============================================================\n")
}

# Add arguments
parser <- ArgumentParser()
parser$add_argument("--input_path", required=TRUE, help = "Path to annotation file (GENE, SNP_ID, ANNOTATION)")
parser$add_argument("--markers_path", required=TRUE, help = "Path to dominance-encoded marker list")
parser$add_argument("--ac_path", required=TRUE, help = "Path to allele count file (.frqx)")
parser$add_argument("--popmax_exclude", required=TRUE, help = "Path to popmax exclusion list")
parser$add_argument("--output_path", required=TRUE, help = "Output path for SAIGE group file")
parser$add_argument("--max_af", required=TRUE, help = "Maximum allele frequency threshold")
parser$add_argument("--min_mac", default=10, help = "Minimum minor allele count (default: 10)")
parser$add_argument("--min_hom_carriers", default=5, help = "Minimum number of homozygous minor allele carriers (default: 5)")
args <- parser$parse_args()

main(args)
