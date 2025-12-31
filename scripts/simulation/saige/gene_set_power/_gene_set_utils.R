#!/usr/bin/env Rscript
# Gene set utilities for simulation framework
# Handles SAIGE weight file parsing, gene/variant filtering, and genotype loading

library(data.table)
library(BEDMatrix)


#' Parse SAIGE weight file format
#'
#' SAIGE weight files have 3 rows per gene:
#'   Row 1: GENE_ID var variant1 variant2 ...
#'   Row 2: GENE_ID anno anno1 anno2 ...
#'   Row 3: GENE_ID weight weight1 weight2 ...
#'
#' @param saige_weight_file Path to SAIGE weight file (can be gzipped)
#' @param annotations_to_keep Vector of annotations to keep (e.g., c("pLoF", "damaging_missense"))
#' @return data.table with columns: gene, variant, annotation, weight
#' @examples
#' dt <- parse_saige_weights("weights.txt.gz", c("pLoF", "damaging_missense"))
parse_saige_weights <- function(saige_weight_file, annotations_to_keep = NULL) {
    cat(paste0("Parsing SAIGE weight file: ", saige_weight_file, "\n"))

    # Read all lines
    if (grepl("\\.gz$", saige_weight_file)) {
        lines <- readLines(gzfile(saige_weight_file))
    } else {
        lines <- readLines(saige_weight_file)
    }

    cat(paste0("  Total lines: ", length(lines), "\n"))

    # Process in groups of 3
    n_genes <- length(lines) / 3
    if (length(lines) %% 3 != 0) {
        stop("SAIGE weight file does not have a multiple of 3 lines")
    }

    cat(paste0("  Genes in file: ", n_genes, "\n"))

    # Parse each gene
    gene_list <- list()
    for (i in seq(1, length(lines), by = 3)) {
        var_line <- strsplit(lines[i], split = " ")[[1]]
        anno_line <- strsplit(lines[i+1], split = " ")[[1]]
        weight_line <- strsplit(lines[i+2], split = " ")[[1]]

        gene_id <- var_line[1]
        variants <- var_line[-c(1, 2)]  # Remove gene ID and "var"
        annotations <- anno_line[-c(1, 2)]  # Remove gene ID and "anno"
        weights <- as.numeric(weight_line[-c(1, 2)])  # Remove gene ID and "weight"

        # Sanity check
        if (length(variants) != length(annotations) || length(variants) != length(weights)) {
            stop(sprintf("Mismatch in gene %s: %d variants, %d annotations, %d weights",
                        gene_id, length(variants), length(annotations), length(weights)))
        }

        gene_list[[gene_id]] <- data.table(
            gene = gene_id,
            variant = variants,
            annotation = annotations,
            weight = weights
        )
    }

    # Combine all genes
    dt <- rbindlist(gene_list)

    cat(paste0("  Total variants before filtering: ", nrow(dt), "\n"))

    # Filter by annotation if specified
    if (!is.null(annotations_to_keep)) {
        dt <- dt[annotation %in% annotations_to_keep]
        cat(paste0("  Total variants after annotation filter: ", nrow(dt), "\n"))
        cat(paste0("  Annotations kept: ", paste(unique(dt$annotation), collapse = ", "), "\n"))
    }

    cat(paste0("  Final genes: ", length(unique(dt$gene)), "\n"))

    return(dt)
}


#' Load allele counts and calculate MAF
#'
#' Reads PLINK .frqx file and calculates MAC, MAF, carrier counts
#'
#' @param frqx_file Path to PLINK .frqx file (can be gzipped)
#' @return data.table with columns: variant (CHR:POS:REF:ALT), MAC, MAF, N_HET, N_HOM_ALT
load_allele_counts <- function(frqx_file) {
    cat(paste0("Loading allele counts from: ", frqx_file, "\n"))

    dt <- fread(frqx_file)

    cat(paste0("  Variants in AC file: ", nrow(dt), "\n"))

    # Calculate allele counts
    dt[, AC_A1 := 2*`C(HOM A1)` + `C(HET)`]
    dt[, AC_A2 := 2*`C(HOM A2)` + `C(HET)`]
    dt[, MAC := pmin(AC_A1, AC_A2)]
    dt[, MAF := MAC/(AC_A1 + AC_A2)]

    # Carrier counts
    dt[, N_HET := `C(HET)`]
    dt[, N_HOM_ALT := pmin(`C(HOM A1)`, `C(HOM A2)`)]  # Homozygous for minor allele

    # Create variant ID in format chr:pos:ref:alt
    # SNP column format is typically chr:pos_ref/alt or similar
    # We'll keep the original SNP column as variant identifier
    setnames(dt, "SNP", "variant")

    # Select relevant columns
    dt <- dt[, .(variant, MAC, MAF, N_HET, N_HOM_ALT)]

    cat(paste0("  MAF range: [", sprintf("%.6f", min(dt$MAF)), ", ", sprintf("%.6f", max(dt$MAF)), "]\n"))
    cat(paste0("  MAC range: [", min(dt$MAC), ", ", max(dt$MAC), "]\n"))

    return(dt)
}


#' Filter genes by criteria
#'
#' Applies multiple filters:
#'   - Annotation types
#'   - MAF range
#'   - Minimum variants per gene
#'   - Note: ≥5 homozygous carriers already enforced in VCF generation
#'
#' @param gene_data data.table from parse_saige_weights()
#' @param allele_counts data.table from load_allele_counts()
#' @param min_variants Minimum number of variants per gene (default: 1)
#' @param maf_range Numeric vector of length 2: c(min_maf, max_maf)
#' @param allowed_annotations Vector of annotation types to keep
#' @return Filtered data.table with merged allele count information
filter_genes <- function(gene_data, allele_counts,
                         min_variants = 1,
                         maf_range = c(0.001, 0.05),
                         allowed_annotations = c("pLoF", "damaging_missense")) {

    cat("\n=== Filtering genes ===\n")
    cat(paste0("Starting variants: ", nrow(gene_data), "\n"))
    cat(paste0("Starting genes: ", length(unique(gene_data$gene)), "\n\n"))

    # 1. Filter by annotation (already done in parse, but apply again for safety)
    dt <- gene_data[annotation %in% allowed_annotations]
    cat(paste0("After annotation filter: ", nrow(dt), " variants, ",
               length(unique(dt$gene)), " genes\n"))
    cat(paste0("  Annotations: ", paste(allowed_annotations, collapse = ", "), "\n"))

    # 2. Merge with allele counts to get MAF
    dt <- merge(dt, allele_counts, by = "variant", all.x = FALSE, all.y = FALSE)
    cat(paste0("After merging with AC: ", nrow(dt), " variants, ",
               length(unique(dt$gene)), " genes\n"))

    # 3. Filter by MAF range
    dt <- dt[MAF >= maf_range[1] & MAF < maf_range[2]]
    cat(paste0("After MAF filter [", maf_range[1], ", ", maf_range[2], "): ",
               nrow(dt), " variants, ", length(unique(dt$gene)), " genes\n"))

    # 4. Count variants per gene
    gene_variant_counts <- dt[, .(n_variants = .N), by = gene]

    # 5. Filter genes with sufficient variants
    eligible_genes <- gene_variant_counts[n_variants >= min_variants]$gene
    dt <- dt[gene %in% eligible_genes]

    cat(paste0("After min_variants filter (>= ", min_variants, "): ",
               nrow(dt), " variants, ", length(unique(dt$gene)), " genes\n"))

    # Summary statistics
    cat("\n=== Summary ===\n")
    cat(paste0("Final genes: ", length(unique(dt$gene)), "\n"))
    cat(paste0("Final variants: ", nrow(dt), "\n"))
    cat(paste0("Variants per gene - Mean: ", sprintf("%.1f", mean(gene_variant_counts[gene %in% eligible_genes]$n_variants)),
               ", Median: ", median(gene_variant_counts[gene %in% eligible_genes]$n_variants),
               ", Range: [", min(gene_variant_counts[gene %in% eligible_genes]$n_variants),
               ", ", max(gene_variant_counts[gene %in% eligible_genes]$n_variants), "]\n"))

    return(dt)
}


#' Select K causal genes randomly
#'
#' @param gene_data data.table with gene information
#' @param K Number of causal genes to select
#' @param seed Random seed
#' @return Character vector of gene IDs
select_causal_genes <- function(gene_data, K, seed = NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }

    available_genes <- unique(gene_data$gene)
    n_available <- length(available_genes)

    if (K > n_available) {
        stop(sprintf("Requested %d causal genes, but only %d available", K, n_available))
    }

    selected_genes <- sample(available_genes, K, replace = FALSE)

    cat(paste0("\nSelected ", K, " causal genes from ", n_available, " available\n"))

    return(selected_genes)
}


#' Select N causal variants per gene
#'
#' Supports multiple selection strategies:
#'   - "random": Uniform random selection
#'   - "lowest_maf": Select N rarest variants
#'   - "highest_maf": Select N most common variants
#'   - "weighted": Probabilistic selection based on beta weights
#'
#' @param gene_data data.table with variant information (must have MAF and weight columns)
#' @param causal_genes Character vector of gene IDs
#' @param N_per_gene Number of causal variants to select per gene
#' @param selection_strategy Strategy for selection (default: "random")
#' @param seed Random seed
#' @return data.table with selected causal variants
select_causal_variants_per_gene <- function(gene_data,
                                            causal_genes,
                                            N_per_gene = 1,
                                            selection_strategy = "random",
                                            seed = NULL) {

    if (!is.null(seed)) {
        set.seed(seed)
    }

    cat(paste0("\n=== Selecting causal variants ===\n"))
    cat(paste0("Strategy: ", selection_strategy, "\n"))
    cat(paste0("N per gene: ", N_per_gene, "\n"))

    # Filter to causal genes only
    gene_data_causal <- gene_data[gene %in% causal_genes]

    # Select variants per gene
    causal_list <- lapply(causal_genes, function(gene_id) {
        variants_in_gene <- gene_data_causal[gene == gene_id]
        n_available <- nrow(variants_in_gene)

        if (N_per_gene > n_available) {
            warning(sprintf("Gene %s has only %d variants, using all (requested %d)",
                          gene_id, n_available, N_per_gene))
            N_select <- n_available
        } else {
            N_select <- N_per_gene
        }

        # Apply selection strategy
        if (selection_strategy == "random") {
            idx <- sample(1:n_available, N_select, replace = FALSE)
        } else if (selection_strategy == "lowest_maf") {
            idx <- order(variants_in_gene$MAF)[1:N_select]
        } else if (selection_strategy == "highest_maf") {
            idx <- order(variants_in_gene$MAF, decreasing = TRUE)[1:N_select]
        } else if (selection_strategy == "weighted") {
            prob <- variants_in_gene$weight / sum(variants_in_gene$weight)
            idx <- sample(1:n_available, N_select, prob = prob, replace = FALSE)
        } else {
            stop(sprintf("Unknown selection strategy: %s", selection_strategy))
        }

        return(variants_in_gene[idx, ])
    })

    # Combine all selected variants
    causal_variants <- rbindlist(causal_list)

    cat(paste0("Total causal variants: ", nrow(causal_variants), "\n"))
    cat(paste0("MAF range: [", sprintf("%.6f", min(causal_variants$MAF)),
               ", ", sprintf("%.6f", max(causal_variants$MAF)), "]\n"))
    cat(paste0("Mean MAF: ", sprintf("%.6f", mean(causal_variants$MAF)), "\n"))

    return(causal_variants)
}


#' Load genotypes from PLINK files
#'
#' Reads PLINK .bed/.bim/.fam files for specific variants
#' Uses BEDMatrix for memory-efficient reading
#'
#' @param bed_file Path to PLINK .bed file (without extension)
#' @param variant_list Character vector of variant IDs (chr:pos:ref:alt format)
#' @return Matrix of genotypes [samples x variants]
load_genotypes_from_plink <- function(bed_file, variant_list) {

    cat(paste0("\n=== Loading genotypes from PLINK ===\n"))
    cat(paste0("PLINK prefix: ", bed_file, "\n"))
    cat(paste0("Variants requested: ", length(variant_list), "\n"))

    # Read PLINK files using BEDMatrix (memory-efficient)
    bed <- BEDMatrix(paste0(bed_file, ".bed"))

    # Read .fam file for sample IDs
    fam <- fread(paste0(bed_file, ".fam"), header = FALSE)
    colnames(fam) <- c("FID", "IID", "father", "mother", "sex", "pheno")
    sample_ids <- fam$IID

    # Read .bim file for variant IDs
    bim <- fread(paste0(bed_file, ".bim"), header = FALSE)
    colnames(bim) <- c("chr", "ID", "cm", "pos", "A1", "A2")

    cat(paste0("PLINK file contains: ", nrow(bed), " samples x ", ncol(bed), " variants\n"))

    # Match requested variants to PLINK variants
    # PLINK IDs should be in format: chr:pos:ref:alt (from --set-all-var-ids '@:#:$r:$a')
    variant_indices <- match(variant_list, bim$ID)

    # Check for missing variants
    n_missing <- sum(is.na(variant_indices))
    if (n_missing > 0) {
        cat(sprintf("WARNING: %d/%d requested variants not found in PLINK file\n",
                   n_missing, length(variant_list)))

        # Show a few missing variants for debugging
        missing_vars <- variant_list[is.na(variant_indices)]
        cat("First few missing variants:\n")
        print(head(missing_vars, 5))
    }

    # Filter to variants that exist
    valid_indices <- variant_indices[!is.na(variant_indices)]
    valid_variant_list <- variant_list[!is.na(variant_indices)]

    if (length(valid_indices) == 0) {
        stop("No matching variants found between request and PLINK file")
    }

    cat(paste0("Loading ", length(valid_indices), " variants...\n"))

    # Extract genotype matrix for requested variants
    # BEDMatrix is [samples x variants], so we extract columns
    genotype_matrix <- as.matrix(bed[, valid_indices, drop = FALSE])

    # Set column names to variant IDs
    colnames(genotype_matrix) <- valid_variant_list
    rownames(genotype_matrix) <- sample_ids

    cat(paste0("Genotype matrix: ", nrow(genotype_matrix), " samples x ",
               ncol(genotype_matrix), " variants\n"))

    # Check for missing data
    n_missing_geno <- sum(is.na(genotype_matrix))
    if (n_missing_geno > 0) {
        cat(sprintf("WARNING: %d missing genotypes (%.2f%%)\n",
                   n_missing_geno, 100 * n_missing_geno / length(genotype_matrix)))
    }

    return(genotype_matrix)
}


#' Apply MAF bin filtering (like power simulation)
#'
#' @param gene_data data.table with variant information
#' @param maf_bins String specifying MAF bins (e.g., "0.001-0.01,0.01-0.05" or "all")
#' @return Filtered data.table
apply_maf_bins <- function(gene_data, maf_bins) {

    if (maf_bins == "all" || is.null(maf_bins) || maf_bins == "") {
        cat("MAF bins: Using all variants (no additional filtering)\n")
        return(gene_data)
    }

    cat(paste0("Applying MAF bin filtering: ", maf_bins, "\n"))

    # Parse bins
    bins <- unlist(strsplit(maf_bins, ","))

    filtered_list <- list()
    for (bin in bins) {
        bounds <- as.numeric(unlist(strsplit(bin, "-")))
        if (length(bounds) != 2) {
            warning(sprintf("Invalid MAF bin format: %s. Skipping.", bin))
            next
        }

        lower <- bounds[1]
        upper <- bounds[2]

        bin_data <- gene_data[MAF >= lower & MAF < upper]
        filtered_list[[bin]] <- bin_data

        cat(sprintf("  MAF bin [%.4f, %.4f): %d variants\n", lower, upper, nrow(bin_data)))
    }

    # Combine all bins
    maf_filtered_data <- rbindlist(filtered_list)

    # Remove duplicates (in case bins overlap)
    maf_filtered_data <- unique(maf_filtered_data)

    cat(paste0("Total variants after MAF filtering: ", nrow(maf_filtered_data), "\n"))

    return(maf_filtered_data)
}


# Example usage and testing
if (!interactive()) {
    cat("Gene set utilities loaded successfully.\n")
}
