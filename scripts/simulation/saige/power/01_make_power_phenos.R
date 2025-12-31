#!/usr/bin/env Rscript
# Simulate quantitative phenotypes under diverse genetic architectures
# for power analysis of SAIGE genetic association tests

options(scipen = 999)

# Required packages
packages <- c('data.table', 'dplyr', 'optparse', 'BEDMatrix')

for (p in packages) {
    if (!require(p, character.only = TRUE, quietly = TRUE)) {
        cat(paste0("Installing package: ", p, "\n"))
        install.packages(p, repos = "http://cran.us.r-project.org")
        library(p, character.only = TRUE)
    }
}

# Source utility functions
# Get script directory (works with Rscript)
get_script_path <- function() {
    cmd_args <- commandArgs(trailingOnly = FALSE)
    file_arg <- cmd_args[grep("^--file=", cmd_args)]
    if (length(file_arg) > 0) {
        script_path <- sub("^--file=", "", file_arg)
        return(dirname(script_path))
    }
    # Fallback: assume current directory
    return(".")
}

script_dir <- get_script_path()
source(file.path(script_dir, "_encoding_utils.R"))
source(file.path(script_dir, "_effect_size_utils.R"))


#' Read genotypes from PLINK bed file
#'
#' @param bed_file Path to .bed file (without extension)
#' @return List with genotype_matrix, sample_ids, variant_ids
read_plink_genotypes <- function(bed_file) {
    cat(paste0("Reading PLINK genotypes from: ", bed_file, "\n"))

    # Read using BEDMatrix (memory-efficient)
    bed <- BEDMatrix(paste0(bed_file, ".bed"))

    # Read .fam file for sample IDs
    fam <- fread(paste0(bed_file, ".fam"), header = FALSE)
    colnames(fam) <- c("FID", "IID", "father", "mother", "sex", "pheno")
    sample_ids <- fam$IID

    # Read .bim file for variant IDs
    bim <- fread(paste0(bed_file, ".bim"), header = FALSE)
    colnames(bim) <- c("chr", "ID", "cm", "pos", "A1", "A2")
    variant_ids <- bim$ID

    cat(paste0("  Samples: ", length(sample_ids), "\n"))
    cat(paste0("  Variants: ", length(variant_ids), "\n"))

    return(list(
        genotypes = bed,  # BEDMatrix object (can be indexed like a matrix)
        sample_ids = sample_ids,
        variant_ids = variant_ids,
        variant_info = bim
    ))
}


#' Compute variant statistics (MAF, carrier counts)
#'
#' Pre-computes MAF and carrier counts for all variants (called once)
#'
#' @param genotype_data List from read_plink_genotypes()
#' @return List with variant_mafs, variant_carriers, variant_hets
compute_variant_statistics <- function(genotype_data) {
    n_variants <- length(genotype_data$variant_ids)

    cat("\nComputing MAF and carrier counts for all variants...\n")
    variant_info <- lapply(1:n_variants, function(idx) {
        geno <- as.numeric(genotype_data$genotypes[, idx])
        geno <- geno[!is.na(geno)]
        props <- compute_genotype_proportions(geno)
        n_hom_alt <- sum(geno == 2)
        n_het <- sum(geno == 1)
        maf <- calculate_maf(props['r'], props['h'], props['a'])
        return(list(maf = maf, n_hom_alt = n_hom_alt, n_het = n_het))
    })

    variant_mafs <- sapply(variant_info, function(x) x$maf)
    variant_carriers <- sapply(variant_info, function(x) x$n_hom_alt)
    variant_hets <- sapply(variant_info, function(x) x$n_het)

    return(list(
        variant_mafs = variant_mafs,
        variant_carriers = variant_carriers,
        variant_hets = variant_hets
    ))
}


#' Select causal variants from genotype data with MAF bin filtering
#'
#' Randomly selects M variants, with optional MAF bin filtering
#' Filters to variants with ≥5 homozygous alternate AND ≥5 heterozygotes
#' (heterozygotes required for stable dominance encoding)
#'
#' @param genotype_data List from read_plink_genotypes()
#' @param variant_stats Pre-computed variant statistics from compute_variant_statistics()
#' @param M Number of causal variants to select
#' @param maf_bins Character string with MAF bins (e.g., "0.001-0.01,0.01-0.05")
#' @param seed Random seed
#' @return List with indices, IDs, and genotype proportions
select_causal_variants <- function(genotype_data, variant_stats, M, maf_bins = NULL, seed = NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }

    variant_mafs <- variant_stats$variant_mafs
    variant_carriers <- variant_stats$variant_carriers
    variant_hets <- variant_stats$variant_hets

    # Filter to variants with ≥5 biallelic carriers AND ≥5 heterozygotes
    # (heterozygotes needed for stable dominance encoding)
    carrier_filter <- variant_carriers >= 5 & variant_hets >= 5
    n_pass_carrier <- sum(carrier_filter)
    n_fail_carrier <- sum(!carrier_filter)
    n_fail_hom_only <- sum(variant_carriers >= 5 & variant_hets < 5)

    cat(sprintf("Carrier count filtering (≥5 homozygous alternate AND ≥5 heterozygotes):\n"))
    cat(sprintf("  Pass: %d variants\n", n_pass_carrier))
    cat(sprintf("  Fail: %d variants\n", n_fail_carrier))
    if (n_fail_hom_only > 0) {
        cat(sprintf("    (including %d with ≥5 hom alt but <5 het)\n", n_fail_hom_only))
    }

    if (n_pass_carrier == 0) {
        stop("No variants have ≥5 homozygous alternate carriers AND ≥5 heterozygotes")
    }

    # Get indices of variants passing carrier filter
    eligible_indices <- which(carrier_filter)

    # Filter by MAF bins if specified (applied after carrier filtering)
    if (!is.null(maf_bins) && maf_bins != "" && maf_bins != "all") {
        cat(paste0("\nApplying MAF bin filtering: ", maf_bins, "\n"))

        # Parse MAF bins (e.g., "0.001-0.01,0.01-0.05")
        bins <- unlist(strsplit(maf_bins, split = ','))
        maf_filtered_indices <- c()

        for (bin in bins) {
            bounds <- as.numeric(unlist(strsplit(bin, split = '-')))
            if (length(bounds) != 2) {
                warning(sprintf("Invalid MAF bin format: %s. Skipping.", bin))
                next
            }

            lower <- bounds[1]
            upper <- bounds[2]

            # Apply MAF filter only to variants already passing carrier filter
            bin_indices <- eligible_indices[variant_mafs[eligible_indices] >= lower &
                                           variant_mafs[eligible_indices] < upper]
            maf_filtered_indices <- c(maf_filtered_indices, bin_indices)

            cat(sprintf("  MAF bin [%.4f, %.4f): %d variants (with ≥5 het AND ≥5 hom alt)\n",
                       lower, upper, length(bin_indices)))
        }

        eligible_indices <- unique(maf_filtered_indices)
        cat(sprintf("Total eligible variants after MAF filtering: %d\n", length(eligible_indices)))

        if (length(eligible_indices) < M) {
            stop(sprintf("Only %d variants pass both carrier and MAF filters, but %d requested",
                        length(eligible_indices), M))
        }

        # Select from eligible variants only
        selected_indices <- sample(eligible_indices, M, replace = FALSE)

    } else {
        # No MAF filtering - select from carrier-filtered variants
        cat("\nNo MAF filtering specified. Selecting from carrier-filtered variants.\n")

        if (M > length(eligible_indices)) {
            stop(sprintf("Requested %d causal variants, but only %d pass carrier filter (≥5)",
                        M, length(eligible_indices)))
        }

        selected_indices <- sample(eligible_indices, M, replace = FALSE)
    }

    selected_ids <- genotype_data$variant_ids[selected_indices]
    selected_mafs <- variant_mafs[selected_indices]

    cat(paste0("\nSelected ", M, " causal variants:\n"))
    cat(sprintf("MAF range: [%.4f, %.4f]\n", min(selected_mafs), max(selected_mafs)))
    cat(sprintf("Mean MAF: %.4f\n", mean(selected_mafs)))

    # Compute genotype proportions for each selected variant
    genotype_props <- lapply(seq_along(selected_indices), function(i) {
        idx <- selected_indices[i]
        geno <- as.numeric(genotype_data$genotypes[, idx])
        geno <- geno[!is.na(geno)]  # Remove missing

        props <- compute_genotype_proportions(geno)

        # Count carriers
        n_hom_alt <- sum(geno == 2)
        n_het <- sum(geno == 1)

        if (i <= 10 || i == length(selected_indices)) {  # Show first 10 and last
            cat(sprintf("  %s: MAF=%.4f, r=%.3f, h=%.3f, a=%.3f (n_het=%d, n_hom_alt=%d)\n",
                       selected_ids[i],
                       selected_mafs[i],
                       props['r'], props['h'], props['a'], n_het, n_hom_alt))
        } else if (i == 11) {
            cat(sprintf("  ... (%d more variants)\n", length(selected_indices) - 11))
        }

        return(props)
    })

    return(list(
        indices = selected_indices,
        ids = selected_ids,
        genotype_props = genotype_props,
        mafs = selected_mafs
    ))
}


#' Simulate phenotypes for one configuration
#'
#' @param genotype_data Genotype data object
#' @param causal_variants Output from select_causal_variants()
#' @param architecture Genetic architecture
#' @param h2_total Total heritability
#' @param M Number of causal variants
#' @param seed Random seed
#' @return List with phenotypes and metadata
simulate_one_phenotype <- function(genotype_data, causal_variants,
                                   architecture, h2_total, M, seed = NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }

    # Per-variant target heritability
    h2_per_variant <- h2_total / M

    # Compute effect sizes for each causal variant
    effect_list <- lapply(1:M, function(j) {
        props <- causal_variants$genotype_props[[j]]
        compute_variant_effects(architecture, props['r'], props['h'], props['a'], h2_per_variant)
    })

    # Extract genotypes for causal variants
    n_samples <- nrow(genotype_data$genotypes)
    genotype_matrix <- matrix(0, nrow = n_samples, ncol = M)

    for (j in 1:M) {
        idx <- causal_variants$indices[j]
        genotype_matrix[, j] <- as.numeric(genotype_data$genotypes[, idx])
    }

    # Handle missing data (set to most common genotype)
    for (j in 1:M) {
        missing_idx <- is.na(genotype_matrix[, j])
        if (sum(missing_idx) > 0) {
            mode_geno <- as.numeric(names(which.max(table(genotype_matrix[, j]))))
            genotype_matrix[missing_idx, j] <- mode_geno
        }
    }

    # Simulate phenotypes
    y <- simulate_phenotypes(genotype_matrix, effect_list, h2_total, seed = seed)

    # Collect metadata
    metadata <- data.table(
        variant_id = causal_variants$ids,
        variant_index = causal_variants$indices,
        maf = causal_variants$mafs,
        architecture = architecture,
        h2_total = h2_total,
        M = M,
        h2_per_variant = h2_per_variant,
        beta_A = sapply(effect_list, function(x) x$beta_A),
        beta_D = sapply(effect_list, function(x) x$beta_D),
        h2_A = sapply(effect_list, function(x) x$h2_A),
        h2_D = sapply(effect_list, function(x) x$h2_D),
        r = sapply(effect_list, function(x) x$r),
        h = sapply(effect_list, function(x) x$h),
        a = sapply(effect_list, function(x) x$a)
    )

    return(list(
        phenotype = y,
        metadata = metadata
    ))
}


main <- function(opt) {
    # Check required inputs
    if (is.null(opt$bed_file) || is.null(opt$ancestry) || is.null(opt$out)) {
        print_help(opt_parser)
        stop("Required parameters: --bed_file, --ancestry, --out", call. = FALSE)
    }

    cat("=============================================================\n")
    cat("Simulating phenotypes for power analysis\n")
    cat("=============================================================\n\n")

    # Parse parameters
    architectures <- unlist(strsplit(opt$architectures, split = ','))
    h2_values <- as.numeric(unlist(strsplit(opt$h2_total, split = ',')))
    M_values <- as.integer(unlist(strsplit(opt$polygenicity, split = ',')))
    n_reps <- opt$n_reps

    cat(paste0("Ancestry: ", opt$ancestry, "\n"))
    cat(paste0("Architectures: ", paste(architectures, collapse = ", "), "\n"))
    cat(paste0("Heritabilities: ", paste(h2_values, collapse = ", "), "\n"))
    cat(paste0("Polygenicities (M): ", paste(M_values, collapse = ", "), "\n"))
    cat(paste0("Replicates: ", n_reps, "\n"))
    cat(paste0("Random seed: ", opt$seed, "\n\n"))

    # Read genotype data
    genotype_data <- read_plink_genotypes(opt$bed_file)

    # Compute variant statistics ONCE (MAF, carrier counts)
    variant_stats <- compute_variant_statistics(genotype_data)

    # Initialize output data frame with sample IDs
    output_dt <- data.table(eid = genotype_data$sample_ids)

    # Store all metadata
    all_metadata <- list()
    metadata_counter <- 1

    # Loop through all combinations
    for (M in M_values) {
        cat(paste0("\n==================== M = ", M, " ====================\n"))

        # Select causal variants once per M value
        set.seed(opt$seed + M)  # Different seed for each M
        causal_variants <- select_causal_variants(
            genotype_data, variant_stats, M,
            maf_bins = opt$maf_bins,
            seed = opt$seed + M
        )

        for (architecture in architectures) {
            cat(paste0("\n--- Architecture: ", architecture, " ---\n"))

            for (h2 in h2_values) {
                cat(paste0("  h² = ", h2, "\n"))

                for (rep in 1:n_reps) {
                    # Unique seed for each replicate
                    rep_seed <- opt$seed + M * 1000 + which(architectures == architecture) * 100 +
                               which(h2_values == h2) * 10 + rep

                    # Simulate phenotype
                    result <- simulate_one_phenotype(
                        genotype_data, causal_variants,
                        architecture, h2, M, seed = rep_seed
                    )

                    # Create phenotype name
                    pheno_name <- sprintf("arch_%s_h2_%s_M_%d_rep_%d",
                                         architecture,
                                         gsub("\\.", "", sprintf("%.6f", h2)),
                                         M, rep)

                    # Add phenotype to output
                    output_dt[[pheno_name]] <- result$phenotype

                    # Store metadata
                    result$metadata$phenotype_name <- pheno_name
                    result$metadata$replicate <- rep
                    result$metadata$seed <- rep_seed
                    all_metadata[[metadata_counter]] <- result$metadata
                    metadata_counter <- metadata_counter + 1

                    if (rep %% 10 == 0 || rep == n_reps) {
                        cat(sprintf("    Completed %d/%d replicates\n", rep, n_reps))
                    }
                }
            }
        }
    }

    # Merge with covariates if provided
    if (!is.null(opt$covars)) {
        cat(paste0("\nMerging with covariates from: ", opt$covars, "\n"))
        covars_dt <- fread(opt$covars)

        sample_id_col <- opt$sample_id_col
        covar_id_col <- opt$covar_id_col

        # Rename ID column if needed
        if (!(sample_id_col %in% colnames(output_dt))) {
            if ("eid" %in% colnames(output_dt)) {
                setnames(output_dt, "eid", sample_id_col)
            }
        }

        cat(paste0("  Phenotypes before merge: ", nrow(output_dt), " samples\n"))
        cat(paste0("  Covariates: ", nrow(covars_dt), " samples\n"))

        output_dt <- merge(output_dt, covars_dt,
                          by.x = sample_id_col,
                          by.y = covar_id_col,
                          all.x = TRUE)

        cat(paste0("  After merge: ", nrow(output_dt), " samples\n"))
    }

    # Ensure final ID column is named 'eid'
    if ("IID" %in% colnames(output_dt)) {
        setnames(output_dt, "IID", "eid")
    }

    # Save phenotypes
    pheno_file <- paste0(opt$out, ".tsv.gz")
    cat(paste0("\nSaving phenotypes to: ", pheno_file, "\n"))
    fwrite(output_dt, file = pheno_file, sep = "\t")

    # Save metadata
    metadata_dt <- rbindlist(all_metadata, fill = TRUE)
    metadata_file <- paste0(opt$out, ".metadata.tsv.gz")
    cat(paste0("Saving metadata to: ", metadata_file, "\n"))
    fwrite(metadata_dt, file = metadata_file, sep = "\t")

    cat("\n=============================================================\n")
    cat("Simulation completed successfully!\n")
    cat(paste0("Total phenotypes generated: ", ncol(output_dt) - 1, "\n"))  # -1 for eid column
    cat(paste0("Total samples: ", nrow(output_dt), "\n"))
    cat("=============================================================\n")
}


# Command-line options
option_list <- list(
    make_option(c("-b", "--bed_file"), type = "character", default = NULL,
        help = "Path to PLINK .bed file (without extension)", metavar = "character"),
    make_option(c("-c", "--covars"), type = "character", default = NULL,
        help = "Path to covariates file", metavar = "character"),
    make_option(c("--sample_id_col"), type = "character", default = "IID",
        help = "Column name for sample ID in phenotype data", metavar = "character"),
    make_option(c("--covar_id_col"), type = "character", default = "IID",
        help = "Column name for sample ID in covariate data", metavar = "character"),
    make_option(c("-a", "--ancestry"), type = "character", default = NULL,
        help = "Ancestry label (eur, sas, afr, eas)", metavar = "character"),
    make_option(c("--architectures"), type = "character",
        default = "additive,partially_recessive,recessive,overdominant",
        help = "Comma-separated genetic architectures", metavar = "character"),
    make_option(c("--h2_total"), type = "character", default = "0.01,0.05,0.1,0.2",
        help = "Comma-separated total heritabilities", metavar = "character"),
    make_option(c("--polygenicity"), type = "character", default = "1,10,50",
        help = "Comma-separated M values (number of causal variants)", metavar = "character"),
    make_option(c("--maf_bins"), type = "character", default = "all",
        help = "MAF bins for variant selection (e.g., '0.001-0.01,0.01-0.05' or 'all')", metavar = "character"),
    make_option(c("-r", "--n_reps"), type = "integer", default = 100,
        help = "Number of replicates per configuration", metavar = "integer"),
    make_option(c("-s", "--seed"), type = "integer", default = 42,
        help = "Random seed", metavar = "integer"),
    make_option(c("-o", "--out"), type = "character", default = NULL,
        help = "Output file prefix", metavar = "character")
)

opt_parser <- OptionParser(add_help_option = TRUE, option_list = option_list)
opt <- parse_args(opt_parser)

main(opt)
