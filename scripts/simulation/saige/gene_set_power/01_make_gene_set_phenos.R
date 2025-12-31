#!/usr/bin/env Rscript
# Simulate quantitative phenotypes from real gene sets
# Uses variants from SAIGE weight files and dominance-encoded VCFs

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
get_script_path <- function() {
    cmd_args <- commandArgs(trailingOnly = FALSE)
    file_arg <- cmd_args[grep("^--file=", cmd_args)]
    if (length(file_arg) > 0) {
        script_path <- sub("^--file=", "", file_arg)
        return(dirname(script_path))
    }
    return(".")
}

script_dir <- get_script_path()
source(file.path(script_dir, "_gene_set_utils.R"))
source(file.path(script_dir, "_encoding_utils.R"))
source(file.path(script_dir, "_effect_size_utils.R"))


#' Compute genotype proportions from genotype vector
#'
#' For dominance-encoded genotypes, we need to map back to original {0,1,2}
#' However, dominance VCFs store the encoded values directly
#' We'll compute proportions from the raw genotypes
#'
#' @param genotypes Numeric vector (can be encoded or raw)
#' @return Named vector with r, h, a proportions
compute_genotype_props_from_vector <- function(genotypes) {
    # Remove NAs
    geno <- genotypes[!is.na(genotypes)]

    # For dominance-encoded data, we need original genotype counts
    # The VCF should store original GT in a separate field, or we derive from dosage
    # For now, assume we have access to raw genotypes {0, 1, 2}

    # Round to nearest integer if needed
    geno_rounded <- round(geno)

    n <- length(geno_rounded)
    r <- sum(geno_rounded == 0) / n
    h <- sum(geno_rounded == 1) / n
    a <- sum(geno_rounded == 2) / n

    return(c(r = r, h = h, a = a))
}


# Note: Genotype loading function is now in _gene_set_utils.R
# load_genotypes_from_plink() handles reading PLINK files


#' Simulate phenotype for one configuration
#'
#' @param raw_genotypes Matrix of raw genotypes [samples x variants] for computing r,h,a
#' @param causal_variants data.table with causal variant info
#' @param architecture Genetic architecture
#' @param h2_total Total heritability
#' @param seed Random seed
#' @return List with phenotype vector and metadata
simulate_one_phenotype <- function(raw_genotypes, causal_variants,
                                   architecture, h2_total, seed = NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }

    n_samples <- nrow(raw_genotypes)
    n_causal <- ncol(raw_genotypes)

    # Per-variant heritability
    h2_per_variant <- h2_total / n_causal

    cat(sprintf("\nSimulating phenotype: arch=%s, h2=%.6f, variants=%d, h2_per_variant=%.6f\n",
                architecture, h2_total, n_causal, h2_per_variant))

    # Compute effect sizes for each variant
    effect_list <- list()
    genotype_props_list <- list()

    for (j in 1:n_causal) {
        # Get raw genotypes for this variant
        geno <- raw_genotypes[, j]

        # Compute genotype proportions
        props <- compute_genotype_props_from_vector(geno)
        genotype_props_list[[j]] <- props

        # Compute effect sizes based on architecture
        effects <- compute_variant_effects(
            architecture,
            r = props['r'],
            h = props['h'],
            a = props['a'],
            h2_target = h2_per_variant
        )

        effect_list[[j]] <- effects

        if (j <= 3 || j == n_causal) {
            cat(sprintf("  Variant %d: r=%.3f, h=%.3f, a=%.3f, MAF=%.4f, beta_A=%.4f, beta_D=%.4f\n",
                       j, props['r'], props['h'], props['a'],
                       causal_variants$MAF[j],
                       effects$beta_A, effects$beta_D))
        } else if (j == 4) {
            cat(sprintf("  ... (%d more variants)\n", n_causal - 4))
        }
    }

    # Compute genetic values
    # For each individual, sum contributions from all causal variants
    genetic_values <- rep(0, n_samples)

    for (j in 1:n_causal) {
        geno <- raw_genotypes[, j]
        props <- genotype_props_list[[j]]
        effects <- effect_list[[j]]

        # Compute additive and dominance encodings
        X_A <- compute_additive_encoding(props['r'], props['h'], props['a'])
        X_D <- compute_dominance_encoding(props['r'], props['h'], props['a'], X_A)

        # Map genotypes to encoded values
        # Handle missing data
        geno_clean <- geno
        geno_clean[is.na(geno_clean)] <- round(props['r'] * 0 + props['h'] * 1 + props['a'] * 2)

        X_A_individual <- X_A[geno_clean + 1]  # +1 for R indexing
        X_D_individual <- X_D[geno_clean + 1]

        # Add contribution
        genetic_values <- genetic_values + effects$beta_A * X_A_individual + effects$beta_D * X_D_individual
    }

    # Verify variance of genetic component
    var_genetic <- var(genetic_values)
    cat(sprintf("Genetic variance: %.6f (target h2: %.6f)\n", var_genetic, h2_total))

    # Add environmental noise
    var_environmental <- (1 - h2_total) / h2_total * var_genetic
    environmental_noise <- rnorm(n_samples, mean = 0, sd = sqrt(var_environmental))

    # Final phenotype
    phenotype <- genetic_values + environmental_noise

    # Standardize to mean 0, sd 1
    phenotype <- scale(phenotype)[,1]

    # Verify heritability
    realized_h2 <- var_genetic / (var_genetic + var_environmental)
    cat(sprintf("Realized h2: %.6f\n", realized_h2))

    # Create metadata
    metadata <- data.table(
        variant_id = causal_variants$variant,
        gene_id = causal_variants$gene,
        MAF = causal_variants$MAF,
        annotation = causal_variants$annotation,
        weight = causal_variants$weight,
        architecture = architecture,
        h2_total = h2_total,
        h2_per_variant = h2_per_variant,
        beta_A = sapply(effect_list, function(x) x$beta_A),
        beta_D = sapply(effect_list, function(x) x$beta_D),
        h2_A = sapply(effect_list, function(x) x$h2_A),
        h2_D = sapply(effect_list, function(x) x$h2_D),
        r = sapply(genotype_props_list, function(x) x['r']),
        h = sapply(genotype_props_list, function(x) x['h']),
        a = sapply(genotype_props_list, function(x) x['a'])
    )

    return(list(
        phenotype = phenotype,
        metadata = metadata
    ))
}


main <- function(opt) {
    # Check required inputs
    required <- c("saige_weights", "bed_file", "frqx_file", "ancestry", "out")
    missing <- required[!required %in% names(opt) | sapply(opt[required], is.null)]
    if (length(missing) > 0) {
        stop(sprintf("Missing required parameters: %s", paste(missing, collapse = ", ")))
    }

    cat("=============================================================\n")
    cat("Gene Set-Based Phenotype Simulation\n")
    cat("=============================================================\n\n")

    # Parse parameters
    architectures <- unlist(strsplit(opt$architectures, ","))
    h2_values <- as.numeric(unlist(strsplit(opt$h2_total, ",")))
    K_values <- as.integer(unlist(strsplit(opt$K_genes, ",")))
    N_values <- as.integer(unlist(strsplit(opt$N_per_gene, ",")))
    allowed_annotations <- unlist(strsplit(opt$annotations, ","))
    n_reps <- opt$n_reps

    cat(paste0("Ancestry: ", opt$ancestry, "\n"))
    cat(paste0("Architectures: ", paste(architectures, collapse = ", "), "\n"))
    cat(paste0("Heritabilities: ", paste(h2_values, collapse = ", "), "\n"))
    cat(paste0("K (causal genes): ", paste(K_values, collapse = ", "), "\n"))
    cat(paste0("N (variants per gene): ", paste(N_values, collapse = ", "), "\n"))
    cat(paste0("Annotations: ", paste(allowed_annotations, collapse = ", "), "\n"))
    cat(paste0("MAF bins: ", opt$maf_bins, "\n"))
    cat(paste0("Replicates: ", n_reps, "\n"))
    cat(paste0("Random seed: ", opt$seed, "\n\n"))

    # Load SAIGE weight file
    gene_data <- parse_saige_weights(opt$saige_weights, allowed_annotations)

    # Load allele counts
    allele_counts <- load_allele_counts(opt$frqx_file)

    # Apply MAF bin filtering
    gene_data <- apply_maf_bins(gene_data, opt$maf_bins)

    # Merge and filter genes
    gene_data_filtered <- filter_genes(
        gene_data,
        allele_counts,
        min_variants = opt$min_variants_per_gene,
        maf_range = c(0.001, 0.05),  # Standard range
        allowed_annotations = allowed_annotations
    )

    # Initialize output
    all_metadata <- list()
    metadata_counter <- 1
    phenotype_list <- list()

    # Loop through K values (number of causal genes)
    for (K in K_values) {
        cat(paste0("\n", paste(rep("=", 60), collapse = ""), "\n"))
        cat(paste0("K = ", K, " causal genes\n"))
        cat(paste0(paste(rep("=", 60), collapse = ""), "\n"))

        # Select causal genes
        causal_gene_seed <- opt$seed + K * 10000
        causal_genes <- select_causal_genes(gene_data_filtered, K, seed = causal_gene_seed)

        cat("\nCausal genes selected:\n")
        for (i in 1:min(10, length(causal_genes))) {
            n_var <- nrow(gene_data_filtered[gene == causal_genes[i]])
            cat(sprintf("  %d. %s (%d variants)\n", i, causal_genes[i], n_var))
        }
        if (length(causal_genes) > 10) {
            cat(sprintf("  ... and %d more\n", length(causal_genes) - 10))
        }

        # Loop through N values (variants per gene)
        for (N_per_gene in N_values) {
            cat(paste0("\n", paste(rep("-", 60), collapse = ""), "\n"))
            cat(paste0("N = ", N_per_gene, " variant(s) per gene\n"))
            cat(paste0(paste(rep("-", 60), collapse = ""), "\n"))

            # Select causal variants
            variant_seed <- causal_gene_seed + N_per_gene * 100
            causal_variants <- select_causal_variants_per_gene(
                gene_data_filtered,
                causal_genes,
                N_per_gene = N_per_gene,
                selection_strategy = opt$variant_selection,
                seed = variant_seed
            )

            total_causal_variants <- nrow(causal_variants)

            # Load genotypes for causal variants from PLINK
            raw_genotypes <- load_genotypes_from_plink(
                opt$bed_file,
                causal_variants$variant
            )

            sample_ids <- rownames(raw_genotypes)

            # Loop through architectures
            for (architecture in architectures) {
                cat(paste0("\nArchitecture: ", architecture, "\n"))

                # Loop through h2 values
                for (h2 in h2_values) {

                    # Loop through replicates
                    for (rep in 1:n_reps) {
                        # Unique seed for each replicate
                        rep_seed <- opt$seed + K * 1000000 + N_per_gene * 100000 +
                                   which(architectures == architecture) * 10000 +
                                   which(h2_values == h2) * 1000 + rep

                        # Simulate phenotype
                        result <- simulate_one_phenotype(
                            raw_genotypes,
                            causal_variants,
                            architecture,
                            h2,
                            seed = rep_seed
                        )

                        # Create phenotype name
                        pheno_name <- sprintf("arch_%s_h2_%s_K_%d_N_%d_rep_%d",
                                             architecture,
                                             gsub("\\.", "", sprintf("%.6f", h2)),
                                             K, N_per_gene, rep)

                        # Store phenotype
                        phenotype_list[[pheno_name]] <- result$phenotype

                        # Store metadata
                        result$metadata$phenotype_name <- pheno_name
                        result$metadata$replicate <- rep
                        result$metadata$seed <- rep_seed
                        result$metadata$K_genes <- K
                        result$metadata$N_per_gene <- N_per_gene
                        result$metadata$total_causal_variants <- total_causal_variants
                        result$metadata$variant_selection_strategy <- opt$variant_selection

                        all_metadata[[metadata_counter]] <- result$metadata
                        metadata_counter <- metadata_counter + 1

                        if (rep %% 10 == 0 || rep == n_reps) {
                            cat(sprintf("      Completed %d/%d replicates for h2=%.4f\n", rep, n_reps, h2))
                        }
                    }
                }
            }
        }
    }

    # Create output data frame
    output_dt <- data.table(eid = sample_ids)
    for (pheno_name in names(phenotype_list)) {
        output_dt[[pheno_name]] <- phenotype_list[[pheno_name]]
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
    cat(paste0("Total phenotypes generated: ", length(phenotype_list), "\n"))
    cat(paste0("Total samples: ", nrow(output_dt), "\n"))
    cat(paste0("Total causal gene-variant combinations: ", nrow(metadata_dt), "\n"))
    cat("=============================================================\n")
}


# Command-line options
option_list <- list(
    make_option("--saige_weights", type = "character", default = NULL,
        help = "Path to SAIGE weight file (gzipped)", metavar = "character"),
    make_option("--bed_file", type = "character", default = NULL,
        help = "Path to PLINK .bed file (without extension)", metavar = "character"),
    make_option("--frqx_file", type = "character", default = NULL,
        help = "Path to allele count file (.frqx)", metavar = "character"),
    make_option("--covars", type = "character", default = NULL,
        help = "Path to covariates file", metavar = "character"),
    make_option("--sample_id_col", type = "character", default = "IID",
        help = "Column name for sample ID in phenotype data", metavar = "character"),
    make_option("--covar_id_col", type = "character", default = "participant_eid",
        help = "Column name for sample ID in covariate data", metavar = "character"),
    make_option("-a", "--ancestry", type = "character", default = NULL,
        help = "Ancestry label (eur, sas, afr, eas)", metavar = "character"),
    make_option("--annotations", type = "character", default = "pLoF,damaging_missense",
        help = "Comma-separated annotation types to include", metavar = "character"),
    make_option("--architectures", type = "character",
        default = "additive,recessive,partially_recessive_0.2",
        help = "Comma-separated genetic architectures", metavar = "character"),
    make_option("--h2_total", type = "character", default = "0.005,0.01",
        help = "Comma-separated total heritabilities", metavar = "character"),
    make_option("--K_genes", type = "character", default = "10",
        help = "Comma-separated K values (number of causal genes)", metavar = "character"),
    make_option("--N_per_gene", type = "character", default = "1",
        help = "Comma-separated N values (causal variants per gene)", metavar = "character"),
    make_option("--variant_selection", type = "character", default = "random",
        help = "Variant selection strategy: random/lowest_maf/highest_maf/weighted", metavar = "character"),
    make_option("--maf_bins", type = "character", default = "all",
        help = "MAF bins for variant selection (e.g., '0.001-0.01,0.01-0.05' or 'all')", metavar = "character"),
    make_option("--min_variants_per_gene", type = "integer", default = 1,
        help = "Minimum variants per gene", metavar = "integer"),
    make_option("-r", "--n_reps", type = "integer", default = 100,
        help = "Number of replicates per configuration", metavar = "integer"),
    make_option("-s", "--seed", type = "integer", default = 42,
        help = "Random seed", metavar = "integer"),
    make_option("-o", "--out", type = "character", default = NULL,
        help = "Output file prefix", metavar = "character")
)

opt_parser <- OptionParser(add_help_option = TRUE, option_list = option_list)
opt <- parse_args(opt_parser)

main(opt)
