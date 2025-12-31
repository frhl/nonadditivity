#!/usr/bin/env Rscript

options(scipen = 999)

# Required packages
packages <- c('data.table', 'dplyr', 'optparse', 'tidyverse', 'ggplot2', 'stringr', 'Matrix')

for (p in packages) {
    if (!require(p, character.only = TRUE, quietly = TRUE)) {
        cat(paste0("Installing package: ", p, "\n"))
        install.packages(p, repos = "http://cran.us.r-project.org")
        library(p, character.only = TRUE)
    }
}

#' Read sparse GRM in Matrix Market format
#'
#' @param grm_file Path to sparse GRM file (.mtx format)
#' @return Sparse symmetric matrix
read_sparse_grm <- function(grm_file) {
    cat(paste0("Reading sparse GRM from: ", grm_file, "\n"))

    # Read the Matrix Market file
    grm <- readMM(grm_file)

    # Make symmetric (GRM files often store only lower triangle)
    grm <- as(grm, "symmetricMatrix")

    cat(paste0("  GRM dimensions: ", nrow(grm), " x ", ncol(grm), "\n"))
    cat(paste0("  Non-zero elements: ", length(grm@x), "\n"))

    return(grm)
}

#' Generate GRM-correlated random noise
#'
#' Uses Cholesky decomposition to generate correlated noise from GRM
#' noise ~ MVN(0, sigma^2 * GRM)
#'
#' @param grm Sparse GRM matrix
#' @param sigma Standard deviation of noise
#' @param seed Random seed
#' @return Vector of correlated noise
generate_grm_correlated_noise <- function(grm, sigma, seed = NULL) {
    if (!is.null(seed)) {
        set.seed(seed)
    }

    n <- nrow(grm)

    # Generate independent standard normals
    z <- rnorm(n)

    # Add small diagonal regularization for numerical stability
    grm_reg <- grm + Diagonal(n, 1e-4)

    # Compute Cholesky decomposition: GRM = L * L^T
    cat("  Computing Cholesky decomposition of GRM...\n")
    L <- tryCatch({
        chol(grm_reg, pivot = FALSE)
    }, error = function(e) {
        cat("  Warning: Standard Cholesky failed, using pivoted Cholesky\n")
        chol(grm_reg, pivot = TRUE)
    })

    # Generate correlated noise: sigma * L^T * z
    # where z ~ N(0, I)
    # This gives noise ~ N(0, sigma^2 * GRM)
    noise <- sigma * as.vector(t(L) %*% z)

    return(noise)
}

main <- function(opt)
{
    if (is.null(opt$samples)) {
        print_help(opt_parser)
        stop("At least one argument must be supplied (samples file).",
            call.=FALSE)
    }

    # Read sample IDs
    n <- nrow(fread(opt$samples))
    n_phenos_sim <- opt$n
    h2_vec <- as.numeric(unlist(strsplit(opt$heritabilities, split=',')))
    n_h2 <- length(h2_vec)
    binary_info <- as.numeric(unlist(strsplit(opt$prevalence, split=',')))
    output_dt <- fread(opt$samples, col.names=c("userID"))

    # Read sparse GRM if provided
    grm <- NULL
    if (!is.null(opt$grm)) {
        grm <- read_sparse_grm(opt$grm)
        if (nrow(grm) != n) {
            stop(paste0("GRM dimensions (", nrow(grm), ") do not match number of samples (", n, ")"))
        }
        cat(paste0("Using GRM-correlated noise for phenotype simulation\n\n"))
    } else {
        cat(paste0("No GRM provided - using uncorrelated noise\n\n"))
    }

    if (mean(binary_info) > 1) {
        n_cases <- binary_info
        prevalences <- n_cases / n
    } else {
        prevalences <- binary_info
        n_cases <- prevalences * n
    }

    n_prevalences <- length(prevalences)
    qnorms <- qnorm(prevalences)
    n_draws <- opt$n_reps
    n_phenos <- ((n_prevalences + 1) * n_draws)
    n_phenos_sim_required <- n_phenos * n_h2

    if ((!opt$sanity_check) |
        (opt$sanity_check & (n_phenos_sim_required > n_phenos_sim))) {
        n_phenos_sim <- n_phenos_sim_required
    }

    cat(paste0(
        'N samples: ', n, "\n",
        'Number of phenotypes to simulate: ', n_phenos_sim, "\n",
        'Number of trait heritabilities to simulate: ', n_h2, "\n",
        'Number of prevalences: ', n_prevalences, "\n")
    )

    for (j in 1:n_h2)
    {
        h2 <- h2_vec[j]
        cat(paste("Current h2:", h2, "\n"))

        # Generate genetic component for null phenotypes (uncorrelated across samples)
        # For null simulations, genetic variance = 0, so this contributes nothing when h2 = 0
        set.seed(opt$seed + j * 1000)
        genetic <- matrix(rnorm(n * n_phenos), nrow=n, ncol=n_phenos)

        # Generate environmental noise (GRM-correlated if GRM provided)
        if (!is.null(grm)) {
            cat("  Generating GRM-correlated noise...\n")
            # Generate n_phenos vectors of GRM-correlated noise
            noise <- sapply(1:n_phenos, function(pheno_idx) {
                # Use different seed for each phenotype
                pheno_seed <- opt$seed + j * 10000 + pheno_idx
                generate_grm_correlated_noise(grm, sigma = sqrt(1 - h2), seed = pheno_seed)
            })
        } else {
            # Uncorrelated noise
            noise <- matrix(rnorm(length(use_phenos)), nrow=nrow(use_phenos))
            noise <- noise * sqrt(1 - h2)
        }

        # Create n_phenos phenotypes with heritability h2
        # For GRM-correlated noise: genetic * sqrt(h2) + noise (already scaled)
        # For uncorrelated noise: genetic * sqrt(h2) + noise * sqrt(1 - h2)
        use_phenos <- genetic * sqrt(h2) + noise

        renaming_binary <- function(x) {
            if (n_draws == 1) {
                paste0('p_', h2, '_', binary_info[i], '_1')
            } else {
                paste0('p_', h2, '_', binary_info[i], '_', gsub('X', '', x))
            }
        }

        renaming_cts <- function(x) {
            if (n_draws == 1) {
                paste0('p_', h2, '_continuous_1')
            } else {
                paste0('p_', h2, '_continuous_', gsub('X', '', x))
            }
        }

        # Create n_prevalence binary phenotypes with the current heritability
        if (opt$trait_type %in% c('binary', 'both')) {
            for (i in 1:n_prevalences) {
                start <- (i - 1) * n_draws + 1
                end <- i * n_draws
                out <- data.frame(use_phenos[,start:end] < qnorms[i])
                output_dt <- output_dt %>%
                    bind_cols(out %>% rename_with(renaming_binary)) %>%
                    mutate_all(as.numeric)
            }
        }

        # Create continuous phenotypes
        if (opt$trait_type %in% c('continuous', 'both')) {
            continuous_df <- data.frame(
                use_phenos[, (n_prevalences*n_draws + 1):n_phenos])
            output_dt <- output_dt %>%
                bind_cols(continuous_df %>% rename_with(renaming_cts))
        }
    }

    # Merge with covariates if provided
    if (!is.null(opt$covars)) {
        cat(paste0("Reading covariates from: ", opt$covars, "\n"))
        covars_dt <- fread(opt$covars)

        # Rename ID columns for merging
        sample_id_col <- opt$sample_id_col
        covar_id_col <- opt$covar_id_col

        # Ensure the sample ID column exists in output_dt
        if (!(sample_id_col %in% colnames(output_dt))) {
            if ("userID" %in% colnames(output_dt)) {
                colnames(output_dt)[colnames(output_dt) == "userID"] <- sample_id_col
            }
        }

        cat(paste0("Merging on phenotype column '", sample_id_col,
                   "' and covariate column '", covar_id_col, "'\n"))
        cat(paste0("Phenotypes before merge: ", nrow(output_dt), " samples\n"))
        cat(paste0("Covariates: ", nrow(covars_dt), " samples\n"))

        # Perform merge
        output_dt <- merge(output_dt, covars_dt,
                          by.x = sample_id_col,
                          by.y = covar_id_col,
                          all.x = TRUE)

        cat(paste0("After merge: ", nrow(output_dt), " samples\n"))
    }

    # Ensure final ID column is named 'eid'
    if ("IID" %in% colnames(output_dt)) {
        colnames(output_dt)[colnames(output_dt) == "IID"] <- "eid"
    } else if ("userID" %in% colnames(output_dt)) {
        colnames(output_dt)[colnames(output_dt) == "userID"] <- "eid"
    }

    print(output_dt %>% summary)
    fwrite(data.table(output_dt), file=paste0(opt$out, ".tsv.gz"), sep="\t")
}

# Add arguments
option_list = list(
    make_option(c("-s", "--samples"), type="character", default=NULL,
        help="path to the sample ID file", metavar="character"),
    make_option(c("-g", "--grm"), type="character", default=NULL,
        help="path to the sparse GRM file (.mtx format)", metavar="character"),
    make_option(c("-c", "--covars"), type="character", default=NULL,
        help="path to the covariates file", metavar="character"),
    make_option(c("--sample_id_col"), type="character", default="IID",
        help="column name for sample ID in phenotype data", metavar="character"),
    make_option(c("--covar_id_col"), type="character", default="IID",
        help="column name for sample ID in covariate data", metavar="character"),
    make_option(c("--heritabilities"), type="character", default='0.5',
        help="heritability of the phenotypes (a list of comma split values)",
        metavar="character"),
    make_option(c("-n", "--n"), type="integer", default=1000,
        help="number of raw phenotypes needed", metavar="integer"),
    make_option(c("-p", "--prevalence"), type="character",
        default='0.001,0.01,0.1,0.2,0.5',
        help=paste0("number of cases (a list of comma split integers) or ",
            "prevalences (a list of comma split value between 0 and 1)"),
            metavar="character"),
    make_option(c("-r", "--n_reps"), type="integer", default=10,
        help="number of draws to run", metavar="integer"),
    make_option(c("-i", "--seed"), type="double", default=663,
        help="set seed for random sampling", metavar="double"),
    make_option(c("-t", "--trait_type"), type="character", default='both',
        help=paste0("the type of phenotypes to simulate, (choose from ",
            "`continuous`, `binary` or `both`)"),
            metavar="character"),
    make_option(c("-o", "--out"), type="character", default="random_phenos",
        help="output file name [default= %default]", metavar="character"),
    make_option(c("--sanity_check"), type="logical", default=FALSE,
        action="store_true",
        help="sanity check option (currently disabled)", metavar="logical")
)
opt_parser <- OptionParser(add_help_option=TRUE, option_list=option_list)
opt <- parse_args(opt_parser)

main(opt)
