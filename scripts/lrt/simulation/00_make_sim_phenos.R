#!/usr/bin/env Rscript

options(scipen = 999)

# Required packages
packages <- c('data.table', 'dplyr', 'optparse', 'tidyverse', 'ggplot2', 'stringr')

for (p in packages) {
    if (!require(p, character.only = TRUE)) {
        install.packages(p, repos = "http://cran.us.r-project.org")
    }
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

    # Generate uncorrelated random normal phenotypes
    set.seed(opt$seed)
    phenos <- matrix(rnorm(n * n_phenos_sim), nrow=n_phenos_sim, ncol=n)

    for (j in 1:n_h2)
    {
        h2 <- h2_vec[j]
        cat(paste("Current h2:", h2, "\n"))

        # Subset to the n_phenos unused raw phenotypes for this round
        start <- (j-1) * n_phenos + 1
        end <- j * n_phenos
        use_phenos <- t(phenos[start:end,])

        # Generate genetic component and environmental noise
        genetic <- matrix(rnorm(length(use_phenos)), nrow=nrow(use_phenos))
        noise <- matrix(rnorm(length(use_phenos)), nrow=nrow(use_phenos))

        # Create n_phenos phenotypes with heritability h2
        use_phenos <- genetic * sqrt(h2) + noise * sqrt(1 - h2)

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
        if (opt$trait_type %in% c('binary', 'both') | !is.null(n_cases)) {
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

    colnames(output_dt)[colnames(output_dt)=="IID"] <- 'eid'
    print(output_dt %>% summary)
    fwrite(data.table(output_dt), file=paste0(opt$out, ".tsv.gz"), sep="\t")
}

# Add arguments
option_list = list(
    make_option(c("-s", "--samples"), type="character", default=NULL,
        help="path to the sample ID file", metavar="character"),
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
