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

# Check if file exists on DNAnexus
dx_file_exists <- function(file_path) {
    cmd <- paste0("dx ls '", file_path, "' 2>/dev/null | wc -l")
    result <- system(cmd, intern = TRUE)
    count <- as.integer(result)
    return(count > 0)
}

create_step1_sim_cmd <- function(
    phenotype, trait_type, covariates, categorical_covariates,
    ancestry, sim_phenos_file, grm_prefix, step0_location, sample_ids_file,
    instance_type="mem3_ssd1_v2_x4", priority="low",
    destination)
{
    name <- paste0("step1_", phenotype, "_", ancestry)
    output_prefix <- paste0("out/", phenotype, "_", ancestry)

    grm_file <- paste0(step0_location, grm_prefix, ".sparseGRM.mtx")
    grm_samples <- paste0(step0_location, grm_prefix, ".sparseGRM.mtx.sampleIDs.txt")

    cmd <- paste0("dx run saige-universal-step-1",
        " -i output_prefix=", output_prefix,
        " -i sample_ids=", sample_ids_file,
        " -i genotype_bed=", step0_location, "UKB.array.", ancestry, ".plink_for_var_ratio.bed",
        " -i genotype_bim=", step0_location, "UKB.array.", ancestry, ".plink_for_var_ratio.bim",
        " -i genotype_fam=", step0_location, "UKB.array.", ancestry, ".plink_for_var_ratio.fam",
        " -i pheno_list=", sim_phenos_file,
        " -i pheno=", phenotype,
        " -i GRM=", grm_file,
        " -i GRM_samples=", grm_samples,
        " -i covariates=", covariates,
        " -i categorical_covariates=", categorical_covariates,
        " -i trait_type=", trait_type,
        " --instance-type ", instance_type,
        " --priority ", priority,
        " --destination .", destination,
        " -y",
        " --name ", name)
    return(cmd)
}

main <- function(opt)
{
    if (is.null(opt$pheno_file)) {
        print_help(opt_parser)
        stop("Phenotype file must be supplied.", call.=FALSE)
    }

    # Generate phenotype names based on parameters
    cat(paste0("Generating phenotype names from parameters\n"))
    cat(paste0("Heritabilities: ", opt$heritabilities, "\n"))
    cat(paste0("N reps: ", opt$n_reps, "\n"))

    h2_vec <- as.numeric(unlist(strsplit(opt$heritabilities, split=',')))

    cts_phenos <- c()
    for (h2 in h2_vec) {
        for (rep in 1:opt$n_reps) {
            pheno_name <- paste0("p_", h2, "_continuous_", rep)
            cts_phenos <- c(cts_phenos, pheno_name)
        }
    }

    binary_phenos <- c()  # No binary phenotypes for now

    cat(paste0("Generated ", length(cts_phenos), " continuous phenotypes\n"))
    cat(paste0("Phenotypes: ", paste(cts_phenos, collapse=", "), "\n"))

    # Set up covariates
    covariates <- opt$covariates
    categorical_covariates <- opt$categorical_covariates

    cat(paste0("Covariates: ", covariates, "\n"))
    cat(paste0("Categorical covariates: ", categorical_covariates, "\n"))

    # GRM prefix (without file extension)
    grm_prefix <- paste0("UKB.array.", opt$ancestry,
                        "_relatednessCutoff_0.05_5000_randomMarkersUsed")

    # Sample IDs file
    sample_ids_file <- paste0(opt$sample_dir, "/UKB.wes.qced.", opt$ancestry, ".samples")

    # Submit jobs for continuous phenotypes
    if (length(cts_phenos) > 0 && opt$trait_type %in% c("quantitative", "both")) {
        cat(paste0("\nChecking ", length(cts_phenos), " continuous phenotype jobs\n"))
        n_submitted <- 0
        n_skipped <- 0

        for (pheno in cts_phenos) {
            # Check if output files already exist
            output_rda <- paste0(opt$destination, "/", pheno, "_", opt$ancestry, ".rda")
            output_vr <- paste0(opt$destination, "/", pheno, "_", opt$ancestry, ".varianceRatio.txt")

            if (dx_file_exists(output_rda) && dx_file_exists(output_vr)) {
                cat(paste0("Output already exists for ", pheno, ". Skipping.\n"))
                n_skipped <- n_skipped + 1
            } else {
                cmd <- create_step1_sim_cmd(
                    phenotype=pheno,
                    trait_type="quantitative",
                    covariates=covariates,
                    categorical_covariates=categorical_covariates,
                    ancestry=opt$ancestry,
                    sim_phenos_file=opt$pheno_file,
                    grm_prefix=grm_prefix,
                    step0_location=opt$step0_location,
                    sample_ids_file=sample_ids_file,
                    instance_type=opt$instance_type,
                    priority=opt$priority,
                    destination=opt$destination
                )
                cat(paste0("\nSubmitting: ", pheno, "\n"))
                cat(paste0(cmd, "\n"))
                if (!opt$dry_run) {
                    system(cmd)
                }
                n_submitted <- n_submitted + 1
            }
        }
        cat(paste0("\nSummary: ", n_submitted, " jobs submitted, ", n_skipped, " skipped\n"))
    }

    # Submit jobs for binary phenotypes
    if (length(binary_phenos) > 0 && opt$trait_type %in% c("binary", "both")) {
        cat(paste0("\nChecking ", length(binary_phenos), " binary phenotype jobs\n"))
        n_submitted <- 0
        n_skipped <- 0

        for (pheno in binary_phenos) {
            # Check if output files already exist
            output_rda <- paste0(opt$destination, "/", pheno, "_", opt$ancestry, ".rda")
            output_vr <- paste0(opt$destination, "/", pheno, "_", opt$ancestry, ".varianceRatio.txt")

            if (dx_file_exists(output_rda) && dx_file_exists(output_vr)) {
                cat(paste0("Output already exists for ", pheno, ". Skipping.\n"))
                n_skipped <- n_skipped + 1
            } else {
                cmd <- create_step1_sim_cmd(
                    phenotype=pheno,
                    trait_type="binary",
                    covariates=covariates,
                    categorical_covariates=categorical_covariates,
                    ancestry=opt$ancestry,
                    sim_phenos_file=opt$pheno_file,
                    grm_prefix=grm_prefix,
                    step0_location=opt$step0_location,
                    sample_ids_file=sample_ids_file,
                    instance_type=opt$instance_type,
                    priority=opt$priority,
                    destination=opt$destination
                )
                cat(paste0("\nSubmitting: ", pheno, "\n"))
                cat(paste0(cmd, "\n"))
                if (!opt$dry_run) {
                    system(cmd)
                }
                n_submitted <- n_submitted + 1
            }
        }
        cat(paste0("\nSummary: ", n_submitted, " jobs submitted, ", n_skipped, " skipped\n"))
    }

    cat("\nDone!\n")
}

# Add arguments
option_list = list(
    make_option(c("-p", "--pheno_file"), type="character", default=NULL,
        help="path to the phenotype file on DNAnexus", metavar="character"),
    make_option(c("-a", "--ancestry"), type="character", default="eur",
        help="ancestry (eur, sas, afr, eas)", metavar="character"),
    make_option(c("--heritabilities"), type="character", default="0.01,0.10,0.2,0.3,0.5",
        help="comma-separated heritabilities matching phenotype file", metavar="character"),
    make_option(c("--n_reps"), type="integer", default=4,
        help="number of replicates per heritability", metavar="integer"),
    make_option(c("--step0_location"), type="character",
        default="/mnt/project/wes_ko_ukbb/data/saige/step0/vr_20k/",
        help="path to step0 files", metavar="character"),
    make_option(c("--sample_dir"), type="character",
        default="/wes_ko_ukbb/data/samples",
        help="path to sample ID files", metavar="character"),
    make_option(c("-c", "--covariates"), type="character",
        default="age,age2,age_sex,age2_sex,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10",
        help="comma-separated list of covariates", metavar="character"),
    make_option(c("--categorical_covariates"), type="character",
        default="sex",
        help="comma-separated list of categorical covariates", metavar="character"),
    make_option(c("-t", "--trait_type"), type="character", default="both",
        help="trait type to run: quantitative, binary, or both", metavar="character"),
    make_option(c("-d", "--destination"), type="character",
        default="/wes_ko_ukbb/data/saige/step1/simulated/",
        help="destination folder on DNAnexus", metavar="character"),
    make_option(c("--instance_type"), type="character",
        default="mem3_ssd1_v2_x16",
        help="instance type for DNAnexus jobs", metavar="character"),
    make_option(c("--priority"), type="character", default="low",
        help="job priority", metavar="character"),
    make_option(c("--dry_run"), action="store_true", default=FALSE,
        help="print commands without executing")
)

opt_parser <- OptionParser(add_help_option=TRUE, option_list=option_list)
opt <- parse_args(opt_parser)

main(opt)
