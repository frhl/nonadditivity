#!/usr/bin/env bash
# author: frederik lassen
# note: this script submits SAIGE step1 jobs for simulated phenotypes

set -euo pipefail

# Define paths
readonly pheno_dir="/wes_ko_ukbb/data/phenotypes/simulated/revisions/051025"
readonly step0_location="/wes_ko_ukbb/data/saige/step0/vr_20k/"
readonly destination="/wes_ko_ukbb/data/saige/step1/simulated/revisions/051025/"

# Create output directory
dx mkdir -p ${destination}

# Define parameters
#readonly ancestries=("eur")
readonly ancestries=("sas" "afr" "eas")
readonly heritabilities="0.01,0.1,0.2,0.3,0.5"
#readonly heritabilities="0.1"
readonly n_reps="4"
readonly covariates="age_at_recruitment,age2,age_sex,age2_sex,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
readonly categorical_covariates="sex"
readonly instance_type="mem3_ssd1_v2_x4"
readonly priority="high"
readonly trait_type="quantitative"  # or "binary" or "both"
readonly dry_run=""  # Set to "--dry_run" to print commands without executing

for ancestry in "${ancestries[@]}"; do
  pheno_file="${pheno_dir}/051025_${ancestry}_null_phenos.txt.gz.tsv.gz"

  echo "Submitting step1 jobs for ancestry: ${ancestry}"
  echo "Using phenotype file: ${pheno_file}"

  # Run the R script locally to submit jobs
  Rscript 01_sim_step1.R \
    --pheno_file ${pheno_file} \
    --ancestry ${ancestry} \
    --heritabilities ${heritabilities} \
    --n_reps ${n_reps} \
    --step0_location ${step0_location} \
    --covariates ${covariates} \
    --categorical_covariates ${categorical_covariates} \
    --trait_type ${trait_type} \
    --destination ${destination} \
    --instance_type ${instance_type} \
    --priority ${priority} \
    ${dry_run}
done

echo "Done submitting step1 jobs!"
