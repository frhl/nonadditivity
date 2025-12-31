#!/usr/bin/env bash
# author: frederik lassen
# note: this script submits REGENIE step1 jobs for simulated phenotypes

set -exuo pipefail

readonly worker_script_local="02_regenie_step1.sh"
readonly worker_script_remote="wes_ko_ukbb/scripts/simulation/regenie/02_regenie_step1.sh"

# Upload worker script to remote
source "/Users/flassen/Projects/00_dx_utils/dx_utils/dx_utils.sh"
dx mkdir -p wes_ko_ukbb/scripts/simulation/regenie
dx_update_remote ${worker_script_remote} ${worker_script_local}

# Define paths
readonly pheno_dir="/wes_ko_ukbb/data/phenotypes/simulated/revisions/051025/regenie_format"
readonly geno_dir="/nbaya/regenie/data/genotypes"
readonly sample_dir="/wes_ko_ukbb/data/samples"
readonly out_dir="/wes_ko_ukbb/data/regenie/step1/simulated/revisions/051025"

# Create output directory
dx mkdir -p ${out_dir}

# Define parameters
readonly ancestries=("eur")
#readonly ancestries=("eur" "sas" "afr" "eas")
readonly heritabilities=("0.01")
#readonly heritabilities=("0.01" "0.10" "0.2" "0.3" "0.5")
readonly n_reps=1
#readonly n_reps=4
readonly covariates="age_at_recruitment,age2,age_sex,age2_sex,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
readonly categorical_covariates="sex"
readonly instance_type="mem3_ssd1_v2_x8"
readonly priority="low"
readonly n_threads=8
readonly bsize=1000

# Loop through ancestries
for anc in "${ancestries[@]}"; do

  # Convert ancestry to uppercase for genotype files
  anc_upper=$(echo ${anc} | tr '[:lower:]' '[:upper:]')

  pheno_file="${pheno_dir}/051025_${anc}_null_phenos.regenie.txt.gz"

  # REGENIE step1 uses genotype array data
  bfile="${geno_dir}/ukb22418_b0_v2.autosomes"
  keep_file="${sample_dir}/UKB.wes.qced.${anc}.samples"
  extract_file="${geno_dir}/ukb22418_b0_v2.autosomes.qced.${anc_upper}.snplist"

  echo "Processing ancestry: ${anc}"

  # Loop through heritabilities
  for h2 in "${heritabilities[@]}"; do
    # Loop through replicates
    for rep in $(seq 1 ${n_reps}); do
      pheno="p_${h2}_continuous_${rep}"

      out_prefix="regenie_step1_${anc}_${pheno}"

      echo "Submitting REGENIE step1: ${out_prefix}"

      # Check if output already exists
      if [[ $(dx ls "${out_dir}/${out_prefix}_pred.list" 2>/dev/null | wc -l) -eq 0 ]]; then

        dx run app-swiss-army-knife \
          -iin="/${worker_script_remote}" \
          -iin="${bfile}.bed" \
          -iin="${bfile}.bim" \
          -iin="${bfile}.fam" \
          -iin="${pheno_file}" \
          -iin="${keep_file}" \
          -iin="${extract_file}" \
          -icmd="
            bash /mnt/project/${worker_script_remote} \
              ukb22418_b0_v2.autosomes \
              051025_${anc}_null_phenos.regenie.txt.gz \
              ${pheno} \
              UKB.wes.qced.${anc}.samples \
              ukb22418_b0_v2.autosomes.qced.${anc_upper}.snplist \
              ${covariates} \
              ${categorical_covariates} \
              ${n_threads} \
              ${bsize} \
              ${out_prefix} &&
            echo '!!!!$(date)'
          " \
          --instance-type ${instance_type} \
          --priority ${priority} \
          --destination "${out_dir}" \
          -y \
          --name "regenie_step1_${anc}_${pheno}"

      else
        >&2 echo "Output already exists: ${out_prefix}_pred.list. Skipping."
      fi

    done
  done
done

echo "Done submitting REGENIE step1 jobs!"
