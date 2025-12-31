#!/usr/bin/env bash
# author: frederik lassen
# note: this script submits REGENIE step2 jobs for simulated phenotypes (variant tests)

set -euo pipefail

readonly worker_script_local="_regenie_step2.sh"
readonly worker_script_remote="wes_ko_ukbb/scripts/simulation/regenie/_regenie_step2.sh"

# Upload worker script to remote
source "/Users/flassen/Projects/00_dx_utils/dx_utils/dx_utils.sh"
dx mkdir -p wes_ko_ukbb/scripts/simulation/regenie
dx_update_remote ${worker_script_remote} ${worker_script_local}

# Define paths
readonly pheno_dir="/wes_ko_ukbb/data/phenotypes/simulated/revisions/051025/regenie_format"
readonly step1_dir="/wes_ko_ukbb/data/regenie/step1/simulated/revisions/051025"
readonly out_dir="/wes_ko_ukbb/data/regenie/step2/simulated/revisions/051025"
readonly sample_dir="/wes_ko_ukbb/data/samples"

# Create output directory
dx mkdir -p ${out_dir}

# Define parameters
readonly ancestries=("eur")
#readonly ancestries=("eur" "sas" "afr" "eas")
readonly heritabilities=("0.01" "0.1" "0.2" "0.3" "0.5")
readonly n_reps=4
readonly covariates="age_at_recruitment,age2,age_sex,age2_sex,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
readonly categorical_covariates="sex"
readonly priority="low"
readonly instance_type="mem3_ssd1_v2_x8"

# Exome file parameters
readonly group="spliceai=0.50_cadd=28.1_revel=0.773"
readonly annotations=("pLoF_damaging_missense")
readonly modes=("additive" "recessive" "dominance")
readonly af="05"
readonly pp="0.90"

# Loop through ancestries
for anc in "${ancestries[@]}"; do

  # Convert ancestry to uppercase for sample files
  anc_upper=$(echo ${anc} | tr '[:lower:]' '[:upper:]')

  pheno_file="${pheno_dir}/051025_${anc}_null_phenos.regenie.txt.gz"
  keep_file="${sample_dir}/UKB.wes.qced.${anc}.samples"
  exome_dir="/wes_ko_ukbb/data/phased/encode_alt_qced/${anc}/${group}/vcf_plus_plink"
  bgen_dir="/wes_ko_ukbb/data/phased/encode_alt_qced_canonical_bgen/${anc}/${group}"

  echo "Processing ancestry: ${anc}"

  # Step1 output files (multi-phenotype mode)
  pred_file="${step1_dir}/regenie_step1_${anc}_all_phenos_pred.list"

  # Check if pred file exists
  if [[ $(dx ls "${pred_file}" 2>/dev/null | wc -l) -eq 0 ]]; then
    >&2 echo "Step1 pred file not found: ${pred_file}. Skipping ancestry ${anc}."
    continue
  fi

  # Download pred file to determine LOCO file mapping for each phenotype
  echo "Reading pred file to map phenotypes to LOCO files..."
  pred_local=$(mktemp)
  dx cat "${pred_file}" > "${pred_local}"

  # Loop through heritabilities
  for h2 in "${heritabilities[@]}"; do
    # Loop through replicates
    for rep in $(seq 1 ${n_reps}); do
      pheno="p_${h2}_continuous_${rep}"

      # Look up which LOCO file this phenotype uses from the pred file
      # Note: grep returns non-zero if no match, so we need to handle that
      loco_filename=$(grep "^${pheno} " "${pred_local}" | awk '{print $2}' | xargs basename || true)

      if [[ -z "${loco_filename}" ]]; then
        echo "WARNING: Phenotype ${pheno} not found in pred file. Skipping this phenotype."
        continue
      fi

      loco_file="${step1_dir}/${loco_filename}"
      echo "Phenotype ${pheno} uses LOCO file: ${loco_filename}"

      # Loop through annotations
      for anno in "${annotations[@]}"; do
        # Loop through modes
        for mode in "${modes[@]}"; do

          out_prefix="regenie_step2.${anc}.${pheno}.${anno}.${mode}"

          echo "Submitting REGENIE step2: ${out_prefix}"

          # Check if output already exists
          if [[ $(dx ls "${out_dir}/${out_prefix}.regenie.gz" 2>/dev/null | wc -l) -eq 0 ]]; then

            # Use BGEN for all modes (dominance, additive, recessive)
            bgen_prefix="UKB.wes.merged.phased.qc_final.${anc}.af${af}.popmax.pp${pp}.${group}.${anno}.${mode}.auto"
            bgen_file="${bgen_dir}/${bgen_prefix}.bgen"
            sample_file="${bgen_dir}/${bgen_prefix}.sample"

            # Check if BGEN file exists
            if [[ $(dx ls "${bgen_file}" 2>/dev/null | wc -l) -gt 0 ]]; then

              dx run app-swiss-army-knife \
                -iin="/${worker_script_remote}" \
                -iin="${bgen_file}" \
                -iin="${sample_file}" \
                -iin="${pheno_file}" \
                -iin="${keep_file}" \
                -iin="${pred_file}" \
                -iin="${loco_file}" \
                -icmd="
                  bash /mnt/project/${worker_script_remote} \
                    ${bgen_prefix}.bgen \
                    regenie_step1_${anc}_all_phenos_pred.list \
                    051025_${anc}_null_phenos.regenie.txt.gz \
                    ${pheno} \
                    UKB.wes.qced.${anc}.samples \
                    ${covariates} \
                    ${categorical_covariates} \
                    ${out_prefix} &&
                  echo '!!!!$(date)'
                " \
                --instance-type ${instance_type} \
                --priority ${priority} \
                --destination "${out_dir}" \
                -y \
                --name "regenie_step2_${anc}_${pheno}_${anno}_${mode}"

            else
              >&2 echo "BGEN file not found: ${bgen_file}. Skipping."
            fi

          else
            >&2 echo "Output already exists: ${out_prefix}.regenie.gz. Skipping."
          fi

        done
      done
    done
  done
done

echo "Done submitting REGENIE step2 jobs!"
