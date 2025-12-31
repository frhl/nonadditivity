#!/usr/bin/env bash
# author: frederik lassen
# note: this script submits REGENIE step2 jobs for simulated phenotypes (set-based tests)

set -euo

readonly worker_script_local="_regenie_step2_set.sh"
readonly worker_script_remote="wes_ko_ukbb/scripts/simulation/regenie/set_based/_regenie_step2_set.sh"
readonly mask_def_local="mask_definitions.txt"
readonly mask_def_remote="wes_ko_ukbb/scripts/simulation/regenie/set_based/mask_definitions.txt"

# Upload worker script and mask definitions to remote
source "/Users/flassen/Projects/00_dx_utils/dx_utils/dx_utils.sh"
dx mkdir -p wes_ko_ukbb/scripts/simulation/regenie/set_based
dx_update_remote ${worker_script_remote} ${worker_script_local}
dx_update_remote ${mask_def_remote} ${mask_def_local}

# Define paths
readonly pheno_dir="/wes_ko_ukbb/data/phenotypes/simulated/revisions/051025/regenie_format"
readonly step1_dir="/wes_ko_ukbb/data/regenie/step1/simulated/revisions/051025"
readonly out_dir="/wes_ko_ukbb/data/regenie/step2_set/simulated/revisions/051025"
readonly sample_dir="/wes_ko_ukbb/data/samples"

# Create output directory
dx mkdir -p ${out_dir}

# Define parameters
readonly ancestries=("eur")
#readonly ancestries=("eur" "sas" "afr" "eas")
#readonly heritabilities=("0.01" "0.1" "0.2" "0.3" "0.5")
readonly heritabilities=("0.01" "0.2")
#readonly n_reps=4
readonly n_reps=1
readonly covariates="age_at_recruitment,age2,age_sex,age2_sex,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
readonly categorical_covariates="sex"
readonly priority="low"
readonly instance_type="mem3_ssd1_v2_x8"

# Gene/set-based test parameters
readonly group="spliceai=0.50_cadd=28.1_revel=0.773"
readonly af="05"

# Per-chromosome paths (chr20-22)
readonly bgen_dir="/wes_ko_ukbb/data/phased/encoded_dominance_gt_bgen/group_dominance/gene_map"
readonly set_list_dir="/wes_ko_ukbb/data/genesets/regenie_format"  # Directory containing annotation, set list, and AAF files

# Loop through ancestries
for anc in "${ancestries[@]}"; do

  # Convert ancestry to uppercase for sample files
  anc_upper=$(echo ${anc} | tr '[:lower:]' '[:upper:]')

  pheno_file="${pheno_dir}/051025_${anc}_null_phenos.regenie.txt.gz"
  keep_file="${sample_dir}/UKB.wes.qced.${anc}.samples"

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

  # Loop through chromosomes (20-22)
  #for chr in {1..22}; do
  for chr in 22; do

    echo "  Processing chromosome: ${chr}"

    # Loop through heritabilities
    for h2 in "${heritabilities[@]}"; do
      # Loop through replicates
      for rep in $(seq 1 ${n_reps}); do
        pheno="p_${h2}_continuous_${rep}"

        # Look up which LOCO file this phenotype uses from the pred file
        loco_filename=$(grep "^${pheno} " "${pred_local}" | awk '{print $2}' | xargs basename || true)

        if [[ -z "${loco_filename}" ]]; then
          echo "WARNING: Phenotype ${pheno} not found in pred file. Skipping this phenotype."
          continue
        fi

        loco_file="${step1_dir}/${loco_filename}"
        echo "Phenotype ${pheno} uses LOCO file: ${loco_filename}"

        # Per-chromosome set-based test
        out_prefix="regenie_step2_set.${anc}.chr${chr}.${pheno}"

        echo "Submitting REGENIE step2 set test: ${out_prefix}"

        # Check if output already exists
        if [[ $(dx ls "${out_dir}/${out_prefix}.regenie.gz" 2>/dev/null | wc -l) -eq 0 ]]; then

          # Use per-chromosome dominance BGEN files (from script 05)
          bgen_prefix="UKB.wes.chr${chr}.phased.qc_final.${anc}.af${af}.popmax.variants.group_dominance_scaling.gene_map"
          bgen_file="${bgen_dir}/${bgen_prefix}.bgen"
          sample_file="${bgen_dir}/${bgen_prefix}.sample"

          # Files from script 08 (annotation, set list, and AAF)
          anno_prefix="UKB.wes.chr${chr}.phased.qc_final.${anc}.af${af}.popmax.variants.${group}"
          anno_file="${set_list_dir}/${anno_prefix}.anno.txt.gz"
          set_list_file="${set_list_dir}/${anno_prefix}.setlist.txt.gz"
          aaf_file="${set_list_dir}/${anno_prefix}.aaf.txt.gz"

          # Check if all required files exist
          bgen_exists=$(dx ls "${bgen_file}" 2>/dev/null | wc -l)
          anno_exists=$(dx ls "${anno_file}" 2>/dev/null | wc -l)
          set_exists=$(dx ls "${set_list_file}" 2>/dev/null | wc -l)
          aaf_exists=$(dx ls "${aaf_file}" 2>/dev/null | wc -l)

          if [[ ${bgen_exists} -gt 0 && ${anno_exists} -gt 0 && ${set_exists} -gt 0 && ${aaf_exists} -gt 0 ]]; then

            dx run app-swiss-army-knife \
              -iin="/${worker_script_remote}" \
              -iin="/${mask_def_remote}" \
              -iin="${bgen_file}" \
              -iin="${sample_file}" \
              -iin="${pheno_file}" \
              -iin="${keep_file}" \
              -iin="${pred_file}" \
              -iin="${loco_file}" \
              -iin="${set_list_file}" \
              -iin="${anno_file}" \
              -iin="${aaf_file}" \
              -icmd="
                bash /mnt/project/${worker_script_remote} \
                  ${bgen_prefix}.bgen \
                  regenie_step1_${anc}_all_phenos_pred.list \
                  051025_${anc}_null_phenos.regenie.txt.gz \
                  ${pheno} \
                  UKB.wes.qced.${anc}.samples \
                  ${covariates} \
                  ${categorical_covariates} \
                  $(basename ${set_list_file}) \
                  $(basename ${anno_file}) \
                  $(basename ${aaf_file}) \
                  mask_definitions.txt \
                  ${out_prefix} &&
                echo '!!!!$(date)'
              " \
              --instance-type ${instance_type} \
              --priority ${priority} \
              --destination "${out_dir}" \
              -y \
              --name "regenie_step2_set_${anc}_chr${chr}_${pheno}"

          else
            echo "  WARNING: Not all required files exist for ${out_prefix}"
            echo "    BGEN: ${bgen_exists}, Annotation: ${anno_exists}, Set list: ${set_exists}, AAF: ${aaf_exists}"
          fi

        else
          echo "  Output already exists: ${out_prefix}.regenie.gz. Skipping."
        fi

      done
    done
  done
done

echo "Done submitting REGENIE step2 set-based jobs!"
