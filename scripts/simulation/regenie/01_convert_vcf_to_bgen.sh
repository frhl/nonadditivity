#!/usr/bin/env bash
# author: frederik lassen
# note: this script converts VCF files (dominance, additive, recessive) to BGEN format for REGENIE

set -euo pipefail

readonly worker_script_local="_convert_vcf_to_bgen.sh"
readonly worker_script_remote="wes_ko_ukbb/scripts/simulation/regenie/_convert_vcf_to_bgen.sh"

# Upload worker script to remote
source "/Users/flassen/Projects/00_dx_utils/dx_utils/dx_utils.sh"
dx mkdir -p wes_ko_ukbb/scripts/simulation/regenie
dx_update_remote ${worker_script_remote} ${worker_script_local}

# Define paths
readonly in_dir="/wes_ko_ukbb/data/phased/encode_alt_qced_canonical"
readonly out_dir="/wes_ko_ukbb/data/phased/encode_alt_qced_canonical_bgen"
readonly group="spliceai=0.50_cadd=28.1_revel=0.773"
readonly annotations=("pLoF_damaging_missense")
readonly modes=("dominance" "additive" "recessive")
readonly af="05"
readonly pp="0.90"

# Create output directory
dx mkdir -p ${out_dir}

# Define parameters
readonly ancestries=("eur")
#readonly ancestries=("eur" "sas" "afr" "eas")
readonly priority="low"
readonly instance_type="mem3_ssd1_v2_x8"

# Loop through ancestries
for anc in "${ancestries[@]}"; do

  echo "Processing ancestry: ${anc}"

  # Create ancestry-specific output directory
  dx mkdir -p ${out_dir}/${anc}/${group}

  # Loop through annotations
  for anno in "${annotations[@]}"; do

    echo "  Processing annotation: ${anno}"

    # Loop through modes
    for mode in "${modes[@]}"; do

      echo "    Processing mode: ${mode}"

      # VCF path
      vcf_dir="${in_dir}/${anc}/${group}/vcf_plus_plink/force_chr_name"
      vcf_prefix="UKB.wes.merged.phased.qc_final.${anc}.af${af}.popmax.pp${pp}.${group}.${anno}.${mode}.auto"
      vcf_file="${vcf_dir}/${vcf_prefix}.vcf.gz"
      vcf_index="${vcf_dir}/${vcf_prefix}.vcf.gz.csi"

      # Output prefix
      out_prefix="${vcf_prefix}"
      out_full_path="${out_dir}/${anc}/${group}/${out_prefix}"

      echo "      Input VCF: ${vcf_file}"
      echo "      Output prefix: ${out_full_path}"

      # Check if output already exists
      if [[ $(dx ls "${out_full_path}.bgen" 2>/dev/null | wc -l) -eq 0 ]]; then

        # Check if input VCF exists
        if [[ $(dx ls "${vcf_file}" 2>/dev/null | wc -l) -gt 0 ]]; then

          dx run app-swiss-army-knife \
            -iin="/${worker_script_remote}" \
            -iin="${vcf_file}" \
            -iin="${vcf_index}" \
            -icmd="
              bash /mnt/project/${worker_script_remote} \
                ${vcf_prefix}.vcf.gz \
                ${out_prefix} &&
              echo '!!!!$(date)'
            " \
            --instance-type ${instance_type} \
            --priority ${priority} \
            --destination "${out_dir}/${anc}/${group}" \
            -y \
            --name "vcf_to_bgen_${anc}_${anno}_${mode}"

        else
          >&2 echo "VCF file not found: ${vcf_file}. Skipping."
        fi

      else
        >&2 echo "Output already exists: ${out_full_path}.bgen. Skipping."
      fi

    done

  done
done

echo "Done submitting VCF to BGEN conversion jobs!"
