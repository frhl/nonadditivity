#!/usr/bin/env bash
# author: frederik lassen
# note: this script converts per-chromosome dominance VCF files (chr20-22) to BGEN format for REGENIE

set -euo pipefail

readonly worker_script_local="_convert_vcf_to_bgen.sh"
readonly worker_script_remote="wes_ko_ukbb/scripts/simulation/regenie/set_based/_convert_vcf_to_bgen.sh"

# Upload worker script to remote
source "/Users/flassen/Projects/00_dx_utils/dx_utils/dx_utils.sh"
dx mkdir -p wes_ko_ukbb/scripts/simulation/regenie/set_based
dx_update_remote ${worker_script_remote} ${worker_script_local}

# Define paths
readonly in_dir="/wes_ko_ukbb/data/phased/encoded_dominance_gt/group_dominance/gene_map"
readonly out_dir="/wes_ko_ukbb/data/phased/encoded_dominance_gt_bgen/group_dominance/gene_map"
readonly group="spliceai=0.50_cadd=28.1_revel=0.773"

# Create output directory
dx mkdir -p ${out_dir}

# Define parameters
readonly ancestries=("eur")
#readonly ancestries=("eur" "sas" "afr" "eas")
readonly priority="low"
readonly instance_type="mem1_ssd1_v2_x2"

# Loop through ancestries
for anc in "${ancestries[@]}"; do

  echo "Processing ancestry: ${anc}"

  # Loop through chromosomes (20-22 only)
  for chr in {1..22}; do

    echo "  Processing chromosome: ${chr}"

    # Per-chromosome VCF path
    vcf_prefix="UKB.wes.chr${chr}.phased.qc_final.${anc}.af05.popmax.variants.group_dominance_scaling.gene_map"
    vcf_file="${in_dir}/${vcf_prefix}.vcf.gz"
    vcf_index="${in_dir}/${vcf_prefix}.vcf.gz.csi"

    # Output prefix
    out_prefix="${vcf_prefix}"
    out_full_path="${out_dir}/${out_prefix}"

    echo "    Input VCF: ${vcf_file}"
    echo "    Output prefix: ${out_full_path}"

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
          --destination "${out_dir}" \
          -y \
          --name "vcf_to_bgen_${anc}_chr${chr}"

      else
        >&2 echo "VCF file not found: ${vcf_file}. Skipping."
      fi

    else
      >&2 echo "Output already exists: ${out_full_path}.bgen. Skipping."
    fi

  done
done

echo "Done submitting VCF to BGEN conversion jobs for chr20-22!"
