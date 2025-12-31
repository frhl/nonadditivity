#!/usr/bin/env bash
# author: frederik lassen
# note: this script compares SAIGE and REGENIE results for the same genes/phenotypes

readonly rscript_local="_compare_saige_regenie.R"
readonly rscript_remote="wes_ko_ukbb/scripts/simulation/_compare_saige_regenie.R"

# Upload R script to remote
source "/Users/flassen/Projects/00_dx_utils/dx_utils/dx_utils.sh"
dx mkdir -p wes_ko_ukbb/scripts/simulation
dx_update_remote ${rscript_remote} ${rscript_local}

# Define paths
readonly saige_dir="/wes_ko_ukbb/data/saige/step2/simulated/revisions/051025/combined"
readonly regenie_dir="/wes_ko_ukbb/data/regenie/step2/simulated/revisions/051025/combined"
readonly out_dir="/wes_ko_ukbb/data/comparison/saige_vs_regenie/simulated/revisions/051025"

# Create output directory
dx mkdir -p ${out_dir}

# Define parameters
readonly ancestries=("eur")
readonly annotations=("pLoF_damaging_missense")
readonly priority="low"
readonly instance_type="mem1_ssd1_v2_x8"

# Loop through ancestries and annotations
for anc in "${ancestries[@]}"; do
  for anno in "${annotations[@]}"; do

    # Define output file
    out_prefix="saige_vs_regenie_${anc}_${anno}"
    out_file="${out_dir}/${out_prefix}.pdf"

    echo "Submitting comparison job for: ${out_prefix}"

    # Check if output already exists
    echo "  Checking if output exists: ${out_file}"
    out_exists=$(dx ls "${out_file}" 2>/dev/null | wc -l)
    echo "    Output exists: ${out_exists}"

    if [[ ${out_exists} -eq 0 ]]; then

      # Submit job
      dx run app-swiss-army-knife \
        -iimage_file="/docker/rsuite.tar.gz" \
        -iin="/${rscript_remote}" \
        -icmd="
          Rscript /mnt/project/${rscript_remote} \
            --saige_dir /mnt/project${saige_dir} \
            --regenie_dir /mnt/project${regenie_dir} \
            --ancestry ${anc} \
            --annotation ${anno} \
            --output ${out_prefix}.pdf &&
          echo '!!!!$(date)'
        " \
        --instance-type ${instance_type} \
        --priority ${priority} \
        --destination "${out_dir}" \
        -y \
        --name "compare_${anc}_${anno}"

    else
      echo "  Output already exists: ${out_file}. Skipping."
    fi

  done
done

echo "Done submitting comparison jobs!"
