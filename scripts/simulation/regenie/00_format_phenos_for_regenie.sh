#!/usr/bin/env bash
# author: frederik lassen
# note: this script formats phenotype files for REGENIE (requires FID IID format)

set -euo pipefail

readonly rscript_local="00_format_phenos_for_regenie.R"
readonly rscript_remote="wes_ko_ukbb/scripts/simulation/regenie/00_format_phenos_for_regenie.R"

# Upload R script to remote
source "/Users/flassen/Projects/00_dx_utils/dx_utils/dx_utils.sh"
dx mkdir -p wes_ko_ukbb/scripts/simulation/regenie
dx_update_remote ${rscript_remote} ${rscript_local}

# Define paths
readonly pheno_dir="/wes_ko_ukbb/data/phenotypes/simulated/revisions/051025"
readonly out_dir="/wes_ko_ukbb/data/phenotypes/simulated/revisions/051025/regenie_format"

# Create output directory
dx mkdir -p ${out_dir}

# Define parameters
readonly ancestries=("eur")
#readonly ancestries=("eur" "sas" "afr" "eas")

for anc in "${ancestries[@]}"; do
  input_file="${pheno_dir}/051025_${anc}_null_phenos.txt.gz.tsv.gz"
  output_file="051025_${anc}_null_phenos.regenie.txt.gz"

  echo "Submitting formatting job for ancestry: ${anc}"
  echo "Input: ${input_file}"
  echo "Output: ${out_dir}/${output_file}"

  dx run app-swiss-army-knife \
    -iimage_file="/docker/rsuite.tar.gz" \
    -icmd="
      Rscript /mnt/project/${rscript_remote} \
        --input /mnt/project/${input_file} \
        --output ${output_file%.gz} &&
      gzip ${output_file%.gz} &&
      echo '!!!!$(date)'
    " \
    --instance-type mem1_ssd1_v2_x2 \
    --folder=".${out_dir}" \
    --priority low \
    --name format_phenos_regenie_${anc} -y
done

echo "Done submitting formatting jobs!"
