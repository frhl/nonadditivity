#!/usr/bin/env bash
# author: frederik lassen
# note: Encode dominance (nonadditive) effects in VCF files for set-based testing
# Uses the transform tool to create X^D encoding with gene mapping
# Uses original QC'd BCF files - filtering is done at the weight creation step

source "/Users/flassen/Projects/00_dx_utils/dx_utils/dx_utils.sh"

set -euo pipefail

readonly threads=24
readonly out_dir="/wes_ko_ukbb/data/phased/simulated/set_based/dominance_encoded"
readonly in_dir="/wes_ko_ukbb/data/phased/final_qc"
readonly gene_map_dir="/wes_ko_ukbb/data/genesets/dominance_weights"
readonly docker="/docker/call_chets_0.1.14.tar.gz"

# Create output directory
dx mkdir -p ${out_dir}

# Define parameters
readonly ancestries=("eur" "afr")
#readonly ancestries=("eur" "sas" "afr" "eas")
readonly priority="high"
readonly instance_type="mem3_ssd1_v2_x8"

echo "=========================================="
echo "Encoding dominance for set-based testing"
echo "=========================================="
echo ""

for pop in "${ancestries[@]}"; do
  echo "Processing ancestry: ${pop}"

  for CHR in {20..22}; do

    # Gene map file (from existing pipeline)
    gene_map="${gene_map_dir}/UKB.chr${CHR}.phased_sites.spliceai=0.50_cadd=28.1_revel=0.773.canonical.gene_map.txt"

    # Input and output files
    input_vcf="${in_dir}/UKB.wes.chr${CHR}.phased.qc_final.${pop}.bcf"
    out_prefix="UKB.wes.chr${CHR}.phased.qc_final.${pop}.group_dominance_scaling.gene_map"

    echo "  Chr${CHR}: ${out_prefix}"

    # Check if output already exists
    if [[ $(dx_file_exists "${out_dir}/${out_prefix}.vcf.gz") -eq 0 ]]; then

      # Check if input and gene map exist
      input_exists=$(dx ls "${input_vcf}" 2>/dev/null | wc -l)
      gene_map_exists=$(dx ls "${gene_map}" 2>/dev/null | wc -l)

      if [[ ${input_exists} -gt 0 && ${gene_map_exists} -gt 0 ]]; then

        dx run app-swiss-army-knife \
          -iimage_file="${docker}" \
          -icmd="
            transform --input /mnt/project${input_vcf} --all-info --gene-map /mnt/project${gene_map} --set-variant-id --scale-dosages | bgzip > ${out_prefix}.vcf.gz &&
            tabix -C ${out_prefix}.vcf.gz &&
            zcat ${out_prefix}.vcf.gz | grep -v '##' | cut -f3 | tail -n+2 > ${out_prefix}.sites.txt &&
            echo '!!!!$(date)'
          " \
          --instance-type ${instance_type} \
          --folder=".${out_dir}" \
          --priority ${priority} \
          --name "encode_dom_chr${CHR}_${pop}" \
          -y

      else
        >&2 echo "  ERROR: Missing files - Input: ${input_exists}, Gene map: ${gene_map_exists}"
      fi

    else
      echo "  Output already exists. Skipping."
    fi

  done

  echo ""

done

echo "=========================================="
echo "Done submitting dominance encoding jobs!"
echo "=========================================="
