#!/usr/bin/env bash
# author: frederik lassen
# note: Create SAIGE group/weight files for set-based testing with homozygous carrier filtering
# Ensures all variants have at least N homozygous carriers to make dominance encoding meaningful
# Uses pre-calculated allele frequencies from QC'd BCF files and filters variants in R script

source "/Users/flassen/Projects/00_dx_utils/dx_utils/dx_utils.sh"

readonly rscript_local="_create_saige_weights.R"
readonly rscript_remote="wes_ko_ukbb/scripts/simulation/saige/set_based/_create_saige_weights.R"

# Upload R script to remote
dx mkdir -p wes_ko_ukbb/scripts/simulation/saige/set_based
dx_update_remote ${rscript_remote} ${rscript_local}

# Define paths
readonly annotation="spliceai=0.50_cadd=28.1_revel=0.773"
readonly anno_dir="/wes_ko_ukbb/data/vep_loftee/csqs/${annotation}"
readonly markers_dir="/wes_ko_ukbb/data/phased/simulated/set_based/dominance_encoded"
readonly ac_dir="/wes_ko_ukbb/data/phased/final_qc"
readonly popmax_dir="/wes_ko_ukbb/data/variants/popmax_exclude"
readonly out_dir="/wes_ko_ukbb/data/genesets/simulated/set_based/saige_weights"

# Create output directory
dx mkdir -p ${out_dir}

# Define parameters
#readonly ancestries=("eur")
readonly ancestries=("eur" "sas" "afr" "eas")
readonly af="05"              # Filter to AF < 0.05
readonly min_mac=10           # Minimum minor allele count
readonly min_hom_carriers=2   # Minimum homozygous carriers
readonly priority="high"
readonly instance_type="mem3_ssd1_v2_x4"

echo "=========================================="
echo "Creating SAIGE group/weight files"
echo "=========================================="
echo "Minimum MAC: ${min_mac}"
echo "Minimum homozygous carriers: ${min_hom_carriers}"
echo "Max AF: 0.${af}"
echo ""

for pop in "${ancestries[@]}"; do
  echo "Processing ancestry: ${pop}"

  for CHR in {20..22}; do

    # Input files
    ac_path="${ac_dir}/UKB.wes.chr${CHR}.phased.qc_final.${pop}.frqx.gz"
    anno_path="${anno_dir}/UKB.chr${CHR}.phased_sites.${annotation}.canonical.txt.gz"
    markers_path="${markers_dir}/UKB.wes.chr${CHR}.phased.qc_final.${pop}.group_dominance_scaling.gene_map.sites.txt"
    popmax_exclude="${popmax_dir}/gnomad.exomes.r2.1.1.grch38.popmax${af}.tsv"
  
    echo "${anno_path}"

    # Output file
    out_prefix="UKB.wes.chr${CHR}.phased.qc_final.${pop}.af${af}.popmax.${annotation}.min_mac_${min_mac}.min_hom_${min_hom_carriers}.saige_weights"

    echo "  Chr${CHR}: ${out_prefix}"

    # Check if output already exists
    if [[ $(dx_file_exists "${out_dir}/${out_prefix}.txt.gz") -eq 0 ]]; then

      # Check if all input files exist
      anno_exists=$(dx ls "${anno_path}" 2>/dev/null | wc -l)
      markers_exists=$(dx ls "${markers_path}" 2>/dev/null | wc -l)
      ac_exists=$(dx ls "${ac_path}" 2>/dev/null | wc -l)
      popmax_exists=$(dx ls "${popmax_exclude}" 2>/dev/null | wc -l)

      if [[ ${anno_exists} -gt 0 && ${markers_exists} -gt 0 && ${ac_exists} -gt 0 && ${popmax_exists} -gt 0 ]]; then

        dx run app-swiss-army-knife \
          -iimage_file="/docker/rsuite.tar.gz" \
          -iin="/${rscript_remote}" \
          -iin="${anno_path}" \
          -iin="${markers_path}" \
          -iin="${ac_path}" \
          -iin="${popmax_exclude}" \
          -icmd="
            Rscript /mnt/project/${rscript_remote} \
              --input_path $(basename ${anno_path}) \
              --markers_path $(basename ${markers_path}) \
              --ac_path $(basename ${ac_path}) \
              --popmax_exclude $(basename ${popmax_exclude}) \
              --max_af 0.${af} \
              --min_mac ${min_mac} \
              --min_hom_carriers ${min_hom_carriers} \
              --output_path ${out_prefix}.txt &&
            gzip ${out_prefix}.txt &&
            echo '!!!!$(date)'
          " \
          --instance-type ${instance_type} \
          --folder=".${out_dir}" \
          --priority ${priority} \
          --name "saige_weights_chr${CHR}_${pop}" \
          -y

      else
        >&2 echo "  ERROR: Missing files - Anno: ${anno_exists}, Markers: ${markers_exists}, AC: ${ac_exists}, Popmax: ${popmax_exists}"
      fi

    else
      echo "  Output already exists. Skipping."
    fi

  done

  echo ""

done

echo "=========================================="
echo "Done submitting SAIGE weight jobs!"
echo "=========================================="
