#!/usr/bin/env bash
# author: frederik lassen
# note: this script converts VEP canonical annotation files to REGENIE annotation and set list format

set -euo

readonly rscript_local="_create_regenie_anno.R"
readonly rscript_remote="wes_ko_ukbb/scripts/simulation/regenie/set_based/_create_regenie_anno.R"

# Upload R script to remote
source "/Users/flassen/Projects/00_dx_utils/dx_utils/dx_utils.sh"
dx mkdir -p wes_ko_ukbb/scripts/simulation/regenie/set_based
dx_update_remote ${rscript_remote} ${rscript_local}

# Define paths
readonly group="spliceai=0.50_cadd=28.1_revel=0.773"
readonly vep_dir="/wes_ko_ukbb/data/vep_loftee/csqs/${group}"
readonly markers_dir="/wes_ko_ukbb/data/phased/encoded_dominance_gt/group_dominance/gene_map"
readonly out_dir="/wes_ko_ukbb/data/genesets/regenie_format"

# Create output directory
dx mkdir -p ${out_dir}

# Define parameters
readonly ancestries=("eur")
#readonly ancestries=("eur" "sas" "afr" "eas")
readonly af="05"
readonly priority="low"
readonly instance_type="mem1_ssd1_v2_x2"

# Loop through ancestries
for anc in "${ancestries[@]}"; do

  echo "Processing ancestry: ${anc}"

  # Loop through chromosomes (20-22)
  for chr in {1..22}; do

    echo "  Processing chromosome: ${chr}"

    # Input VEP canonical annotation file
    vep_file="${vep_dir}/UKB.chr${chr}.phased_sites.${group}.canonical.txt.gz"

    # Markers file (variants to keep from dominance encoding)
    markers_prefix="UKB.wes.chr${chr}.phased.qc_final.${anc}.af${af}.popmax.variants.group_dominance_scaling.gene_map"
    markers_file="${markers_dir}/${markers_prefix}.sites.txt"

    # Allele count file (for AAF calculation and weights)
    ac_file="/wes_ko_ukbb/data/phased/final_qc/UKB.wes.chr${chr}.phased.qc_final.${anc}.frqx.gz"

    # Output files
    out_prefix="UKB.wes.chr${chr}.phased.qc_final.${anc}.af${af}.popmax.variants.${group}"
    out_anno="${out_dir}/${out_prefix}.anno.txt"
    out_setlist="${out_dir}/${out_prefix}.setlist.txt"
    out_aaf="${out_dir}/${out_prefix}.aaf.txt"

    echo "    VEP file: ${vep_file}"
    echo "    Markers: ${markers_file}"
    echo "    AC file: ${ac_file}"
    echo "    Output anno: ${out_anno}.gz"
    echo "    Output setlist: ${out_setlist}.gz"
    echo "    Output AAF: ${out_aaf}.gz"

    # Check if output already exists
    anno_exists=$(dx ls "${out_anno}.gz" 2>/dev/null | wc -l)
    setlist_exists=$(dx ls "${out_setlist}.gz" 2>/dev/null | wc -l)
    aaf_exists=$(dx ls "${out_aaf}.gz" 2>/dev/null | wc -l)

    if [[ ${anno_exists} -eq 0 || ${setlist_exists} -eq 0 || ${aaf_exists} -eq 0 ]]; then

      # Check if input files exist
      vep_exists=$(dx ls "${vep_file}" 2>/dev/null | wc -l)
      markers_exists=$(dx ls "${markers_file}" 2>/dev/null | wc -l)
      ac_exists=$(dx ls "${ac_file}" 2>/dev/null | wc -l)

      if [[ ${vep_exists} -gt 0 && ${markers_exists} -gt 0 && ${ac_exists} -gt 0 ]]; then

        dx run app-swiss-army-knife \
          -iimage_file="/docker/rsuite.tar.gz" \
          -iin="/${rscript_remote}" \
          -iin="${vep_file}" \
          -iin="${markers_file}" \
          -iin="${ac_file}" \
          -icmd="
            Rscript /mnt/project/${rscript_remote} \
              --vep_file UKB.chr${chr}.phased_sites.${group}.canonical.txt.gz \
              --markers_file ${markers_prefix}.sites.txt \
              --ac_file UKB.wes.chr${chr}.phased.qc_final.${anc}.frqx.gz \
              --out_anno ${out_prefix}.anno.txt \
              --out_setlist ${out_prefix}.setlist.txt \
              --out_aaf ${out_prefix}.aaf.txt \
              --chr chr${chr} &&
            gzip ${out_prefix}.anno.txt &&
            gzip ${out_prefix}.setlist.txt &&
            gzip ${out_prefix}.aaf.txt &&
            echo '!!!!$(date)'
          " \
          --instance-type ${instance_type} \
          --priority ${priority} \
          --destination "${out_dir}" \
          -y \
          --name "regenie_anno_${anc}_chr${chr}"

      else
        echo "  WARNING: Input files not found"
        echo "    VEP exists: ${vep_exists}"
        echo "    Markers exists: ${markers_exists}"
        echo "    AC exists: ${ac_exists}"
      fi

    else
      echo "  Output files already exist. Skipping."
    fi

  done
done

echo "Done submitting REGENIE annotation file creation jobs!"
