#!/usr/bin/env bash
# author: frederik lassen
# note: Generate FRQX file from existing PLINK files

source "/Users/flassen/Projects/00_dx_utils/dx_utils/dx_utils.sh"

set -o errexit
set -o nounset

# PLINK files location
readonly plink_dir="/wes_ko_ukbb/data/phased/gene_set_simulation/plink"

# Parameters
readonly ancestries=("eur")
readonly af="05"
readonly group="spliceai=0.50_cadd=28.1_revel=0.773"
readonly min_mac="10"

# Instance configuration
readonly instance_type="mem1_ssd1_v2_x2"
readonly priority="high"

for anc in "${ancestries[@]}"; do

    echo "=========================================="
    echo "Processing ancestry: ${anc}"
    echo "=========================================="

    for chr in 21; do

        echo ""
        echo "Chromosome: ${chr}"

        # PLINK file prefix
        plink_prefix="${plink_dir}/UKB.wes.chr${chr}.phased.qc_final.${anc}.af${af}.popmax.${group}.min_mac_${min_mac}"

        echo "PLINK files: ${plink_prefix}.{bed,bim,fam}"
        echo "Output: ${plink_prefix}.frqx.gz"
        echo ""

        # Check if FRQX already exists
        if [[ $(dx_file_exists "${plink_prefix}.frqx.gz") -eq "0" ]]; then

            # Check if PLINK files exist
            bed_exists=$(dx ls "${plink_prefix}.bed" 2>/dev/null | wc -l)
            bim_exists=$(dx ls "${plink_prefix}.bim" 2>/dev/null | wc -l)
            fam_exists=$(dx ls "${plink_prefix}.fam" 2>/dev/null | wc -l)

            if [[ ${bed_exists} -gt 0 && ${bim_exists} -gt 0 && ${fam_exists} -gt 0 ]]; then

                echo "Submitting FRQX generation job..."

                dx run app-swiss-army-knife \
                    -iin="${plink_prefix}.bed" \
                    -iin="${plink_prefix}.bim" \
                    -iin="${plink_prefix}.fam" \
                    -icmd="
                        echo 'Generating FRQX file from PLINK files'
                        echo 'Date: \$(date)'
                        echo ''

                        # Get base name
                        base=\$(basename ${plink_prefix})

                        echo 'PLINK files:'
                        ls -lh \${base}.*
                        echo ''

                        # Generate frequency counts
                        echo 'Running plink2 --freq counts...'
                        plink2 \
                            --bfile \${base} \
                            --freq counts \
                            --out \${base}

                        echo ''
                        echo 'Compressing FRQX file...'
                        gzip -f \${base}.frqx

                        echo ''
                        echo 'FRQX file created:'
                        ls -lh \${base}.frqx.gz

                        echo ''
                        echo 'First 5 lines:'
                        zcat \${base}.frqx.gz | head -5

                        echo ''
                        echo 'Complete! \$(date)'
                    " \
                    --instance-type ${instance_type} \
                    --folder=".${plink_dir}" \
                    --priority ${priority} \
                    --name "create_frqx_chr${chr}_${anc}" \
                    -y

                echo "Job submitted!"

            else
                >&2 echo "ERROR: PLINK files not found"
                >&2 echo "  .bed exists: ${bed_exists}"
                >&2 echo "  .bim exists: ${bim_exists}"
                >&2 echo "  .fam exists: ${fam_exists}"
            fi

        else
            echo "FRQX file already exists. Skipping..."
        fi

        echo ""

    done

done

echo "=========================================="
echo "All jobs submitted!"
echo "=========================================="
