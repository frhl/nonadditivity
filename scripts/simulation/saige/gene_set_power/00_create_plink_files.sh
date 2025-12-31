#!/usr/bin/env bash
# author: frederik lassen
# note: Convert raw BCF files to PLINK format for gene set simulation
# Filters to variants present in SAIGE weight files (pLoF + damaging_missense)

source "/Users/flassen/Projects/00_dx_utils/dx_utils/dx_utils.sh"

set -o errexit
set -o nounset

# =============================================================================
# CONFIGURATION
# =============================================================================

# Input: Raw BCF files (additive genotypes)
readonly in_dir="/wes_ko_ukbb/data/phased/export_alt_qced"

# SAIGE weight files (to get variant list)
readonly saige_dir="/wes_ko_ukbb/data/genesets/dominance_weights/af05"

# Output: PLINK files for simulation
readonly out_dir="/wes_ko_ukbb/data/phased/gene_set_simulation/plink"

# Create output directory
dx mkdir -p ${out_dir}

# Parameters
readonly ancestries=("eur")
readonly af="05"
readonly group="spliceai=0.50_cadd=28.1_revel=0.773"
readonly min_mac="10"

# Instance configuration
readonly instance_type="mem3_ssd1_v2_x4"
readonly priority="high"

# =============================================================================
# MAIN LOOP - START WITH CHR21 FOR DEBUGGING
# =============================================================================

for anc in "${ancestries[@]}"; do

    echo "=========================================="
    echo "Processing ancestry: ${anc}"
    echo "=========================================="

    # Start with chr21 only
    for chr in 21; do

        echo ""
        echo "Chromosome: ${chr}"

        # Input BCF (raw additive genotypes)
        input_bcf="${in_dir}/UKB.wes.chr${chr}.phased.qc_final.${anc}.af${af}.popmax.variants.bcf"

        # SAIGE weight file (to get variant list)
        saige_weights="${saige_dir}/UKB.wes.chr${chr}.phased.full_qc.${anc}.af${af}.popmax.variants.${group}.min_mac_${min_mac}.saige_with_${anc}_weights.txt.gz"

        # Output prefix
        out_prefix="UKB.wes.chr${chr}.phased.qc_final.${anc}.af${af}.popmax.${group}.min_mac_${min_mac}"

        echo "Input BCF: ${input_bcf}"
        echo "SAIGE weights: ${saige_weights}"
        echo "Output: ${out_dir}/${out_prefix}.{bed,bim,fam}"
        echo ""

        # Check if output already exists
        if [[ $(dx_file_exists "${out_dir}/${out_prefix}.bed") -eq "0" ]]; then

            # Check if input files exist
            bcf_exists=$(dx ls "${input_bcf}" 2>/dev/null | wc -l)
            weights_exists=$(dx ls "${saige_weights}" 2>/dev/null | wc -l)

            if [[ ${bcf_exists} -gt 0 && ${weights_exists} -gt 0 ]]; then

                echo "Submitting conversion job..."

                dx run app-swiss-army-knife \
                    -iin="${input_bcf}" \
                    -iin="${saige_weights}" \
                    -icmd="
                        echo 'Starting BCF to PLINK conversion for chr${chr}'
                        echo 'Date: \$(date)'
                        echo ''

                        # Extract variant list from SAIGE weight file
                        echo 'Extracting variant IDs from SAIGE weight file...'
                        zcat $(basename ${saige_weights}) | \
                            grep '^ENSG' | \
                            awk 'NR % 3 == 1' | \
                            cut -d' ' -f3- | \
                            tr ' ' '\n' | \
                            grep -v '^$' > variant_list.txt

                        n_variants=\$(wc -l < variant_list.txt)
                        echo \"Found \${n_variants} variants in SAIGE weight file\"
                        echo ''

                        # Show first few variants
                        echo 'First 5 variants:'
                        head -n5 variant_list.txt
                        echo ''

                        # Convert directly from BCF to PLINK with variant extraction
                        # PLINK2 can read BCF directly and extract specific variants
                        # Note: Do NOT use --set-all-var-ids because some indels are too long
                        # We'll keep the variant IDs as they are in the BCF (chr:pos:ref:alt)
                        echo 'Converting BCF to PLINK format (with variant extraction)...'
                        echo 'This may take a few minutes...'
                        plink2 \
                            --bcf $(basename ${input_bcf}) \
                            --extract variant_list.txt \
                            --make-bed \
                            --out ${out_prefix} \
                            --allow-extra-chr

                        echo ''
                        echo 'PLINK conversion complete!'

                        # Show file sizes
                        echo ''
                        echo 'Output files:'
                        ls -lh ${out_prefix}.*

                        # Show variant count
                        n_plink=\$(wc -l < ${out_prefix}.bim)
                        echo ''
                        echo \"PLINK files contain \${n_plink} variants\"

                        # Show sample count
                        n_samples=\$(wc -l < ${out_prefix}.fam)
                        echo \"PLINK files contain \${n_samples} samples\"

                        # Show first few variants from .bim file
                        echo ''
                        echo 'First 5 variants in PLINK .bim file:'
                        head -n5 ${out_prefix}.bim

                        # Generate FRQX file from PLINK files
                        echo ''
                        echo 'Generating allele frequency file (.frqx)...'
                        plink2 \
                            --bfile ${out_prefix} \
                            --freq counts \
                            --out ${out_prefix}

                        echo ''
                        echo 'Compressing FRQX file...'
                        gzip -f ${out_prefix}.frqx

                        echo ''
                        echo 'FRQX file created:'
                        ls -lh ${out_prefix}.frqx.gz

                        echo ''
                        echo 'Complete! \$(date)'
                    " \
                    --instance-type ${instance_type} \
                    --folder=".${out_dir}" \
                    --priority ${priority} \
                    --name "bcf_to_plink_chr${chr}_${anc}" \
                    -y

                echo "Job submitted!"

            else
                >&2 echo "ERROR: Input files not found"
                >&2 echo "  BCF exists: ${bcf_exists}"
                >&2 echo "  SAIGE weights exists: ${weights_exists}"
            fi

        else
            echo "Output already exists. Skipping..."
        fi

        echo ""

    done

done

echo "=========================================="
echo "All conversion jobs submitted!"
echo "=========================================="
