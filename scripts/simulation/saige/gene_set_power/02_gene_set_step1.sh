#!/usr/bin/env bash
# author: frederik lassen
# note: Run SAIGE Step1 (fit null models) for gene set-based simulated phenotypes

source "/Users/flassen/Projects/00_dx_utils/dx_utils/dx_utils.sh"

set -o errexit
set -o nounset

# =============================================================================
# CONFIGURATION
# =============================================================================
readonly DATE="251228"

# Input/output directories
readonly pheno_dir="/wes_ko_ukbb/data/phenotypes/simulated/gene_set_power/${DATE}"
readonly step0_dir="/wes_ko_ukbb/data/saige/step0/vr_20k"
readonly step1_dir="/wes_ko_ukbb/data/saige/step1/simulated/gene_set_power/${DATE}"

# Create output directory
dx mkdir -p ${step1_dir}

# Instance configuration
readonly instance_type="mem1_ssd1_v2_x4"
readonly priority="high"

# =============================================================================
# PARAMETERS
# =============================================================================

readonly ancestries=("eur")

# =============================================================================
# MAIN LOOP
# =============================================================================

for anc in "${ancestries[@]}"; do

    echo "=========================================="
    echo "Processing ancestry: ${anc}"
    echo "=========================================="

    # Start with chr21 for debugging
    for chr in 21; do

        # Phenotype file from step 0
        pheno_prefix="${DATE}_chr${chr}_${anc}_gene_set_power"
        pheno_file="${pheno_dir}/${pheno_prefix}.tsv.gz"

        echo "Phenotype file: ${pheno_file}"

        # Check if phenotype file exists
        if [[ $(dx_file_exists "${pheno_file}") -eq "0" ]]; then
            >&2 echo "ERROR: Phenotype file not found: ${pheno_file}"
            continue
        fi

        # Get list of phenotypes from file
        echo "Extracting phenotype column names..."
        pheno_cols=$(dx cat "${pheno_file}" | gunzip | head -n1 | tr '\t' '\n' | grep "^arch_" || true)

        if [[ -z "${pheno_cols}" ]]; then
            >&2 echo "ERROR: No phenotypes found in file"
            continue
        fi

        # Count phenotypes
        n_phenos=$(echo "${pheno_cols}" | wc -l | tr -d ' ')
        echo "Found ${n_phenos} phenotypes to process"

        # Process each phenotype
        idx=0
        for pheno in ${pheno_cols}; do
            idx=$((idx + 1))

            echo ""
            echo "[${idx}/${n_phenos}] Processing: ${pheno}"

            # Check if output already exists
            if [[ $(dx_file_exists "${step1_dir}/${pheno}_${anc}.rda") -eq "0" ]]; then

                echo "  Submitting SAIGE Step1 job..."

                dx run saige-universal-step-1 \
                    -i phenotype_file="${pheno_file}" \
                    -i phenotype_name="${pheno}" \
                    -i phenotype_id_col="eid" \
                    -i sample_file="${step0_dir}/UKB.array.${anc}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt" \
                    -i GRM="${step0_dir}/UKB.array.${anc}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx" \
                    -i output_prefix="${pheno}_${anc}" \
                    -i trait_type="quantitative" \
                    -i inv_normalize="TRUE" \
                    --instance-type ${instance_type} \
                    --priority ${priority} \
                    --destination "${step1_dir}" \
                    --name "step1_${pheno}_${anc}" \
                    -y

            else
                echo "  Output already exists. Skipping."
            fi

            # Progress update every 10 phenotypes
            if [[ $((idx % 10)) -eq 0 ]]; then
                echo "  Progress: ${idx}/${n_phenos} phenotypes submitted"
            fi

        done

        echo ""
        echo "Finished submitting Step1 jobs for chr${chr}"

    done

done

echo ""
echo "=========================================="
echo "All SAIGE Step1 jobs submitted!"
echo "=========================================="
