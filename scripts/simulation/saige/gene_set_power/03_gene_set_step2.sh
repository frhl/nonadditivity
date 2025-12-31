#!/usr/bin/env bash
# author: frederik lassen
# note: Run SAIGE Step2 (gene set burden tests) for simulated phenotypes

source "/Users/flassen/Projects/00_dx_utils/dx_utils/dx_utils.sh"

set -o errexit
set -o nounset

# =============================================================================
# CONFIGURATION
# =============================================================================
readonly DATE="251228"

# Directories
readonly step0_dir="/wes_ko_ukbb/data/saige/step0/vr_20k"
readonly step1_dir="/wes_ko_ukbb/data/saige/step1/simulated/gene_set_power/${DATE}"
readonly step2_dir="/wes_ko_ukbb/data/saige/step2/simulated/gene_set_power/${DATE}"

# Create output directory
dx mkdir -p ${step2_dir}

# Instance configuration
readonly instance_type="mem3_ssd1_v2_x2"
readonly priority="high"
readonly test_type="group"

# =============================================================================
# DATA PARAMETERS
# =============================================================================

readonly ancestries=("eur")
readonly group="spliceai=0.50_cadd=28.1_revel=0.773"
readonly af="05"
readonly min_mac="10"
readonly annotations="pLoF:damaging_missense"

# =============================================================================
# MAIN LOOP
# =============================================================================

for anc in "${ancestries[@]}"; do

    echo "=========================================="
    echo "Processing ancestry: ${anc}"
    echo "=========================================="

    # Start with chr21 for debugging
    for chr in 21; do

        echo ""
        echo "Chromosome: ${chr}"

        # SAIGE group file
        saige_group_dir="/wes_ko_ukbb/data/genesets/dominance_weights/af05"
        group_file="${saige_group_dir}/UKB.wes.chr${chr}.phased.full_qc.${anc}.af${af}.popmax.variants.${group}.min_mac_${min_mac}.saige_with_${anc}_weights.txt.gz"

        # Dominance-encoded VCF
        vcf_dir="/wes_ko_ukbb/data/phased/encoded_dominance_gt/group_dominance/gene_map"
        vcf_prefix="UKB.wes.chr${chr}.phased.full_qc.${anc}.af${af}.popmax.variants.group_dominance_scaling.gene_map"
        vcf_file="${vcf_dir}/${vcf_prefix}.vcf.gz"
        vcf_index="${vcf_dir}/${vcf_prefix}.vcf.gz.csi"

        echo "Group file: ${group_file}"
        echo "VCF file: ${vcf_file}"

        # Check if files exist
        group_exists=$(dx ls "${group_file}" 2>/dev/null | wc -l)
        vcf_exists=$(dx ls "${vcf_file}" 2>/dev/null | wc -l)

        if [[ ${group_exists} -eq 0 || ${vcf_exists} -eq 0 ]]; then
            >&2 echo "ERROR: Required files not found"
            >&2 echo "  Group file exists: ${group_exists}"
            >&2 echo "  VCF exists: ${vcf_exists}"
            continue
        fi

        # Get list of phenotypes from step1 directory
        echo "Finding phenotypes in Step1 directory..."
        step1_files=$(dx ls "${step1_dir}" | grep "\.rda$" | grep "_${anc}\.rda$" || true)

        if [[ -z "${step1_files}" ]]; then
            >&2 echo "ERROR: No Step1 output files found in ${step1_dir}"
            continue
        fi

        # Count phenotypes
        n_phenos=$(echo "${step1_files}" | wc -l | tr -d ' ')
        echo "Found ${n_phenos} phenotypes to process"

        # Process each phenotype
        idx=0
        for step1_file in ${step1_files}; do
            idx=$((idx + 1))

            # Extract phenotype name
            # Format: arch_X_h2_Y_K_Z_N_W_rep_R_eur.rda
            pheno=$(echo "${step1_file}" | sed "s/_${anc}\.rda$//")

            echo ""
            echo "[${idx}/${n_phenos}] Processing: ${pheno}"

            # Step1 files
            model_file="${step1_dir}/${pheno}_${anc}.rda"
            variance_ratio="${step1_dir}/${pheno}_${anc}.varianceRatio.txt"

            # Output prefix
            out_prefix="UKB.auto.chr${chr}.${pheno}.${anc}"

            # Check if output already exists
            if [[ $(dx_file_exists "${step2_dir}/${out_prefix}.txt.gz") -eq "0" ]]; then

                # Check if model file exists
                if [[ $(dx_file_exists "${model_file}") -gt 0 ]]; then

                    echo "  Submitting SAIGE Step2 job..."

                    dx run saige-universal-step-2-group \
                        -i chrom="chr${chr}" \
                        -i output_prefix="${out_prefix}" \
                        -i model_file="${model_file}" \
                        -i variance_ratio="${variance_ratio}" \
                        -i group_file="${group_file}" \
                        -i annotations="${annotations}" \
                        -i test_type=${test_type} \
                        -i vcf_file="${vcf_file}" \
                        -i vcf_index_file="${vcf_index}" \
                        -i GRM="${step0_dir}/UKB.array.${anc}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx" \
                        -i GRM_samples="${step0_dir}/UKB.array.${anc}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt" \
                        --instance-type ${instance_type} \
                        --priority ${priority} \
                        --destination "${step2_dir}" \
                        -y \
                        --name "step2_${pheno}_chr${chr}"

                else
                    >&2 echo "  WARNING: Model file not found: ${model_file}"
                fi

            else
                echo "  Output already exists. Skipping."
            fi

            # Progress update every 10 phenotypes
            if [[ $((idx % 10)) -eq 0 ]]; then
                echo "  Progress: ${idx}/${n_phenos} phenotypes submitted"
            fi

        done

        echo ""
        echo "Finished submitting Step2 jobs for chr${chr}"

    done

done

echo ""
echo "=========================================="
echo "All SAIGE Step2 jobs submitted!"
echo "=========================================="
