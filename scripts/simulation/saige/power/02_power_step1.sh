#!/usr/bin/env bash
# author: frederik lassen
# note: Submit SAIGE step1 jobs for power phenotypes

set -euo pipefail

# =============================================================================
# CONFIGURATION - Must match 01_make_power_phenos.sh
# =============================================================================
readonly DATE="251227b"  # Must match the date in 01_make_power_phenos.sh

# Define paths
readonly pheno_dir="/wes_ko_ukbb/data/phenotypes/simulated/power/${DATE}"
readonly step0_location="/wes_ko_ukbb/data/saige/step0/vr_20k/"
readonly destination="/wes_ko_ukbb/data/saige/step1/simulated/power/${DATE}/"

# Create output directory
dx mkdir -p ${destination}

# Define parameters (must match simulation script settings)
readonly ancestries=("afr")
#readonly ancestries=("eur" "sas" "afr" "eas")
readonly architectures=("additive" "partially_recessive_0.05" "partially_recessive_0.1" "partially_recessive_0.2" "partially_recessive_0.3" "recessive")
# AFR heritabilities
readonly h2_values=("0.0005" "0.001" "0.01" "0.05" "0.1" "0.2" "0.3")
# EUR heritabilities
#readonly h2_values=("0.0005" "0.001" "0.002" "0.003" "0.004" "0.005" "0.0075" "0.01" "0.015" "0.02")
readonly M_values=("100")
readonly n_reps=1

readonly covariates="age_at_recruitment,age2,age_sex,age2_sex,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
readonly categorical_covariates="sex"
readonly instance_type="mem3_ssd1_v2_x8"
readonly priority="low"  # Changed from "low" to speed up execution
readonly trait_type="quantitative"

# Loop through ancestries
for anc in "${ancestries[@]}"; do

    pheno_file="${pheno_dir}/${DATE}_${anc}_power_phenos.tsv.gz"
    grm_file="${step0_location}/UKB.array.${anc}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
    grm_samples="${step0_location}/UKB.array.${anc}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
    sample_ids_file="/wes_ko_ukbb/data/samples/UKB.wes.qced.${anc}.samples"

    # Variance ratio PLINK files
    plink_bed="${step0_location}/UKB.array.${anc}.plink_for_var_ratio.bed"
    plink_bim="${step0_location}/UKB.array.${anc}.plink_for_var_ratio.bim"
    plink_fam="${step0_location}/UKB.array.${anc}.plink_for_var_ratio.fam"

    # Read pheno job ID if available (for dependency chaining)
    pheno_depends_arg=""
    pheno_job_running=false

    if [[ -f ".pheno_job_${anc}.txt" ]]; then
        pheno_job_id=$(cat ".pheno_job_${anc}.txt")
        echo "Found phenotype job ID: ${pheno_job_id}"

        # Check if phenotype file already exists
        if [[ $(dx ls "${pheno_file}" 2>/dev/null | wc -l) -gt 0 ]]; then
            echo "Phenotype file already exists"
        else
            # Check job status
            job_state=$(dx describe "${pheno_job_id}" --json 2>/dev/null | jq -r '.state' || echo "unknown")

            if [[ "${job_state}" == "done" ]]; then
                echo "Phenotype job completed, file should be available"
            elif [[ "${job_state}" == "running" || "${job_state}" == "runnable" ]]; then
                echo "WARNING: Phenotype job is still ${job_state}"
                echo "You must wait for it to complete before running this script."
                echo "Run: dx wait ${pheno_job_id}"
                echo ""
                echo "Exiting..."
                exit 1
            elif [[ "${job_state}" == "failed" ]]; then
                echo "ERROR: Phenotype job failed!"
                exit 1
            else
                # Job state is unknown, file doesn't exist, try to proceed with dependency
                pheno_depends_arg="--depends-on ${pheno_job_id}"
                echo "Step1 jobs will wait for phenotype generation (job state: ${job_state})"
            fi
        fi
    fi

    # Initialize step1 job IDs file
    > ".step1_jobs_${anc}.txt"

    echo "Processing ancestry: ${anc}"
    echo "Using phenotype file: ${pheno_file}"

    # Loop through all phenotype combinations
    for arch in "${architectures[@]}"; do
        for h2 in "${h2_values[@]}"; do
            for M in "${M_values[@]}"; do
                for rep in $(seq 1 ${n_reps}); do

                    # Create phenotype name (matching simulation script naming)
                    # Must match R's sprintf("%.6f", h2) format (always 6 decimal places)
                    h2_str=$(printf "%.6f" ${h2} | sed 's/\.//g')  # 0.00001 -> 0.000010 -> 000010
                    pheno="arch_${arch}_h2_${h2_str}_M_${M}_rep_${rep}"

                    # Output prefix
                    out_prefix="${pheno}_${anc}"

                    # Check if output already exists
                    if [[ $(dx ls "${destination}/${out_prefix}.rda" 2>/dev/null | wc -l) -eq 0 ]]; then

                        echo "  Submitting: ${pheno}"

                        step1_job_id=$(dx run saige-universal-step-1 \
                            -i output_prefix="out/${out_prefix}" \
                            -i sample_ids="${sample_ids_file}" \
                            -i genotype_bed="${plink_bed}" \
                            -i genotype_bim="${plink_bim}" \
                            -i genotype_fam="${plink_fam}" \
                            -i pheno_list="${pheno_file}" \
                            -i pheno="${pheno}" \
                            -i GRM="${grm_file}" \
                            -i GRM_samples="${grm_samples}" \
                            -i covariates="${covariates}" \
                            -i categorical_covariates="${categorical_covariates}" \
                            -i trait_type="${trait_type}" \
                            --instance-type ${instance_type} \
                            --priority ${priority} \
                            --destination "${destination}" \
                            ${pheno_depends_arg} \
                            --brief \
                            -y \
                            --name "step1_${pheno}_${anc}")

                        # Save job ID for step2
                        echo "${step1_job_id}" >> ".step1_jobs_${anc}.txt"

                    else
                        echo "  Output already exists: ${out_prefix}.rda. Skipping."
                    fi

                done
            done
        done
    done

done

echo "Done submitting step1 jobs!"
