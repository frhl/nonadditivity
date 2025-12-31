#!/usr/bin/env bash
# author: frederik lassen
# note: Submit SAIGE step2 jobs for power phenotypes
# Tests multiple encodings: additive, recessive, dominance (nonadditive)

# =============================================================================
# CONFIGURATION - Must match 01_make_power_phenos.sh
# =============================================================================
readonly DATE="251227"  # Must match the date in 01_make_power_phenos.sh

# Define paths
readonly step0_dir="/wes_ko_ukbb/data/saige/step0/vr_20k"
readonly step1_dir="/wes_ko_ukbb/data/saige/step1/simulated/power/${DATE}"
readonly out_dir="/wes_ko_ukbb/data/saige/step2/simulated/power/${DATE}"

# Create output directory
dx mkdir -p ${out_dir}

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

readonly priority="low"
readonly instance_type="mem3_ssd1_v2_x4"

# Exome file parameters
readonly group="spliceai=0.50_cadd=28.1_revel=0.773"
readonly annotation="pLoF_damaging_missense"
readonly af="05"
readonly pp="0.90"

# Test encodings to run
# - additive: standard dosage [0, 1, 2]
# - recessive: carrier [0, 0, 1]
# - dominance: nonadditive (X^D from VCF)
readonly test_encodings=("additive" "recessive" "dominance")

# Loop through ancestries
for anc in "${ancestries[@]}"; do

    # Genotype file directory
    exome_dir="/wes_ko_ukbb/data/phased/encode_alt_qced_canonical/${anc}/${group}/vcf_plus_plink/force_chr_name"

    # Read step1 job IDs if available (for dependency chaining)
    step1_depends_arg=""
    if [[ -f ".step1_jobs_${anc}.txt" ]]; then
        # Check if any step1 jobs are still running
        any_running=false
        while IFS= read -r job_id; do
            if [[ -n "${job_id}" ]]; then
                job_state=$(dx describe "${job_id}" --json 2>/dev/null | jq -r '.state' || echo "unknown")
                if [[ "${job_state}" == "running" || "${job_state}" == "runnable" ]]; then
                    any_running=true
                    break
                fi
            fi
        done < ".step1_jobs_${anc}.txt"

        if [[ "${any_running}" == "true" ]]; then
            echo "WARNING: Some step1 jobs are still running"
            echo "You must wait for them to complete before running this script."
            echo "Run: dx wait .step1_jobs_${anc}.txt"
            echo "Or check status with: cat .step1_jobs_${anc}.txt | xargs -I {} dx describe {}"
            echo ""
            echo "Exiting..."
            exit 1
        fi

        # Build --depends-on arguments (one flag per job ID)
        while IFS= read -r job_id; do
            if [[ -n "${job_id}" ]]; then
                step1_depends_arg="${step1_depends_arg} --depends-on ${job_id}"
            fi
        done < ".step1_jobs_${anc}.txt"

        if [[ -n "${step1_depends_arg}" ]]; then
            echo "Found step1 job IDs from: .step1_jobs_${anc}.txt"
            echo "All step1 jobs completed, ready to run step2"
        fi
    fi

    echo "Processing ancestry: ${anc}"
    echo ""

    # Loop through test encodings
    for test_enc in "${test_encodings[@]}"; do

        echo "==================== Test Encoding: ${test_enc} ===================="
        echo ""

        # Determine file type and test type based on encoding
        if [[ "${test_enc}" == "dominance" ]]; then
            # Dominance uses VCF with nonadditive encoding
            use_vcf=true
            test_type="variant"
            file_mode="dominance"
            chrom_arg="chr1"  # VCF uses "chr1" prefix
        else
            # Additive and recessive use PLINK
            use_vcf=false
            test_type="variant"
            file_mode="${test_enc}"
            chrom_arg="1"  # PLINK uses no prefix
        fi

        # Genotype file prefix (merged, all chromosomes)
        exome_prefix="UKB.wes.merged.phased.qc_final.${anc}.af${af}.popmax.pp${pp}.${group}.${annotation}.${file_mode}.auto"

        # Loop through all phenotype combinations
        for arch in "${architectures[@]}"; do
            for h2 in "${h2_values[@]}"; do
                for M in "${M_values[@]}"; do
                    for rep in $(seq 1 ${n_reps}); do

                        # Create phenotype name (must match R's sprintf("%.6f", h2) format)
                        h2_str=$(printf "%.6f" ${h2} | sed 's/\.//g')  # 0.00001 -> 0.000010 -> 000010
                        pheno="arch_${arch}_h2_${h2_str}_M_${M}_rep_${rep}"

                        # Check if step1 model file exists
                        model_file="${step1_dir}/${pheno}_${anc}.rda"
                        variance_ratio="${step1_dir}/${pheno}_${anc}.varianceRatio.txt"

                        # Output prefix (no chromosome since merged file)
                        out_prefix="UKB.auto.${pheno}.${test_enc}.${anc}"

                        # Only show every 10th phenotype to reduce output
                        if [[ $((rep % 10)) -eq 1 || ${rep} -eq ${n_reps} ]]; then
                            echo "  Processing: ${pheno} (rep ${rep}/${n_reps})"
                        fi

                        # Check if output already exists
                        if [[ $(dx ls "${out_dir}/${out_prefix}.txt.gz" 2>/dev/null | wc -l) -eq 0 ]]; then

                            # Check if model file exists (skip check if we have step1 dependencies)
                            model_exists=0
                            if [[ -n "${step1_depends_arg}" ]]; then
                                # Trust that step1 will create the files
                                model_exists=1
                            else
                                # No dependencies, check if files exist now
                                model_exists=$(dx ls "${model_file}" 2>/dev/null | wc -l)
                            fi

                            if [[ ${model_exists} -gt 0 ]]; then

                                if [[ "${use_vcf}" == "true" ]]; then
                                    # Use VCF for dominance test
                                    vcf_file="${exome_dir}/${exome_prefix}.vcf.gz"
                                    vcf_index="${exome_dir}/${exome_prefix}.vcf.gz.csi"

                                    # Check if VCF exists
                                    if [[ $(dx ls "${vcf_file}" 2>/dev/null | wc -l) -gt 0 ]]; then

                                        dx run saige-universal-step-2 \
                                            -i chrom="${chrom_arg}" \
                                            -i output_prefix="${out_prefix}" \
                                            -i model_file="${model_file}" \
                                            -i variance_ratio="${variance_ratio}" \
                                            -i test_type="${test_type}" \
                                            -i vcf_file="${vcf_file}" \
                                            -i vcf_index_file="${vcf_index}" \
                                            -i GRM="${step0_dir}/UKB.array.${anc}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx" \
                                            -i GRM_samples="${step0_dir}/UKB.array.${anc}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt" \
                                            --instance-type ${instance_type} \
                                            --priority ${priority} \
                                            --destination "${out_dir}" \
                                            ${step1_depends_arg} \
                                            -y \
                                            --name "step2_${pheno}_${test_enc}_${anc}" \
                                            --brief

                                    else
                                        if [[ ${rep} -eq 1 ]]; then
                                            >&2 echo "    WARNING: VCF not found: ${vcf_file}"
                                        fi
                                    fi

                                else
                                    # Use PLINK for additive/recessive test
                                    plink_bed="${exome_dir}/${exome_prefix}.bed"
                                    plink_bim="${exome_dir}/${exome_prefix}.bim"
                                    plink_fam="${exome_dir}/${exome_prefix}.fam"

                                    # Check if PLINK files exist
                                    if [[ $(dx ls "${plink_bed}" 2>/dev/null | wc -l) -gt 0 ]]; then

                                        dx run saige-universal-step-2 \
                                            -i chrom="${chrom_arg}" \
                                            -i output_prefix="${out_prefix}" \
                                            -i model_file="${model_file}" \
                                            -i variance_ratio="${variance_ratio}" \
                                            -i test_type="${test_type}" \
                                            -i exome_bed="${plink_bed}" \
                                            -i exome_bim="${plink_bim}" \
                                            -i exome_fam="${plink_fam}" \
                                            -i GRM="${step0_dir}/UKB.array.${anc}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx" \
                                            -i GRM_samples="${step0_dir}/UKB.array.${anc}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt" \
                                            --instance-type ${instance_type} \
                                            --priority ${priority} \
                                            --destination "${out_dir}" \
                                            ${step1_depends_arg} \
                                            -y \
                                            --name "step2_${pheno}_${test_enc}_${anc}" \
                                            --brief

                                    else
                                        if [[ ${rep} -eq 1 ]]; then
                                            >&2 echo "    WARNING: PLINK files not found: ${plink_bed}"
                                        fi
                                    fi
                                fi

                            else
                                if [[ ${rep} -eq 1 ]]; then
                                    >&2 echo "    WARNING: Model file not found: ${model_file}"
                                fi
                            fi

                        fi

                    done  # reps
                done  # M values
            done  # h2 values
        done  # architectures

    done  # test encodings

done  # ancestries

echo ""
echo "Done submitting step2 jobs!"
