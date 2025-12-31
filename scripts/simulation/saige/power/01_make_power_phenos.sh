#!/usr/bin/env bash
# author: frederik lassen
# note: Submit phenotype simulation jobs for power analysis

# =============================================================================
# CONFIGURATION - Set date and parameters here
# =============================================================================
readonly DATE="251227b"  # Change this to set a consistent date across all files/folders

# Source dx utilities
source "/Users/flassen/Projects/00_dx_utils/dx_utils/dx_utils.sh"

# Define remote directory and scripts
readonly remote_dir="wes_ko_ukbb/scripts/simulation/saige/power"
readonly rscript_local="01_make_power_phenos.R"
readonly encoding_utils_local="_encoding_utils.R"
readonly effect_utils_local="_effect_size_utils.R"

# Upload R scripts to remote
dx mkdir -p ${remote_dir}
dx_update_remote "${remote_dir}/${rscript_local}" ${rscript_local}
dx_update_remote "${remote_dir}/${encoding_utils_local}" ${encoding_utils_local}
dx_update_remote "${remote_dir}/${effect_utils_local}" ${effect_utils_local}

# Define paths
readonly covars_dir="/wes_ko_ukbb/data/phenotypes/covariates"
readonly covars_file="${covars_dir}/core_covariates_with_PCs.tsv.gz"
readonly out_dir="/wes_ko_ukbb/data/phenotypes/simulated/power/${DATE}"

# Create output directory
dx mkdir -p ${out_dir}

# Define parameters
readonly ancestries=("afr")
#readonly ancestries=("eur" "sas" "afr" "eas")
# Note: For partially recessive, specify heterozygote value after underscore (e.g., partially_recessive_0.1)
#readonly architectures="additive,partially_recessive_0.1,partially_recessive_0.2,partially_recessive_0.3,recessive,overdominant"
readonly architectures="additive,partially_recessive_0.05,partially_recessive_0.1,partially_recessive_0.2,partially_recessive_0.3,recessive"
#readonly h2_values="0.00001,0.00005,0.000075,0.0001,0.0002"
#readonly h2_values="0.0005,0.001,0.002,0.003,0.004,0.005,0.0075,0.01,0.015,0.02"
readonly h2_values="0.0005,0.001,0.01,0.05,0.1,0.2,0.3"
#readonly h2_values="0.0005,0.001,0.005,0.0075,0.01,0.02"
readonly polygenicity="25"  # M values
readonly n_reps=1
readonly seed=42
readonly instance_type="mem1_ssd1_v2_x16"
readonly priority="high"

# Genotype file parameters (merged, all chromosomes, ~18k variants)
# Note: Simulation script will filter to variants with ≥5 hom alt AND ≥5 heterozygotes
readonly group="spliceai=0.50_cadd=28.1_revel=0.773"
readonly annotation="pLoF_damaging_missense"
readonly af="05"
readonly pp="0.90"

# MAF bins for variant selection (optional)
# Ultra-rare variants (MAF < 0.001) often have no heterozygotes - use MAF ≥ 0.001
readonly maf_bins="0.01-0.05"  # Options: "all", "0.001-0.01", "0.01-0.05", etc.

# Loop through ancestries
for anc in "${ancestries[@]}"; do

    # Path to merged PLINK files (all chromosomes combined)
    # Use additive encoding to preserve raw genotypes [0, 1, 2] for simulation
    bed_prefix="/wes_ko_ukbb/data/phased/encode_alt_qced_canonical/${anc}/${group}/vcf_plus_plink/force_chr_name/UKB.wes.merged.phased.qc_final.${anc}.af${af}.popmax.pp${pp}.${group}.${annotation}.additive.auto"

    # Output prefix
    out_prefix="${DATE}_${anc}_power_phenos"

    echo "=========================================="
    echo "Submitting phenotype simulation for: ${anc}"
    echo "=========================================="
    echo "Using genotypes from: ${bed_prefix}"
    echo "Output: ${out_dir}/${out_prefix}.tsv.gz"
    echo ""

    # Check if output already exists
    if [[ $(dx_file_exists "${out_dir}/${out_prefix}.tsv.gz") -eq "0" ]]; then

        # Check if PLINK files exist
        bed_exists=$(dx ls "${bed_prefix}.bed" 2>/dev/null | wc -l)
        bim_exists=$(dx ls "${bed_prefix}.bim" 2>/dev/null | wc -l)
        fam_exists=$(dx ls "${bed_prefix}.fam" 2>/dev/null | wc -l)

        if [[ ${bed_exists} -gt 0 && ${bim_exists} -gt 0 && ${fam_exists} -gt 0 ]]; then

            echo "Submitting job..."

            pheno_job_id=$(dx run app-swiss-army-knife \
                -iimage_file="/docker/rsuite.tar.gz" \
                -iin="${bed_prefix}.bed" \
                -iin="${bed_prefix}.bim" \
                -iin="${bed_prefix}.fam" \
                -iin="${covars_file}" \
                -iin="/${remote_dir}/${rscript_local}" \
                -iin="/${remote_dir}/${encoding_utils_local}" \
                -iin="/${remote_dir}/${effect_utils_local}" \
                -icmd="
                    Rscript /mnt/project/${remote_dir}/${rscript_local} \
                        --bed_file $(basename ${bed_prefix}) \
                        --covars $(basename ${covars_file}) \
                        --covar_id_col participant_eid \
                        --ancestry ${anc} \
                        --architectures ${architectures} \
                        --h2_total ${h2_values} \
                        --polygenicity ${polygenicity} \
                        --maf_bins ${maf_bins} \
                        --n_reps ${n_reps} \
                        --seed ${seed} \
                        --out ${out_prefix} &&
                    echo '!!!!$(date)'
                " \
                --instance-type ${instance_type} \
                --folder=".${out_dir}" \
                --priority ${priority} \
                --name "power_phenos_${anc}" \
                --brief \
                -y)

            echo "Job submitted successfully! Job ID: ${pheno_job_id}"

            # Save job ID for step1 to use
            echo "${pheno_job_id}" > ".pheno_job_${anc}.txt"

        else
            >&2 echo "ERROR: PLINK files not found:"
            >&2 echo "  .bed exists: ${bed_exists}"
            >&2 echo "  .bim exists: ${bim_exists}"
            >&2 echo "  .fam exists: ${fam_exists}"
        fi

    else
        echo "Output already exists: ${out_dir}/${out_prefix}.tsv.gz"
        echo "Skipping..."
    fi

    echo ""

done

echo "=========================================="
echo "All submission jobs completed!"
echo "=========================================="
