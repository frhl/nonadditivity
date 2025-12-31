#!/usr/bin/env bash
# author: frederik lassen
# note: Submit gene set-based phenotype simulation jobs

# =============================================================================
# CONFIGURATION
# =============================================================================
readonly DATE="251228"  # Today's date for output folder naming

# Source dx utilities
source "/Users/flassen/Projects/00_dx_utils/dx_utils/dx_utils.sh"

set -o errexit
set -o nounset

# Define remote directory and scripts
readonly remote_dir="wes_ko_ukbb/scripts/simulation/saige/gene_set_power"
readonly rscript_local="01_make_gene_set_phenos.R"
readonly gene_set_utils_local="_gene_set_utils.R"
readonly encoding_utils_local="_encoding_utils.R"
readonly effect_utils_local="_effect_size_utils.R"

# Upload R scripts to remote
dx mkdir -p ${remote_dir}
dx_update_remote "${remote_dir}/${rscript_local}" ${rscript_local}
dx_update_remote "${remote_dir}/${gene_set_utils_local}" ${gene_set_utils_local}
dx_update_remote "${remote_dir}/${encoding_utils_local}" ${encoding_utils_local}
dx_update_remote "${remote_dir}/${effect_utils_local}" ${effect_utils_local}

# =============================================================================
# DATA PATHS
# =============================================================================

# Covariates
readonly covars_dir="/wes_ko_ukbb/data/phenotypes/covariates"
readonly covars_file="${covars_dir}/core_covariates_with_PCs.tsv.gz"

# Output directory
readonly out_dir="/wes_ko_ukbb/data/phenotypes/simulated/gene_set_power/${DATE}"
dx mkdir -p ${out_dir}

# =============================================================================
# SIMULATION PARAMETERS
# =============================================================================

# Ancestries to process
readonly ancestries=("eur")

# Genetic architectures
readonly architectures="additive,recessive,partially_recessive_0.1,partially_recessive_0.2"

# Heritabilities (h² per variant)
readonly h2_values="0.001,0.005,0.01,0.02"

# K: Number of causal genes
readonly K_genes="10"

# N: Number of causal variants per gene
readonly N_per_gene="1"

# Variant selection strategy
readonly variant_selection="random"

# Annotations to include
readonly annotations="pLoF,damaging_missense"

# MAF bins for variant selection
readonly maf_bins="0.001-0.05"

# Minimum variants per gene
readonly min_variants_per_gene=1

# Number of replicates per configuration
readonly n_reps=2

# Random seed
readonly seed=42

# Instance configuration
readonly instance_type="mem3_ssd1_v2_x16"
readonly priority="high"

# =============================================================================
# INPUT DATA CONFIGURATION
# =============================================================================

# SAIGE weight file parameters
readonly group="spliceai=0.50_cadd=28.1_revel=0.773"
readonly af="05"
readonly min_mac="10"

# =============================================================================
# MAIN LOOP - START WITH CHR21 FOR DEBUGGING
# =============================================================================

for anc in "${ancestries[@]}"; do

    echo "=========================================="
    echo "Processing ancestry: ${anc}"
    echo "=========================================="

    # Start with chr21 for debugging
    for chr in 21; do

        echo ""
        echo "Chromosome: ${chr}"

        # SAIGE weight file
        saige_weights="/wes_ko_ukbb/data/genesets/dominance_weights/af05/UKB.wes.chr${chr}.phased.full_qc.${anc}.af${af}.popmax.variants.${group}.min_mac_${min_mac}.saige_with_${anc}_weights.txt.gz"

        # PLINK files (raw additive genotypes)
        plink_dir="/wes_ko_ukbb/data/phased/gene_set_simulation/plink"
        bed_prefix="${plink_dir}/UKB.wes.chr${chr}.phased.qc_final.${anc}.af${af}.popmax.${group}.min_mac_${min_mac}"

        # Allele count file (generated from PLINK files in same directory)
        frqx_file="${bed_prefix}.frqx.gz"

        # Output prefix
        out_prefix="${DATE}_chr${chr}_${anc}_gene_set_power"

        echo "SAIGE weights: ${saige_weights}"
        echo "PLINK files: ${bed_prefix}.{bed,bim,fam}"
        echo "FRQX file: ${frqx_file}"
        echo "Output: ${out_dir}/${out_prefix}.tsv.gz"
        echo ""

        # Check if output already exists
        if [[ $(dx_file_exists "${out_dir}/${out_prefix}.tsv.gz") -eq "0" ]]; then

            # Check if input files exist
            weights_exists=$(dx ls "${saige_weights}" 2>/dev/null | wc -l)
            bed_exists=$(dx ls "${bed_prefix}.bed" 2>/dev/null | wc -l)
            bim_exists=$(dx ls "${bed_prefix}.bim" 2>/dev/null | wc -l)
            fam_exists=$(dx ls "${bed_prefix}.fam" 2>/dev/null | wc -l)
            frqx_exists=$(dx ls "${frqx_file}" 2>/dev/null | wc -l)

            if [[ ${weights_exists} -gt 0 && ${bed_exists} -gt 0 && ${bim_exists} -gt 0 && ${fam_exists} -gt 0 && ${frqx_exists} -gt 0 ]]; then

                echo "Submitting job..."

                pheno_job_id=$(dx run app-swiss-army-knife \
                    -iimage_file="/docker/rsuite.tar.gz" \
                    -iin="${saige_weights}" \
                    -iin="${bed_prefix}.bed" \
                    -iin="${bed_prefix}.bim" \
                    -iin="${bed_prefix}.fam" \
                    -iin="${frqx_file}" \
                    -iin="${covars_file}" \
                    -iin="/${remote_dir}/${rscript_local}" \
                    -iin="/${remote_dir}/${gene_set_utils_local}" \
                    -iin="/${remote_dir}/${encoding_utils_local}" \
                    -iin="/${remote_dir}/${effect_utils_local}" \
                    -icmd="
                        Rscript /mnt/project/${remote_dir}/${rscript_local} \
                            --saige_weights $(basename ${saige_weights}) \
                            --bed_file $(basename ${bed_prefix}) \
                            --frqx_file $(basename ${frqx_file}) \
                            --covars $(basename ${covars_file}) \
                            --covar_id_col participant_eid \
                            --ancestry ${anc} \
                            --annotations ${annotations} \
                            --architectures ${architectures} \
                            --h2_total ${h2_values} \
                            --K_genes ${K_genes} \
                            --N_per_gene ${N_per_gene} \
                            --variant_selection ${variant_selection} \
                            --maf_bins ${maf_bins} \
                            --min_variants_per_gene ${min_variants_per_gene} \
                            --n_reps ${n_reps} \
                            --seed ${seed} \
                            --out ${out_prefix} &&
                        echo '!!!!$(date)'
                    " \
                    --instance-type ${instance_type} \
                    --folder=".${out_dir}" \
                    --priority ${priority} \
                    --name "gene_set_power_chr${chr}_${anc}" \
                    --brief \
                    -y)

                echo "Job submitted! Job ID: ${pheno_job_id}"

                # Save job ID for downstream use
                echo "${pheno_job_id}" > ".pheno_job_chr${chr}_${anc}.txt"

            else
                >&2 echo "ERROR: Input files not found:"
                >&2 echo "  SAIGE weights exists: ${weights_exists}"
                >&2 echo "  PLINK .bed exists: ${bed_exists}"
                >&2 echo "  PLINK .bim exists: ${bim_exists}"
                >&2 echo "  PLINK .fam exists: ${fam_exists}"
                >&2 echo "  FRQX exists: ${frqx_exists}"
            fi

        else
            echo "Output already exists. Skipping..."
        fi

        echo ""

    done

done

echo "=========================================="
echo "All jobs submitted!"
echo "=========================================="
