#!/usr/bin/env bash
# author: frederik lassen
# note: this script submits SAIGE step2 set-based tests for simulated phenotypes
# This matches the regenie set-based test configuration using the same VCF files and step1 outputs

# Define paths
readonly step0_dir="/wes_ko_ukbb/data/saige/step0/vr_20k"
readonly step1_dir="/wes_ko_ukbb/data/saige/step1/simulated/revisions/051025"
readonly out_dir="/wes_ko_ukbb/data/saige/step2_set/simulated/revisions/051025"

# Create output directory
dx mkdir -p ${out_dir}

# Define parameters
readonly ancestries=("eur")
#readonly ancestries=("eur" "sas" "afr" "eas")
readonly heritabilities=("0.01" "0.1" "0.2" "0.3" "0.5")
#readonly heritabilities=("0.01")
readonly n_reps=4
#readonly n_reps=1
readonly test_type="group"
readonly priority="high"
readonly instance_type="mem3_ssd1_v2_x4"

# Set-based test parameters
readonly group="spliceai=0.50_cadd=28.1_revel=0.773"
readonly af="05"
readonly pp="0.90"

# SAIGE group file directory (NEW: with min homozygous carrier filter)
readonly min_mac="10"
readonly min_hom="2"
readonly saige_group_dir="/wes_ko_ukbb/data/genesets/simulated/set_based/saige_weights"

# VCF directory (NEW: from simulation workflow with dominance encoding)
readonly vcf_dir="/wes_ko_ukbb/data/phased/simulated/set_based/dominance_encoded"

# Annotation specification (matches regenie configuration)
readonly annotations="pLoF:damaging_missense"

# Loop through ancestries
for anc in "${ancestries[@]}"; do

  echo "Processing ancestry: ${anc}"

  # Loop through chromosomes
  #for chr in {1..22}; do
  for chr in {20..22}; do

    echo "  Processing chromosome: ${chr}"

    # Paths to VCF files (NEW: from simulation workflow)
    vcf_prefix="UKB.wes.chr${chr}.phased.qc_final.${anc}.group_dominance_scaling.gene_map"
    vcf_file="${vcf_dir}/${vcf_prefix}.vcf.gz"
    vcf_index="${vcf_dir}/${vcf_prefix}.vcf.gz.csi"

    # SAIGE group file (NEW: with min homozygous carrier filter)
    group_file="${saige_group_dir}/UKB.wes.chr${chr}.phased.qc_final.${anc}.af${af}.popmax.${group}.min_mac_${min_mac}.min_hom_${min_hom}.saige_weights.txt.gz"

    # Loop through heritabilities
    for h2 in "${heritabilities[@]}"; do
      # Loop through replicates
      for rep in $(seq 1 ${n_reps}); do
        pheno="p_${h2}_continuous_${rep}"

        # Check if step1 model file exists
        model_file="${step1_dir}/${pheno}_${anc}.rda"
        variance_ratio="${step1_dir}/${pheno}_${anc}.varianceRatio.txt"

        # Output prefix
        out_prefix="UKB.auto.chr${chr}.${pheno}.${anc}"

        echo "    Submitting: ${out_prefix}"

        # Check if output already exists
        if [[ $(dx ls "${out_dir}/${out_prefix}.txt.gz" 2>/dev/null | wc -l) -eq 0 ]]; then
          # Check if model file exists
          if [[ $(dx ls "${model_file}" 2>/dev/null | wc -l) -gt 0 ]]; then

            # Check if VCF and group files exist
            vcf_exists=$(dx ls "${vcf_file}" 2>/dev/null | wc -l)
            group_exists=$(dx ls "${group_file}" 2>/dev/null | wc -l)

            if [[ ${vcf_exists} -gt 0 && ${group_exists} -gt 0 ]]; then

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
                --destination "${out_dir}" \
                -y \
                --name "saige_step2_set_${pheno}_${anc}_chr${chr}"

            else
              >&2 echo "    WARNING: Not all required files exist for ${out_prefix}"
              >&2 echo "      VCF exists: ${vcf_exists}, Group file exists: ${group_exists}"
            fi

          else
            >&2 echo "    Model file not found: ${model_file}. Skipping."
          fi
        else
          >&2 echo "    Output already exists: ${out_prefix}.txt.gz. Skipping."
        fi

      done
    done
  done
done

echo "Done submitting SAIGE step2 set-based test jobs!"
