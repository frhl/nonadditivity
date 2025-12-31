#!/usr/bin/env bash
# author: frederik lassen
# note: this script submits SAIGE step2 jobs for simulated phenotypes

set -o errexit
set -o nounset

source "/Users/flassen/Projects/00_dx_utils/dx_utils/dx_utils.sh"

# Define paths
readonly step0_dir="/wes_ko_ukbb/data/saige/step0/vr_20k"
readonly step1_dir="/wes_ko_ukbb/data/saige/step1/simulated/revisions/051025"
readonly out_dir="/wes_ko_ukbb/data/saige/step2/simulated/revisions/051025"

# Create output directory
dx mkdir -p ${out_dir}

# Define parameters
#readonly ancestries=("eur")
#readonly ancestries=("eur" "sas" "afr" "eas")
readonly ancestries=("sas" "afr" "eas")
#readonly heritabilities=("0.01")
readonly heritabilities=("0.01" "0.1" "0.2" "0.3" "0.5")
readonly n_reps=4
#readonly n_reps=4
readonly test_type="variant"
readonly priority="high"
readonly instance_type="mem3_ssd1_v2_x4"

# Exome file parameters
readonly group="spliceai=0.50_cadd=28.1_revel=0.773"
readonly annotations=("pLoF_damaging_missense")
#readonly annotations=("pLoF" "damaging_missense" "synonymous")
readonly modes=("dominance" "additive" "recessive")
#readonly modes=("recessive" "additive" "dominance")
readonly af="05"
readonly pp="0.90"

# Loop through ancestries
for anc in "${ancestries[@]}"; do
  exome_dir="/wes_ko_ukbb/data/phased/encode_alt_qced_canonical/${anc}/${group}/vcf_plus_plink/force_chr_name"

  echo "=========================================="
  echo "Processing ancestry: ${anc}"
  echo "=========================================="

  # Loop through heritabilities
  for h2 in "${heritabilities[@]}"; do
    # Loop through replicates
    for rep in $(seq 1 ${n_reps}); do
      pheno="p_${h2}_continuous_${rep}"

      # Check if step1 model file exists
      model_file="${step1_dir}/${pheno}_${anc}.rda"
      variance_ratio="${step1_dir}/${pheno}_${anc}.varianceRatio.txt"

      # Validate step1 files exist
      model_exists=$(dx_file_exists "${model_file}")
      variance_exists=$(dx_file_exists "${variance_ratio}")

      if [[ ${model_exists} -eq 0 ]]; then
        >&2 echo "ERROR: Model file not found: ${model_file}. Skipping ${pheno}."
        continue
      fi

      if [[ ${variance_exists} -eq 0 ]]; then
        >&2 echo "ERROR: Variance ratio file not found: ${variance_ratio}. Skipping ${pheno}."
        continue
      fi

      # Loop through annotations
      for anno in "${annotations[@]}"; do
        # Loop through modes
        for mode in "${modes[@]}"; do

          # Construct exome prefix and output prefix
          exome_prefix="${exome_dir}/UKB.wes.merged.phased.qc_final.${anc}.af${af}.popmax.pp${pp}.${group}.${anno}.${mode}.auto"
          out_prefix="UKB.auto.qc_final.${pheno}.${mode}.${anc}.af${af}.pp${pp}.${anno}"
          output_file="${out_dir}/${out_prefix}.txt.gz"

          # Check if output already exists and is non-empty
          output_exists=$(dx_file_exists "${output_file}")

          if [[ ${output_exists} -eq 1 ]]; then
            # File exists, check if it's empty (could be from failed job)
            output_is_empty=$(dx_is_empty "${output_file}")
            if [[ ${output_is_empty} -eq 1 ]]; then
              echo "  Removing empty/incomplete file: ${output_file}"
              dx rm "${output_file}"
            else
              echo "Output already exists (non-empty): ${out_prefix}.txt.gz. Skipping."
              continue
            fi
          fi

          # Validate input files exist before submitting
          if [[ "${mode}" == "dominance" ]]; then
            # Check VCF files for dominance mode
            vcf_file="${exome_prefix}.vcf.gz"
            vcf_index="${exome_prefix}.vcf.gz.csi"

            vcf_exists=$(dx_file_exists "${vcf_file}")
            vcf_index_exists=$(dx_file_exists "${vcf_index}")

            if [[ ${vcf_exists} -eq 0 ]]; then
              >&2 echo "ERROR: VCF file not found: ${vcf_file}"
              >&2 echo "       Skipping ${out_prefix}"
              continue
            fi

            if [[ ${vcf_index_exists} -eq 0 ]]; then
              >&2 echo "ERROR: VCF index not found: ${vcf_index}"
              >&2 echo "       Skipping ${out_prefix}"
              continue
            fi
          else
            # Check PLINK files for additive/recessive modes
            bed_file="${exome_prefix}.bed"
            bim_file="${exome_prefix}.bim"
            fam_file="${exome_prefix}.fam"

            bed_exists=$(dx_file_exists "${bed_file}")
            bim_exists=$(dx_file_exists "${bim_file}")
            fam_exists=$(dx_file_exists "${fam_file}")

            if [[ ${bed_exists} -eq 0 ]]; then
              >&2 echo "ERROR: PLINK bed file not found: ${bed_file}"
              >&2 echo "       Skipping ${out_prefix}"
              continue
            fi

            if [[ ${bim_exists} -eq 0 ]]; then
              >&2 echo "ERROR: PLINK bim file not found: ${bim_file}"
              >&2 echo "       Skipping ${out_prefix}"
              continue
            fi

            if [[ ${fam_exists} -eq 0 ]]; then
              >&2 echo "ERROR: PLINK fam file not found: ${fam_file}"
              >&2 echo "       Skipping ${out_prefix}"
              continue
            fi
          fi

          # All checks passed, submit the job
          echo "Submitting: ${out_prefix}"

          if [[ "${mode}" == "dominance" ]]; then
            # Use VCF for dominance mode
            dx run saige-universal-step-2 \
              -i chrom="chr1" \
              -i output_prefix="${out_prefix}" \
              -i model_file="${model_file}" \
              -i variance_ratio="${variance_ratio}" \
              -i test_type=${test_type} \
              -i vcf_file="${exome_prefix}.vcf.gz" \
              -i vcf_index_file="${exome_prefix}.vcf.gz.csi" \
              -i GRM="${step0_dir}/UKB.array.${anc}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx" \
              -i GRM_samples="${step0_dir}/UKB.array.${anc}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt" \
              --instance-type ${instance_type} \
              --priority ${priority} \
              --destination "${out_dir}" \
              -y \
              --name "step2_${pheno}_${anc}_${mode}_${anno}"
          else
            # Use PLINK for additive/recessive modes
            dx run saige-universal-step-2 \
              -i chrom="1" \
              -i output_prefix="${out_prefix}" \
              -i model_file="${model_file}" \
              -i variance_ratio="${variance_ratio}" \
              -i test_type=${test_type} \
              -i exome_bed="${exome_prefix}.bed" \
              -i exome_bim="${exome_prefix}.bim" \
              -i exome_fam="${exome_prefix}.fam" \
              -i GRM="${step0_dir}/UKB.array.${anc}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx" \
              -i GRM_samples="${step0_dir}/UKB.array.${anc}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt" \
              --instance-type ${instance_type} \
              --priority ${priority} \
              --destination "${out_dir}" \
              -y \
              --name "step2_${pheno}_${anc}_${mode}_${anno}"
          fi

        done
      done
    done
  done
done

echo "Done submitting step2 jobs!"
