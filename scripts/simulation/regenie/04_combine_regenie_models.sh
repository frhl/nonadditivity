#!/usr/bin/env bash
# author: frederik lassen
# note: this script combines additive, recessive, and dominance REGENIE results into a single file

readonly rscript_local="_combine_regenie_models.R"
readonly rscript_remote="wes_ko_ukbb/scripts/simulation/regenie/_combine_regenie_models.R"

# Upload R script to remote
source "/Users/flassen/Projects/00_dx_utils/dx_utils/dx_utils.sh"
dx mkdir -p wes_ko_ukbb/scripts/simulation/regenie
dx_update_remote ${rscript_remote} ${rscript_local}

# Define paths
readonly step2_dir="/wes_ko_ukbb/data/regenie/step2/simulated/revisions/051025"
readonly out_dir="/wes_ko_ukbb/data/regenie/step2/simulated/revisions/051025/combined"

# Create output directory
dx mkdir -p ${out_dir}

# Define parameters
readonly ancestries=("eur")
#readonly ancestries=("eur" "sas" "afr" "eas")
readonly heritabilities=("0.01" "0.1" "0.2" "0.3" "0.5")
readonly n_reps=4
readonly annotations=("pLoF_damaging_missense")
readonly priority="low"
readonly instance_type="mem1_ssd1_v2_x2"

# Loop through ancestries
for anc in "${ancestries[@]}"; do

  echo "Processing ancestry: ${anc}"

  # Loop through heritabilities
  for h2 in "${heritabilities[@]}"; do
    # Loop through replicates
    for rep in $(seq 1 ${n_reps}); do
      pheno="p_${h2}_continuous_${rep}"

      # Loop through annotations
      for anno in "${annotations[@]}"; do

        # Define input files
        add_file="${step2_dir}/regenie_step2.${anc}.${pheno}.${anno}.additive.regenie.gz"
        rec_file="${step2_dir}/regenie_step2.${anc}.${pheno}.${anno}.recessive.regenie.gz"
        dom_file="${step2_dir}/regenie_step2.${anc}.${pheno}.${anno}.dominance.regenie.gz"

        # Output file
        out_prefix="regenie_step2.${anc}.${pheno}.${anno}.combined"
        out_file="${out_dir}/${out_prefix}.regenie.tsv.gz"

        echo "Submitting combine job for: ${out_prefix}"

        # Check if output already exists
        echo "  Checking if output exists: ${out_file}"
        out_exists=$(dx ls "${out_file}" 2>/dev/null | wc -l)
        echo "    Output exists: ${out_exists}"

        if [[ ${out_exists} -eq 0 ]]; then

          # Check if all input files exist
          echo "  Checking input files..."
          add_exists=$(dx ls "${add_file}" 2>/dev/null | wc -l)
          rec_exists=$(dx ls "${rec_file}" 2>/dev/null | wc -l)
          dom_exists=$(dx ls "${dom_file}" 2>/dev/null | wc -l)

          echo "    Additive: ${add_exists} (${add_file})"
          echo "    Recessive: ${rec_exists} (${rec_file})"
          echo "    Dominance: ${dom_exists} (${dom_file})"

          if [[ ${add_exists} -gt 0 && ${rec_exists} -gt 0 && ${dom_exists} -gt 0 ]]; then

            dx run app-swiss-army-knife \
              -iimage_file="/docker/rsuite.tar.gz" \
              -iin="/${rscript_remote}" \
              -iin="${add_file}" \
              -iin="${rec_file}" \
              -iin="${dom_file}" \
              -icmd="
                Rscript /mnt/project/${rscript_remote} \
                  --additive regenie_step2.${anc}.${pheno}.${anno}.additive.regenie.gz \
                  --recessive regenie_step2.${anc}.${pheno}.${anno}.recessive.regenie.gz \
                  --dominance regenie_step2.${anc}.${pheno}.${anno}.dominance.regenie.gz \
                  --output ${out_prefix}.regenie.tsv &&
                gzip ${out_prefix}.regenie.tsv &&
                echo '!!!!$(date)'
              " \
              --instance-type ${instance_type} \
              --priority ${priority} \
              --destination "${out_dir}" \
              -y \
              --name "combine_${anc}_${pheno}_${anno}"

          else
            echo "  WARNING: Not all input files exist for ${out_prefix}"
            echo "    Additive: ${add_exists}, Recessive: ${rec_exists}, Dominance: ${dom_exists}"
            echo "  Skipping."
          fi

        else
          echo "  Output already exists: ${out_file}. Skipping."
        fi

      done
    done
  done
done

echo "Done submitting combine jobs!"
