# author: frederik lassen
# note: this script simulates phenotypes for different ancestries
source "/Users/flassen/Projects/00_dx_utils/dx_utils/dx_utils.sh"

set -euo pipefail


readonly remote_dir="wes_ko_ukbb/scripts"
readonly rscript_local="00_make_sim_phenos.R"
readonly rscript_remote="${remote_dir}/00_make_sim_phenos.R"

# Upload R script to remote
dx_update_remote ${rscript_remote} ${rscript_local}

# Define paths
readonly RAP_location="/mnt/project/wes_ko_ukbb/data/saige/step0/vr_20k/"
readonly RAP_covars_dir="/mnt/project/wes_ko_ukbb/data/phenotypes/covariates/"
readonly RAP_covars="${RAP_covars_dir}/core_covariates_with_PCs.tsv.gz"
readonly out_dir="/wes_ko_ukbb/data/phenotypes/simulated/revisions/051025"

# Create output directory
dx mkdir -p ${out_dir}

# Define parameters
readonly ancestries=("eur" "sas" "afr" "eas")
readonly heritabilities="0.01,0.1,0.2,0.3,0.5"
readonly prevalence="0.5"
readonly seed="1236"
readonly n_reps="4"
readonly trait_type="continuous"

for ancestry in "${ancestries[@]}"; do
  sampleIDs="UKB.array.${ancestry}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx.sampleIDs.txt"
  sparseGRM="UKB.array.${ancestry}_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"
  out_prefix="051025_${ancestry}_null_phenos.txt.gz"
  
  # Check if output already exists
  if [[ $(dx_file_exists "${out_dir}/${out_prefix}") -eq "0" ]]; then
    echo "Submitting job for ancestry: ${ancestry}"
    echo "Using samples from: ${RAP_location}${sampleIDs}"
    echo "Using GRM from: ${RAP_location}${sparseGRM}"
    
    dx run app-swiss-army-knife \
      -iimage_file="/docker/rsuite.tar.gz" \
      -icmd="
        Rscript /mnt/project/${rscript_remote} \
          --samples ${RAP_location}${sampleIDs} \
          --grm ${RAP_location}${sparseGRM} \
          --covars ${RAP_covars} \
          --covar_id_col participant_eid \
          --prevalence ${prevalence} \
          --heritabilities ${heritabilities} \
          --seed ${seed} \
          --out ${out_prefix} \
          --n_reps ${n_reps} \
          --trait_type ${trait_type} &&
        echo '!!!!$(date)'
      " \
      --instance-type mem1_ssd1_v2_x16 \
      --folder=".${out_dir}" \
      --priority low \
      --name sim_phenos_${ancestry} -y
  else
    >&2 echo "Output already exists: ${out_dir}/${out_prefix}"
  fi
done
