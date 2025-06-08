# author: frederik lassen

# note: this script uses a github repo to get the genesets to be processed https://github.com/frhl/genesets_public

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

readonly remote_dir="wes_ko_ukbb/scripts"
readonly rscript_local="01_lrt_null_traits.R"
readonly rscript_remote="${remote_dir}/01_lrt_null_traits.R"
dx_update_remote ${rscript_remote} ${rscript_local}

readonly covariates="age,age2,age_sex,age2_sex,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
readonly categorical_covariates="sex"

readonly samples="/mnt/project/wes_ko_ukbb/data/samples/UKB.wes.qced.eur.unrelated.samples"

readonly group="spliceai=0.50_cadd=28.1_revel=0.773"
readonly anno="pLoF_damaging_missense"

# Number of simulations to run
readonly num_simulations=50
# Starting seed
readonly base_seed=12345

for pop in "eur"; do
  out_dir="/wes_ko_ukbb/data/lrt_testing/lrt_null_simulations"
  dx mkdir -p ${out_dir}
  for pp in "0.90"; do
    for af in "05"; do
      for ((sim=1; sim<=${num_simulations}; sim++)); do
        # Calculate unique seed for this simulation
        sim_seed=$((base_seed + sim))
        
        # Create a unique name for this simulation
        sim_name="sim_${sim}"
        
        in_dir="/mnt/project/wes_ko_ukbb/data/phased/encode_alt_qced_canonical/${pop}/${group}"
        in_file="${in_dir}/UKB.wes.merged.phased.qc_final.eur.af${af}.popmax.pp${pp}.${group}.${anno}.txt.gz"
        out_file="UKB.wes.merged.phased.qc_final.eur.af${af}.popmax.pp${pp}.${group}.${anno}.${sim_name}.null_sim.txt.gz"
        
        if [[ $(dx_file_exists "${out_dir}/${out_file}") -eq "0" ]]; then 
            dx run app-swiss-army-knife \
              -iimage_file="/docker/rsuite.tar.gz"\
              -icmd="
                 Rscript /mnt/project/${rscript_remote} \
                    --samples_path ${samples} \
                    --covariates ${covariates} \
                    --categorical_covariates ${categorical_covariates} \
                    --dosages_path ${in_file} \
                    --out_path ${out_file} \
                    --sim_seed ${sim_seed}
                   echo '!!!!$(date)'
                "\
              --instance-type mem2_ssd1_v2_x2 \
              --folder=".${out_dir}" \
              --priority low \
              --name lrt_null_${sim_name} -y
        else
           >&2 echo "Skipping existing '${out_file}'.."
        fi
      done
    done
  done
done
