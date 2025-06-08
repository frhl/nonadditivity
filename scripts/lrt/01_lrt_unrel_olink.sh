# author: frederik lassen

# note: this script uses a github repo to get the genests to be processed https://github.com/frhl/genesets_public

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

readonly remote_dir="wes_ko_ukbb/scripts"
readonly rscript_local="01_lrt_unrel.R"
readonly rscript_remote="${remote_dir}/01_lrt_unrel.R"
dx_update_remote ${rscript_remote} ${rscript_local}

readonly pheno_dir="/mnt/project/wes_ko_ukbb/data/phenotypes"
readonly phenos="${pheno_dir}/olink_3k.npx_adj_lod.with_covariates.phenos"
readonly phenos_path="${pheno_dir}/olink_3k.npx_adj_lod.with_covariates.txt.gz"

readonly covariates="age,age2,age_sex,age2_sex,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
readonly categorical_covariates="sex"

readonly samples="/mnt/project/wes_ko_ukbb/data/samples/UKB.wes.qced.eur.unrelated.samples"

readonly group="spliceai=0.50_cadd=28.1_revel=0.773"

readonly anno="pLoF_damaging_missense"

for pop in "eur"; do
  out_dir="/wes_ko_ukbb/data/lrt_testing/lrt_unrel/olink/2df_chi2_trick_final"
  dx mkdir -p ${out_dir}
  for pp in  "0.90"; do
    for af in "05"; do
      echo "${phenos}"
      for pheno in $(dx cat $(wo_mnt_project $phenos)); do
        echo "Phenotype='${pheno}'"
        in_dir="/mnt/project/wes_ko_ukbb/data/phased/encode_alt_qced_canonical/${pop}/${group}"
        in_file="${in_dir}/UKB.wes.merged.phased.qc_final.eur.af${af}.popmax.pp${pp}.${group}.${anno}.txt.gz"
        out_file="UKB.wes.merged.phased.qc_final.eur.af${af}.popmax.pp${pp}.${group}.${anno}.${pheno}.unrel.chi2_trick.txt.gz"
        if [[ $(dx_file_exists "${out_dir}/${out_file}") -eq "0" ]]; then 
            dx run app-swiss-army-knife \
              -iimage_file="/docker/rsuite.tar.gz"\
              -icmd="
                 Rscript /mnt/project/${rscript_remote} \
                    --phenotypes_path ${phenos_path} \
                    --covariates ${covariates} \
                    --categorical_covariates ${categorical_covariates} \
                    --dosages_path ${in_file} \
                    --out_path ${out_file} \
                    --samples_path ${samples} \
                    --pheno_col ${pheno}
                   echo '!!!!$(date)'
                "\
              --instance-type mem2_ssd1_v2_x2 \
              --folder=".${out_dir}" \
              --priority low \
              --name lrt_lm_cts_${pheno} -y
        else
           >&2 echo "Skipping existing '${out_file}'.."
        fi
      done
    done
  done
done





