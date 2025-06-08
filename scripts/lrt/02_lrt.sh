# author: frederik lassen

# note: this script uses a github repo to get the genests to be processed https://github.com/frhl/genesets_public

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

readonly remote_dir="wes_ko_ukbb/scripts"
readonly rscript_local="02_lrt.R"
readonly rscript_remote="${remote_dir}/02_lrt.R"
dx_update_remote ${rscript_remote} ${rscript_local}

readonly pheno_dir="/mnt/project/wes_ko_ukbb/data/phenotypes/quantitative/with_whr"
readonly phenos_path="${pheno_dir}/ukb_cts_phenotypes.inst0.with_covariates.txt.gz"

readonly pheno_path_dir="/wes_ko_ukbb/data/phenotypes/brava/"
readonly phenos="${pheno_path_dir}/brava_cts_phenos.txt"

readonly grm_path="/mnt/project/wes_ko_ukbb/data/saige/step0/vr_20k/UKB.array.eur_relatednessCutoff_0.05_5000_randomMarkersUsed.sparseGRM.mtx"

readonly covariates="age,age2,age_sex,age2_sex,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
readonly categorical_covariates="sex"

#readonly samples="/mnt/project/wes_ko_ukbb/data/samples/UKB.wes.qced.eur.samples"
readonly group="spliceai=0.50_cadd=28.1_revel=0.773"
readonly anno="pLoF_damaging_missense"

for pop in "eur"; do
  out_dir="/wes_ko_ukbb/lrt_grm/${pop}"
  dx mkdir -p ${out_dir}
  for pp in  "0.90"; do
    for af in "05"; do
      for pheno in $(dx cat $phenos | grep choles | head -n1); do
        for mm in "add"; do
          in_dir="/mnt/project/wes_ko_ukbb/data/phased/encode_alt_qced_canonical/${pop}/${group}"
          in_file="${in_dir}/UKB.wes.merged.phased.full_qc.eur.af${af}.popmax.hwe.pp${pp}.${group}.${anno}.txt.gz"
          out_file="UKB.wes.merged.phased.full_qc.eur.af${af}.popmax.hwe.pp${pp}.${group}.${anno}.${pheno}.${mm}.with_grm.txt.gz"
          if [[ $(dx_file_exists "${out_dir}/${out_file}") -eq "0" ]]; then 
              dx run app-swiss-army-knife \
                -iimage_file="/docker/rsuite.tar.gz"\
                -icmd="
                   Rscript /mnt/project/${rscript_remote} \
                      --phenotypes_path ${phenos_path} \
                      --sparse_grm_path ${grm_path} \
                      --covariates ${covariates} \
                      --categorical_covariates ${categorical_covariates} \
                      --dosages_path ${in_file} \
                      --out_path ${out_file} \
                      --pheno_col ${pheno} \
                      --mixed_model ${mm}
                     echo '!!!!$(date)'
                  "\
                --instance-type mem3_ssd1_v2_x4 \
                --folder=".${out_dir}" \
                --priority normal \
                --name wald_mixed_${anno}_${pheno} -y
          else
             >&2 echo "Error! '${in_vcf_auto}' does not exists."
          fi
        done
      done
    done
  done
done





