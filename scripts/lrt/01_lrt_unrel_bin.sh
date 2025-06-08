# author: frederik lassen

# note: this script uses a github repo to get the genests to be processed https://github.com/frhl/genesets_public

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

readonly remote_dir="wes_ko_ukbb/scripts"
readonly rscript_local="01_lrt_unrel_bin.R"
readonly rscript_remote="${remote_dir}/01_lrt_unrel_bin.R"
dx_update_remote ${rscript_remote} ${rscript_local}

# call correct phenotypes
readonly pheno_dir="/wes_ko_ukbb/data/phenotypes/brava"
readonly phenos_path="/mnt/project${pheno_dir}/brava_icd10_icd9_primary_secondary_cancer.txt.gz"
readonly phenos_bin="${pheno_dir}/brava_icd10_icd9_primary_secondary_cancer.eur.counts"
readonly icd_sex_map="${pheno_dir}/two_digit_icd_codes_sex_specificity_v2.txt"

# get sex for corresponding icd code
readonly tmp_sex_map=$(mktemp)
dx cat ${icd_sex_map} |  grep -wE "(F|M)" > ${tmp_sex_map}
echo "** showing header of sex-specific ICD codes:"
cat ${tmp_sex_map} | head

readonly covariates="age,age2,age_sex,age2_sex,sex,PC1,PC2,PC3,PC4,PC5,PC6,PC7,PC8,PC9,PC10"
readonly categorical_covariates="sex"

readonly samples="/mnt/project/wes_ko_ukbb/data/samples/UKB.wes.qced.eur.unrelated.samples"

readonly group="spliceai=0.50_cadd=28.1_revel=0.773"

readonly anno="pLoF_damaging_missense"

for pop in "eur"; do
  out_dir="/wes_ko_ukbb/lrt_unrel/oct/${pop}"
  dx mkdir -p ${out_dir}
  for pp in  "0.90"; do
    for af in "05"; do
       for pheno in "J45"; do
       #for pheno in 'M16' 'E53' 'J30' 'J45' 'L20' 'L30' 'H90' 'H91' 'C18' 'K50'; do
      #for pheno in $(dx cat $phenos); do
        in_dir="/mnt/project/wes_ko_ukbb/data/phased/encode_alt_qced_canonical/${pop}/${group}"
        in_file="${in_dir}/UKB.wes.merged.phased.full_qc.eur.af${af}.popmax.hwe.pp${pp}.${group}.${anno}.txt.gz"
        out_file="UKB.wes.merged.phased.full_qc.eur.af${af}.popmax.hwe.pp${pp}.${group}.${anno}.${pheno}.unrel.txt"
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
                    --samples ${samples} \
                    --pheno_col ${pheno}
                   echo '!!!!$(date)'
                "\
              --instance-type mem3_ssd1_v2_x4 \
              --folder=".${out_dir}" \
              --priority normal \
              --name lrt_glm_binary_${pheno} -y
        else
           >&2 echo "Skipping existing '${out_file}'.."
        fi
      done
    done
  done
done





