# author: frederik lassen

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

set -eu
readonly threads=24
readonly in_dir="/mnt/project/wes_ko_ukbb/data/conditional/imputed/recoded/v2"
readonly out_dir="/wes_ko_ukbb/data/conditional/imputed/recoded/v2"
readonly maf="01"

dx mkdir -p ${out_dir}
for CHR in {1..22}; do
  for pop in "eur"; do
    in_prefix="UKB.gel_imputed.sorted.maf${maf}.chr${CHR}.${pop}.final.v2.recessive"
    out_prefix="UKB.gel_imputed.sorted.maf${maf}.chr${CHR}.${pop}.final.v2.recessive"
    in="${in_dir}/${in_prefix}.vcf.gz"
    dx ls $(wo_mnt_project ${in}) 
    # run recessive recoding
    if [[ $(dx_file_exists "${out_dir}/${out_prefix}.bed") -eq "0" ]]; then
      dx run app-swiss-army-knife \
        -icmd="
             plink2 --vcf ${in} dosage=DS --make-bed --out ${out_prefix} &&
             echo '$(date)'
              "\
        --instance-type mem2_ssd1_v2_x2 \
        --folder=".${out_dir}" \
        --priority low \
        --name c${CHR}_recessive_to_plink -y
    else
      >&2 echo "${out_prefix}.recessive.vcf.gz already exists!"
    fi

  done
done

