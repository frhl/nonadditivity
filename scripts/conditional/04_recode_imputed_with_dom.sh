# author: frederik lassen

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

set -eu
readonly threads=24
readonly in_dir="/mnt/project/wes_ko_ukbb/data/conditional/imputed/v2"
readonly out_dir="/wes_ko_ukbb/data/conditional/imputed/recoded/v2"
readonly docker="/docker/call_chets_0.3.2.tar.gz"
readonly maf="01"

dx mkdir -p ${out_dir}
for CHR in {1..22}; do
  for pop in "eur"; do
    in_prefix="UKB.gel_imputed.sorted.maf${maf}.chr${CHR}.${pop}.final.v2"
    out_prefix="UKB.gel_imputed.sorted.maf${maf}.chr${CHR}.${pop}.final.v2"
    in="${in_dir}/${in_prefix}.bcf"
    if [[ $(dx_file_exists "${out_dir}/${out_prefix}.dominance.vcf.gz") -eq "0" ]]; then
      dx run app-swiss-army-knife \
        -iimage_file=${docker}\
        -icmd="
             recode_vcf --input ${in} --scale-dosages --mode dominance | bgzip > ${out_prefix}.dominance.vcf.gz &&
             tabix -C ${out_prefix}.dominance.vcf.gz &&
             echo '$(date)'
              "\
        --instance-type mem2_ssd1_v2_x2 \
        --folder=".${out_dir}" \
        --priority low \
        --name c${CHR}_recode_imputed_dominance -y
    else
      >&2 echo "${out_prefix}.vcf.gz already exists!"
    fi

    # run recessive recoding
    if [[ $(dx_file_exists "${out_dir}/${out_prefix}.recessive.vcf.gz") -eq "0" ]]; then
      dx run app-swiss-army-knife \
        -iimage_file=${docker}\
        -icmd="
             recode_vcf --input ${in} --mode recessive | bgzip > ${out_prefix}.recessive.vcf.gz &&
             tabix -C ${out_prefix}.recessive.vcf.gz &&
             echo '$(date)'
              "\
        --instance-type mem2_ssd1_v2_x2 \
        --folder=".${out_dir}" \
        --priority low \
        --name c${CHR}_recode_imputed_recessive -y
    else
      >&2 echo "${out_prefix}.recessive.vcf.gz already exists!"
    fi

  done
done

