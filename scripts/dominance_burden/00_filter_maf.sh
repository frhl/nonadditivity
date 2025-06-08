# author: frederik lassen

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"


set -eu
readonly threads=12
readonly out_dir="/wes_ko_ukbb/data/phased/export_alt_qced"
readonly in_dir="/mnt/project/wes_ko_ukbb/data/phased/final_qc"
readonly pop_max_dir="/mnt/project/wes_ko_ukbb/data/variants/popmax_exclude"

set -x
dx mkdir -p ${out_dir}
#for pop in "eur" "afr" "eas" "sas"; do
#for pop in "afr" "eas" "sas"; do
for pop in "eur"; do
  for af in "05"; do
    for CHR in {1..22} X; do
      variants_to_exclude="${pop_max_dir}/gnomad.exomes.r2.1.1.grch38.popmax${af}.tsv"
      input_vcf="${in_dir}/UKB.wes.chr${CHR}.phased.qc_final.${pop}.bcf"
      out_prefix="UKB.wes.chr${CHR}.phased.qc_final.${pop}.af${af}.popmax.variants"
      dx ls $(wo_mnt_project ${input_vcf})
      if [[ $(dx_file_exists "${out_dir}/${out_prefix}.txt.gz") -eq 0 ]]; then
        dx run app-swiss-army-knife \
          -icmd="
               bcftools view ${input_vcf} --threads ${threads} --max-af 0.${af} -Ou | bcftools view -e 'ID=@${variants_to_exclude}' -Ob -o ${out_prefix}.bcf &&
               bcftools index ${out_prefix}.bcf
        "\
          --instance-type mem1_ssd1_v2_x2 \
          --folder=".${out_dir}" \
          --priority low \
          --name export_alt_chr${CHR} -y
        else
          >&2 echo "${out_prefix}.tar.gz already exists. Skipping.."
        fi
    done
  done
done



