# author: frederik lassen

set -o errexit
set -o nounset

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

threads=24

in_dir="/mnt/project/wes_ko_ukbb/data/phased/qced_samples_and_sites"
out_dir="/wes_ko_ukbb/data/phased/qced_samples_and_sites"

dx mkdir -p ${out_dir}
for CHR in {1..22} X; do
  for pop in "afr" "eas" "sas"; do
    in_bcf="${in_dir}/UKB.wes.chr${CHR}.phased.full_qc.${pop}.bcf"
    out_prefix="UKB.wes.chr${CHR}.phased.full_qc.${pop}"
    if [[ $(dx_file_exists "${out_dir}/${out_prefix}.frqx.gz") -eq 0 ]]; then
      dx run app-swiss-army-knife \
        -icmd="
          plink2 --bcf ${in_bcf} --threads ${threads} --output-chr chrMT --make-bed --out tmp &&
          plink --bfile tmp --threads ${threads} --freqx --out ${out_prefix} && gzip ${out_prefix}.frqx &&
          rm ${out_prefix}.log &&
          rm -f ${out_prefix}.nosex &&
          rm -f tmp*
          "\
        --instance-type mem2_ssd1_v2_x8 \
        --folder=".${out_dir}" \
        --priority high \
        --depends-on job-Gk7QFX0Jg8JQVZX2GbKQJf8k \
        --name c${CHR}_${pop}_frqx -y
      else
        >&2 echo "${out_prefix}.frqx.gz already exists. Skipping.."
      fi
   done
done


