# author: frederik lassen

set -o errexit
set -o nounset

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

threads=24

readonly in_dir="/mnt/project/wes_ko_ukbb/data/phased/qced_samples_and_sites"
readonly out_dir="/wes_ko_ukbb/data/phased/qced_samples_and_sites"
readonly samples_dir="/mnt/project/wes_ko_ukbb/data/samples"

dx mkdir -p ${out_dir}
for CHR in {1..22} X; do
  for pop in "afr" "eas" "sas"; do
#  for pop in "eur"; do
    in_bcf="${in_dir}/UKB.wes.chr${CHR}.phased.full_qc.${pop}.bcf"
    out_prefix="UKB.wes.chr${CHR}.phased.full_qc.${pop}.unrelated"
    samples="${samples_dir}/UKB.wes.qced.${pop}.unrelated.samples"
    if [[ $(dx_file_exists "${out_dir}/${out_prefix}.frqx.gz") -eq 0 ]]; then
      dx run app-swiss-army-knife \
        -icmd="
          plink2 --bcf ${in_bcf} --keep ${samples} --threads ${threads} --output-chr chrMT --make-bed --out tmp &&
          plink --bfile tmp --threads ${threads} --freqx --out ${out_prefix} && gzip ${out_prefix}.frqx &&
          rm ${out_prefix}.log &&
          rm -f ${out_prefix}.nosex &&
          rm -f tmp*
          "\
        --instance-type mem1_ssd1_v2_x2 \
        --folder=".${out_dir}" \
        --priority high \
        --name c${CHR}_${pop}_unrel_frqx -y
      else
        >&2 echo "${out_prefix}.frqx.gz already exists. Skipping.."
      fi
   done
done


