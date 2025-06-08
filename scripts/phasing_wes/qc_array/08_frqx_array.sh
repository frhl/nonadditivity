
# author: frederik lassen

set -eu

readonly threads=16
readonly out_dir="/Phasing/PhasingSNParray/step4_liftover"
readonly qced_eur="/mnt/project/wes_ko_ukbb/data/samples/UKB.wes.qced.eur.samples"
readonly in_dir="/mnt/project/Phasing/PhasingSNParray/step4_liftover"

dx mkdir -p ${out_dir}



for CHR in {1..20} 22 X; do
  array="${in_dir}/full_c$CHR\_b0_v2.b38.annotated.vcf.gz"
  out_prefix="full_c${CHR}_b0_v2.b38.annotated"
  dx run app-swiss-army-knife -icmd="
    plink --vcf ${array} --freqx --keep  <( cat ${qced_eur} | awk '{print \$0 \"\t\" \$0}') --out ${out_prefix}.eur &&
    plink --vcf ${array} --freqx --out ${out_prefix} &&
    rm *.log &&
    rm *.nosex
  "\
  --tag chr${CHR}\
  --instance-type mem1_ssd1_v2_x8\
  --folder=".${out_dir}"\
  --name frqx_array_chr${CHR}\
  --priority normal -y

done

