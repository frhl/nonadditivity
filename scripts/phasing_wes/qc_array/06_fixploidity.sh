
# author: frederik lassen

set -eu

readonly threads=16
readonly out_dir="/Phasing/PhasingSNParray/step4_liftover"
dx mkdir -p ${out_dir}

for CHR in X; do
  array="/mnt/project/Phasing/PhasingSNParray/step4_liftover/full_c$CHR\_b0_v2.b38.sorted.vcf.gz"
  out_prefix="full_c${CHR}_b0_v2.b38.fixploidity"
  dx run app-swiss-army-knife -icmd="
    bcftools index ${out_prefix}.vcf.gz
  "\
  --tag chr${CHR}\
  --instance-type mem1_ssd1_v2_x8\
  --folder=".${out_dir}"\
  --name fix_ploidity_chr${CHR}\
  --priority normal -y
done
#bcftools +fixploidy ${array} -Oz -o ${out_prefix}.vcf.gz 

