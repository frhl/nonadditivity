
# author: frederik lassen

set -eu

readonly threads=16
readonly in_dir="/mnt/project/Phasing/PhasingSNParray/step4_liftover"
readonly out_dir="/Phasing/PhasingSNParray/step4_liftover"
readonly qced_snps="/mnt/project/Phasing/PhasingSNParray/step1_dataqc/SNPlist.filtered.QC.txt"
dx mkdir -p ${out_dir}

for CHR in {1..22} X; do
  array="${in_dir}/full_c$CHR\_b0_v2.b38.sorted.vcf.gz"
  if [[ "${CHR}" == "X" ]]; then
    array="${in_dir}/full_c$CHR\_b0_v2.b38.fixploidity.vcf.gz"
  fi
  out_prefix="full_c${CHR}_b0_v2.b38.annotated"
  dx run app-swiss-army-knife -icmd="
    bcftools view -i 'ID=@${qced_snps}' ${array} | bcftools annotate --threads ${threads} --set-id '%CHROM\:%POS\:%REF\:%ALT' -Oz -o ${out_prefix}.vcf.gz &&
    bcftools index ${out_prefix}.vcf.gz
  "\
  --tag chr${CHR}\
  --instance-type mem1_ssd1_v2_x8\
  --folder=".${out_dir}"\
  --name filter_qc_chr${CHR}\
  --priority normal -y
done

