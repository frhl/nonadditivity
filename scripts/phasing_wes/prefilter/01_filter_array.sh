# author: frederik lassen

set -eu

readonly threads=16
readonly out_dir="/Phasing/PhasingWES/step0_merge_qc/support"
dx mkdir -p ${out_dir}

for CHR in {1..20} 22 X; do
  array="/mnt/project/Phasing/PhasingSNParray/step4_liftover/full_c$CHR\_b0_v2.b38.annotated.filtered.bcf"
  array_samples="/mnt/project/Phasing/PhasingWES/step0_merge_qc/array.samples.txt"
  keep_samples="/mnt/project/Phasing/PhasingWES/step0_merge_qc/overlapping_samples.qced.txt"
  out_prefix="UKB.chr${CHR}.qc.array"
  out_expt="${out_dir}/${out_prefix}.bed"
  if [[ "$(dx ls ${out_expt} | wc -l)" -eq "0" ]]; then
    echo "* Submitting chr${CHR}\n"
    dx run app-swiss-army-knife -icmd="
      bcftools reheader --threads ${threads} --samples ${array_samples} ${array} | bcftools view --threads ${threads} --samples-file ${keep_samples} -Oz -o tmp0.bcf &&
      bcftools index --threads ${threads} tmp0.bcf &&
      bcftools query -f '%ID\n' tmp0.bcf > chr${CHR}.array_snps_kept.txt &&
      plink2 --bcf tmp0.bcf --make-bed --out ${out_prefix} &&
      rm ${out_prefix}.log &&
      rm tmp0.bcf*
    "\
    --tag chr${CHR}\
    --instance-type mem1_ssd1_v2_x16\
    --folder=".${out_dir}"\
    --name filter_SNParray_chr${CHR}\
    --priority normal -y
  else
    >&2 echo "${out_expt} already exists! Skipping."
  fi

done



