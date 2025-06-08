# author: Frederik Lassen

set -eu
readonly out_dir="/Phasing/PhasingWES/step0_merge_qc"
readonly instance_type="mem2_ssd1_v2_x16"
readonly missing="0.10"
readonly threads="96"
readonly priority="high"
for CHR in 17; do
    exome_dir="/Phasing/PhasingWES/step0_merge_qc/support"
    array_dir="/Phasing/PhasingWES/step0_merge_qc/support"
    mnt_exome_dir="/mnt/project${exome_dir}"
    mnt_array_dir="/mnt/project${array_dir}"
    keep_samples="/mnt/project/wes_ko_ukbb/data/samples/UKB.wes.qced.combined.overlapping_exome_array" 
    array_prefix="${mnt_array_dir}/UKB.chr${CHR}.qc.array"
    exome_prefix="${mnt_exome_dir}/UKB.chr${CHR}.full_qc.exome"
    out_prefix=UKB.chr${CHR}.exome_array.full_qc
    out_expt="${out_dir}/${out_prefix}.bed"
    if [[ "$(dx ls ${out_expt} | wc -l)" -eq "0" ]]; then 
      dx run app-swiss-army-knife -icmd="
        plink --threads ${threads} --bfile ${exome_prefix} --bmerge ${array_prefix} --geno ${missing} --keep-allele-order --make-bed --out ${out_prefix} &&
        plink2 --threads ${threads} --bfile ${out_prefix} --keep ${keep_samples} --output-chr chrMT --recode vcf bgz --out ${out_prefix} &&
        bcftools index --threads ${threads} ${out_prefix}.vcf.gz &&
        bcftools view --threads ${threads} ${out_prefix}.vcf.gz -Ob -o ${out_prefix}.bcf &&
        bcftools index --threads ${threads} ${out_prefix}.bcf &&
        rm ${out_prefix}.vcf.gz* &&
        rm ${out_prefix}.nosex &&
        rm ${out_prefix}.log &&
        echo 'ciao'
      " \
      --instance-type ${instance_type} \
      --folder=".${out_dir}" \
      --name WESarray_merge_chr${CHR} \
      --tag "chr${CHR}" \
      --priority ${priority} -y
   else
     >&2 echo "${out_expt} already exists! Skipping."
   fi
  
done



