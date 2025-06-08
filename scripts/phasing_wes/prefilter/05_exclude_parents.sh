# author: Frederik Lassen

readonly threads=48
readonly in_dir="/mnt/project/Phasing/PhasingWES/step0_merge_qc"
readonly out_dir="/Phasing/PhasingWES/step0_merge_qc"
readonly file_parents="${in_dir}/ukb_eur_random_trios_seed104_parents.txt" # testing with qced data

for CHR in 17; do
  input_vcf="${in_dir}/UKB.chr${CHR}.exome_array.full_qc.bcf"
  out_prefix_no_parents="UKB.chr${CHR}.exome_array.full_qc.no_parents"
  out_expt="${out_dir}/${out_prefix_no_parents}.bcf"
  if [[ "$(dx ls ${out_expt} | wc -l)" -eq "0" ]]; then
    dx run app-swiss-army-knife -icmd="
          bcftools view ${input_vcf} --threads ${threads} -S ^${file_parents} --force-samples -Oz -o ${out_prefix_no_parents}.no_tags.vcf.gz &&
          bcftools index --threads ${threads} ${out_prefix_no_parents}.no_tags.vcf.gz &&
          bcftools +fill-tags --threads ${threads} ${out_prefix_no_parents}.no_tags.vcf.gz -Ob -o ${out_prefix_no_parents}.bcf -- -t AN,AC &&
          bcftools index --threads ${threads} ${out_prefix_no_parents}.bcf &&
          rm ${out_prefix_no_parents}.no_tags.vcf.gz* &&
          echo 'ok'
    " \
    --instance-type mem2_ssd1_v2_x16 \
    --depends-on job-GjQ43b8Jg8JXqvPv2v05F8Jj \
    --folder=".${out_dir}" \
    --name WESarray_exclude_parents_chr${CHR} \
    --priority normal -y
  else
     >&2 echo "${out_expt} already exists! Skipping."
  fi
done

#plink2 --threads ${threads} --bfile ${in_prefix} --remove ${file_parents} --keep ${keep_samples} --output-chr chrMT --recode vcf bgz --out ${out_prefix_no_parents}.no_tags &&
  # rm ${out_prefix_no_parents}.no_tags.log
