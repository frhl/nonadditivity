# author: Frederik Lassen

in_dir="/mnt/project/Phasing/PhasingWES/step0_merge_qc"
out_dir="/Phasing/PhasingWES/step0_merge_qc"
threads=8

for CHR in 17; do
    in_file="${in_dir}/UKB.chr${CHR}.exome_array.full_qc.no_parents.bcf"
    out_file="UKB.chr${CHR}.exome_array.full_qc.no_parents.variants_only.vcf.gz"
    dx rm  ${out_dir}/${out_file}
    dx run app-swiss-army-knife -icmd="
       bcftools view --drop-genotypes --threads ${threads} ${in_file} -Oz -o ${out_file}
    " \
  --instance-type mem2_ssd1_v2_x4 \
  --folder=".${out_dir}" \
  --name export_variants_chr${CHR} \
  --priority normal -y
done

