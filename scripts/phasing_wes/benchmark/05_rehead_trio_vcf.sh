# author: frederik lassen

set -eu

in_dir="/mnt/project/Phasing/PhasingWES/step3_benchmark/polished"
out_dir="/Phasing/PhasingWES/step3_benchmark/trio_merge/polished"

reference="/mnt/project/Phasing/Reference/Homo_sapiens_assembly38.fasta"
reference_index="${reference}.fai"

dx mkdir -p ${out_dir}
for CHR in 1; do
  in_file="${in_dir}/UKB.chr${CHR}.exome_array.shapeit5.ligated.polished.trio.vcf.gz"
  out_file="UKB.chr${CHR}.exome_array.shapeit5.ligated.polished.trio.reheader.vcf.gz"
  dx run app-swiss-army-knife \
    -icmd="
       echo -e '##fileformat=VCFv4.2' > tmp.vcf &&
       bcftools view ${in_file} | grep '#CHROM' >> tmp.vcf &&
       bcftools reheader -f ${reference_index} tmp.vcf > fai.vcf &&
       bcftools concat --no-version -Oz -o ${out_file} fai.vcf ${in_file} && 
       rm tmp.vcf &&
       rm fai.vcf
      "\
    --instance-type mem1_ssd1_v2_x8 \
    --folder=".${out_dir}" \
    --priority normal \
    --name rehead_chr${CHR} -y
done

