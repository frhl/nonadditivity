# author: frederik lassen

set -eu

in_dir="/mnt/project/Phasing/PhasingWES/step3_benchmark/trio_merge/polished"
out_dir="/Phasing/PhasingWES/step3_benchmark/trio_merge/polished"

dx mkdir -p ${out_dir}
in_files="${in_dir}/UKB.chr*.exome_array.shapeit5.ligated.polished.trio.reheader.vcf.gz"
out_file="UKB.autosomes.exome_array.shapeit5.ligated.polished.trio.vcf.gz"
dx run app-swiss-army-knife \
  -icmd="
    bcftools concat ${in_files} -Oz -o tmp.vcf.gz &&
    bcftools sort tmp.vcf.gz -Oz -o ${out_file}
  "\
  --instance-type mem1_ssd1_v2_x8 \
  --folder=".${out_dir}" \
  --priority normal \
  --name concat_trios_by_chrom -y


