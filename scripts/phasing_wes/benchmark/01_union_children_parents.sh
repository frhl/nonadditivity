# author: frederik lassen

set -eu

out_dir="/Phasing/PhasingWES/step3_benchmark_qc"
in_dir_unphased="/mnt/project/Phasing/PhasingWES/step0_merge_qc"
in_dir_phased="/mnt/project/Phasing/PhasingWES/step2_phase_rare_qc/ligated"
#in_dir_phased="/mnt/project/Phasing/PhasingWES/step4_polish_phase/rephased"
in_dir_pedigree="/mnt/project/Phasing/PhasingWES/step0_merge_qc"

dx mkdir -p ${out_dir}
for CHR in 17; do
  
  # use polished or ligated VCFs  
  phased="${in_dir_phased}/UKB_chr${CHR}.exome_array.shapeit5.ligated.bcf"
  
  # should not change
  unphased="${in_dir_unphased}/UKB.chr${CHR}.exome_array.full_qc.bcf"
  pedigree_list="${in_dir_pedigree}/ukb_eur_random_trios_seed104.ped"
  trios_list="${in_dir_pedigree}/ukb_eur_random_trios_seed104_all.txt"
  parents_list="${in_dir_pedigree}/ukb_eur_random_trios_seed104_parents.txt"
  
  # setup out_prefix
  out_prefix="UKB.chr${CHR}.exome_array.full_qc.shapeit5.ligated.trio"
  out_file="${out_prefix}.vcf.gz"
  out_ser="${out_prefix}.switches"
  dx run app-swiss-army-knife \
    -icmd="
      bcftools view --samples-file <(sort ${parents_list}) ${unphased} --force-samples -Ob -o unphased.parents.bcf &&
      bcftools index unphased.parents.bcf &&
      bcftools query -l unphased.parents.bcf | wc -l &&
      bcftools view --samples-file <(sort ${trios_list}) ${phased} --force-samples -Ob -o phased.children.bcf &&
      bcftools index phased.children.bcf &&
      bcftools query -l phased.children.bcf | wc -l &&
      bcftools merge phased.children.bcf unphased.parents.bcf -Oz -o ${out_file} &&
      bcftools index ${out_file} &&
      bcftools +trio-switch-rate ${out_file} -- -p ${pedigree_list} > ${out_ser} &&
      rm unphased.parents.bcf* &&
      rm phased.children.bcf* &&
      echo 'ciao'
      "\
    --instance-type mem1_ssd1_v2_x8 \
    --folder=".${out_dir}" \
    --priority normal \
    --name trio_ser_chr${CHR} -y
done


