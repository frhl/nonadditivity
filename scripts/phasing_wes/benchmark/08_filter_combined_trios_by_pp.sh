
source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"
set -eu

in_dir="/mnt/project/Phasing/PhasingWES/step3_benchmark/trio_merge/polished"
out_dir="/Phasing/PhasingWES/step3_benchmark/trio_merge/polished"
in_dir_pedigree="/mnt/project/Phasing/PhasingWES/step0_merge"
pedigree_list="${in_dir_pedigree}/ukb_wes_200k_trios.ped"

dx mkdir -p ${out_dir}
in_file="${in_dir}/UKB.autosomes.exome_array.shapeit5.ligated.polished.trio.vcf.gz"
out_file="UKB.autosomes.exome_array.shapeit5.ligated.polished.pp_lt_90.trio.vcf.gz"
out_ser="UKB.autosomes.exome_array.shapeit5.ligated.polished.pp_lt_90.switches"
if [[ $(dx_file_exists "${out_dir}/${out_file}") -eq "0" ]]; then
  dx run app-swiss-army-knife \
    -icmd="
      bcftools +setGT ${in_file}  -- -t q -n . -e '(FMT/PP > 0.90)|(FMT/PP=\".\")' | bcftools view -Oz -o ${out_file}
      bcftools index ${out_file} &&
      bcftools +trio-switch-rate ${out_file} -- -p ${pedigree_list} > ${out_ser}
      "\
    --instance-type mem1_ssd1_v2_x8 \
    --folder=".${out_dir}" \
    --priority normal \
    --name filter_trio_autosome_by_pp -y
fi

