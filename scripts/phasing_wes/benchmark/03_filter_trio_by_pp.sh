# author: frederik lassen


source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"


set -eu

readonly in_dir="/mnt/project/Phasing/PhasingWES/step3_benchmark_qc"
readonly out_dir="/Phasing/PhasingWES/step3_benchmark_qc/pp_gt_90"
readonly in_dir_pedigree="/mnt/project/Phasing/PhasingWES/step0_merge_qc"
readonly pedigree_list="${in_dir_pedigree}/ukb_eur_random_trios_seed104.ped"

dx mkdir -p ${out_dir}
for CHR in {1..22} X; do
  in_file="${in_dir}/UKB.chr${CHR}.exome_array.full_qc.shapeit5.ligated.trio.vcf.gz"
  out_file="UKB.chr${CHR}.exome_array.shapeit5.full_qc.ligated.pp_lt_90.trio.vcf.gz"
  out_ser="UKB.chr${CHR}.exome_array.full_qc.shapeit5.ligated.pp_lt_90.trio.switches"
  echo "VCF = ${in_file} ($(dx_size_in_bytes $(wo_mnt_project ${in_file}))B)"
  
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
      --name trio_ser_chr${CHR} -y
  fi
done


