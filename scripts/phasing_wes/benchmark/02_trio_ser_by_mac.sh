# author: frederik lassen

# need to make sure that the PED file overlaps with the merged file!
source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

set -eu
readonly out_dir="/Phasing/PhasingWES/step3_benchmark/ligated"
readonly in_dir_info="/mnt/project/Phasing/PhasingWES/step2_phase_rare_qc/ligated/info"
readonly in_dir_merged="/mnt/project/Phasing/PhasingWES/step3_benchmark_qc"
readonly in_dir_pedigree="/mnt/project/Phasing/PhasingWES/step0_merge_qc"
dx mkdir -p ${out_dir}
for CHR in 17; do
  info="${in_dir_info}/UKB_chr${CHR}.exome_array.shapeit5.ligated.all.tags.txt.gz"
  merged="${in_dir_merged}/UKB.chr${CHR}.exome_array.full_qc.shapeit5.ligated.trio.vcf.gz"
  #pedigree="${in_dir_pedigree}/ukb_wes_200k_trios.ped"
  pedigree="${in_dir_pedigree}/ukb_eur_random_trios_seed104.ped"
  out="UKB.chr${CHR}.exome_array.full_qc.shapeit5.ligated.trio.switches"
  
  echo "## INPUT FILE STATUS" 
  echo "INFO = ${info} ($(dx_size_in_bytes $(wo_mnt_project ${info}))B)"
  echo "MERGED = ${merged} ($(dx_size_in_bytes $(wo_mnt_project ${merged}))B)"
  echo "PED = ${pedigree} ($(dx_size_in_bytes $(wo_mnt_project ${pedigree}))B)"

  dx run app-swiss-army-knife \
    -iimage_file="/docker/bcfhacks.tar.gz"\
    -icmd="
      zcat ${info} | cut -f1-6 > info.txt &&
      switch_error_by_site ${merged} ${pedigree} info.txt ${out} &&
      rm info.txt
      "\
    --instance-type mem1_ssd1_v2_x8 \
    --folder=".${out_dir}" \
    --priority normal \
    --name trio_ser_by_mac_chr${CHR} -y
done


