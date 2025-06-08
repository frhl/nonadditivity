# author: frederik lassen

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

set -eu
readonly out_dir="/Phasing/PhasingWES/step3_benchmark/pp_gt_90"
readonly in_dir_info="/mnt/project/Phasing/PhasingWES/step2_phase_rare/ligated/info"
readonly in_dir_merged="/mnt/project/Phasing/PhasingWES/step3_benchmark/pp_gt_90"
readonly in_dir_pedigree="/mnt/project/Phasing/PhasingWES/step0_merge"
dx mkdir -p ${out_dir}

for CHR in 1; do
  
  info="${in_dir_info}/UKB_chr${CHR}.exome_array.shapeit5.ligated.eur.tags.txt.gz"
  merged="${in_dir_merged}/UKB.chr${CHR}.exome_array.shapeit5.ligated.polished.pp_lt_90.trio.vcf.gz"
  pedigree="${in_dir_pedigree}/ukb_wes_200k_trios.ped"
  out="UKB.chr${CHR}.exome_array.shapeit5.ligated.polished.pp_lt_90.trio.ser_by_site"

  echo "## INPUT FILE STATUS"
  echo "INFO = ${info} ($(dx_size_in_bytes $(wo_mnt_project ${info}))B)"
  echo "MERGED = ${merged} ($(dx_size_in_bytes $(wo_mnt_project ${merged}))B)"
  echo "PED = ${pedigree} ($(dx_size_in_bytes $(wo_mnt_project ${pedigree}))B)"

  if [[ $(dx_file_exists "${out_dir}/${out}") -eq "0" ]]; then
    dx run app-swiss-army-knife \
      -iimage_file="/docker/bcfhacks.tar.gz"\
      -icmd="
        zcat ${info} | cut -f1-6 > info.txt &&
        head info.txt &&
        switch_error_by_site ${merged} ${pedigree} info.txt ${out} &&
        rm info.txt
        "\
      --instance-type mem1_ssd1_v2_x4 \
      --folder=".${out_dir}" \
      --priority normal \
      --name trio_ser_by_mac_after_pp_chr${CHR} -y
  fi
done


