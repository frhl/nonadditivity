# author: frederik lassen

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

set -eu
out_dir="/Phasing/PhasingWES/step3_benchmark/trio_merge/polished/"
in_dir_info="/mnt/project/Phasing/PhasingWES/step2_phase_rare/ligated/info"
in_dir_merged="/mnt/project/Phasing/PhasingWES/step3_benchmark/trio_merge/polished"
in_dir_pedigree="/mnt/project/Phasing/PhasingWES/step0_merge"

dx mkdir -p ${out_dir}

info="${in_dir_info}/UKB_chrCHR.exome_array.shapeit5.ligated.eur.tags.txt.gz"
merged="${in_dir_merged}/UKB.autosomes.exome_array.shapeit5.ligated.polished.pp_lt_90.trio.vcf.gz"
pedigree="${in_dir_pedigree}/ukb_wes_200k_trios.ped"
out="UKB.autosomes.exome_array.shapeit5.polished.ligated.pp_lt_90.trio_ser_by_site"

if [[ $(dx_file_exists "${out_dir}/${out}") -eq "0" ]]; then
  dx run app-swiss-army-knife \
    -iimage_file="/docker/bcfhacks.tar.gz"\
    -icmd="
      zcat '${info/CHR/1}' | cut -f1-6 > info.txt &&
      for idx in {1..22}; do zcat ${info/CHR/\$idx}| cut -f1-6 | tail -n+2 >> info.txt; done &&
      switch_error_by_site ${merged} ${pedigree} info.txt ${out} &&
      rm info.txt && 
      echo 'ok!'
      "\
    --instance-type mem1_ssd1_v2_x16 \
    --folder=".${out_dir}" \
    --priority normal \
    --name trio_ser_by_mac_after_pp -y
fi


