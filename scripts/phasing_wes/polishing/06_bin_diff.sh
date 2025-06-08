
module purge && module load jq

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

set -o errexit
set -o nounset


for chr in {1..22} X; do

  # get sample IDs
  in_dir_pp_extract="/mnt/project/Phasing/PhasingWES/step4_polish_phase/pp_extract"
  samples_id="${in_dir_pp_extract}/UKB_chr${chr}.exome_array.shapeit5.ligated.bcf.samples.txt"
  vcf_id="${in_dir_pp_extract}/UKB_chr${chr}.exome_array.shapeit5.ligated.bcf_vars.bcf"

  # original
  in_dir_bin1="/mnt/project/Phasing/PhasingWES/step4_polish_phase/pp_extract"
  bin1_id="${in_dir_bin1}/UKB_chr${chr}.exome_array.shapeit5.ligated.bcf_hets.bin"  

  # get re-phased
  in_dir_bin2="/mnt/project/Phasing/PhasingWES/step4_polish_phase/phase_calling"
  bin2_id="${in_dir_bin2}/UKB_chr${chr}.exome_array.shapeit5.ligated.bcf_hets.bin_sub_merged.bin"

  out_dir="/Phasing/PhasingWES/step4_polish_phase/bin_diff"
  out_file="UKB_chr${chr}.exome_array.shapeit5.ligated.bcf_hets.bin_diff.csv"
  docker_image="docker/pp_toolkit_v1.4.tar.gz"

  echo "VCF/BCF: ${vcf_id}"
  echo "samples: ${samples_id}"
  echo "Binary 1: ${bin1_id}"
  echo "Binary 2: ${bin2_id}"

  instance="mem2_ssd1_v2_x8"
  threads="24" # 3xinstance

  docker_script_dir="/usr/src/pp/bin_tools"
  script="${docker_script_dir}/bin_diff"

  dx mkdir -p ${out_dir}
  dx run app-swiss-army-knife \
    -iimage_file="/docker/pp_toolkit_v1.4.tar.gz"\
    -icmd="
      ${script} --vcf-file ${vcf_id} --samples ${samples_id} --binary1 ${bin1_id} --binary2 ${bin2_id} --more --extra-info --ac  > ${out_file} &&
       gzip ${out_file} 
     " \
    --instance-type mem1_ssd1_v2_x8 \
    --folder=".${out_dir}" \
    --priority high \
    --name c${chr}_bin_diff -y
done



