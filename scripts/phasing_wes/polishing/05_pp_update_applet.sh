
source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

set -o errexit
set -o nounset

in_dir_original_vcf="/Phasing/PhasingWES/step2_phase_rare/ligated"
in_dir_rephased_bin="/Phasing/PhasingWES/step4_polish_phase/phase_calling"
out_dir="/Phasing/PhasingWES/step4_polish_phase/rephased"
dx mkdir -p ${out_dir}
for chr in 18; do
  original_vcf_file="${in_dir_original_vcf}/UKB_chr${chr}.exome_array.shapeit5.ligated.bcf"
  rephased_bin_file="${in_dir_rephased_bin}/UKB_chr${chr}.exome_array.shapeit5.ligated.bcf_hets.bin_sub_merged.bin"
  out_vcf_file="UKB_chr${chr}.exome_array.shapeit5.ligated.polished.bcf"
  if [[ $(dx_file_exists "${out_dir}/${out_vcf_file}") -eq 0 ]]; then
    dx run pp-update-applet \
      -i original_vcf_file="${original_vcf_file}" \
      -i rephased_binary_file="${rephased_bin_file}" \
      -i output_vcf_filename="${out_vcf_file}" \
      --instance-type "mem2_ssd1_v2_x8" --priority high --destination ".${out_dir}" -y --name "c${chr}_update_pp"
  else
    >&2 echo "${out_vcf_file} already exists. Skipping."
  fi
done

