# https://github.com/rwk-unil/pp/blob/main/dnanexus_app/bin-splitter-applet/dxapp.json

#dx build -f saige-universal-step-2

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

set -o errexit
set -o nounset

in_dir="/Phasing/PhasingWES/step4_polish_phase/rephased"
out_dir="/Phasing/PhasingWES/step5_polish_phase_wgs/pp_extract"

dx mkdir -p ${out_dir}

for chr in 18; do
  in_bcf="${in_dir}/UKB_chr${chr}.exome_array.shapeit5.ligated.polished.bcf"
  dx run pp-extract-applet \
    -i vcf_bcf_file="${in_bcf}" \
    --instance-type "mem2_ssd1_v2_x4" --priority high --destination ".${out_dir}" -y --name "c${chr}_extract_pp_wgs"
done


