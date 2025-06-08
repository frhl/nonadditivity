#dx build -f saige-universal-step-2

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

set -o errexit
set -o nounset

in_dir="/Phasing/PhasingWES/step4_polish_phase/pp_extract"


for chr in X; do
  out_dir="/Phasing/PhasingWES/step4_polish_phase/bin_splitter/new/chr${chr}"
  in_bin="${in_dir}/UKB_chr${chr}.exome_array.shapeit5.ligated.bcf_hets.bin"
  split_size=10000
  dx mkdir -p ${out_dir}
  echo ${in_bin}
  dx run bin-splitter-applet \
    -i binary_file_to_split="${in_bin}" \
    -i split_size="${split_size}" \
    --instance-type "mem2_ssd1_v2_x4" --priority high --destination ".${out_dir}" -y --name "c${chr}_bin_splitter"
done

