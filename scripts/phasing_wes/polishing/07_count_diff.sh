# author: frederik lassen

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"
set -eu

remote_dir="wes_ko_ukbb/scripts"
rscript_local="07_count_diff.R"
rscript_remote="${remote_dir}/07_count_diff.R"
dx_update_remote ${rscript_remote} ${rscript_local}


# get i/o
in_dir="/mnt/project/Phasing/PhasingWES/step4_polish_phase/bin_diff"
out_dir="Phasing/PhasingWES/step4_polish_phase/benchmark"
out_prefix="UKB.exome_array.shapeit5.sapphire.overview"
dx mkdir -p ${out_dir}
dx run app-swiss-army-knife \
  -iimage_file="/docker/rsuite.tar.gz"\
  -icmd="Rscript /mnt/project/${rscript_remote} --in_dir ${in_dir} --out_prefix ${out_prefix} && echo 'done'" \
  --instance-type mem1_ssd1_v2_x8 \
  --folder=".${out_dir}" \
  --priority normal \
  --name summarize_bin_diff -y





