#!/bin/bash


source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"
set -eu

# keep track of updates to rscript
remote_dir="wes_ko_ukbb/scripts"
rscript_local="11_plot_indiv_chrom.R"
rscript_remote="${remote_dir}/11_plot_indiv_chrom.R"
dx_update_remote ${rscript_remote} ${rscript_local}

in_dir="/mnt/project/Phasing/PhasingWES/step3_benchmark/trio_merge/new"
out_dir="/Phasing/PhasingWES/step3_benchmark/plots"

switch_low_pp="/mnt/project/Phasing/PhasingWES/step3_benchmark/new/"
switch_high_pp="/mnt/project/Phasing/PhasingWES/step3_benchmark/pp_lt_90/new/"
pattern_low_pp="switches$"
pattern_high_pp="ser_by_site$"

out_prefix="UKB.indiv_chrom.exome_array.shapeit5.ligated"

dx mkdir -p ${out_dir}

dx run app-swiss-army-knife \
  -iimage_file="/docker/rsuite.tar.gz"\
  -icmd="
     Rscript /mnt/project/${rscript_remote} \
       --dir_switch_low ${switch_low_pp} \
       --pattern_switch_low ${pattern_low_pp} \
       --dir_switch_high ${switch_high_pp} \
       --pattern_switch_high ${pattern_high_pp} \
       --out_prefix ${out_prefix}
    "\
  --instance-type mem2_ssd1_v2_x32 \
  --folder=".${out_dir}" \
  --priority normal \
  --name plot_indiv_chrom -y




