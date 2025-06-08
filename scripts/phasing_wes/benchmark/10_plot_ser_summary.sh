#!/bin/bash


source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"
set -eu

# keep track of updates to rscript
remote_dir="wes_ko_ukbb/scripts"
rscript_local="10_plot_ser_summary.R"
rscript_remote="${remote_dir}/10_plot_ser_summary.R"
dx_update_remote ${rscript_remote} ${rscript_local}

in_dir="/mnt/project/Phasing/PhasingWES/step3_benchmark/trio_merge/new"
out_dir="/Phasing/PhasingWES/step3_benchmark/plots"

switch_low_pp="${in_dir}/UKB.autosomes.exome_array.shapeit5.ligated.trio_ser_by_site.by.mac"
switch_high_pp="${in_dir}/UKB.autosomes.exome_array.shapeit5.ligated.pp_lt_90.trio_ser_by_site.by.mac"
out_prefix="UKB.autosomes.exome_array.shapeit5.ligated.trio_ser_by_mac"

dx mkdir -p ${out_dir}

dx run app-swiss-army-knife \
  -iimage_file="/docker/rsuite.tar.gz"\
  -icmd="
     Rscript /mnt/project/${rscript_remote} \
       --input_switch_low ${switch_low_pp} \
       --input_switch_high ${switch_high_pp} \
       --out_prefix ${out_prefix}
    "\
  --instance-type mem1_ssd1_v2_x8 \
  --folder=".${out_dir}" \
  --priority normal \
  --name plot_ser_summary -y




