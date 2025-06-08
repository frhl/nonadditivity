# author: frederik lassen

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"
set -eu

remote_dir="wes_ko_ukbb/scripts"
rscript_local="03_get_ps_conf.R"
rscript_remote="${remote_dir}/03_get_ps_conf.R"
dx_update_remote ${rscript_remote} ${rscript_local}
out_dir="/wes_ko_ukbb/data/phased/calls/sites/polished/pp"
in_dir="/mnt/project/wes_ko_ukbb/data/phased/calls/sites/polished"

dx mkdir -p ${out_dir}
for CHR in {1..22} X; do
  for af in "01" "05"; do
     geno="${in_dir}/UKB.chr${CHR}.exome_array.shapeit5.sapphire.sites.af${af}.txt.gz"
     out_prefix="UKB.chr${CHR}.af${af}.shapeit5.sapphire.sites_by_pp"
     dx run app-swiss-army-knife \
          -iimage_file="/docker/rsuite.tar.gz"\
          -icmd="
            Rscript /mnt/project/${rscript_remote} \
             --input_file ${geno} \
             --n_breaks 100 \
             --ac_max 1 \
             --out_prefix ${out_prefix} &&
            echo 'done' 
          " \
          --instance-type mem1_ssd1_v2_x8 \
          --folder=".${out_dir}" \
          --priority normal \
          --name get_ps_conf_chr${CHR} -y
   done
done





