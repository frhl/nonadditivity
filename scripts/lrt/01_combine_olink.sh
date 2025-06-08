# author: frederik lassen

# note: this script uses a github repo to get the genests to be processed https://github.com/frhl/genesets_public

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

readonly remote_dir="wes_ko_ukbb/scripts"
readonly rscript_local="01_combine_olink.R"
readonly rscript_remote="${remote_dir}/01_combine_olink.R"
dx_update_remote ${rscript_remote} ${rscript_local}

readonly group="spliceai=0.50_cadd=28.1_revel=0.773"
readonly pop="eur"
readonly pp="0.90"
readonly af="0.50"

readonly in_dir="/mnt/project/wes_ko_ukbb/data/lrt_testing/lrt_unrel/olink/2df_chi2_trick_final"
readonly out_dir="/wes_ko_ukbb/data/lrt_testing/lrt_unrel/olink/"
readonly out_file="UKB.wes.merged.phased.qc_final.eur.af${af}.popmax.pp${pp}.${group}.combined.unrel.chi2_trick.txt.gz"
readonly pattern="UKB.wes.merged.phased.qc_final.eur"

dx run app-swiss-army-knife \
  -iimage_file="/docker/rsuite.tar.gz"\
  -icmd="
     Rscript /mnt/project/${rscript_remote} \
        --in_dir ${in_dir} \
        --pattern ${pattern} \
        --out_path ${out_file} &&
       echo '!!!!$(date)'
    "\
  --instance-type mem2_ssd1_v2_x2 \
  --folder=".${out_dir}" \
  --priority high \
  --name olink_combine -y





