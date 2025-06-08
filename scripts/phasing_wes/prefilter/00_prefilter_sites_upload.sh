


set +eu
bmrc_bim_dir="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/data/unphased/variants/qced_variants"
rap_bim_dir="/Phasing/PhasingWES/step0_merge/support"
dx mkdir -p ${rap_bim_dir}
for CHR in  X; do
  bmrc_path="${bmrc_bim_dir}/ukb_wes_450k.qced.prefilter.chr${CHR}.bim"
  rap_path="${rap_bim_dir}/ukb_wes_450k.qced.prefilter.chr${CHR}.bim"
  dx upload ${bmrc_path} --path ${rap_path}
done

