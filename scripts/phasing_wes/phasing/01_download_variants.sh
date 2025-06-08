# download variants only for interval creation

target_dir="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/data/unphased/variants/full_qc_new"
mkdir -p ${target_dir}
for CHR in 17; do
  download_path="Phasing/PhasingWES/step0_merge_qc/UKB.chr${CHR}.exome_array.full_qc.no_parents.variants_only.vcf.gz"
  file=$( basename ${download_path} )
  dx download ${download_path}
  #dx download ${download_path}.csi
  mv ${file} ${target_dir}/.
  #mv ${file}.csi ${target_dir}/.
done

