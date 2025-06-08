
readonly in_dir="/Phasing/PhasingSNParray/step4_liftover"
readonly out_dir="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/data/unphased/array"

mkdir -p ${out_dir}
set -x
for CHR in {1..22} X; do
  fname="full_c${CHR}_b0_v2.b38.annotated.eur.frqx"
  echo "Getting ${fname}.."
  dx cat ${in_dir}/${fname} > ${out_dir}/${fname}
done


