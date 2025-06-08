# author: Frederik Lassen

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

in_dir="/Phasing/PhasingWES/step2_phase_rare/ligated/info"
out_dir="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/data/phased/wes_union_calls/450k/info/download"

threads=8

for anc in "eas" "eur" "sas"; do
  for CHR in {1..22} X; do
      prefix="UKB_chr${CHR}.exome_array.shapeit5.ligated.${anc}.tags.frqx.gz"
      in_file="${in_dir}/${prefix}"
      out_file="${out_dir}/${prefix}"
      if [[ ! -f ${out_file} ]]; then
        dx download ${in_file}
        mv ${prefix} ${out_file}
        echo "File downloaded and moved to '${out_file}'"
      else
        >&2 echo "File already exists: ${out_file}"
      fi
  done
done
