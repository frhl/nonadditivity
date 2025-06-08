


source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

# get dirs
readonly variant_list="/Barney/qc2/08_0_final_variant_qc/08_final_qc.pop.keep.variant_list"
readonly local_dir="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/data/unphased/variants/qced_variants"
readonly remote_dir="/wes_ko_ukbb/data/variants/qced"

readonly summary="${local_dir}/summary.txt"
echo "QC Summary Report" > ${summary}

for POP in "EUR" "AFR" "EAS" "SAS"; do

  # keep track of filenames
  pop="$(echo ${POP} | tr '[:upper:]' '[:lower:]')"
  file_name="UKB.wes.qced_variants.${pop}.txt"
  local_file="${local_dir}/${file_name}"
  remote_file="${remote_dir}/${file_name}"
  
  # keep variants locally
  echo "Saving variants to ${local_file}.."
  dx cat ${variant_list} | awk -F'\t' -v pop="$POP" '$1 == pop && NF==3 { gsub(/\[|\]|\"/, "", $3); split($3,a,","); print $2":"a[1]":"a[2] }' > ${local_file}
  
  # quality check the file
  for chr in {1..22} X; do
    pop_n_variants=$(cat ${local_file} | grep "chr${chr}:" | wc -l)
    echo "chr${chr} (${pop}) - ${pop_n_variants}" >> ${summary}
    if [[ "${pop_n_variants}" -lt 100 ]]; then
      >&2 echo "Error: chr${chr} (${pop} does not have all variants accounted for!)"
      rm -f ${local_file}
      exit 1 
    fi
  done
  
  # update remote file if it changed
  dx_update_remote ${remote_file} ${local_file}
done


