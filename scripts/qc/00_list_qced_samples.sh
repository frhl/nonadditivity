
set -o errexit
set -o nounset

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"
readonly pop_fname="superpopulation_labels.tsv"
readonly samples_pop_list="/Barney/qc2/05_estimate_superpopulation/${pop_fname}"
readonly samples_list="/Barney/qc2/09_0_final_sample_qc/09_final_qc.keep.BRaVa.sample_list"
readonly females_dir="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/data/phenotypes"
readonly females="${females_dir}/bin_matrix_eur.txt"
readonly local_dir="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/data/samples"
readonly remote_dir="/wes_ko_ukbb/data/samples"

# These are the samples used prior to phasing. Note, that
# this step does not take into account parents removed for benchmarking
readonly prephase_qc_samples_dir="/Phasing/PhasingWES/step0_merge"
readonly prephase_qc_samples_fname="overlapping_samples.qced.txt"

# these are the samples remaining after phasing:
# Subsetting pre-qced samples, those that are overlapping
# between array and WES and removing parents for benchmarking
#readonly final_phased_samples="${local_dir}/post_phasing_samples.txt"
readonly final_phased_samples="${local_dir}/UKB.wes.full_qc.post_phasing.samples"

readonly local_summary="${local_dir}/summary.txt"
readonly qced_samples_fname="09_final_qc.keep.BRaVa.sample_list"

mkdir -p ${local_dir}
dx mkdir -p ${remote_dir}
if [[ ! -f ${local_dir}/${pop_fname} ]]; then
  echo "Downloading ${samples_pop_list}.."
  dx download ${samples_pop_list} 
  mv "${pop_fname}" "${local_dir}/${pop_fname}"
fi

if [[ ! -f ${local_dir}/${qced_samples_fname} ]]; then
  echo "Downloading ${samples_list}.."
  dx download ${samples_list} 
  mv "${qced_samples_fname}" "${local_dir}/${qced_samples_fname}"
fi

if [[ ! -f ${local_dir}/${prephase_qc_samples_fname} ]]; then
  echo "Downloading ${prephase_qc_samples_fname} (pre-phased samples).."
  dx download "${prephase_qc_samples_dir}/${prephase_qc_samples_fname}" 
  mv "${prephase_qc_samples_fname}" "${local_dir}/${prephase_qc_samples_fname}"
fi

# Initialize summary file
echo "QC Summary Report" > ${local_summary}

for POP in "EUR" "AFR" "EAS" "SAS"; do
#for POP in "EUR"; do

  # subset to ancestry specific file
  pop="$(echo ${POP} | tr '[:upper:]' '[:lower:]')"
  remote_file_pop="${remote_dir}/UKB.wes.qced.${pop}.samples"
  remote_file_pop_females="${remote_dir}/UKB.wes.qced.${pop}.females.samples"
  remote_file_pop_males="${remote_dir}/UKB.wes.qced.${pop}.males.samples"
  local_file_pop="${local_dir}/UKB.wes.qced.${pop}.samples"
  local_file_pop_females="${local_dir}/UKB.wes.qced.${pop}.females.samples"
  local_file_pop_males="${local_dir}/UKB.wes.qced.${pop}.males.samples"

  # save file locally and then upload to RAP
  echo "Syncing ${pop}.. saving to ${local_summary}"

  # get samples before QC  
  cat ${local_dir}/${pop_fname} | awk -F'\t' -v pop_subset="${POP}" '$2==pop_subset' | cut -f1 | grep -vE "^-" | sort > ${local_file_pop}.tmp
  before_qc_count=$(cat ${local_file_pop}.tmp | wc -l) 
  echo "Before QC (${pop}): ${before_qc_count}" >> ${local_summary}

  # Filter 500K directly using new QC
  comm -12 <(sort "${local_dir}/${qced_samples_fname}" ) <(sort ${local_file_pop}.tmp) > ${local_file_pop}.qced
  after_qc_count=$(cat ${local_file_pop}.qced | wc -l) 
  echo "After QC (${pop}): ${after_qc_count}" >> ${local_summary}

  # filter to overlapping samples and those that a
  comm -12 <(sort ${local_file_pop}.qced) <(sort ${local_dir}/${prephase_qc_samples_fname}) > ${local_file_pop}.overlap
  after_prephase_qc_count=$(cat ${local_file_pop}.overlap | wc -l)
  echo "After QC + prephase QC + overlapping WES+Array samples (${pop}): ${after_prephase_qc_count}"  >> ${local_summary}
  
  # filter by final phased samples
  comm -12 <(sort ${local_file_pop}.overlap) <(sort ${final_phased_samples}) > ${local_file_pop}
  after_prephase_qc_count=$(cat ${local_file_pop} | wc -l)
  echo "After QC + prephase QC + overlapping WES+Array + Final phased samples (${pop}): ${after_prephase_qc_count}"  >> ${local_summary}
  dx_update_remote ${remote_file_pop} ${local_file_pop}

  # get females
  comm -12 <(sort "${local_file_pop}" ) <( cat ${females}  | cut -f1,2 | awk '$2=="0"' | cut -f1 | sort) > ${local_file_pop_females}
  echo "After QC (${pop}) females: $( cat ${local_file_pop_females} | wc -l )" >> ${local_summary}
  dx_update_remote ${remote_file_pop_females} ${local_file_pop_females}

  # get males
  comm -12 <(sort "${local_file_pop}" ) <( cat ${females}  | cut -f1,2 | awk '$2=="1"' | cut -f1 | sort) > ${local_file_pop_males}
  echo "After QC (${pop}) males: $( cat ${local_file_pop_males} | wc -l )" >> ${local_summary}
  dx_update_remote ${remote_file_pop_males} ${local_file_pop_males}

  # clean up
  rm "${local_file_pop}.tmp"

done


