
module purge && module load jq

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

set -o errexit
set -o nounset

script_dir="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/scripts/phasing_wes/polishing/pp/dnanexus_scripts"
script="${script_dir}/run_phase_caller_batch.sh"

# get project ID
project_id="23143"

echo "sleeping 2h.."
sleep 2h

for chr in X; do

  # get sample IDs
  in_dir_pp_extract="Phasing/PhasingWES/step4_polish_phase/pp_extract"
  samples_id="${in_dir_pp_extract}/UKB_chr${chr}.exome_array.shapeit5.ligated.bcf.samples.txt"
  vcf_id="${in_dir_pp_extract}/UKB_chr${chr}.exome_array.shapeit5.ligated.bcf_vars.bcf"

  # get qualifying example for binary ID
  in_dir_bin_split="Phasing/PhasingWES/step4_polish_phase/bin_splitter/new/chr${chr}"
  bin_id="${in_dir_bin_split}/UKB_chr${chr}.exome_array.shapeit5.ligated.bcf_hets.bin_sub_1"

  # check CRAM dir
  # note, that i mofified ${script} to take the exome path directly to avoid annoying 
  # related to having spaces in the filename.
  #cram_dir="/mnt/project/Bulk/Exome sequences/Exome OQFE CRAM files"

  out_dir="/Phasing/PhasingWES/step4_polish_phase/phase_calling/chunks/chr${chr}"
  docker_image="docker/pp_toolkit_v1.4.tar.gz"

  instance="mem2_ssd1_v2_x8"
  threads="24" # 3xinstance

  dx mkdir -p ${out_dir}

  files_found="$(dx ls ${out_dir} | wc -l)"
  files_expected="42"
  if [[ ! "${files_found}" -eq "${files_expected}" ]]; then
    bash $script \
      --project-id ${project_id} \
      --destination ${out_dir} \
      --samples-id ${samples_id} \
      --vcf-id ${vcf_id} \
      --bin-id ${bin_id} \
      --docker-image ${docker_image} \
      --instance ${instance} \
      --threads ${threads} \
      --tag chr${chr}
  else
    >&2 echo "All ${files_expected} output files accounted for (${out_dir})! Skipping chrom."
  fi
done



