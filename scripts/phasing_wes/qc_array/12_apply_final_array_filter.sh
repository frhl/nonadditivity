# author: frederik lassen

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

set -eu

readonly threads=16
readonly in_dir="/mnt/project/Phasing/PhasingSNParray/step4_liftover"
readonly out_dir="/Phasing/PhasingSNParray/step4_liftover"

# for sycning
readonly local_dir="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/data/unphased/array/gnomad_exact_test"
readonly remote_dir="Phasing/PhasingSNParray/step1_dataqc"

dx mkdir -p ${out_dir}

for CHR in {1..20} 22 X; do
  
  # update remote with local
  file="full_c${CHR}_b0_v2.b38.eur.exact_test_failed.txt"
  file_local="${local_dir}/${file}"
  file_remote="${remote_dir}/${file}"
  dx_update_remote ${file_remote} ${file_local}

  # setup paths
  snps_to_remove="/mnt/project/${file_remote}"
  array="${in_dir}/full_c$CHR\_b0_v2.b38.annotated.vcf.gz"
  out_prefix="full_c${CHR}_b0_v2.b38.annotated.filtered"
  dx run app-swiss-army-knife -icmd="
    bcftools view -e 'ID=@${snps_to_remove}' ${array} -Ob -o ${out_prefix}.bcf &&
    bcftools index ${out_prefix}.bcf
  "\
  --instance-type mem1_ssd1_v2_x4\
  --folder=".${out_dir}"\
  --name filter_gnomad_fail_chr${CHR}\
  --priority normal -y
done

