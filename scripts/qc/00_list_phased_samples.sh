# author: frederik lassen

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"

readonly out_dir="/wes_ko_ukbb/data/samples"
readonly in_dir="/mnt/project/Phasing/PhasingWES/step2_phase_rare_qc/ligated"
readonly pop="eur"
readonly input_vcf="${in_dir}/UKB_chr21.exome_array.shapeit5.ligated.bcf"
readonly out_prefix="UKB.wes.full_qc.post_phasing.samples"
dx run app-swiss-army-knife \
  -icmd="
  bcftools query -l ${input_vcf} > ${out_prefix}
"\
  --instance-type mem1_ssd1_v2_x8 \
  --folder=".${out_dir}" \
  --priority high \
  --name phased_samples -y

