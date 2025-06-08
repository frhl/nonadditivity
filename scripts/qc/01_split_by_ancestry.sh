# author: frederik lassen

source "/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/utils/dx_utils.sh"


set -eu
readonly threads=24
readonly out_dir="/wes_ko_ukbb/data/phased/qced_samples_and_sites"
readonly in_dir="/mnt/project/Phasing/PhasingWES/step2_phase_rare_qc/ligated"

dx mkdir -p ${out_dir}
#for CHR in {1..22} X; do
for CHR in {1..22} X; do
  for pop in "eur"; do
    sample_qc="/mnt/project/wes_ko_ukbb/data/samples/UKB.wes.qced.${pop}.samples"
    variant_qc="/mnt/project/wes_ko_ukbb/"
    input_vcf="${in_dir}/UKB_chr${CHR}.exome_array.shapeit5.ligated.bcf"
    out_prefix="UKB.wes.chr${CHR}.phased.full_qc.${pop}"
      dx run app-swiss-army-knife \
        -icmd="
             bcftools view ${input_vcf} --threads ${threads} -S ${sample_qc} --force-samples -Ob -o ${out_prefix}.bcf &&
             bcftools index ${out_prefix}.bcf
      "\
        --instance-type mem1_ssd1_v2_x8 \
        --folder=".${out_dir}" \
        --priority high \
        --name c${CHR}_${pop}_apply_qc -y
  done
done

