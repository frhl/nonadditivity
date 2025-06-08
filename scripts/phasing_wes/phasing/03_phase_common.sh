# author: Frederik Lassen

set -eu

readonly MAF=0.001
readonly THREADS=180

readonly in_dir="Phasing/PhasingWES/step0_merge_qc"
readonly in_dir_mnt="/mnt/project/${in_dir}"
readonly out_dir="/Phasing/PhasingWES/step1_phase_common_qc"

readonly samples_dir="/mnt/project/wes_ko_ukbb/data/samples"
readonly males="${samples_dir}/UKB.wes.qced.combined.males.samples"
readonly pedigree="${samples_dir}/ukb_all_trios.formatted.ped"

dx mkdir -p ${out_dir}
for CHR in 17; do
  MAP=/mnt/project/data/shapeit_maps/chr${CHR}.b38.gmap.gz
  BCF="${in_dir_mnt}/UKB.chr${CHR}.exome_array.full_qc.no_parents.bcf"
  OUT=UKB_chr${CHR}.exome_array.full_qc.no_parents.shapeit5.common.bcf
  LOG=UKB_chr${CHR}.exome_array.full_qc.no_parents.shapeit5.common.log
  TIM=UKB_chr${CHR}.exome_array.full_qc.no_parents.shapeit5.common.time
  dx ls ${out_dir}/${OUT}
  haploid_expr=""
  if [[ "${CHR}" == "X" ]]; then
    haploid_expr="--haploids ${males}"
  fi
  dx run app-swiss-army-knife \
    -iimage_file="/docker/shapeit5_v5.1.1.docker.tar.gz" \
    -icmd="/usr/bin/time -vo $TIM phase_common_static --input $BCF --map $MAP --pedigre ${pedigree} ${haploid_expr} --output $OUT --thread $THREADS --log $LOG --filter-maf $MAF --region chr${CHR} && bcftools index -f $OUT --threads $THREADS && echo 'done'" \
    --folder=".${out_dir}" \
    --instance-type mem2_ssd1_v2_x64 \
    --depends-on job-GjQ44G0Jg8JQ6V4XF9p2j0vK \
    --priority normal \
    --name shapeit5_full_qc_common_chr${CHR} -y 
done



