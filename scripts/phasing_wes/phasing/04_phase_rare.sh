# author: Frederik Lassen

set +eu

# setup directories
#bmrc_regions_dir="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/data/phased/wes_union_calls/450k/shapeit5/intervals_qced"
readonly bmrc_regions_dir="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/data/phased/wes_union_calls/450k/shapeit5/intervals_full_qc_new"
readonly scaffold_dir="/mnt/project/Phasing/PhasingWES/step1_phase_common_qc"
readonly bcf_dir="/mnt/project/Phasing/PhasingWES/step0_merge_qc"
readonly out_dir="/Phasing/PhasingWES/step2_phase_rare_qc/chunks"

readonly samples_dir="/mnt/project/wes_ko_ukbb/data/samples"
readonly males="${samples_dir}/UKB.wes.qced.combined.males.samples"
readonly pedigree="${samples_dir}/ukb_all_trios.formatted.ped"

#readonly THREADS=192
readonly THREADS=96

dx mkdir -p ${out_dir}

# rare variant phasing for chunks
for CHR in 17; do
  MAP=/mnt/project/data/shapeit_maps/chr${CHR}.b38.gmap.gz
  path_phasing_region="${bmrc_regions_dir}/phasing_intervals_u1000_s500000_o25000_chr${CHR}.tsv"
  path_scaffold_region="${bmrc_regions_dir}/scaffold_intervals_u1000_s500000_o25000_chr${CHR}.tsv"
  path_trimming_region="${bmrc_regions_dir}/scaffold_intervals_u1000_s500000_o25000_chr${CHR}.tsv"
  #path_phasing_region="${bmrc_regions_dir}/phasing_intervals_u1000_s300000_o25000_chr${CHR}.tsv"
  #path_scaffold_region="${bmrc_regions_dir}/scaffold_intervals_u1000_s300000_o25000_chr${CHR}.tsv"
  #path_trimming_region="${bmrc_regions_dir}/scaffold_intervals_u1000_s300000_o25000_chr${CHR}.tsv"
  TOTAL_CHUNKS=$( wc -l ${path_phasing_region} | awk '{print $1}' )

  # get haploids
  haploid_expr=""
  if [[ "${CHR}" == "X" ]]; then
    haploid_expr="--haploids ${males}"
  fi

  for CURRENT_CHUNK in $(seq 1 ${TOTAL_CHUNKS}); do
     
      # Need to add 1 to line count because of header 
      trim_line=$(sed -n $((${CURRENT_CHUNK}+1))p ${path_trimming_region})
      scaffold_line=$(sed -n $((${CURRENT_CHUNK}+1))p ${path_scaffold_region})
      phasing_line=$(sed -n $((${CURRENT_CHUNK}))p ${path_phasing_region})

      # setup phasing chunks
      CHUNK_NAME="${CURRENT_CHUNK}of${TOTAL_CHUNKS}"
      INPUT_REG=$(echo $phasing_line)
      INPUT_REG_START=$(echo ${INPUT_REG} | cut -d":" -f 2 | cut -d"-" -f1)
      INPUT_REG_END=$(echo ${INPUT_REG} | cut -d":" -f 2 | cut -d"-" -f2)
      INPUT_REG_NAME=${INPUT_REG_START}_${INPUT_REG_END}
      
      # setup scaffolds
      SCAFFOLD_REG=$(echo $scaffold_line | cut -d" " -f5)
      SCAFFOLD_REG_START=$(echo ${SCAFFOLD_REG} | cut -d":" -f 2 | cut -d"-" -f1)
      SCAFFOLD_REG_END=$(echo ${SCAFFOLD_REG} | cut -d":" -f 2 | cut -d"-" -f2)
      SCAFFOLD_REG_NAME=${SCAFFOLD_REG_START}_${SCAFFOLD_REG_END}
 
      # setup trimming of edges
      TRIM_REG=$(echo $trim_line | cut -d" " -f5)
      TRIM_REG_START=$(echo ${TRIM_REG} | cut -d":" -f 2 | cut -d"-" -f1)
      TRIM_REG_END=$(echo ${TRIM_REG} | cut -d":" -f 2 | cut -d"-" -f2)
      TRIM_REG_NAME=${TRIM_REG_START}_${TRIM_REG_END}
      
      # i/o
      BCF_TO_PHASE=${bcf_dir}/UKB.chr${CHR}.exome_array.full_qc.no_parents.bcf
      SCAFFOLD=${scaffold_dir}/UKB_chr${CHR}.exome_array.full_qc.no_parents.shapeit5.common.bcf
      TMP=UKB.chr${CHR}.phased.exome_array.full_qc.no_parents.shapeit5.common.rare.b${CHUNK_NAME}.pr${INPUT_REG_NAME}.sr${SCAFFOLD_REG_NAME}.bcf
      LOG=UKB.chr${CHR}.phased.exome_array.full_qc.no_parents.shapeit5.common.rare.b${CHUNK_NAME}.pr${INPUT_REG_NAME}.sr${SCAFFOLD_REG_NAME}.log
      TIM=UKB.chr${CHR}.phased.exome_array.full_qc.no_parents.shapeit5.common.rare.b${CHUNK_NAME}.pr${INPUT_REG_NAME}.sr${SCAFFOLD_REG_NAME}.time
      OUT=UKB.chr${CHR}.phased.exome_array.full_qc.no_parents.shapeit5.common.rare.b${CHUNK_NAME}.pr${INPUT_REG_NAME}.sr${SCAFFOLD_REG_NAME}.tr${TRIM_REG_NAME}.bcf
      found=$( dx ls ${out_dir} | grep ${OUT} | wc -l | cut -d" " -f1)
      echo "###"
      echo ${OUT}
      echo ${BCF_TO_PHASE}
      if [[ "${found}" -eq "0" ]]; then 
       echo "pr: $INPUT_REG"
       echo "sr: $SCAFFOLD_REG"
       echo "tr: $TRIM_REG"
       dx run app-swiss-army-knife -iimage_file="/docker/shapeit5_v5.1.1.docker.tar.gz" --folder=".${out_dir}" -icmd="/usr/bin/time -vo $TIM phase_rare_static --input $BCF_TO_PHASE --scaffold $SCAFFOLD --map $MAP --pedigre ${pedigree} ${haploid_expr} --output $TMP --log $LOG --scaffold-region $SCAFFOLD_REG --input-region $INPUT_REG --pbwt-depth-rare 4 --pbwt-depth-common 4 --thread $THREADS && bcftools index -f $TMP --threads $THREADS && bcftools view -r $TRIM_REG $TMP -Ob -o $OUT && bcftools index --threads $THREADS $OUT " \
              --instance-type mem2_ssd1_v2_x32 \
              --depends-on job-GjPyX70Jg8JXJb0yPPJ2K6qk \
              --priority normal \
              --name shapeit5_rare_c${CHR}_${CHUNK_NAME} -y
    else
      echo "${OUT} already exists. Skipping!"
    fi
    done
done



