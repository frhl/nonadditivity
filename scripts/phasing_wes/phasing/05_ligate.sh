# author: Frederik Lassen

 set -o nounset


# to get paths in correct order
module load R
rscript="/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus/scripts/phasing_wes/phasing/05_ligate.R"

# get directories
in_dir="/Phasing/PhasingWES/step2_phase_rare_qc/chunks"
#in_dir="/Phasing/PhasingWES/step2_phase_rare/chunks"
out_dir="/Phasing/PhasingWES/step2_phase_rare_qc/ligated"
#out_dir="/Phasing/PhasingWES/step2_phase_rare/ligated"
threads=16

dx mkdir -p ${out_dir}

set -eu
set -x
for CHR in 17; do
    #file_prefix="${in_dir}/UKB.chr${CHR}.phased.exome_array" 
    file_prefix="${in_dir}/UKB.chr${CHR}.phased.exome_array.full_qc" 
    files=$(Rscript ${rscript} --in_prefix ${file_prefix} --out_prefix "/mnt/project")
    out="UKB_chr${CHR}.exome_array.shapeit5.ligated.bcf"
    dx run app-swiss-army-knife -icmd="
      bcftools concat --threads ${threads} --ligate ${files} -Ob -o ${out} &&
      bcftools index --threads ${threads} ${out} && echo 'done!'
      " \
      --instance-type mem2_ssd1_v2_x64 \
      --folder=".${out_dir}" \
      --name ligate_chr${CHR} \
      --priority normal -y
done



