#!/usr/bin/env bash
#
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=saige_weights
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus
#SBATCH --output=logs/saige_weights.log
#SBATCH --error=logs/saige_weights.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-23

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/dominance_burden/03_saige_weights.R"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly revel_cutoff="0.773"
readonly cadd_cutoff="28.1"
readonly spliceai_cutoff="0.50"
readonly annotation="spliceai=${spliceai_cutoff}_cadd=${cadd_cutoff}_revel=${revel_cutoff}"

readonly in_dir="data/vep/vep_loftee/csqs/${annotation}"
readonly in_path="${in_dir}/UKB.chr${chr}.phased_sites.${annotation}.canonical.txt.gz"

readonly ac_pop="eur"
readonly markers_dir="data/vep/dominance_scalings"
readonly markers="${markers_dir}/UKB.wes.chr${chr}.phased.full_qc.${ac_pop}.af05.popmax.variants.group_dominance_scaling.sites.txt"

readonly ac_dir="data/phased/wes_union_calls/450k/frqx"
readonly ac_path="${ac_dir}/UKB.wes.chr${chr}.phased.full_qc.${ac_pop}.frqx.gz"

readonly min_mac=5
readonly out_dir="data/vep/saige_rescaled_weights/af05"
#readonly out_saige="${out_dir}/UKB.chr${chr}.phased_sites.${annotation}.saige.txt"
readonly out_saige="${out_dir}/UKB.wes.chr${chr}.phased.full_qc.${ac_pop}.af05.popmax.variants.${annotation}.min_mac_${min_mac}.saige.txt"

#readonly out_saige_weights="${out_dir}/UKB.chr${chr}.phased_sites.${annotation}.saige_with_${ac_pop}_weights.txt"
readonly out_saige_weights="${out_dir}/UKB.wes.chr${chr}.phased.full_qc.${ac_pop}.af05.popmax.variants.${annotation}.min_mac_${min_mac}.saige_with_${ac_pop}_weights.txt"

#readonly scaling_dir="data/vep/dominance_scalings"
#readonly scaling="${scaling_dir}/UKB.wes.chr${chr}.phased.full_qc.${ac_pop}.gt.dom.scaled.scalings.txt"
#readonly dom_scaling_max="$(cat ${scaling} | head -n1)"
#readonly dom_scaling_min="$(cat ${scaling} | tail -n1)"


mkdir -p ${out_dir}

set_up_rpy

# Run with weights. Weights are created from a beta distribution
# based on MAF specifically dbeta(MAF, 1, 25) 
Rscript ${rscript} \
    --input_path "${in_path}" \
    --markers_path "${markers}" \
    --min_mac ${min_mac} \
    --output_path "${out_saige_weights}" \
    --weights_by_ac_path "${ac_path}" \
    --delimiter " "
gzip -f ${out_saige_weights}

# without weights
Rscript ${rscript} \
    --input_path "${in_path}" \
    --markers_path "${markers}" \
    --output_path "${out_saige}" \
    --min_mac "${min_mac}" \
    --delimiter " "
gzip -f ${out_saige}







