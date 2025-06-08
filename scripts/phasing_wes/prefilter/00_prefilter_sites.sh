#!/usr/bin/env bash
#
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=prefilter_sites
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus
#SBATCH --output=logs/prefilter_sites.log
#SBATCH --error=logs/prefilter_sites.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-22

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly cluster=$( get_current_cluster )
readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly spark_dir="data/help/tmp/spark"
readonly hail_script="scripts/phasing_wes/prefilter/00_prefilter_sites.py"

readonly input_dir="data/unphased/variants/qced_variants"
readonly input_path="${input_dir}/UKB.wes.qced_variants.combined.txt"

readonly target_path_dir="data/support/target_regions"
readonly target_path="${target_path_dir}/xgen_plus_spikein.GRCh38.bed"

readonly lcr_path_dir="data/support/lcr"
readonly lcr_path="${lcr_path_dir}/btu356_LCR-hs38.bed"

readonly out_dir="data/unphased/variants/qced_variants"
readonly out_prefix="${out_dir}/UKB.wes.qced_variants.combined.prefilter.chr${chr}"

readonly target_padding=50

mkdir -p ${spark_dir}

set_up_hail
set_up_pythonpath_legacy
set -x && python3 "${hail_script}" \
   --input_path ${input_path} \
   --lcr_path ${lcr_path} \
   --target_path ${target_path} \
   --target_padding ${target_padding} \
   --out_prefix "${out_prefix}"

