#!/usr/bin/env bash
#
# get list of parents (but not their offspring) to exclude downstream for SER calculation.
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=get_fam
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus
#SBATCH --output=logs/get_fam.log
#SBATCH --error=logs/get_fam.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh
source utils/vcf_utils.sh

readonly spark_dir="data/help/tmp/spark"
readonly hail_script="scripts/phasing_wes/phasing/00_get_fam.py"

readonly samples_dir="data/samples"
readonly overlapping_path="${samples_dir}/UKB.wes.qced.eur.overlapping_exome_array"

readonly out_dir="data/samples/trios"
readonly out_prefix="${out_dir}/ukb_450k_fam"

mkdir -p ${spark_dir}
mkdir -p ${samples_dir}

set_up_hail
set_up_pythonpath_legacy
set -x && python3 "${hail_script}" \
   --extract_samples ${overlapping_path} \
   --out_prefix "${out_prefix}"

