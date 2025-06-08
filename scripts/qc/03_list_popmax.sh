#!/usr/bin/env bash

#SBATCH --account=lindgren.prj
#SBATCH --job-name=list_popmaxs
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus
#SBATCH --output=logs/list_popmaxs.log
#SBATCH --error=logs/list_popmaxs.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1

set -o errexit
set -o nounset

source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly spark_dir="data/tmp/spark"

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly in_dir="/well/lindgren/UKBIOBANK/dpalmer"
readonly in="${in_dir}/gnomad.exomes.r2.1.1.sites.liftover_grch38.ht"

readonly cutoff="0.10"

readonly out_dir="data/vep/popmax"
readonly out_prefix="${out_dir}/gnomad.exomes.r2.1.1.grch38.popmax${cutoff}"

readonly hail_script="scripts/qc/03_list_popmax.py"

mkdir -p ${out_dir}
mkdir -p ${spark_dir}

set_up_hail 0.2.97
set_up_pythonpath_legacy
python3 ${hail_script} \
     --input_path "${in}" \
     --cutoff ${cutoff} \
     --out_prefix "${out_prefix}"
