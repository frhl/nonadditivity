#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=final_qced_array_variants
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus
#SBATCH --output=logs/final_qced_array_variants.log
#SBATCH --error=logs/final_qced_array_variants.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-23

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly in_dir="data/unphased/array/gnomad_exact_test"
readonly in_file="${in_dir}/full_c${chr}_b0_v2.b38.eur.exact_test.txt"

readonly out_dir="data/unphased/array/gnomad_exact_test"
readonly out_file_fail="${out_dir}/full_c${chr}_b0_v2.b38.eur.exact_test_failed.txt"
readonly out_file_pass="${out_dir}/ukb_c${chr}_array_variants_kept.txt"

cut -f1,9 ${in_file} | grep TRUE | cut -f1 > ${out_file_fail}
cut -f1,9 ${in_file} | grep FALSE | cut -f1 > ${out_file_pass}




