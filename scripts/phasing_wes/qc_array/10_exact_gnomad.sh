#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=exact_gnomad
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus
#SBATCH --output=logs/exact_gnomad.log
#SBATCH --error=logs/exact_gnomad.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-23

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/phasing_wes/qc_array/10_exact_gnomad.R"

readonly task_id=$( get_array_task_id )
readonly chr=$( get_chr ${task_id} )

readonly gnomad_dir="data/gnomad/genomes"
readonly gnomad="${gnomad_dir}/gnomad.genomes.v4.0.sites.chr${chr}.counts.txt.gz"

readonly array_dir="data/unphased/array"
readonly array="${array_dir}/full_c${chr}_b0_v2.b38.annotated.eur.frqx"

readonly out_dir="data/unphased/array/gnomad_exact_test"
readonly out_prefix="${out_dir}/full_c${chr}_b0_v2.b38.eur.exact_test"

mkdir -p ${out_dir}

set_up_rpy
Rscript ${rscript} \
  --gnomad_path ${gnomad} \
  --array_path ${array} \
  --out_prefix ${out_prefix}



