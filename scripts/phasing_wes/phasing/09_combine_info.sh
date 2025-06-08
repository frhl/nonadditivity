#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=combine_info_text
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus
#SBATCH --output=logs/combine_info_text.log
#SBATCH --error=logs/combine_info_text.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 2

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly rscript="scripts/phasing_wes/phasing/10_combine_info_text.R"

readonly input_pattern="tags"
readonly input_dir="data/phased/wes_union_calls/450k/info/download"

readonly out_dir="data/phased/wes_union_calls/450k/info"
readonly out_prefix="${out_dir}/UKB.exome_array.superpop"

readonly cols="ID,MAF,AN,AC,AC_Het,AC_Hom"
readonly cols_maf="ID,MAF"

set_up_rpy
Rscript ${rscript} \
  --input_dir ${input_dir} \
  --input_pattern ${input_pattern} \
  --columns ${cols_maf} \
  --out_prefix ${out_prefix}.only_maf

Rscript ${rscript} \
  --input_dir ${input_dir} \
  --input_pattern ${input_pattern} \
  --columns ${cols} \
  --out_prefix ${out_prefix}.full







