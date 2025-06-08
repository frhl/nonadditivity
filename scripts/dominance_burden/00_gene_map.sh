#!/usr/bin/env bash
#
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=gene_map
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus
#SBATCH --output=logs/gene_map.log
#SBATCH --error=logs/gene_map.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=1-23

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh

readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

readonly revel_cutoff="0.773"
readonly cadd_cutoff="28.1"
readonly spliceai_cutoff="0.50"
readonly annotation="spliceai=${spliceai_cutoff}_cadd=${cadd_cutoff}_revel=${revel_cutoff}"

readonly in_dir="data/vep/vep_loftee/csqs/${annotation}"
readonly in_path="${in_dir}/UKB.chr${chr}.phased_sites.${annotation}.canonical.txt.gz"

readonly out_dir="data/vep/saige_rescaled_weights"
readonly out_file="${out_dir}/UKB.chr${chr}.phased_sites.${annotation}.canonical.gene_map.txt"

mkdir -p ${out_dir}
zcat ${in_path} | cut -f1,2 | tail -n+2 > ${out_file}




