#!/usr/bin/env bash
#
#SBATCH --account=lindgren.prj
#SBATCH --job-name=write_intervals
#SBATCH --chdir=/well/lindgren-ukbb/projects/ukbb-11867/flassen/projects/KO/wes_ko_ukbb_nexus
#SBATCH --output=logs/write_intervals.log
#SBATCH --error=logs/write_intervals.errors.log
#SBATCH --partition=short
#SBATCH --cpus-per-task 1
#SBATCH --array=17

set -o errexit
set -o nounset

source utils/bash_utils.sh
source utils/qsub_utils.sh
source utils/hail_utils.sh

readonly hail_script="scripts/phasing_wes/phasing/02_write_intervals.py"
readonly rscript_trimming="scripts/phasing_wes/phasing/02_write_trimming_intervals.R"
readonly rscript_scaffolds="scripts/phasing_wes/phasing/02_write_scaffolds.R"
readonly spark_dir="data/tmp/spark"

mkdir -p ${spark_dir}

set -eu

readonly cluster=$( get_current_cluster )
readonly array_idx=$( get_array_task_id )
readonly chr=$( get_chr ${array_idx} )

# Number of variants within each interval
readonly min_interval_unit=1000
# Default size of phasing window in terms of variant count (should be a multiple of min_interval_unit)
#readonly phasing_region_size=300000
readonly phasing_region_size=500000
# Minimum overlap between adjacent phasing windows
readonly phasing_region_overlap=25000 #$(( ${phasing_region_size}/8 ))
# Maximum size of phasing window allowed, only used at the end of a chromosome
# Must be larger than phasing_region_size
#readonly max_phasing_region_size=350000
readonly max_phasing_region_size=1000000

# intervals for phasing and for later edge trimming
readonly interval_dir="data/phased/wes_union_calls/450k/shapeit5/intervals_full_qc_new"
readonly interval_path="${interval_dir}/intervals_min_${min_interval_unit}_chr${chr}.tsv"
readonly region_path="${interval_dir}/phasing_intervals_u${min_interval_unit}_s${phasing_region_size}_o${phasing_region_overlap}_chr${chr}.tsv"
readonly region_trimming_path="${interval_dir}/trimming_intervals_u${min_interval_unit}_s${phasing_region_size}_o${phasing_region_overlap}_chr${chr}.tsv"
readonly region_scaffold_path="${interval_dir}/scaffold_intervals_u${min_interval_unit}_s${phasing_region_size}_o${phasing_region_overlap}_chr${chr}.tsv"
readonly phasing_interval_flags="--chrom ${chr} --min_interval_unit ${min_interval_unit}"
mkdir -p ${interval_dir}

# VCF with variants only
readonly vcf_dir="data/unphased/variants/full_qc_new"
readonly vcf_to_phase="${vcf_dir}/UKB.chr${chr}.exome_array.full_qc.no_parents.variants_only.vcf.gz"

if [ -z "${interval_path}" ]; then
  raise_error "Getting intervals path failed"
fi

echo "starting hail.."
module purge
set_up_pythonpath_legacy
set_up_hail
#
# Write phasing (minumum unit) intervals to slice later for phasing intervals
if [ ! -f ${interval_path} ]; then
  mkdir -p $( dirname ${interval_path} )
  SECONDS=0
  set -x
  python3 ${hail_script} \
    ${phasing_interval_flags} \
    --write_intervals \
    --interval_path ${interval_path} \
    --target_vcf ${vcf_to_phase} \
    && print_update "Finished writing intervals for chr${chr}" ${SECONDS} \
    || raise_error "Writing intervals for chr${chr} failed"
  set +x
else
  print_update "${interval_path} already exists!"
fi


# get max phasing index
readonly min_phasing_idx=1
readonly max_phasing_idx=$( python3 ${hail_script} ${phasing_interval_flags} \
    --phasing_region_size ${phasing_region_size} \
    --phasing_region_overlap ${phasing_region_overlap} \
    --max_phasing_region_size ${max_phasing_region_size} \
    --get_max_phasing_idx \
    --interval_path ${interval_path} )

# setup flags
readonly interval_flags="--chrom ${chr}
  --min_interval_unit ${min_interval_unit}
  --phasing_region_size ${phasing_region_size}
  --phasing_region_overlap ${phasing_region_overlap}
  --max_phasing_region_size ${max_phasing_region_size}"

# create phasing intervals to be submitted
echo ${interval_flags} && rm -f ${region_path}
for phasing_idx in $(seq ${min_phasing_idx}  ${max_phasing_idx}); do
  cur_region=$( python3 ${hail_script} ${interval_flags} --get_interval --phasing_idx ${phasing_idx} --interval_path ${interval_path} )
  echo "chr${cur_region}" >> ${region_path}
done

# create trimming intervals to be submitted
module purge
set_up_rpy
Rscript ${rscript_trimming} \
  --interval_path "${interval_path}" \
  --current_intervals_path "${region_path}" \
  --out_path "${region_trimming_path}" \
  --overlap_desired 5 # 5000 variants overlap

#create scaffold intervals to be submitted
Rscript ${rscript_scaffolds} \
  --interval_path "${interval_path}" \
  --current_intervals_path "${region_path}" \
  --out_path "${region_scaffold_path}" \
  --padding 50 # 50,000 BP padding




