#!/usr/bin/env bash
# author: frederik lassen
# note: Wrapper script to run local REGENIE Type I error assessment

set -euo pipefail

# Default parameters
dx_dir="/wes_ko_ukbb/data/regenie/step2/simulated/revisions/051025/combined"
pattern="combined.regenie.tsv.gz"
output="regenie_type1_error_summary.tsv"

# Allow overriding via command line args
if [[ $# -ge 1 ]]; then dx_dir="$1"; fi
if [[ $# -ge 2 ]]; then pattern="$2"; fi
if [[ $# -ge 3 ]]; then output="$3"; fi

echo "Running REGENIE Type I error assessment..."
echo "  DX Directory: ${dx_dir}"
echo "  Pattern: ${pattern}"
echo "  Output: ${output}"

Rscript 06_assess_type1_error_local.R \
    --dx_dir "${dx_dir}" \
    --pattern "${pattern}" \
    --output "${output}"

echo "Done!"
