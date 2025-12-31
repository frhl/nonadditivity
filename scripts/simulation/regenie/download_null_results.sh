#!/usr/bin/env bash
# author: frederik lassen
# note: this script downloads REGENIE null simulation combined results from DNAnexus

set -euo pipefail

# DNAnexus path
readonly dx_null_dir="/wes_ko_ukbb/data/regenie/step2/simulated/revisions/051025/combined"

# Local output path (simplified to 2 levels)
readonly local_base="/Users/flassen/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/data"
readonly local_null_dir="${local_base}/null_simulation_regenie"

# Create local directory if it doesn't exist
mkdir -p "${local_null_dir}"

echo "====================================="
echo "Downloading REGENIE Null Simulation Results"
echo "====================================="
echo ""
echo "From: ${dx_null_dir}"
echo "To: ${local_null_dir}"
echo ""

# List .gz files in the directory
null_files=$(dx ls "${dx_null_dir}" | grep '\.gz$' || true)

if [[ -z "${null_files}" ]]; then
  echo "No .gz files found in ${dx_null_dir}"
else
  null_count=0
  while IFS= read -r file; do
    echo "Downloading: ${file}"
    dx download "${dx_null_dir}/${file}" -o "${local_null_dir}/${file}" --no-progress || true
    ((null_count++))
  done <<< "${null_files}"
  echo ""
  echo "Downloaded ${null_count} file(s)"
fi

echo ""
echo "====================================="
echo "Download complete!"
echo "====================================="
echo ""
echo "Results saved to: ${local_null_dir}"
echo ""
