#!/usr/bin/env bash
# author: frederik lassen
# note: this script downloads SAIGE power simulation results from DNAnexus

set -euo pipefail

# =============================================================================
# CONFIGURATION - Must match 01_make_power_phenos.sh
# =============================================================================
readonly DATE="251227"  # Must match the date in 01_make_power_phenos.sh

# DNAnexus path
readonly dx_power_dir="/wes_ko_ukbb/data/saige/step2/simulated/power/${DATE}"

# Local output path (simplified to 2 levels)
readonly local_base="/Users/flassen/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/data/${DATE}"
readonly local_power_dir="${local_base}/power_simulation"

# Create local directory if it doesn't exist
mkdir -p "${local_power_dir}"

echo "====================================="
echo "Downloading Power Simulation Results"
echo "====================================="
echo ""
echo "From: ${dx_power_dir}"
echo "To: ${local_power_dir}"
echo ""

# List .gz files in the directory
power_files=$(dx ls "${dx_power_dir}" | grep '\.gz$' || true)

if [[ -z "${power_files}" ]]; then
  echo "No .gz files found in ${dx_power_dir}"
else
  power_count=0
  skipped_count=0
  while IFS= read -r file; do
    local_file="${local_power_dir}/${file}"

    # Check if file already exists locally
    if [[ -f "${local_file}" ]]; then
      echo "Skipping (already exists): ${file}"
      ((skipped_count++))
    else
      echo "Downloading: ${file}"
      dx download "${dx_power_dir}/${file}" -o "${local_file}" --no-progress || true
      ((power_count++))
    fi
  done <<< "${power_files}"
  echo ""
  echo "Downloaded ${power_count} new file(s)"
  echo "Skipped ${skipped_count} existing file(s)"
fi

echo ""
echo "====================================="
echo "Download complete!"
echo "====================================="
echo ""
echo "Results saved to: ${local_power_dir}"
echo ""
