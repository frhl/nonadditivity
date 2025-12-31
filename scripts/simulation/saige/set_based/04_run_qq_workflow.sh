#!/usr/bin/env bash
# author: frederik lassen
# note: this script creates QQ plots from downloaded SAIGE results
# Run 03_download_results.sh first to download the data

set -euo pipefail

# Get the directory where this script is located
script_dir="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"

echo "============================================"
echo "Creating SAIGE QQ Plots"
echo "============================================"
echo ""

Rscript "${script_dir}/04_run_qq_workflow.R"

echo ""
echo "============================================"
echo "Workflow complete!"
echo "============================================"
echo ""
echo "Results are saved in: /Users/flassen/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/derived/saige_simulation_results/plots"
echo ""
echo "Generated plots:"
echo "  - qq_plot_variant_based.png (all variant tests combined - all modes)"
echo "  - qq_plot_gene_based.png (all gene tests combined)"
echo "  - qq_plot_variant_h2_*.png (stratified by heritability)"
echo "  - qq_plot_gene_h2_*.png (stratified by heritability)"
echo "  - qq_plot_combined.png (comparison plot)"
