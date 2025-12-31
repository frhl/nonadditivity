#!/usr/bin/env bash
# author: frederik lassen
# note: this script downloads SAIGE results for both variant-based and gene-based tests

set -euo pipefail

# Define output directory
readonly local_dir="/Users/flassen/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/derived/saige_simulation_results"
mkdir -p ${local_dir}

# Define DNAnexus paths
readonly step2_variant_dir="/wes_ko_ukbb/data/saige/step2/simulated/revisions/051025"
readonly step2_gene_dir="/wes_ko_ukbb/data/saige/step2_set/simulated/revisions/051025"

# Define parameters
readonly ancestries=("eur")  # Only EUR results available
readonly heritabilities=("0.01" "0.1" "0.2" "0.3" "0.5")
readonly modes=("additive" "recessive" "dominance")
readonly n_reps=4

echo "============================================"
echo "Downloading variant-based SAIGE results"
echo "============================================"

# Download variant-based results (single-variant association tests)
for anc in "${ancestries[@]}"; do
  echo "Processing ancestry: ${anc}"

  for mode in "${modes[@]}"; do
    echo "  Processing mode: ${mode}"

    for h2 in "${heritabilities[@]}"; do
      for rep in $(seq 1 ${n_reps}); do
        pheno="p_${h2}_continuous_${rep}"

        # Example file pattern from step2 script:
        # UKB.auto.${pheno}.${mode}.${anc}.af${af}.pp${pp}.${anno}.txt.gz
        pattern="UKB.auto.${pheno}.${mode}.${anc}.af05.pp0.90.pLoF_damaging_missense"

        echo "    Downloading: ${pattern}.txt.gz"
        dx download "${step2_variant_dir}/${pattern}.txt.gz" \
          -o "${local_dir}/" --no-progress 2>/dev/null || \
          echo "      Warning: Could not download ${pattern}.txt.gz"
      done
    done
  done
done

echo ""
echo "============================================"
echo "Downloading gene-based SAIGE results"
echo "============================================"

# Download gene-based results (set-based tests)
for anc in "${ancestries[@]}"; do
  echo "Processing ancestry: ${anc}"

  for chr in 20; do
    echo "  Processing chromosome: ${chr}"

    for h2 in "${heritabilities[@]}"; do
      for rep in $(seq 1 ${n_reps}); do
        pheno="p_${h2}_continuous_${rep}"

        # Example file pattern from step2_set_based script:
        # UKB.auto.chr${chr}.${pheno}.${anc}.txt.gz
        pattern="UKB.auto.chr${chr}.${pheno}.${anc}"

        echo "    Downloading: ${pattern}.txt.gz"
        dx download "${step2_gene_dir}/${pattern}.txt.gz" \
          -o "${local_dir}/" --no-progress 2>/dev/null || \
          echo "      Warning: Could not download ${pattern}.txt.gz"

        # Also download singleAssoc file (variant-level results from gene-based tests)
        echo "    Downloading: ${pattern}.txt.singleAssoc.txt.gz"
        dx download "${step2_gene_dir}/${pattern}.txt.singleAssoc.txt.gz" \
          -o "${local_dir}/" --no-progress 2>/dev/null || \
          echo "      Warning: Could not download ${pattern}.txt.singleAssoc.txt.gz"
      done
    done
  done
done

echo ""
echo "============================================"
echo "Download complete!"
echo "============================================"
echo "Results saved to: ${local_dir}"
echo ""
echo "File counts:"
echo "  Variant-based results (by mode):"
echo "    - Additive: $(ls ${local_dir}/*additive*.txt.gz 2>/dev/null | wc -l | tr -d ' ')"
echo "    - Recessive: $(ls ${local_dir}/*recessive*.txt.gz 2>/dev/null | wc -l | tr -d ' ')"
echo "    - Dominance: $(ls ${local_dir}/*dominance*.txt.gz 2>/dev/null | wc -l | tr -d ' ')"
echo "  Gene-based results:"
echo "    - Gene-level (burden): $(ls ${local_dir}/UKB.auto.chr*.txt.gz 2>/dev/null | grep -v singleAssoc | wc -l | tr -d ' ')"
echo "    - Variant-level (singleAssoc): $(ls ${local_dir}/*singleAssoc*.txt.gz 2>/dev/null | wc -l | tr -d ' ')"
