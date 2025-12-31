# Gene Set-Based Power Simulation Framework

Simulate quantitative phenotypes using **real variants from gene sets** (not pseudo-variants) to evaluate power of gene-based burden tests under diverse genetic architectures.

## Overview

This framework differs from the single-variant power simulation in several key ways:

| Feature | Single-Variant (`power/`) | Gene Set-Based (this framework) |
|---------|---------------------------|----------------------------------|
| **Variants** | Pseudo-variants from PLINK | Real variants from gene sets |
| **Unit** | Individual variants | Gene sets |
| **Causal selection** | M random variants | K genes, N variants per gene |
| **Encoding** | Additive PLINK | Dominance-encoded VCF |
| **Input** | Merged PLINK files | SAIGE weight files + VCF |

---

## Quick Start

```bash
# Step 0: Simulate phenotypes (chr21 for debugging)
cd /Users/flassen/Projects/11_wes_ko_ukbb_nexus/wes_ko_ukbb_nexus/scripts/simulation/saige/gene_set_power
./01_make_gene_set_phenos.sh

# Step 1: Fit null models
./02_gene_set_step1.sh

# Step 2: Run burden tests
./03_gene_set_step2.sh
```

---

## Input Data

### 1. SAIGE Weight Files
**Location:** `/wes_ko_ukbb/data/genesets/dominance_weights/af05/`

**Format:**
```
ENSG00000155307 var chr21:14482965:C:T chr21:14487001:G:A ...
ENSG00000155307 anno non_coding pLoF damaging_missense ...
ENSG00000155307 weight 13.1058 12.0873 ...
```

**Contains:**
- Gene IDs (Ensembl)
- Variant IDs (chr:pos:ref:alt)
- Annotations (pLoF, damaging_missense, synonymous, etc.)
- Beta weights for burden tests

**Filtering already applied:**
- MAF < 0.05 (af05)
- MAC ≥ 10 (min_mac_10)
- Spliceai, CADD, REVEL thresholds
- Popmax exclusions

### 2. Dominance-Encoded VCF Files
**Location:** `/wes_ko_ukbb/data/phased/encoded_dominance_gt/group_dominance/gene_map/`

**Contains:**
- Dominance-encoded genotypes (X^D)
- Gene mappings
- ≥5 homozygous carriers per variant (enforced during encoding)

### 3. Allele Count Files
**Location:** `/wes_ko_ukbb/data/phased/wes_union_calls/450k/frqx/`

**Used for:**
- MAF calculation
- Carrier count verification

---

## Simulation Parameters

### Core Parameters

```bash
# Genetic architectures
--architectures "additive,recessive,partially_recessive_0.1,partially_recessive_0.2"

# Heritabilities (h² per variant, not per gene)
--h2_total "0.001,0.005,0.01,0.02"

# K: Number of causal genes
--K_genes "5,10,20"

# N: Number of causal variants per gene (default: 1)
--N_per_gene "1,3"

# Variant selection strategy
--variant_selection "random"  # Options: random, lowest_maf, highest_maf, weighted

# Annotations to include
--annotations "pLoF,damaging_missense"

# MAF bins
--maf_bins "0.001-0.05"  # Or "all" for no filtering

# Number of replicates
--n_reps 100
```

### Effect Size Calculation

**Equal h² per variant** (not per gene):

```
Total causal variants = K_genes × N_per_gene
h² per variant = h²_total / Total causal variants
```

**Example:**
- K = 10 genes
- N = 3 variants/gene
- Total variants = 30
- h²_total = 0.01
- h²_per_variant = 0.01 / 30 = 0.000333

---

## Output Files

### Phenotype File
`{DATE}_chr{chr}_{ancestry}_gene_set_power.tsv.gz`

**Format:**
```
eid    arch_recessive_h2_0.005_K_10_N_1_rep_1    arch_recessive_h2_0.005_K_10_N_3_rep_1    ...
1001   -0.234                                     0.512                                      ...
1002   1.123                                      -0.834                                     ...
```

**Phenotype naming:**
```
arch_{architecture}_h2_{h2}_K_{K_genes}_N_{N_per_gene}_rep_{replicate}
```

### Metadata File
`{DATE}_chr{chr}_{ancestry}_gene_set_power.metadata.tsv.gz`

**Columns:**
```
phenotype_name, replicate, seed,
architecture, h2_total, K_genes, N_per_gene, total_causal_variants, h2_per_variant,
gene_id, variant_id, MAF, annotation, weight,
beta_A, beta_D, h2_A, h2_D, r, h, a
```

**Contains:**
- All causal genes and variants
- Effect sizes (additive and dominance)
- Genotype proportions
- Selection parameters

---

## Workflow Details

### Script 1: `01_make_gene_set_phenos.R`

**Main simulation script**

1. Parse SAIGE weight files (gene-variant mappings)
2. Load allele counts and calculate MAF
3. Filter genes by:
   - Annotations (pLoF + damaging_missense)
   - MAF range (e.g., 0.001-0.05)
   - Minimum variants per gene
4. Select K causal genes (random)
5. Select N causal variants per gene (strategy: random/lowest_maf/etc.)
6. Load dominance-encoded genotypes from VCF
7. For each variant:
   - Compute genotype proportions (r, h, a)
   - Compute effect sizes (β_A, β_D) based on architecture
8. Simulate phenotype:
   - Genetic value = Σ(β_A × X^A + β_D × X^D)
   - Add environmental noise
   - Standardize to mean=0, sd=1
9. Output phenotypes + metadata

**Key Functions (_gene_set_utils.R):**
- `parse_saige_weights()`: Parse gene sets
- `load_allele_counts()`: Get MAF/MAC
- `filter_genes()`: Apply filters
- `select_causal_genes()`: Random gene selection
- `select_causal_variants_per_gene()`: Variant selection strategies
- `load_raw_genotypes_from_vcf()`: Extract genotypes

---

### Script 2: `02_gene_set_step1.sh`

**SAIGE Step1: Fit null models**

- Reads phenotype file
- Extracts all phenotype columns (arch_*)
- Submits `saige-universal-step-1` for each phenotype
- Quantitative trait with inverse normal transformation

**Output:**
- `.rda` model files
- `.varianceRatio.txt` files

---

### Script 3: `03_gene_set_step2.sh`

**SAIGE Step2: Gene-based burden tests**

- Uses dominance-encoded VCF
- Uses SAIGE weight files (with gene mappings)
- Tests all genes (not just causal ones) for proper calibration
- Group test with annotations: `pLoF:damaging_missense`

**Output:**
- `.txt.gz` results with gene-level p-values

---

## Comparison with Single-Variant Power

This framework enables direct comparison:

| Gene-Based | Single-Variant Equivalent |
|------------|---------------------------|
| K=10 genes, N=1 var/gene | M=10 causal variants |
| K=10 genes, N=3 var/gene | M=30 causal variants |
| K=50 genes, N=1 var/gene | M=50 causal variants |

**Same architectures and h² values** enable cross-framework comparison.

---

## Variant Selection Strategies

### `random` (default)
Uniform random selection within each gene.

### `lowest_maf`
Select the N rarest variants per gene (smallest MAF).

### `highest_maf`
Select the N most common variants per gene (largest MAF).

### `weighted`
Probabilistic selection using beta weights (dbeta(MAF, 1, 25)).

---

## Gene Filtering Criteria

Applied in order:

1. **Annotation filter:** pLoF + damaging_missense only
2. **MAF range:** 0.001 ≤ MAF < 0.05 (configurable)
3. **Homozygous carriers:** ≥5 (enforced in VCF generation)
4. **Minimum variants per gene:** ≥1 (configurable)

---

## Debugging with Chr21

All scripts default to **chr21** for fast debugging:

```bash
for chr in 21; do
    # Only chr21 processed
done
```

**Benefits:**
- Faster execution
- Smaller datasets
- Quick iteration

**To expand:** Change loop to `{1..22}` or `{1..22} X`

---

## File Locations (DNAnexus)

### Inputs
```
/wes_ko_ukbb/data/genesets/dominance_weights/af05/
/wes_ko_ukbb/data/phased/encoded_dominance_gt/group_dominance/gene_map/
/wes_ko_ukbb/data/phased/wes_union_calls/450k/frqx/
/wes_ko_ukbb/data/phenotypes/covariates/
```

### Outputs
```
/wes_ko_ukbb/data/phenotypes/simulated/gene_set_power/{DATE}/
/wes_ko_ukbb/data/saige/step1/simulated/gene_set_power/{DATE}/
/wes_ko_ukbb/data/saige/step2/simulated/gene_set_power/{DATE}/
```

---

## Troubleshooting

### Error: "No variants pass carrier filter"
- Check that VCF has ≥5 homozygous carriers
- Verify MAF range isn't too restrictive

### Error: "SAIGE weight file format mismatch"
- Ensure file has 3 rows per gene (var, anno, weight)
- Check gzip compression

### Error: "VCF index not found"
- Ensure `.csi` index file exists alongside VCF

### Warning: "Gene X has only Y variants, using all"
- Gene has fewer variants than N_per_gene
- All available variants will be used

---

## Key Design Features

✅ **Real gene sets:** Uses actual variant collections, not pseudo-variants
✅ **VEP annotations:** pLoF and damaging_missense variants
✅ **Flexible polygenicity:** Control K genes and N variants/gene independently
✅ **Equal h² per variant:** Fair comparison across configurations
✅ **MAF control:** Filter like single-variant power simulation
✅ **Chr21 debugging:** Fast iteration for development
✅ **Rich metadata:** Track all genes, variants, and effect sizes
✅ **Reusable code:** Leverages existing effect size utilities

---

## Dependencies

- R packages: `data.table`, `dplyr`, `optparse`, `VariantAnnotation`
- DNAnexus applets: `saige-universal-step-1`, `saige-universal-step-2-group`
- Docker: `rsuite.tar.gz`

---

## Example Command (Local Testing)

```bash
Rscript 01_make_gene_set_phenos.R \
  --saige_weights chr21_weights.txt.gz \
  --vcf_file chr21_dominance.vcf.gz \
  --frqx_file chr21_eur.frqx.gz \
  --ancestry eur \
  --annotations "pLoF,damaging_missense" \
  --architectures "additive,recessive" \
  --h2_total "0.005,0.01" \
  --K_genes "5,10" \
  --N_per_gene "1,3" \
  --variant_selection "random" \
  --maf_bins "0.001-0.05" \
  --n_reps 2 \
  --seed 42 \
  --out test_output
```

---

## Authors

Frederik Lassen

## Date

2024-12-28
