#!/usr/bin/env bash
# Worker script that runs REGENIE step1 inside swiss-army-knife job

set -euo pipefail

# Parse arguments
BFILE="$1"              # PLINK bfile prefix (e.g., "ukb22418_b0_v2.autosomes")
PHENOFILE="$2"          # Phenotype file name
PHENOCOL="$3"           # Phenotype column name (e.g., "p_0.01_continuous_1")
KEEP="$4"               # Keep file (sample IDs to include)
EXTRACT="$5"            # Extract file (variant IDs to include)
COVARCOLLIST="$6"       # Comma-separated covariate list
CATEGCOVARCOLLIST="$7"  # Comma-separated categorical covariate list
N_THREADS="$8"          # Number of threads
BSIZE="$9"              # Block size for REGENIE
OUT="${10}"             # Output prefix

echo "Running REGENIE step1 with:"
echo "  BFILE: ${BFILE}"
echo "  PHENOFILE: ${PHENOFILE}"
echo "  PHENOCOL: ${PHENOCOL}"
echo "  KEEP: ${KEEP}"
echo "  EXTRACT: ${EXTRACT}"
echo "  Covariates: ${COVARCOLLIST}"
echo "  Categorical covariates: ${CATEGCOVARCOLLIST}"
echo "  Threads: ${N_THREADS}"
echo "  Block size: ${BSIZE}"
echo "  Output: ${OUT}"

# Decompress phenotype file if needed
PHENOFILE_LOCAL="phenotypes.txt"
if [[ "${PHENOFILE}" == *.gz ]]; then
  gunzip -c "${PHENOFILE}" > "${PHENOFILE_LOCAL}"
else
  cp "${PHENOFILE}" "${PHENOFILE_LOCAL}"
fi

echo "Phenotype file first 5 lines:"
head -5 "${PHENOFILE_LOCAL}"

echo "Checking header format:"
head -1 "${PHENOFILE_LOCAL}" | od -c | head -20

echo "Keep file first 5 lines:"
head -5 "${KEEP}"

# REGENIE expects keep file to have FID and IID columns
# Create a proper keep file with FID IID format
KEEP_LOCAL="keep_file.txt"
awk '{ print $1,$1 }' "${KEEP}" > "${KEEP_LOCAL}"
KEEP="${KEEP_LOCAL}"

echo "Reformatted keep file first 5 lines:"
head -5 "${KEEP}"

echo "Extract file first 5 lines:"
head -5 "${EXTRACT}"

# Run REGENIE step 1
# Note: covariates are in the phenotype file, so we use --covarFile with the same file
regenie \
  --step 1 \
  --bed "${BFILE}" \
  --phenoFile "${PHENOFILE_LOCAL}" \
  --phenoCol "${PHENOCOL}" \
  --covarFile "${PHENOFILE_LOCAL}" \
  --covarColList "${COVARCOLLIST}" \
  --catCovarList "${CATEGCOVARCOLLIST}" \
  --keep "${KEEP}" \
  --extract "${EXTRACT}" \
  --threads ${N_THREADS} \
  --bsize ${BSIZE} \
  --qt --apply-rint \
  --out "${OUT}"

echo "REGENIE step1 completed successfully!"
echo "Output files:"
ls -lh ${OUT}*
