#!/usr/bin/env bash
# Worker script that runs REGENIE step2 inside swiss-army-knife job

set -euo pipefail

# Parse arguments
GENOTYPES="$1"          # PLINK bed file (e.g., "file.bed") OR BGEN file (e.g., "file.bgen")
PRED="$2"               # Pred file from step1
PHENOFILE="$3"          # Phenotype file
PHENOCOL="$4"           # Phenotype column name
KEEP="$5"               # Keep file (sample IDs)
COVARCOLLIST="$6"       # Comma-separated covariate list
CATEGCOVARCOLLIST="$7"  # Comma-separated categorical covariate list
OUT="$8"                # Output prefix

echo "Running REGENIE step2 with:"
echo "  GENOTYPES: ${GENOTYPES}"
echo "  PRED: ${PRED}"
echo "  PHENOFILE: ${PHENOFILE}"
echo "  PHENOCOL: ${PHENOCOL}"
echo "  KEEP: ${KEEP}"
echo "  Covariates: ${COVARCOLLIST}"
echo "  Categorical covariates: ${CATEGCOVARCOLLIST}"
echo "  Output: ${OUT}"

# Clean up pred file (remove output prefix paths)
PRED_LOCAL="pred_file_cleaned.txt"
cat ${PRED} | sed 's/\/home\/dnanexus\/out\/out\///g' > ${PRED_LOCAL}
PRED=${PRED_LOCAL}

echo "Cleaned pred file:"
head ${PRED}

# Decompress phenotype file if needed
PHENOFILE_LOCAL="phenotypes.txt"
if [[ "${PHENOFILE}" == *.gz ]]; then
  gunzip -c "${PHENOFILE}" > "${PHENOFILE_LOCAL}"
else
  cp "${PHENOFILE}" "${PHENOFILE_LOCAL}"
fi

echo "Phenotype file first 5 lines:"
head -5 "${PHENOFILE_LOCAL}"

# Reformat keep file to FID IID format (REGENIE expects both columns)
echo "Keep file first 5 lines (original):"
head -5 "${KEEP}"

KEEP_LOCAL="keep_file.txt"
awk '{ print $1,$1 }' "${KEEP}" > "${KEEP_LOCAL}"
KEEP="${KEEP_LOCAL}"

echo "Keep file first 5 lines (reformatted):"
head -5 "${KEEP}"

# Handle genotype format (PLINK or BGEN)
if [[ "${GENOTYPES}" == *.bed ]]; then
  # PLINK format
  BFILE=$(echo ${GENOTYPES} | sed 's/.bed$//g')

  # Rename FID column in PLINK fam file to match IID
  awk '{ print $2,$2,$3,$4,$5,$6 }' ${BFILE}.fam > ${BFILE}.fam-tmp
  mv ${BFILE}.fam-tmp ${BFILE}.fam

  echo "FAM file first 5 lines:"
  head -5 ${BFILE}.fam

  genotypes_flag="--bed ${BFILE}"

elif [[ "${GENOTYPES}" == *.bgen ]]; then
  # BGEN format
  BGEN_FILE="${GENOTYPES}"
  SAMPLE_FILE=$(echo ${BGEN_FILE} | sed 's/.bgen$/.sample/g')

  echo "BGEN file: ${BGEN_FILE}"
  echo "Sample file: ${SAMPLE_FILE}"
  echo "Sample file first 5 lines:"
  head -5 ${SAMPLE_FILE}

  # BGEN sample files from plink2 have this structure:
  # Line 1: ID_1 ID_2 missing sex (column names)
  # Line 2: 0 0 0 D (data types - NOT actual data)
  # Line 3+: Actual sample IDs
  # We need to use ID_2 (second column) for both FID and IID since ID_1 is often 0
  echo "Original sample file first 5 lines:"
  head -5 ${SAMPLE_FILE}

  SAMPLE_FILE_LOCAL="sample_file_fixed.sample"
  head -2 ${SAMPLE_FILE} > ${SAMPLE_FILE_LOCAL}  # Keep both header lines
  tail -n +3 ${SAMPLE_FILE} | awk '{print $2, $2, $3, $4}' >> ${SAMPLE_FILE_LOCAL}
  SAMPLE_FILE=${SAMPLE_FILE_LOCAL}

  echo "Fixed sample file first 10 lines:"
  head -10 ${SAMPLE_FILE}

  genotypes_flag="--bgen ${BGEN_FILE} --sample ${SAMPLE_FILE} --ref-first"

else
  echo "ERROR: Unknown genotype format: ${GENOTYPES}"
  exit 1
fi

# Run REGENIE step 2
regenie \
  --step 2 \
  ${genotypes_flag} \
  --phenoFile "${PHENOFILE_LOCAL}" \
  --phenoCol "${PHENOCOL}" \
  --covarFile "${PHENOFILE_LOCAL}" \
  --covarColList "${COVARCOLLIST}" \
  --catCovarList "${CATEGCOVARCOLLIST}" \
  --keep "${KEEP}" \
  --pred "${PRED}" \
  --minMAC 0.5 \
  --bsize 400 \
  --qt --apply-rint \
  --out "${OUT}"

# Rename and compress output
mv *regenie ${OUT}.regenie
gzip ${OUT}.regenie

echo "REGENIE step2 completed successfully!"
echo "Output files:"
ls -lh ${OUT}*
