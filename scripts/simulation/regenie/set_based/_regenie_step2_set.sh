#!/usr/bin/env bash
# Worker script that runs REGENIE step2 set-based tests inside swiss-army-knife job

set -euo pipefail

# Parse arguments
GENOTYPES="$1"          # BGEN file (e.g., "file.bgen")
PRED="$2"               # Pred file from step1
PHENOFILE="$3"          # Phenotype file
PHENOCOL="$4"           # Phenotype column name
KEEP="$5"               # Keep file (sample IDs)
COVARCOLLIST="$6"       # Comma-separated covariate list
CATEGCOVARCOLLIST="$7"  # Comma-separated categorical covariate list
SETLIST="$8"            # Set list file (set definitions)
ANNO="$9"               # Annotation file
AAF="${10}"             # AAF file (allele frequencies)
MASKDEF="${11}"         # Mask definition file
OUT="${12}"             # Output prefix

echo "Running REGENIE step2 set-based test with:"
echo "  GENOTYPES: ${GENOTYPES}"
echo "  PRED: ${PRED}"
echo "  PHENOFILE: ${PHENOFILE}"
echo "  PHENOCOL: ${PHENOCOL}"
echo "  KEEP: ${KEEP}"
echo "  Covariates: ${COVARCOLLIST}"
echo "  Categorical covariates: ${CATEGCOVARCOLLIST}"
echo "  Set list: ${SETLIST}"
echo "  Annotation file: ${ANNO}"
echo "  AAF file: ${AAF}"
echo "  Mask definitions: ${MASKDEF}"
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

# Handle BGEN format
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

# Decompress set list file if needed
SETLIST_LOCAL="setlist.txt"
if [[ "${SETLIST}" == *.gz ]]; then
  gunzip -c "${SETLIST}" > "${SETLIST_LOCAL}"
else
  cp "${SETLIST}" "${SETLIST_LOCAL}"
fi

echo "Set list file first 10 lines:"
head -10 "${SETLIST_LOCAL}"

# Decompress annotation file if needed
ANNO_LOCAL="anno.txt"
if [[ "${ANNO}" == *.gz ]]; then
  gunzip -c "${ANNO}" > "${ANNO_LOCAL}"
else
  cp "${ANNO}" "${ANNO_LOCAL}"
fi

echo "Annotation file first 10 lines (before scaling):"
head -10 "${ANNO_LOCAL}"

# Scale weights in column 4 to [0, 1] to avoid issues with dominance encoding
echo "Scaling weights to [0, 1]..."
MAX_WEIGHT=$(awk '{print $4}' "${ANNO_LOCAL}" | sort -g | tail -1)
echo "Maximum weight found: ${MAX_WEIGHT}"

ANNO_SCALED="anno_scaled.txt"
awk -v max_w="${MAX_WEIGHT}" '{print $1, $2, $3, $4/max_w}' "${ANNO_LOCAL}" > "${ANNO_SCALED}"
ANNO_LOCAL="${ANNO_SCALED}"

echo "Annotation file first 10 lines (after scaling):"
head -10 "${ANNO_LOCAL}"

# Decompress AAF file if needed
AAF_LOCAL="aaf.txt"
if [[ "${AAF}" == *.gz ]]; then
  gunzip -c "${AAF}" > "${AAF_LOCAL}"
else
  cp "${AAF}" "${AAF_LOCAL}"
fi

echo "AAF file first 10 lines:"
head -10 "${AAF_LOCAL}"

# Show mask definitions
echo "Mask definitions:"
cat "${MASKDEF}"

# Run REGENIE step 2 with set-based testing
echo "Running REGENIE step2 set-based test..."
#--ref-first \

regenie \
  --step 2 \
  --bgen ${BGEN_FILE} \
  --sample ${SAMPLE_FILE} \
  --phenoFile "${PHENOFILE_LOCAL}" \
  --phenoCol "${PHENOCOL}" \
  --covarFile "${PHENOFILE_LOCAL}" \
  --covarColList "${COVARCOLLIST}" \
  --catCovarList "${CATEGCOVARCOLLIST}" \
  --aaf-file "${AAF_LOCAL}" \
  --keep "${KEEP}" \
  --pred "${PRED}" \
  --anno-file "${ANNO_LOCAL}" \
  --set-list "${SETLIST_LOCAL}" \
  --mask-def "${MASKDEF}" \
  --minMAC 0.5 \
  --weights-col 4 \
  --aaf-bins 1.0 \
  --bsize 400 \
  --apply-rint \
  --debug \
  --out "${OUT}"

# Clean up intermediate files
echo "Cleaning up intermediate files..."
rm -f pred_file_cleaned.txt phenotypes.txt keep_file.txt sample_file_fixed.sample
rm -f setlist.txt anno.txt anno_scaled.txt aaf.txt

# Compress output
echo "Compressing output..."
gzip ${OUT}_${PHENOCOL}.regenie

echo "REGENIE step2 set-based test completed successfully!"
echo "Final output files:"
ls -lh ${OUT}*
