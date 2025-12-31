#!/usr/bin/env bash
# Worker script that converts VCF to BGEN format

# Parse arguments
VCF_FILE="$1"       # VCF file (e.g., "file.vcf.gz")
OUT_PREFIX="$2"     # Output prefix

echo "Converting VCF to BGEN format:"
echo "  VCF: ${VCF_FILE}"
echo "  Output prefix: ${OUT_PREFIX}"

# Check if plink2 is available
if ! command -v plink2 &> /dev/null; then
  echo "ERROR: plink2 not found. Checking for plink..."
  if command -v plink &> /dev/null; then
    plink --version
    echo "ERROR: Only plink1.9 is available. plink2 is required for dosage data."
    exit 1
  else
    echo "ERROR: Neither plink nor plink2 found."
    exit 1
  fi
fi

echo "plink2 version:"
plink2 --version

# Convert VCF to BGEN format
# --vcf dosage=DS: Import dosage data from VCF DS field (for float dosages)
# --import-dosage-certainty 1: Keep dosage probabilities as-is
# --hard-call-threshold 0: Don't convert to hard calls
# --export bgen-1.3 'bits=16': Export to BGEN v1.3 with 16-bit precision
# --out: Output prefix
echo "Running plink2 conversion to BGEN..."
plink2 \
  --vcf ${VCF_FILE} dosage=DS \
  --import-dosage-certainty 1 \
  --hard-call-threshold 0 \
  --export bgen-1.3 'bits=16' \
  --out ${OUT_PREFIX}

echo "Conversion completed successfully!"

# Extract variant information from VCF (first 8 columns)
echo ""
echo "Extracting variant information from VCF..."
zcat ${VCF_FILE} | cut -f1-8 | grep -v "##" | gzip > ${OUT_PREFIX}.variants.txt.gz
echo "Variant info file created: ${OUT_PREFIX}.variants.txt.gz"
echo "First 10 lines:"
zcat ${OUT_PREFIX}.variants.txt.gz | head

echo ""
echo "Output files:"
ls -lh ${OUT_PREFIX}*
