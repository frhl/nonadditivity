#!/usr/bin/env bash
# Helper script to run all steps with proper waiting
# author: frederik lassen

set -euo pipefail

echo "=========================================="
echo "Power Simulation Pipeline"
echo "=========================================="
echo ""

# Step 1: Generate phenotypes
echo "=== Step 1: Generating Phenotypes ==="
if [[ ! -f ".pheno_job_eur.txt" ]]; then
    echo "Running 01_make_power_phenos.sh..."
    bash 01_make_power_phenos.sh
else
    echo "Phenotype job already submitted"
fi

if [[ -f ".pheno_job_eur.txt" ]]; then
    pheno_job=$(cat .pheno_job_eur.txt)
    echo "Waiting for phenotype job: ${pheno_job}"
    dx wait "${pheno_job}"
    echo "✓ Phenotype generation complete"
else
    echo "ERROR: No phenotype job found"
    exit 1
fi

echo ""

# Step 2: SAIGE step1
echo "=== Step 2: SAIGE Step1 (Null Model) ==="
echo "Running 02_power_step1.sh..."
bash 02_power_step1.sh

if [[ -f ".step1_jobs_eur.txt" ]]; then
    echo "Waiting for step1 jobs..."
    cat .step1_jobs_eur.txt | xargs dx wait
    echo "✓ Step1 complete"
else
    echo "WARNING: No step1 jobs found"
fi

echo ""

# Step 3: SAIGE step2
echo "=== Step 3: SAIGE Step2 (Association Tests) ==="
echo "Running 03_power_step2.sh..."
bash 03_power_step2.sh

echo ""
echo "=========================================="
echo "All jobs submitted!"
echo "=========================================="
echo ""
echo "To monitor all jobs:"
echo "  dx find jobs --created-after=-1h --brief"
echo ""
echo "To download results when complete:"
echo "  bash download_power_results.sh"
echo ""
