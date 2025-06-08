#!/bin/bash

# Conditional genetic association analysis using SAIGE
#
# This script performs iterative conditional analysis to identify independent genetic
# signals in a region around a transcript. It conditions on the strongest associated
# variants (not the transcript) in each iteration until either:
# 1. The transcript's p-value becomes non-significant (>threshold)
# 2. No more significant variants are found to condition on
# 3. Maximum iterations are reached

# all functions being sources from here
source /mnt/project/wes_ko_ukbb/scripts/_conditional_utils.sh

# Exit on error or undefined variables
set -e
set -u

# Help function
usage() {
    echo "Usage: $0
    --pheno <phenotype>
    --model_file <SAIGE model file>
    --var_ratio_file <variance ratio file>
    --grm <sparse GRM matrix file>
    --grm_samples <GRM sample IDs file>
    --plink_file <plink file prefix>
    --out_prefix <output prefix>
    --start_pos <start position>
    --end_pos <end position>
    --transcript_id <ensembl transcript ID>
    --chromosome <chromosome number>
    --pseudovariant_to_variant_file <mapping file> [optional]
    --max_iterations <maximum number of iterations> [default: 10]
    --threshold <p-value threshold> [default: 0.0001]
    --padding <padding size in bp> [default: 500000]"
    exit 1
}

# Initialize variables
PHENO=""
MODEL_FILE=""
VAR_RATIO_FILE=""
GRM=""
GRM_SAMPLES=""
PLINK_FILE=""
OUT_PREFIX=""
START_POS=""
END_POS=""
TRANSCRIPT_ID=""
CHROMOSOME=""
PSEUDOVARIANT_MAP_FILE=""
MAX_ITERATION=10  # Default value
THRESHOLD=0.0001  # Default value
PADDING=500000    # Default value

# Parse command line arguments
while [ $# -gt 0 ]; do
    case "$1" in
        --pheno)
            PHENO="$2"
            shift 2
            ;;
        --model_file)
            MODEL_FILE="$2"
            shift 2
            ;;
        --var_ratio_file)
            VAR_RATIO_FILE="$2"
            shift 2
            ;;
        --grm)
            GRM="$2"
            shift 2
            ;;
        --grm_samples)
            GRM_SAMPLES="$2"
            shift 2
            ;;
        --plink_file)
            PLINK_FILE="$2"
            shift 2
            ;;
        --out_prefix)
            OUT_PREFIX="$2"
            shift 2
            ;;
        --start_pos)
            START_POS="$2"
            shift 2
            ;;
        --end_pos)
            END_POS="$2"
            shift 2
            ;;
        --transcript_id)
            TRANSCRIPT_ID="$2"
            shift 2
            ;;
        --chromosome)
            CHROMOSOME="$2"
            shift 2
            ;;
        --pseudovariant_to_variant_file)
            PSEUDOVARIANT_MAP_FILE="$2"
            shift 2
            ;;
        --max_iterations)
            MAX_ITERATION="$2"
            shift 2
            ;;
        --threshold)
            THRESHOLD="$2"
            shift 2
            ;;
        --padding)
            PADDING="$2"
            shift 2
            ;;
        -h|--help)
            usage
            ;;
        *)
            echo "Unknown parameter: $1"
            usage
            ;;
    esac
done

# Check if required arguments are provided
if [[ -z "$PHENO" ]] || [[ -z "$MODEL_FILE" ]] || [[ -z "$VAR_RATIO_FILE" ]] || \
   [[ -z "$GRM" ]] || [[ -z "$GRM_SAMPLES" ]] || [[ -z "$PLINK_FILE" ]] || \
   [[ -z "$OUT_PREFIX" ]] || [[ -z "$START_POS" ]] || [[ -z "$END_POS" ]] || \
   [[ -z "$TRANSCRIPT_ID" ]] || [[ -z "$CHROMOSOME" ]]; then
    echo "Error: Missing required arguments"
    usage
fi

# Validate pseudovariant mapping file if provided
if [[ -n "$PSEUDOVARIANT_MAP_FILE" ]]; then
    if [[ ! -f "$PSEUDOVARIANT_MAP_FILE" ]]; then
        echo "Error: Pseudovariant mapping file does not exist: $PSEUDOVARIANT_MAP_FILE"
        exit 1
    fi
    echo "Using pseudovariant mapping file: $PSEUDOVARIANT_MAP_FILE"
fi

# Validate numeric positions and parameters
if ! [[ "$START_POS" =~ ^[0-9]+$ ]] || ! [[ "$END_POS" =~ ^[0-9]+$ ]]; then
    echo "Error: Start and end positions must be numeric"
    exit 1
fi

if ! [[ "$MAX_ITERATION" =~ ^[0-9]+$ ]]; then
    echo "Error: Maximum iterations must be a positive integer"
    exit 1
fi

if ! [[ "$THRESHOLD" =~ ^[0-9.]+$ ]]; then
    echo "Error: Threshold must be a positive number"
    exit 1
fi

if ! [[ "$PADDING" =~ ^[0-9]+$ ]]; then
    echo "Error: Padding must be a positive integer"
    exit 1
fi

# Validate position order
if [ "$START_POS" -gt "$END_POS" ]; then
    echo "Error: Start position must be less than or equal to end position"
    exit 1
fi

echo "Starting analysis at $(date)"
echo "Creating data directory..."

# Create data directory and change to it
cd ~
mkdir -p data_dir/
cd data_dir/

# Download required files
echo "Downloading input files..."

cp ${MODEL_FILE} .
cp ${VAR_RATIO_FILE} .
cp ${GRM} .
cp ${GRM_SAMPLES} .

# Download PLINK files
echo "Downloading PLINK files..."
echo "prefix: ${PLINK_FILE}"
cp ${PLINK_FILE}.bed .
cp ${PLINK_FILE}.bim .
cp ${PLINK_FILE}.fam .

# Download pseudo-variant files and exclude variants that overlap with plink file




# Return to home directory
cd ~

# Set up local paths for all files
LOCAL_MODEL_FILE="/data_dir/$(basename ${MODEL_FILE})"
LOCAL_VAR_RATIO_FILE="/data_dir/$(basename ${VAR_RATIO_FILE})"
LOCAL_SPARSE_GRM="/data_dir/$(basename ${GRM})"
LOCAL_SPARSE_GRM_SAMPLES="/data_dir/$(basename ${GRM_SAMPLES})"
LOCAL_PLINK_PREFIX="/data_dir/$(basename ${PLINK_FILE})"

# install saige
cd ~
git clone -b dominance https://github.com/frhl/universal-saige.git
mv universal-saige/* .
bash download_resources.sh --saige-image
mkdir -p out

cd ~
# Load the Docker image from the tar.gz file
docker load -i "resources/saige.tar"
image_id=$(docker images --filter=reference='wzhou88/saige:*' --format "{{.ID}}" | head -n 1)

## get IDs to be tested (this speeds up analysis)
exec > >(tee ~/data_dir/output.log)
echo "Creating ids_to_include.txt..."

# Calculate padded regions
PADDED_START=$((START_POS - PADDING))
PADDED_END=$((END_POS + PADDING))

# Ensure we don't get negative coordinates
if [ $PADDED_START -lt 0 ]; then
    PADDED_START=0
fi

echo "Analyzing region:"
echo "  Transcript: ${TRANSCRIPT_ID}"
echo "  Original region: ${START_POS}-${END_POS}"
echo "  Padded region: ${PADDED_START}-${PADDED_END}"
echo "  Padding size: ${PADDING}"

# Create the markers file using the local BIM file
cd ~/data_dir
awk -v start="$PADDED_START" -v end="$PADDED_END" -v transcript="$TRANSCRIPT_ID" '
    $2 == transcript || ( $4 >= start && $4 <= end) {
        print $0
    }
' $(basename ${PLINK_FILE}).bim > marker_list.bim

# Create ids_to_include.txt
cut -f2 marker_list.bim > ids_to_include.txt
echo "Original number of variants: $(wc -l < ids_to_include.txt)"

# For every pseudo-variant, find its constituent variants, and map them to the rsid.
if [[ -n "${PSEUDOVARIANT_MAP_FILE}" ]]; then
    echo "Processing pseudo-variant mapping file to exclude constituent variants..."
    
    # Extract variants used in pseudo-variant construction
    zcat "${PSEUDOVARIANT_MAP_FILE}" | grep "${TRANSCRIPT_ID}" | grep -E "(chet)|(hom)" | \
        cut -f6 | tr '|' '\n' | tr ';' '\n' | sort | uniq > variants_to_exclude.txt
    
    # For each variant, find corresponding rsid/marker from bim file
    echo "Mapping variants.."
    echo "Note: Variants are not mapped if they are not the in imputed .bim file (e.g. too rare) which results in NA."
    echo "Found variants (chr:pos:ref:alt) to exclude from pseudo-variant file:"
    > rsids_to_exclude.txt
    while IFS= read -r variant; do
        if [[ ! -z "$variant" ]]; then
            variant=$(echo "$variant" | tr -d '[:space:]')
            rsid=$(find_rsid_from_bim_file "$variant" "marker_list.bim")
            echo "${variant} -> ${rsid}"
            if [[ "$rsid" != "NA" ]]; then
                echo "$rsid" >> rsids_to_exclude.txt
            fi
        fi
    done < variants_to_exclude.txt

    # if there were no variants mapped, then we don't want to do any filtering
    if [[ -s rsids_to_exclude.txt ]]; then
        echo "Mapped and excluding $(wc -l < rsids_to_exclude.txt) variants from pseudo-variants in imputed data:"
        #echo "----------------------------------------"
        #while IFS= read -r rsid; do
        #    echo "$rsid"
        #done < rsids_to_exclude.txt
        #echo "----------------------------------------"
        
        # Create filtered version of ids_to_include.txt while preserving order
        grep -v -F -f rsids_to_exclude.txt ids_to_include.txt > ids_to_include.txt.filtered
        mv ids_to_include.txt.filtered ids_to_include.txt
        echo "Successfully filtered ids_to_include.txt"
    else
        echo "No variants found to exclude for transcript ${TRANSCRIPT_ID}"
    fi
fi

# Log the number of markers
echo "Final number of markers included: $(wc -l < ids_to_include.txt)"
cd ~

# Log the parameters and paths
echo "Analysis parameters:"
echo "-------------------"
echo "Phenotype: $PHENO"
echo "Model file: $LOCAL_MODEL_FILE"
echo "Variance ratio file: $LOCAL_VAR_RATIO_FILE"
echo "GRM file: $LOCAL_SPARSE_GRM"
echo "GRM samples file: $LOCAL_SPARSE_GRM_SAMPLES"
echo "Plink file prefix: $LOCAL_PLINK_PREFIX"
echo "Output prefix: $OUT_PREFIX"
echo "Start position: $START_POS"
echo "End position: $END_POS"
echo "Transcript ID: $TRANSCRIPT_ID"
echo "Maximum iterations: $MAX_ITERATION"
echo "P-value threshold: $THRESHOLD"
echo "-------------------"

# Set up conditional analysis parameters
#readonly OUT_FILE=~/data_dir/saige.out.txt
COMBINED_OUTPUT=combined_iterations.txt

MARKERS_TO_CONDITION_ON=""
SIGNIFICANT=true
ITERATION=0

while [ "$SIGNIFICANT" = true ]; do
    ITERATION=$((ITERATION + 1))
    echo "Starting iteration $ITERATION..."
    cd ~

	     # Check if we've hit the maximum iterations
    if [ "$ITERATION" -gt "$MAX_ITERATION" ]; then
        echo "Reached maximum number of iterations ($MAX_ITERATION)"
        SIGNIFICANT=false
        break
    fi
    
    # Run SAIGE with current conditioning markers
    echo "Running SAIGE step 2 (iteration $ITERATION)..."
    docker run \
        -e HOME="$HOME" \
        -v ~/data_dir:/data_dir \
        "${image_id}" \
        step2_SPAtests.R \
            --bedFile=$LOCAL_PLINK_PREFIX.bed \
            --bimFile=$LOCAL_PLINK_PREFIX.bim \
            --famFile=$LOCAL_PLINK_PREFIX.fam \
            --chrom="${CHROMOSOME}" \
            --condition="$MARKERS_TO_CONDITION_ON" \
            --idstoIncludeFile /data_dir/ids_to_include.txt \
            --minMAF=0.0001 \
            --minMAC=4 \
            --GMMATmodelFile=${LOCAL_MODEL_FILE} \
            --varianceRatioFile=${LOCAL_VAR_RATIO_FILE} \
            --sparseGRMFile=${LOCAL_SPARSE_GRM} \
            --sparseGRMSampleIDFile=${LOCAL_SPARSE_GRM_SAMPLES} \
            --LOCO=FALSE \
            --is_Firth_beta=TRUE \
            --SPAcutoff=0.5 \
            --pCutoffforFirth=0.10 \
            --is_output_moreDetails=TRUE \
            --is_fastTest=TRUE \
            --maxMAF_in_groupTest=0.5 \
            --SAIGEOutputFile=/data_dir/saige.out.txt \
            --dosage_zerod_MAC_cutoff=0 \
            --dosage_zerod_cutoff=0 \
            --dosage_zerod_MAC_cutoff=0

    echo "SAIGE finalized.."
    cd data_dir

    # For first iteration, include header in combined file
    if [ "$ITERATION" -eq 1 ]; then
        standardize_output "saige.out.txt" "$ITERATION" "NA" > "$COMBINED_OUTPUT"
    else
        standardize_output "saige.out.txt" "$ITERATION" "$MARKERS_TO_CONDITION_ON" | tail -n +2 >> "$COMBINED_OUTPUT"
    fi
    
    # Extract p-values and markers
    HEADER=$(head -n1 "saige.out.txt")
    MAIN_TRANSCRIPT_ROW=$(cat "saige.out.txt" | grep -w "${TRANSCRIPT_ID}")
    MIN_PVALUE_ROW=$(extract_min_pvalue_row "saige.out.txt" "${TRANSCRIPT_ID}" "$MARKERS_TO_CONDITION_ON")
    
    # Get the p-values using the helper function
    TRANSCRIPT_PVAL=$(get_pvalue_from_row "$MAIN_TRANSCRIPT_ROW" "$HEADER")
    MIN_PVAL=$(get_pvalue_from_row "$MIN_PVALUE_ROW" "$HEADER")
    MIN_PVAL_MARKER=$(echo "$MIN_PVALUE_ROW" | cut -f3)
    
		echo "Iteration $ITERATION Results:"
    echo "  Current marker being tested: $MIN_PVAL_MARKER"
    echo "  Minimum p-value found: $MIN_PVAL"
    echo "  Transcript p-value: $TRANSCRIPT_PVAL"
  
    # Store first p-value (only on first iteration)
	  if [ "$ITERATION" -eq 1 ]; then
	  		FIRST_PVALUE="$TRANSCRIPT_PVAL"
	  fi

	  # Store the latest p-value in each iteration
	  FINAL_PVALUE="$TRANSCRIPT_PVAL" 

    # First check if the transcript's p-value has become non-significant
    if awk -v pval="$TRANSCRIPT_PVAL" -v thresh="$THRESHOLD" 'BEGIN{exit !(pval >= thresh)}'; then
        SIGNIFICANT=false
        echo "  Success! Transcript p-value ($TRANSCRIPT_PVAL) is now above threshold"
        break
    fi
    
    # If transcript is still significant, check if there are any other significant signals
    if awk -v pval="$MIN_PVAL" -v thresh="$THRESHOLD" 'BEGIN{exit !(pval < thresh)}'; then
        # Add the new marker to conditioning list
        if [ -z "$MARKERS_TO_CONDITION_ON" ]; then
            MARKERS_TO_CONDITION_ON="$MIN_PVAL_MARKER"
        else
            MARKERS_TO_CONDITION_ON="${MARKERS_TO_CONDITION_ON},${MIN_PVAL_MARKER}"
        fi
        echo "  Found significant marker to condition on"
        echo "  Updated conditioning markers: $MARKERS_TO_CONDITION_ON"
    else
        SIGNIFICANT=false
        echo "  No more significant associations found"
        echo "  Final conditioning markers: $MARKERS_TO_CONDITION_ON"
    fi
    echo "----------------------------------------"
done

echo "Conditional analysis complete after $ITERATION iterations"
echo "Final conditioning markers: $MARKERS_TO_CONDITION_ON"
echo "Initial transcript p-value: $FIRST_PVALUE"
echo "Final transcript p-value: $FINAL_PVALUE"
echo "Combined results saved to: $COMBINED_OUTPUT"

# copy final file to output
#cp ${COMBINED_OUTPUT} /home/dnanexus/out/out/${OUT_PREFIX}.txt
gzip -c ${COMBINED_OUTPUT} > /home/dnanexus/out/out/${OUT_PREFIX}.txt.gz
cp ids_to_include.txt /home/dnanexus/out/out/${OUT_PREFIX}.markers.txt
cp ~/data_dir/output.log /home/dnanexus/out/out/${OUT_PREFIX}.log
