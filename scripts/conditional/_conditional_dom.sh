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

# Help function
usage() {
    echo "Usage: $0
    --pheno <phenotype>
    --model_file <SAIGE model file>
    --var_ratio_file <variance ratio file>
    --grm <sparse GRM matrix file>
    --grm_samples <GRM sample IDs file>
    --plink_recessive_file <plink file prefix>
    --vcf_dominance_file <VCF file>
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
PLINK_RECESSIVE_FILE=""
VCF_DOMINANCE_FILE=""
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
        --plink_recessive_file)
            PLINK_RECESSIVE_FILE="$2"
            shift 2
            ;;
        --vcf_dominance_file)
            VCF_DOMINANCE_FILE="$2"
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
   [[ -z "$GRM" ]] || [[ -z "$GRM_SAMPLES" ]] || [[ -z "$PLINK_RECESSIVE_FILE" ]] || \
   [[ -z "$VCF_DOMINANCE_FILE" ]] || [[ -z "$OUT_PREFIX" ]] || [[ -z "$START_POS" ]] || \
   [[ -z "$END_POS" ]] || [[ -z "$TRANSCRIPT_ID" ]] || [[ -z "$CHROMOSOME" ]]; then
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
echo "prefix: ${PLINK_RECESSIVE_FILE}"
cp ${PLINK_RECESSIVE_FILE}.vcf.gz .
cp ${PLINK_RECESSIVE_FILE}.vcf.gz.csi .
# still using .bim file for padding. Need better solution
# but will do for now.
cp ${PLINK_RECESSIVE_FILE}.bim . 

#cp ${PLINK_RECESSIVE_FILE}.bed .
#cp ${PLINK_RECESSIVE_FILE}.fam .

# Download VCF file
echo "Downloading VCF file and index..."
cp ${VCF_DOMINANCE_FILE}.vcf.gz .
cp ${VCF_DOMINANCE_FILE}.vcf.gz.csi .

# Return to home directory
cd ~

# Set up local paths for all files
LOCAL_MODEL_FILE="/data_dir/$(basename ${MODEL_FILE})"
LOCAL_VAR_RATIO_FILE="/data_dir/$(basename ${VAR_RATIO_FILE})"
LOCAL_SPARSE_GRM="/data_dir/$(basename ${GRM})"
LOCAL_SPARSE_GRM_SAMPLES="/data_dir/$(basename ${GRM_SAMPLES})"
LOCAL_PLINK_RECESSIVE_PREFIX="/data_dir/$(basename ${PLINK_RECESSIVE_FILE})"
LOCAL_VCF_RECESSIVE_PREFIX="/data_dir/$(basename ${PLINK_RECESSIVE_FILE})"
LOCAL_VCF_DOMINANCE_PREFIX="/data_dir/$(basename ${VCF_DOMINANCE_FILE})"


# install BCFtools
#cd ~ & mkdir software && cd software
#git clone --recurse-submodules https://github.com/samtools/htslib.git
#git clone https://github.com/samtools/bcftools.git
#cd bcftools && make
#export PATH=~/software/bcftools/:$PATH
#cd ../htslib && make
#export PATH=~/software/htslib/:$PATH
#unset DX_WORKSPACE_ID && dx cd $DX_PROJECT_CONTEXT_ID: &&
#echo "done!" && cd ~ 

# Exit on error or undefined variables
# Needs to be after BCFtools!
set -e
set -u

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
' $(basename ${LOCAL_PLINK_RECESSIVE_PREFIX}).bim > marker_list.bim

# Create ids_to_include.tx
awk '{print "chr"$1":"$4":"$6":"$5}' marker_list.bim > ids_to_include.txt
echo "Found $(wc -l ids_to_include.txt) variants after restricting to padded regions.."

# For every pseudo-variant, find its constituent variants. We do not want to 
# test a variant if it was used to make up a pseudo-variant.
if [[ -n "${PSEUDOVARIANT_MAP_FILE}" ]]; then
    echo "Processing pseudo-variant mapping file to exclude constituent variants..."
    
    # Extract variants used in pseudo-variant construction
    zcat "${PSEUDOVARIANT_MAP_FILE}" | grep "${TRANSCRIPT_ID}" | grep -E "(chet)|(hom)" | \
        cut -f6 | tr '|' '\n' | tr ';' '\n' | sort | uniq > ids_to_exclude.txt
   
    echo "Found $(wc -l ids_to_exclude.txt) variant(s) that make up the pseudo variant."

    # Create filtered version of ids_to_include.txt while preserving order
    grep -v -F -f ids_to_exclude.txt ids_to_include.txt > ids_to_include.txt.filtered
    mv ids_to_include.txt.filtered ids_to_include.txt
    echo "..$(wc -l ids_to_include.txt) variants left after excluding those in pseudo variants."
fi

cd ~

# Log the parameters and paths
echo "Analysis parameters:"
echo "-------------------"
echo "Phenotype: $PHENO"
echo "Model file: $LOCAL_MODEL_FILE"
echo "Variance ratio file: $LOCAL_VAR_RATIO_FILE"
echo "GRM file: $LOCAL_SPARSE_GRM"
echo "GRM samples file: $LOCAL_SPARSE_GRM_SAMPLES"
echo "Plink recessivefile prefix: $LOCAL_PLINK_RECESSIVE_PREFIX"
echo "VCF dominance file prefix: $LOCAL_VCF_DOMINANCE_PREFIX"
echo "Output prefix: $OUT_PREFIX"
echo "Start position: $START_POS"
echo "End position: $END_POS"
echo "Transcript ID: $TRANSCRIPT_ID"
echo "Maximum iterations: $MAX_ITERATION"
echo "P-value threshold: $THRESHOLD"
echo "-------------------"


# first we run SAIGE with the recessivei encoding to 
# simply get a list of markers with at least 5 biallelic
# carriers. 
echo "Running SAIGE to get biallelic carriers.."
docker run \
  -e HOME="$HOME" \
  -v ~/data_dir:/data_dir \
  "${image_id}" \
  step2_SPAtests.R \
      --vcfFile=$LOCAL_VCF_RECESSIVE_PREFIX.vcf.gz \
      --vcfFileIndex=$LOCAL_VCF_RECESSIVE_PREFIX.vcf.gz.csi \
      --idstoIncludeFile /data_dir/ids_to_include.txt \
      --chrom="${CHROMOSOME}" \
      --minMAF=0.0001 \
      --minMAC=10 \
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
      --SAIGEOutputFile=/data_dir/saige.markers.out.txt \
      --dosage_zerod_cutoff=0 \
      --dosage_zerod_MAC_cutoff=0


# now we re-write IDs to include to only the ones with at least X (see above) homozygous 
# carriers. This is because we are testing a file that has been recoded [0,1,2]->[0,0,1] 
cat ~/data_dir/saige.markers.out.txt | tail -n+2 | awk '{print $1":"$2":"$4":"$5}' > ~/data_dir/ids_to_include.txt
echo "Found $(wc -l ~/data_dir/ids_to_include.txt) markers with at least some homozygotes.."

set -x
# now we need to add the pseudo variant we care about
echo "Adding transcript ID to list of ids to include.."
zgrep -m 1 ${TRANSCRIPT_ID} ~/data_dir/$(basename ${LOCAL_VCF_DOMINANCE_PREFIX}).vcf.gz | cut -f1-10 > enst.txt
cat enst.txt | awk '{print $1":"$2":"$4":"$5}' >> ~/data_dir/ids_to_include.txt

# note: we MUST sort the list! Otherwise, saige can not read the markers!
sort -t: -k2,2n ~/data_dir/ids_to_include.txt > ~/data_dir/sorted_ids_to_include.txt
head ~/data_dir/sorted_ids_to_include.txt

# use BCFtools to reduce the VCF to the region and the required marker
readonly SNP_REGION="chr${CHROMOSOME}:${PADDED_START}-${PADDED_END}"
readonly GENE_START="$(cut -f2 enst.txt | awk '{print ($1-1)}')"
readonly GENE_END="$(cut -f2 enst.txt | awk '{print ($1+1)}')"
readonly GENE_REGION="chr${CHROMOSOME}:${GENE_START}-${GENE_END}"
readonly BCF_QUERY="${GENE_REGION},${SNP_REGION}"

# here we extract the specific gene region and the SNPs in the padded area
bcftools view -r "${BCF_QUERY}" ~/data_dir/$(basename ${LOCAL_VCF_DOMINANCE_PREFIX}).vcf.gz -Oz -o ~/data_dir/$(basename ${LOCAL_VCF_DOMINANCE_PREFIX}).filtered.vcf.gz
tabix -C ~/data_dir/$(basename ${LOCAL_VCF_DOMINANCE_PREFIX}).filtered.vcf.gz

echo "Validating the presence of ID '${TRANSCRIPT_ID}' in VCF:"
zcat ~/data_dir/$(basename ${LOCAL_VCF_DOMINANCE_PREFIX}).filtered.vcf.gz | grep ${TRANSCRIPT_ID} | cut -f1-12

if [[ $(zcat ~/data_dir/$(basename ${LOCAL_VCF_DOMINANCE_PREFIX}).filtered.vcf.gz | grep ${TRANSCRIPT_ID} | wc -l) -eq "0" ]]; then
  >&2 echo "Transcript ID is not available in the VCF! Printing full VCF (col 1-12) and exiting.."
  zcat ~/data_dir/$(basename ${LOCAL_VCF_DOMINANCE_PREFIX}).filtered.vcf.gz | cut -f1-12
  exit 1
fi


readonly n_final_markers_to_test=$(cat ~/data_dir/sorted_ids_to_include.txt | wc -l)
echo "Done! ${n_final_markers_to_test} markers present after adding the pseudo-variant."
echo "----------------------------------------------------------------"
# Set up conditional analysis parameters
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
            --vcfFile=$LOCAL_VCF_DOMINANCE_PREFIX.filtered.vcf.gz \
            --vcfFileIndex=$LOCAL_VCF_DOMINANCE_PREFIX.filtered.vcf.gz.csi \
            --chrom="chr${CHROMOSOME}" \
            --condition="$MARKERS_TO_CONDITION_ON" \
            --minMAF=0 \
            --minMAC=0.5 \
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
            --dosage_zerod_cutoff=0 \
            --dosage_zerod_MAC_cutoff=0

    cd data_dir

    echo "Looking for transcript ID in saige.out.txt:"
    cat saige.out.txt | grep ${TRANSCRIPT_ID}
    echo "done..!"

   	# Save the header from SAIGE output
		head -n 1 saige.out.txt > saige.filtered.out.txt

		# Create formatted lookup from SAIGE output (CHR:POS:A1:A2)
		awk 'NR>1 {print $1":"$2":"$4":"$5}' saige.out.txt > saige_markers.tmp

		# Get line numbers of matching markers (including header line 1)
		(echo "1"; grep -nf sorted_ids_to_include.txt saige_markers.tmp | cut -d: -f1) | sort -n | \
		awk '{print $1+1}' > matching_lines.tmp

		# Extract those lines from original SAIGE file
		awk 'NR==FNR{lines[$1]; next} FNR in lines' matching_lines.tmp saige.out.txt >> saige.filtered.out.txt

    echo "Looking for transcript ID in saige.filtered.out.txt:"
   	if [[ ! $(grep ${TRANSCRIPT_ID} saige.filtered.out.txt) ]]; then
      echo "Error: grep failed" >&2
      exit 1
    fi 
    
    cat saige.filtered.out.txt | grep ${TRANSCRIPT_ID}
    echo "done again..!"

    # how many markers that are in our criteria were tested
    n_markers_tested=$(cat saige.filtered.out.txt | tail -n+2 | wc -l)
    echo "${n_markers_tested} markers were carried forward after filtering.."

    # For first iteration, include header in combined file
    if [ "$ITERATION" -eq 1 ]; then
        standardize_output "saige.filtered.out.txt" "$ITERATION" "NA" > "$COMBINED_OUTPUT"
    else
        standardize_output "saige.filtered.out.txt" "$ITERATION" "$MARKERS_TO_CONDITION_ON" | tail -n +2 >> "$COMBINED_OUTPUT"
    fi
   
    # Extract p-values and markers
    HEADER=$(head -n1 "saige.filtered.out.txt")
    MAIN_TRANSCRIPT_ROW=$(cat "saige.filtered.out.txt" | grep -w "${TRANSCRIPT_ID}")
    MIN_PVALUE_ROW=$(extract_min_pvalue_row "saige.filtered.out.txt" "${TRANSCRIPT_ID}" "$MARKERS_TO_CONDITION_ON")
    
    # Get the p-values using the helper function
    TRANSCRIPT_PVAL=$(get_pvalue_from_row "$MAIN_TRANSCRIPT_ROW" "$HEADER")
    MIN_PVAL=$(get_pvalue_from_row "$MIN_PVALUE_ROW" "$HEADER")
    MIN_PVAL_MARKER=$(echo "$MIN_PVALUE_ROW" | awk '{print $1":"$2":"$4":"$5}')
    
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

echo "Conditional (Dominance) analysis complete after $ITERATION iterations"
echo "Final conditioning markers: $MARKERS_TO_CONDITION_ON"
echo "Initial transcript p-value: $FIRST_PVALUE"
echo "Final transcript p-value: $FINAL_PVALUE"
echo "Combined results saved to: $COMBINED_OUTPUT"

# copy final file to output
#cp ${COMBINED_OUTPUT} /home/dnanexus/out/out/${OUT_PREFIX}.txt
gzip -c ${COMBINED_OUTPUT} > /home/dnanexus/out/out/${OUT_PREFIX}.txt.gz
cp sorted_ids_to_include.txt /home/dnanexus/out/out/${OUT_PREFIX}.markers.txt
cp ~/data_dir/output.log /home/dnanexus/out/out/${OUT_PREFIX}.log
