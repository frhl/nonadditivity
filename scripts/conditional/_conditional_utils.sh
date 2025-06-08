# Add helper functions for conditional analysis
function get_pvalue_from_row() {
    local row="$1"
    local header="$2"

    # Determine which p-value column to use
    local p_value_col
    if [[ $header == *"p.value_c"* ]]; then
        # Condition ON - find p.value_c column
        p_value_col=$(echo "$header" | tr '\t' '\n' | grep -n "p.value_c" | cut -d: -f1)
    else
        # Condition OFF - find p.value column
        p_value_col=$(echo "$header" | tr '\t' '\n' | grep -n "p.value$" | cut -d: -f1)
    fi

    echo "$row" | cut -f"$p_value_col"
}

function extract_min_pvalue_row() {
    local outfile="$1"
    local transcript_id="$2"
    local exclude_markers="$3"  # Optional: comma-separated list of markers to exclude

    # Check inputs
    if [[ ! -f "$outfile" ]]; then
        echo "Error: File $outfile does not exist" >&2
        return 1
    fi

    if [[ -z "$transcript_id" ]]; then
        echo "Error: transcript_id not provided" >&2
        return 1
    fi

    # Get header and check format
    local header=$(head -n1 "$outfile")
    local p_value_col

    if [[ $header == *"p.value_c"* ]]; then
        # Condition ON - find p.value_c column
        p_value_col=$(echo "$header" | tr '\t' '\n' | grep -n "p.value_c" | cut -d: -f1)
    else
        # Condition OFF - find p.value column
        p_value_col=$(echo "$header" | tr '\t' '\n' | grep -n "p.value$" | cut -d: -f1)
    fi

    if [[ -z "$p_value_col" ]]; then
        echo "Error: Could not find p-value column" >&2
        return 1
    fi

    # Prepare exclude pattern for awk
    local exclude_pattern=""
    if [[ ! -z "$exclude_markers" ]]; then
        # Convert comma-separated list to awk pattern
        exclude_pattern=$(echo "$exclude_markers" | tr ',' '|')
    fi

    # Get all rows except the header, transcript row, and excluded markers, then sort by p-value
    local min_pvalue_row
    if [[ ! -z "$exclude_pattern" ]]; then
        min_pvalue_row=$(awk -v transcript="$transcript_id" \
            -v col="$p_value_col" \
            -v exclude="$exclude_pattern" \
            'NR>1 && $0 !~ transcript && $3 !~ "^("exclude")$" {print $0}' "$outfile" | \
            sort -k${p_value_col},${p_value_col} -g | \
            head -n1)
    else
        min_pvalue_row=$(awk -v transcript="$transcript_id" \
            -v col="$p_value_col" \
            'NR>1 && $0 !~ transcript {print $0}' "$outfile" | \
            sort -k${p_value_col},${p_value_col} -g | \
            head -n1)
    fi

    if [[ -z "$min_pvalue_row" ]]; then
        echo "Error: Could not find any rows with p-values" >&2
        return 1
    fi

    echo "$min_pvalue_row"
    return 0
}

function standardize_output() {
    local input_file="$1"
    local iteration="$2"
    local condition_markers="$3"

    # First determine if we're in conditional mode by checking header
    local header=$(head -n1 "$input_file")

    if [[ $header == *"BETA_c"* ]]; then
        # WITH conditioning - need to find and extract _c columns
        awk -v iter="$iteration" -v cond="$condition_markers" '
        BEGIN {
            OFS="\t"
        }
        NR==1 {
            # Find positions of conditional columns
            split($0, cols, /[ \t]+/)
            for(i=1; i<=NF; i++) {
                if(cols[i] == "BETA_c") beta_c_col = i
                if(cols[i] == "SE_c") se_c_col = i
                if(cols[i] == "p.value_c") p_c_col = i
            }

            # Print standardized header
            print "iteration", "CHR", "POS", "MarkerID", "Allele1", "Allele2", "AF", \
                  "BETA", "SE", "p.value", \
                  "BETA_c", "SE_c", "p.value_c", "conditioning_markers"
            next
        }
        {
            # For data rows, extract base columns and conditional columns
            # Note: p.value is column 13, AF_Allele2 is column 7
            print iter, $1, $2, $3, $4, $5, $7, \
                  $9, $10, $13, \
                  $beta_c_col, $se_c_col, $p_c_col, cond
        }' "$input_file"
    else
        # NO conditioning - simpler format
        awk -v iter="$iteration" '
        BEGIN {
            OFS="\t"
        }
        NR==1 {
            print "iteration", "CHR", "POS", "MarkerID", "Allele1", "Allele2", "AF", \
                  "BETA", "SE", "p.value", \
                  "BETA_c", "SE_c", "p.value_c", "conditioning_markers"
            next
        }
        {
            # Note: p.value is column 13, AF_Allele2 is column 7
            print iter, $1, $2, $3, $4, $5, $7, \
                  $9, $10, $13, \
                  "NA", "NA", "NA", "NA"
        }' "$input_file"
    fi
}

find_rsid_from_bim_file() {
    local variant=$1
    local call_chets_file=$2

    local variant=${variant#chr}
    local chr=$(echo $variant | cut -d: -f1)
    local pos=$(echo $variant | cut -d: -f2)
    local ref=$(echo $variant | cut -d: -f3)
    local alt=$(echo $variant | cut -d: -f4)

    # Store result in variable to check if empty
    local result=$(grep -v "ENST" $call_chets_file | awk -v c=$chr -v p=$pos -v r=$ref -v a=$alt \
        '($1==c && $4==p && $5==r && $6==a) || ($1==c && $4==p && $5==a && $6==r) {print $2}')

    if [ -z "$result" ]; then
        echo "Warning: No match found for chr${chr}:${pos}:${ref}:${alt}" >&2
        echo "NA"
    else
        echo "$result"
    fi
}


