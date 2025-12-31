#!/usr/bin/env bash
# author: frederik lassen
# note: Find and remove duplicate files in step1 directory (keep newest)

set -o errexit
set -o nounset

readonly step1_dir="/wes_ko_ukbb/data/saige/step1/simulated/revisions/051025"

echo "=========================================="
echo "Finding duplicate files in ${step1_dir}"
echo "=========================================="
echo ""

# Get all files with full details
temp_file=$(mktemp)
dx ls -l "${step1_dir}" > "${temp_file}"

# Parse and find duplicates
echo "Checking for duplicates..."
echo ""

duplicates_found=0

# Read the file list and extract filenames
while read -r line; do
    # Skip empty lines and directory listings
    if [[ -z "$line" || "$line" == *"/"* ]]; then
        continue
    fi

    # Extract filename (last field before the file-ID)
    # Format: "closed  2025-11-11 14:03:56 2.44 MB   filename (file-ID)"
    filename=$(echo "$line" | awk '{for(i=1;i<=NF;i++) if($i ~ /^\(file-/) print $(i-1)}')

    if [[ -n "$filename" ]]; then
        # Count how many times this filename appears
        count=$(grep -c " ${filename} " "${temp_file}" || true)

        if [[ $count -gt 1 ]]; then
            duplicates_found=$((duplicates_found + 1))

            echo "Found duplicate: ${filename} (${count} copies)"
            echo "----------------------------------------"

            # Show all versions with index numbers
            grep " ${filename} " "${temp_file}" | nl -v 0
            echo ""

            # Get the newest file (first line from grep)
            newest_line=$(grep " ${filename} " "${temp_file}" | head -1)
            newest_id=$(echo "$newest_line" | grep -o 'file-[A-Za-z0-9]*')
            newest_date=$(echo "$newest_line" | awk '{print $2, $3}')

            echo "Keeping newest version:"
            echo "  Date: ${newest_date}"
            echo "  ID: ${newest_id}"
            echo ""

            # Remove older versions (all except first)
            line_num=0
            while IFS= read -r old_line; do
                line_num=$((line_num + 1))
                if [[ $line_num -gt 1 ]]; then
                    old_id=$(echo "$old_line" | grep -o 'file-[A-Za-z0-9]*')
                    old_date=$(echo "$old_line" | awk '{print $2, $3}')

                    echo "Removing older version:"
                    echo "  Date: ${old_date}"
                    echo "  ID: ${old_id}"

                    dx rm "${old_id}"
                    echo "  ✓ Removed"
                fi
            done < <(grep " ${filename} " "${temp_file}")

            echo ""
            echo "=========================================="
            echo ""
        fi
    fi
done < <(sort -u "${temp_file}")

rm "${temp_file}"

if [[ $duplicates_found -eq 0 ]]; then
    echo "No duplicates found!"
else
    echo "Cleaned up ${duplicates_found} duplicate file(s)"
fi

echo ""
echo "Done!"
