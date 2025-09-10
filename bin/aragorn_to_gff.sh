#!/bin/bash

# Script to convert ARAGORN output to GFF format
# Usage: aragorn_to_gff.sh <input_tbl_file> <output_gff_file>

input_file="$1"
output_file="$2"

# Add GFF header
echo "##gff-version 3" > "$output_file"

# Process the ARAGORN output with proper whitespace handling
awk '
BEGIN { 
    OFS="\t"
    feature_count = 0
} 

# Parse header line to get sequence name
/^>/ && !/^>end/ {
    sequence_name = $1
    gsub(/^>/, "", sequence_name)
    next
}

# Skip lines with "genes found" or "end"
/ found/ || /^>end/ { next }

# Process tRNA lines (lines starting with a number)
/^[0-9]+/ {
    feature_count++
    
    # First, normalize multiple spaces/tabs to single spaces for field splitting
    # But preserve the coordinate brackets
    line = $0
    gsub(/[ \t]+/, " ", line)  # Convert multiple spaces/tabs to single space
    
    # Split the normalized line into fields
    split(line, fields, " ")
    
    # Extract tRNA type (should be field 2 after normalization)
    trna_type = fields[2]
    
    # Extract coordinates (should be field 3 after normalization)
    coord_field = fields[3]
    
    # Determine strand and extract coordinates
    strand = "+"
    if (match(coord_field, /^c\[/)) {
        strand = "-"
        gsub(/^c/, "", coord_field)
    }

    gsub(/\[/, "", coord_field)
    gsub(/\]/, "", coord_field)
    split(coord_field, coords, ",")
    start = coords[1]
    end = coords[2]

    # Print GFF3 line
    print sequence_name, "aragorn", "tRNA", start, end, ".", strand, ".", "ID=tRNA" feature_count ";product=" trna_type

}
' "$input_file" >> "$output_file"

