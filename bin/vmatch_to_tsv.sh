#!/bin/bash

# Script to convert vmatch output to GFF format
# Usage: vmatch_to_tsv.sh <input_vmatch_file> <output_tsv_file>

input_file="$1"
output_file="$2"

# Process the vmatch output and append to the GFF file
awk '!/^#/ {
    ice_element = $2
    start_1 = $3 + 1
    end_1 = $3 + $1
    start_2 = $7 + 1
    end_2 = $7 + $5
    print ice_element "\t" start_1 "\t" end_1 "\t" start_2 "\t" end_2
}' "$input_file" > "$output_file"

