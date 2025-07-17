#!/bin/bash

# Simple script to extract contigs from FASTA based on boundaries.tsv
# Usage: ./extract_contigs.sh boundaries.tsv assembly.fasta > output.fasta

if [ $# -ne 2 ]; then
    echo "Usage: $0 <boundaries.tsv> <assembly.fasta>" >&2
    exit 1
fi

BOUNDARIES_FILE="$1"
ASSEMBLY_FASTA="$2"

# Check input files
if [ ! -f "$BOUNDARIES_FILE" ] || [ ! -f "$ASSEMBLY_FASTA" ]; then
    echo "Error: Input file(s) not found!" >&2
    exit 1
fi

# Get unique contig names and create grep pattern using awk instead of paste
CONTIGS=$(tail -n +2 "$BOUNDARIES_FILE" | cut -f2 | sort -u)
PATTERN=$(echo "$CONTIGS" | awk '{printf "%s%s", (NR>1 ? "|" : ""), ">" $0} END {print ""}')

echo "Extracting contigs: $(echo "$CONTIGS" | tr '\n' ' ')" >&2

# Extract sequences
awk -v pattern="$PATTERN" '
/^>/ {
    if (match($1, "^(" pattern ")$")) {
        print_seq = 1
        print $0
        count++
    } else {
        print_seq = 0
    }
    next
}
{
    if (print_seq) print $0
}
END {
    print "Extracted " count " contigs" > "/dev/stderr"
}
' "$ASSEMBLY_FASTA"
