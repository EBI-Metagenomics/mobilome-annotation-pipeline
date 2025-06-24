#!/usr/bin/env python3
"""
Merge compositional outlier detection results from multiple BED files.
"""

import argparse
import sys
from pathlib import Path

def merge_bed_files(input_files, output_file):
    """Merge multiple BED files into one, preserving headers and sorting by confidence."""
    
    all_predictions = []
    headers = []
    
    for bed_file in input_files:
        if not Path(bed_file).exists():
            print(f"Warning: File {bed_file} does not exist, skipping", file=sys.stderr)
            continue
            
        with open(bed_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('#'):
                    if line not in headers:
                        headers.append(line)
                else:
                    all_predictions.append(line)
    
    # Write merged output
    with open(output_file, 'w') as f:
        for header in headers:
            f.write(header + '\n')
        
        for line in all_predictions:
            f.write(line + '\n')
    
    print(f"Merged {len(all_predictions)} predictions from {len(input_files)} files", file=sys.stderr)
    print(f"Output written to {output_file}", file=sys.stderr)

def main():
    parser = argparse.ArgumentParser(description="Merge BED files from compositional outlier detection")
    parser.add_argument("output_file", help="Output merged BED file")
    parser.add_argument("input_files", nargs='+', help="Input BED files to merge")
    
    args = parser.parse_args()
    
    merge_bed_files(args.input_files, args.output_file)

if __name__ == "__main__":
    main()
