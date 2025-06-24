#!/usr/bin/env python3
"""
Split FASTA file into chunks for parallel processing.
"""

import argparse
import sys
from Bio import SeqIO
from pathlib import Path

def split_fasta(input_fasta, output, chunk_size):
    chunk_num = 1
    current_chunk = []
    current_size = 0
    
    for record in SeqIO.parse(input_fasta, "fasta"):
        current_chunk.append(record)
        current_size += 1
        
        if current_size >= chunk_size:
            # Write current chunk
            output_file = f"{output}_chunk_{chunk_num:03d}.fasta"
            SeqIO.write(current_chunk, output_file, "fasta")
            print(f"Written chunk {chunk_num} with {len(current_chunk)} contigs ({current_size:,} bp)")
            
            # Reset for next chunk
            chunk_num += 1
            current_chunk = []
            current_size = 0
    
    # Write remaining sequences
    if current_chunk:
        output_file = f"chunk_{chunk_num:03d}.fasta"
        SeqIO.write(current_chunk, output_file, "fasta")
        print(f"Written chunk {chunk_num} with {len(current_chunk)} contigs ({current_size:,} bp)")
    
    print(f"Total chunks created: {chunk_num}")

def main():
    parser = argparse.ArgumentParser(description="Split FASTA file into chunks")
    parser.add_argument("input_fasta", help="Input FASTA file")
    parser.add_argument("output", help="Output prefix for chunks")
    parser.add_argument("--chunk-size", type=int, default=20, 
                       help="Target chunk size in number of sequences (default: 20 contigs)")
    
    args = parser.parse_args()
    
    split_fasta(args.input_fasta, args.output, args.chunk_size)

if __name__ == "__main__":
    main()
