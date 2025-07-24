#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from Bio import SeqIO
import pandas as pd


def macsyfinder_parser(macsy_input):
    macsy_prots = []

    # Read the file, skipping comment lines that start with #
    with open(self.all_systems_file, "r") as f:
        lines = [line for line in f.readlines() if line.strip()]

    # Find the header line (first non-comment line)
    header_idx = 0
    for i, line in enumerate(lines):
        if not line.startswith("#") and line.strip():
            header_idx = i
            break

    # Read the data starting from the header
    df = pd.read_csv(self.all_systems_file, sep="\t", skiprows=header_idx)

    if "hit_id" in df.columns:
        macsy_prots = df["hit_id"].tolist()
    else:
        print("Warning: 'hit_id' column not found in the macsyfinder output")

    macsy_prots = set(macsy_prots)

    return macsy_prots


def gff_parser(gff):
    proteins_contigs = {}

    with open(gff, "r") as input_table:
        for line in input_table:
            l_line = line.rstrip().split("\t")
            # Annotation lines have exactly 9 columns
            if len(l_line) == 9:
                contig = l_line[0]
                attrib = l_line[8]
                entry_id = attrib.split(";")[0].replace("ID=", "")
                proteins_contigs[entry_id] = contig
    return proteins_contigs


def contigs_finder(macsy_prots, proteins_contigs):
    contigs_list = set()
    for protein in macsy_prots:
        if protein in proteins_contigs:
            contigs_list.add(proteins_contigs[protein])

    return contigs_list


def fasta_parser(contigs_list, assembly, output_file):
    with open(output_file, "w") as output:
        for record in SeqIO.parse(assembly, "fasta"):
            seq_id = str(record.id)
            if seq_id in contigs_list:
                output.write(">" + seq_id + "\n")
                output.write(str(record.seq).upper() + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="Process macsyfinder output for boundary delineation"
    )
    parser.add_argument(
        "--macsy_input",
        required=True,
        help="Path to all_systems.tsv file from macsyfinder",
    )
    parser.add_argument("--assembly", required=True, help="Path to genome FASTA file")
    parser.add_argument("--gff", help="Path to annotation file on GFF3 format")
    parser.add_argument("--output", required=True, help="Output file")
    args = parser.parse_args()

    # Validate input files
    if not os.path.exists(args.macsy_input):
        print(f"Error: all_systems.tsv file not found: {args.macsy_input}")
        sys.exit(1)

    if not os.path.exists(args.assembly):
        print(f"Error: Genome file not found: {args.assembly}")
        sys.exit(1)

    if not os.path.exists(args.gff):
        print(f"Error: Genome file not found: {args.gff}")
        sys.exit(1)

    # Initialize processor
    macsy_prots = macsyfinder_parser(args.macsy_input)
    proteins_contigs = gff_parser(args.gff)
    contigs_list = contigs_finder(macsy_prots, proteins_contigs)
    fasta_parser(contigs_list, args.assembly, args.output)


if __name__ == "__main__":
    main()
