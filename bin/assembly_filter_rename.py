#!/usr/bin/env python

from Bio import SeqIO
import sys
import argparse
import os

##### This script filter by contig size a metagenomic assembly and generates two files to feed the MoMofy pipeline
##### Alejandra Escobar, EMBL-EBI
##### January 19th, 2023

def rename(input_file):
    output_1kb = "1kb_contigs.fasta"
    output_5kb = "5kb_contigs.fasta"
    output_map = "contigID.map"
    counter = 0

    with open(output_1kb, "w") as to_1kb, \
        open(output_5kb, "w") as to_5kb, \
        open(output_map, "w") as to_map:

        for record in SeqIO.parse(input_file, "fasta"):
            counter += 1
            new_id = ">contig_" + str(counter)
            my_chain = str(record.seq).upper()
            to_map.write(new_id + "\t" + str(record.id) + "\n")
            if len(my_chain) > 1000:
                to_1kb.write(new_id + "\n")
                to_1kb.write(my_chain + "\n")
            if len(my_chain) > 5000:
                to_5kb.write(new_id + "\n")
                to_5kb.write(my_chain + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="This script filter by contig size an assembly and generates two files to feed the MoMofy pipeline: contigs > 1kb for insertion sequences prediction and contigs > 5 kb for integron prediction"
    )
    parser.add_argument(
        '--assembly',
        type=str, 
        help="Input fasta file",
        required=True,
    )
    args = parser.parse_args()

    input_file = args.assembly

    rename(input_file)


if __name__ == "__main__":
    main()

