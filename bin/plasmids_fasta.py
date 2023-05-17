#!/usr/bin/env python

from Bio import SeqIO
import argparse
import sys
import os.path
import glob


##### This script generates the fasta file of contigs labelled as plasmids by ppr-meta
##### Alejandra Escobar, EMBL-EBI
##### May 17, 2023


def gff_parser(gff_in):
    ### Saving contig names equivalence
    plasmids_list = {}
    if os.path.exists(gff_in):
        with open(gff_in, "r") as input_gff:
            for line in input_gff:
                line_l = line.rstrip().split("\t")
                if len(line_l) == 9 :
                    (
                        contig,
                        source,
                        seq_type,
                        start,
                        end,
                        score,
                        strand,
                        phase,
                        attributes,
                    ) = line.rstrip().split("\t")
                    if seq_type == 'plasmid':
                        plasmid_id = attributes.split(';')[0].replace('ID=','')
                        plasmids_list[contig]=plasmid_id
    return(plasmids_list)


def fasta_writer(plasmids_list, fasta_in):
    with open('plasmids.fasta', "w") as to_fasta:
        if len(plasmids_list) > 0:
            for record in SeqIO.parse(fasta_in, "fasta"):
                my_chain = str(record.seq).upper()
                my_id = str(record.id)
                if my_id in plasmids_list:
                    plasmid_id = plasmids_list[my_id]
                    to_fasta.write('>'+my_id+'|'+plasmid_id+'\n')
                    to_fasta.write(my_chain+'\n')


def main():
    parser = argparse.ArgumentParser(
        description="This script generates the fasta file of contigs labelled as plasmids by ppr-meta. Please provide the relevant input files"
    )
    parser.add_argument(
        "--fasta",
        type=str,
        help="Original fasta file",
        required=True,
    )
    parser.add_argument(
        "--gff",
        type=str, 
        help="mobilome_nogenes.gff file", 
        required=True,
    )
    args = parser.parse_args()


    ### Setting up variables
    fasta_in = args.fasta
    gff_in = args.gff

    ## Calling functions
    # Saving list of plasmid contigs
    plasmids_list = gff_parser(gff_in)

    # Writting the fasta file
    fasta_writer(plasmids_list, fasta_in)


if __name__ == "__main__":
    main()


