#!/usr/bin/env python

from Bio import SeqIO
import argparse
import sys
import os.path
import glob


##### This script parse the virify output for quality control
##### Alejandra Escobar, EMBL-EBI
##### October 3 2023


def quality_parser(checkv_out):
    checkv_pass = []
    for summ_file in checkv_out:
        if os.path.exists(summ_file):
            head, tail = os.path.split(summ_file)
            current_qual = tail.split("_")[0]
            with open(summ_file, "r") as input_table:
                next(input_table)
                for line in input_table:
                    l_line = line.rstrip().split("\t")
                    phage_id = l_line[0].replace("prophage-0:", "prophage-1:")
                    viral_genes = int(l_line[5])
                    checkv_quality = l_line[7]
                    kmer_freq = float(l_line[12])
                    if current_qual == "high":
                        checkv_pass.append(phage_id)
                    else:
                        if all([viral_genes > 0, checkv_quality != "Not-determined"]):
                            if any(
                                [
                                    checkv_quality == "Low-quality",
                                    checkv_quality == "Medium",
                                ]
                            ):
                                if kmer_freq <= 1.0:
                                    checkv_pass.append(phage_id)
                            else:
                                checkv_pass.append(phage_id)
        else:
            print("No checkV files found\n")

    return checkv_pass


def virify_parser(virify_gff, hq_list):
    with open(virify_gff, "r") as input_table, open("virify_hq.gff", "w") as output_gff:
        for line in input_table:
            line = line.rstrip()
            line_l = line.split("\t")
            # Annotation lines have exactly 9 columns
            if len(line_l) == 9:
                (
                    contig,
                    seq_source,
                    seq_type,
                    start,
                    end,
                    score,
                    strand,
                    phase,
                    attr,
                ) = line.rstrip().split("\t")

                # Saving mobilome predictions
                if seq_source == "VIRify":
                    feature_id = attr.split(";")[0].replace("ID=", "")
                else:
                    id_split = attr.split(";")[0].split("_")
                    id_split.pop(-1)
                    feature_id = "_".join(id_split).replace("ID=", "")

                if feature_id in hq_list:
                    output_gff.write(line + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="This script parse the virify output for quality control."
    )
    parser.add_argument(
        "--virify_gff",
        type=str,
        help="VIRify v2.0 output in GFF format (08-final/gff/)",
    )
    parser.add_argument(
        "--checkv_out",
        nargs="*",
        help="CheckV results generaged by VIRify (07-checkv/*quality_summary.tsv)",
    )
    args = parser.parse_args()

    # Calling functions

    checkv_pass = quality_parser(args.checkv_out)

    virify_parser(args.virify_gff, checkv_pass)


if __name__ == "__main__":
    main()
