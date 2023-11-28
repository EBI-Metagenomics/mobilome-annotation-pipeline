#!/usr/bin/env python

import argparse
import os.path

##### This script integrates the results of amrfinderplus with the mobilome
##### Alejandra Escobar, EMBL-EBI
##### May 31, 2023


def arg_parser(amr_out):
    ### Saving amr predictions
    amr_data = {}
    if os.path.exists(amr_out):
        if os.stat(amr_out).st_size > 0:
            with open(amr_out, "r") as input_csv:
                next(input_csv)
                for line in input_csv:
                    line_l = line.rstrip().split("\t")
                    protein_id = line_l[0]
                    e_type = line_l[4]
                    e_stype = line_l[5]
                    if e_type != e_stype:
                        e_type = e_type + "|" + e_stype
                    e_class = line_l[6]
                    e_sclass = line_l[7]
                    if e_class != e_sclass:
                        e_class = e_class + "|" + e_sclass
                    decriptors = (e_stype, e_class)
                    amr_data[protein_id] = decriptors

    return amr_data


def mob_parser(mobilome):
    ### Saving the proteins in the mobilome
    mob_prots = []
    with open(mobilome, "r") as input_gff:
        for line in input_gff:
            l_line = line.rstrip().split("\t")
            ## Annotation lines have exactly 9 columns
            if len(l_line) == 9:
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

                if seq_type == "CDS":
                    att_fields = attr.split(";")
                    protein_id = att_fields[0].replace("ID=", "")
                    for attribute in att_fields:
                        if attribute == "location=mobilome":
                            mob_prots.append(protein_id)

    return mob_prots


def location_parser(amr_data, mob_prots):
    with open("amr_location.txt", "w") as to_output:
        to_output.write(
            "\t".join(["Gene_id", "pred_type", "pred_class", "location"]) + "\n"
        )
        for gene in amr_data:
            prediction_type = amr_data[gene][0]
            prediction_description = amr_data[gene][1]
            if gene in mob_prots:
                to_output.write(
                    "\t".join(
                        [gene, prediction_type, prediction_description, "mobilome"]
                    )
                    + "\n"
                )
            else:
                to_output.write(
                    "\t".join(
                        [gene, prediction_type, prediction_description, "chromosome"]
                    )
                    + "\n"
                )


def main():
    parser = argparse.ArgumentParser(
        description="This script integrates the results of amrfinderplus with the mobilome. Please provide the relevant input files"
    )
    parser.add_argument(
        "--mobilome",
        type=str,
        help="Clean version of the output of the mobilome annotation pipeline in GFF format (mobilome_prokka.gff)",
    )
    parser.add_argument(
        "--amr_out",
        type=str,
        help="ARG prediction output (amrfinderplus.tsv)",
        required=True,
    )
    args = parser.parse_args()

    ### Calling functions
    amr_data = arg_parser(args.amr_out)
    mob_prots = mob_parser(args.mobilome)
    location_parser(amr_data, mob_prots)


if __name__ == "__main__":
    main()
