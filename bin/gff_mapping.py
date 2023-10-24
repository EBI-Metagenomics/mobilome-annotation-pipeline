#!/usr/bin/env python

import argparse
import sys
import os.path
import glob


##### This script add the mobilome prediction to the user GFF file
##### Alejandra Escobar, EMBL-EBI
##### October 19, 2023


def mobilome_parser(mobilome_prokka_extra):
    # Parsing the mobilome prediction
    proteins_annot, mobilome_annot = {}, {}

    source_tools = [
        "ICEfinder",
        "IntegronFinder",
        "ISEScan",
        "PaliDIS",
        "geNomad",
        "VIRify",
        "geNomad",
        "geNomad_VIRify",
    ]
    extra_annot = [
        "viphog",
        "viphog_taxonomy",
        "mobileOG",
    ]

    if os.stat(mobilome_prokka_extra).st_size > 0:
        with open(mobilome_prokka_extra, "r") as input_table:
            for line in input_table:
                l_line = line.rstrip().split("\t")
                # Annotation lines have exactly 9 columns
                if len(l_line) == 9:
                    contig = l_line[0]
                    annot_source = l_line[1]
                    if annot_source in source_tools:
                        if contig in mobilome_annot:
                            mobilome_annot[contig].append(line.rstrip())
                        else:
                            mobilome_annot[contig] = [line.rstrip()]
                    else:
                        start = l_line[3]
                        end = l_line[4]
                        strand = l_line[6]
                        composite_key = (contig, start, end, strand)
                        attrib = l_line[8]
                        extra_list = []
                        for attr in attrib.split(";"):
                            att_key = attr.split("=")[0]
                            att_val = attr.split("=")[1]
                            if att_key in extra_annot:
                                extra_list.append(attr)
                        if len(extra_list) > 0:
                            extra_val = ";".join(extra_list)
                            proteins_annot[composite_key] = extra_val

    return (proteins_annot, mobilome_annot)


# Adding the mobilome predictions to the user file
def gff_updater(user_gff, proteins_annot, mobilome_annot):
    used_contigs = []
    if os.stat(user_gff).st_size > 0:
        with open(user_gff, "r") as input_table, open(
            "user_mobilome_extra.gff", "w"
        ) as output_table:
            for line in input_table:
                l_line = line.rstrip().split("\t")
                # Annotation lines have exactly 9 columns
                if len(l_line) == 9:
                    contig = l_line[0]
                    start = l_line[3]
                    end = l_line[4]
                    strand = l_line[6]
                    if not contig in used_contigs:
                        used_contigs.append(contig)
                        if contig in mobilome_annot:
                            for mge in mobilome_annot[contig]:
                                output_table.write(mge + "\n")

                    composite_val = (contig, start, end, strand)
                    if composite_val in proteins_annot:
                        extra_annot = proteins_annot[composite_val]
                        output_table.write(line.rstrip() + ";" + extra_annot + "\n")

                else:
                    output_table.write(line.rstrip() + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="This script add the extra annotations to the user GFF file"
    )
    parser.add_argument(
        "--mobilome_prokka_extra",
        type=str,
        help="Mobilome prediction on prokka gff file containing the mobilome proteins with extra annotation only",
        required=True,
    )
    parser.add_argument(
        "--user_gff",
        type=str,
        help="User GFF file",
        required=True,
    )
    args = parser.parse_args()

    ## Calling functions
    # Storing the mobilome predictions
    (proteins_annot, mobilome_annot) = mobilome_parser(args.mobilome_prokka_extra)

    # Adding the mobilome predictions to the user file
    gff_updater(args.user_gff, proteins_annot, mobilome_annot)


if __name__ == "__main__":
    main()
