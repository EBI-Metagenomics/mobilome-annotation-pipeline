#!/usr/bin/env python

import argparse
import os.path


##### This script generates three versions of the prokka+mobilome gff file
##### Alejandra Escobar, EMBL-EBI
##### Feb 13, 2023


def gff_parser(mobilome_prokka_gff):
    ### Parsing the gff file
    mge_coord = {}
    mge_id = {}
    mge_desc = {}
    mges = [
        "insertion_sequence",
        "integron",
        "conjugative_integron",
        "plasmid",
        "viral_sequence",
        "prophage",
        "phage_plasmid",
    ]
    flank = ["terminal_inverted_repeat_element", "attC_site", "direct_repeat_element"]
    valid_attr = ["viphog", "viphog_taxonomy", "mobileOG"]

    if os.path.isfile(mobilome_prokka_gff):
        with open(mobilome_prokka_gff, "r") as input_file, open(
            "mobilome_clean.gff", "w"
        ) as to_clean_gff, open("mobilome_extra.gff", "w") as to_extra_gff, open(
            "mobilome_nogenes.gff", "w"
        ) as to_nogenes_gff:

            for line in input_file:
                l_line = line.rstrip().split("\t")
                # Annotation lines have exactly 9 columns
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
                    if seq_type in flank:
                        to_clean_gff.write(line)
                        to_extra_gff.write(line)
                        to_nogenes_gff.write(line)
                    elif seq_type in mges:
                        to_clean_gff.write(line)
                        to_extra_gff.write(line)
                        to_nogenes_gff.write(line)
                    else:

                        if "location=mobilome" in attr:
                            to_clean_gff.write(line)
                            if any(["viphog=" in attr, "mobileOG=mobileOG" in attr]):
                                att_l = attr.split(";")
                                new_attr = [att_l[0]]
                                for element in att_l:
                                    element_key = element.split("=")[0]
                                    if element_key in valid_attr:
                                        new_attr.append(element)
                                new_attr = ";".join(new_attr)
                                to_extra_gff.write(
                                    "\t".join(
                                        [
                                            contig,
                                            seq_source,
                                            seq_type,
                                            start,
                                            end,
                                            score,
                                            strand,
                                            phase,
                                            new_attr,
                                        ]
                                    )
                                    + "\n"
                                )

                else:
                    to_clean_gff.write(line)
                    to_extra_gff.write(line)
                    to_nogenes_gff.write(line)


def main():
    parser = argparse.ArgumentParser(
        description="This script generates three versions of the prokka+mobilome gff file:"
        + "\n"
        + "1. mobilome_clean.gff: mobilome + associated CDSs"
        + "\n"
        + "3. mobilome_extra.gff: mobilome + viphogs/mobileOG annotation genes"
        + "\n"
        + "2. mobilome_nogenes.gff: mobilome only"
    )
    parser.add_argument(
        "--mobilome_prokka_gff",
        type=str,
        help="Prokka output + mobilome",
        required=True,
    )
    args = parser.parse_args()

    ## Calling functions
    gff_parser(args.mobilome_prokka_gff)


if __name__ == "__main__":
    main()
