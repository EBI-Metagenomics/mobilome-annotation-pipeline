#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2025 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import os.path
import logging

logging.basicConfig(level=logging.INFO)


def gff_parser(mobilome_prokka_gff: str, prefix: str) -> None:
    """Parsing and filtering GFF file containing mobile genetic elements annotations.
    
    :param mobilome_prokka_gff: Path to the input GFF file containing Prokka annotations for mobile elements
    :type mobilome_prokka_gff: str
    :param prefix: Prefix to use for output filenames
    :type prefix: str
    :return: None, creates three output files:
             - {prefix}_mobilome_clean.gff
             - {prefix}_mobilome_extra.gff 
             - {prefix}_mobilome_nogenes.gff
    :rtype: None
    
    The function processes mobile genetic elements including:
    - insertion sequences
    - integrons
    - conjugative integrons 
    - plasmids
    - viral sequences
    - prophages
    - phage plasmids
    
    And flanking elements:
    - terminal inverted repeat elements
    - attC sites
    - direct repeat elements
    """
    mges = [
        "insertion_sequence",
        "integron",
        "conjugative_integron",
        "plasmid",
        "viral_sequence",
        "prophage",
        "phage_plasmid",
        "compositional_outlier",
    ]
    flank = ["attC_site", "direct_repeat_element", "inverted_repeat_element"]
    valid_attr = ["viphog", "viphog_taxonomy", "mobileOG"]

    if not os.path.isfile(mobilome_prokka_gff):
        logging.debug(f"{mobilome_prokka_gff} not found. Skipping")
        return

    with open(mobilome_prokka_gff, "r") as input_file, open(
        f"{prefix}_mobilome_clean.gff", "w"
    ) as to_clean_gff, open(f"{prefix}_mobilome_extra.gff", "w") as to_extra_gff, open(
        f"{prefix}_mobilome_nogenes.gff", "w"
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
        "-m",
        "--mobilome_prokka_gff",
        type=str,
        help="Prokka output + mobilome",
        required=True,
    )
    parser.add_argument(
        "-p",
        "--prefix",
        type=str,
        help="The output prefix",
        required=True
    )
    args = parser.parse_args()

    gff_parser(args.mobilome_prokka_gff, args.prefix)


if __name__ == "__main__":
    main()
