#!/usr/bin/env python

import sys
import os.path


##### This module parse crispr annotation
##### Alejandra Escobar, EMBL-EBI
##### October 13, 2023


def minced_parser(pkka_gff):
    pkka_minced = {}
    # Reading PROKKA crispr annotations
    with open(pkka_gff, "r") as input_table:
        for line in input_table:
            l_line = line.rstrip().split("\t")
            if len(l_line) == 9:
                pred_source = l_line[1].split(":")[0]
                if pred_source == "minced":
                    contig = l_line[0]
                    start = int(l_line[3])
                    end = int(l_line[4])
                    coordinates = (start, end)
                    if contig in pkka_minced:
                        pkka_minced[contig].append(coordinates)
                    else:
                        pkka_minced[contig] = [coordinates]

    return pkka_minced


def crispr_parser(crispr_out, pkka_gff):
    crispr_annot = {}
    # Reading PROKKA crispr annotations
    pkka_minced = minced_parser(pkka_gff)
    counter = 0
    with open(crispr_out, "r") as file_in:
        next(file_in)
        for line in file_in:
            # Ignore empty lines
            if len(line.strip()) == 0:
                continue
            line_l = line.rstrip().split("\t")
            conf_level = line_l[-1]
            if int(conf_level) > 2:
                contig = line_l[1]
                start = int(line_l[5])
                end = int(line_l[6])
                orientation = line_l[8].lower()
                concensus_repeat = line_l[10]
                spacers_nb = int(line_l[14])

                description = ";".join(
                    [
                        "note=CRISPR with " + str(spacers_nb) + " repeat units",
                        "rpt_family=CRISPR",
                        "rpt_type=" + orientation,
                        "rpt_unit_seq=" + concensus_repeat,
                    ]
                )
                coord = (start, end)
                value = (description, coord)

                if contig in pkka_minced:
                    flag = 0
                    curr_range = range(start, end + 1)
                    curr_len = end - start
                    for pkka_prediction in pkka_minced[contig]:
                        pkka_s = pkka_prediction[0]
                        pkka_e = pkka_prediction[1]
                        pkka_range = range(pkka_s, pkka_e + 1)
                        pkka_len = pkka_e - pkka_s
                        intersection = len(list(set(curr_range) & set(pkka_range)))
                        if intersection > 0:
                            p_cov = float(intersection) / float(pkka_len)
                            if p_cov >= 0.1:
                                flag = 1

                    if flag == 0:
                        if contig in crispr_annot:
                            crispr_annot[contig].append(value)
                        else:
                            crispr_annot[contig] = [value]

                else:
                    if contig in crispr_annot:
                        crispr_annot[contig].append(value)
                    else:
                        crispr_annot[contig] = [value]

    return crispr_annot
