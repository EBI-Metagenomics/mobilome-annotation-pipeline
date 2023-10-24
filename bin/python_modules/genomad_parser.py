#!/usr/bin/env python

import sys
import os.path


##### This module parse the geNomad output files
##### Alejandra Escobar, EMBL-EBI
##### October 4, 2023


def genomad_viral(geno_out, mge_data):
    mge_counter = 0
    if os.stat(geno_out).st_size > 0:
        with open(geno_out, "r") as input_table:
            next(input_table)
            for line in input_table:
                line_l = line.rstrip().split("\t")
                pred_id = line_l[0]
                score = float(line_l[6])
                if score > 0.8:
                    mge_counter += 1
                    mge_id = "vir1_" + str(mge_counter)
                    taxonomy = line_l[10].replace(";", "%3B")
                    if "provirus" in pred_id:
                        contig = line_l[0].split("|")[0]
                        description = (
                            "mobile_element_type=prophage;" + "taxonomy=" + taxonomy
                        )
                        start = int(line_l[3].split("-")[0])
                        end = int(line_l[3].split("-")[1])
                    else:
                        contig = line_l[0]
                        description = (
                            "mobile_element_type=viral_sequence;"
                            + "taxonomy="
                            + taxonomy
                        )
                        start = 1
                        end = int(line_l[1])

                    coord = (start, end)
                    value = (contig, description, coord)
                    mge_data[mge_id] = value

    return mge_data


def plasmids_parser(geno_plas, mge_data):
    mge_counter = 0
    if os.stat(geno_plas).st_size > 0:
        with open(geno_plas, "r") as input_table:
            next(input_table)
            for line in input_table:
                line_l = line.rstrip().split("\t")
                score = float(line_l[5])
                if score > 0.8:
                    contig = line_l[0]
                    mge_counter += 1
                    plasmid_id = "plas_" + str(mge_counter)
                    start = 1
                    end = int(line_l[1])
                    coord = (start, end)
                    description = "mobile_element_type=plasmid"
                    composite_value = (contig, description, coord)
                    mge_data[plasmid_id] = composite_value

    return mge_data
