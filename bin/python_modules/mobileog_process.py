#!/usr/bin/env python

import sys
import os.path

##### This module parse the mobileOG annotation
##### Alejandra Escobar, EMBL-EBI
##### October 13, 2023


def mobileog_parser(mog_results):
    ### Parsing the mobileOG annotation file
    mog_annot = {}
    if os.stat(mog_results).st_size > 0:
        with open(mog_results) as input_table:
            for line in input_table:
                line_l = line.rstrip().split("\t")
                target_id = line_l[0]
                query_id = line_l[1].split(" ")[0]
                (
                    mog_id,
                    gene_name,
                    best_hit_id,
                    major,
                    minor,
                    db,
                    evidence,
                ) = target_id.split("|")
                function = mog_id + "|" + major + "|" + minor

                if query_id not in mog_annot:
                    mog_annot[query_id] = function.replace(",", "/")

    return mog_annot
