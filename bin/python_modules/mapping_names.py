#!/usr/bin/env python

import sys
import os.path

##### This module parse the contig names mapping file
##### Alejandra Escobar, EMBL-EBI
##### October 4, 2023


def names_map(map_file):
    names_equiv = {}
    with open(map_file, "r") as input_map:
        for line in input_map:
            new_name, old_name = line.rstrip().split("\t")
            names_equiv[new_name.replace(">", "")] = old_name
    inv_names_equiv = {v: k for k, v in names_equiv.items()}

    return (names_equiv, inv_names_equiv)
