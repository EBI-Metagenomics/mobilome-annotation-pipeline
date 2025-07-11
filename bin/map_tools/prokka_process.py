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

import os.path
from collections import defaultdict


def prokka_parser(prokka_gff: str):
    prots_coord = {}

    contig_prots = defaultdict(list)
    rnas_coord = defaultdict(list)

    if os.stat(prokka_gff).st_size == 0:
        return

    with open(prokka_gff, "r") as input_table:
        for line in input_table:
            l_line = line.rstrip().split("\t")
            # Annotation lines have exactly 9 columns
            if len(l_line) == 9:
                contig = l_line[0]
                source = l_line[1]
                feature_type = l_line[2]
                start = int(l_line[3])
                end = int(l_line[4])
                attrib = l_line[8]
                prot_id = attrib.split(";")[0].replace("ID=", "")
                value = (start, end)
                prots_coord[prot_id] = value
                contig_prots[contig].append(prot_id)

                if "RNA" in feature_type:
                    rnas_coord[contig].append(value)

    return contig_prots, prots_coord, rnas_coord
