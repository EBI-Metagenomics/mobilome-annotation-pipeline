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
