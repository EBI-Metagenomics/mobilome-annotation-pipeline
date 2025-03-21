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


def isescan_parser(mge_data, iss_results):
    ### Parsing ISEScan outputs
    mge_counter = 0
    itr_sites = {}
    if os.stat(iss_results).st_size > 0:
        with open(iss_results, "r") as input_table:
            next(input_table)
            for line in input_table:
                mge_counter += 1
                mge_id = "iss_" + str(mge_counter)
                (
                    contig,
                    family,
                    cluster,
                    isBegin,
                    isEnd,
                    isLen,
                    ncopy4is,
                    start1,
                    end1,
                    start2,
                    end2,
                    score,
                    irId,
                    irLen,
                    nGaps,
                    orfBegin,
                    orfEnd,
                    strand,
                    orfLen,
                    e_value,
                    e_value4copy,
                    iss_type,
                    ov,
                    tir,
                ) = line.rstrip().split("\t")

                ## Keeping complete predictions only
                if iss_type == "c":
                    if int(irLen) == 0:
                        description = "without_TIR"
                    else:
                        description = "with_TIR"
                        ir_1 = (start1, end1)
                        ir_2 = (start2, end2)
                        itr_sites[mge_id] = (ir_1, ir_2)

                    description = "mobile_element_type=" + family + "_" + description
                    coord = (int(isBegin), int(isEnd))
                    value = (contig, description, coord)
                    mge_data[mge_id] = value

    return (mge_data, itr_sites)
