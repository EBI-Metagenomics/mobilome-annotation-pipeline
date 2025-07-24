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


def icf_parser(icf_results):
    ### Saving the ICEfinder sequences
    mge_counter = 0
    mge_data = {}
    icf_dr = {}

    ### Parsing ICEfinder2-lite results
    with open(icf_results, "r") as input_table:
        next(input_table)
        for line in input_table:
            (
                contig,  # contig_1
                ice_id,  # contig_1_ICE2
                ice_type,  # T4SS-type ICE
                ice_location,  # 529362..549932
                ice_length,  # 20571
                gc_content,  # 0.36
                direct_repeats,  # attL:529362..529421(TAGGTTGAGGGCCTAGTGGGTGAATAACCCGTGGAGGTTCAAGTCCTCTCGGCCGCATC),attR:549873..549932(TAGGTTGAGGGCCTAGTGGGTGAATAACCCGTGGAGGTTCAAGTCCTCTCGGCCGCATC)
                relaxase_type,  # MOBT
                mpf_systems,  # typeFA
                close_to_RNA,  # tRNA-Glu(528903..528974)[+]
            ) = line.rstrip().split("\t")

            mge_counter += 1
            mge_id = "icf_" + str(mge_counter)
            start = int(ice_location.split("..")[0])
            end = int(ice_location.split("..")[1])
            mge_coord = (start, end)

            ice_type = ice_type.replace(" ", "_")
            if len(direct_repeats) > 1:
                ice_type = ice_type + "_with_DRs"

                dr1, dr2 = direct_repeats.split(",")
                dr1 = dr1.replace("attL:", "").split("(")[0]
                dr1_start = int(dr1.split("..")[0])
                dr1_end = int(dr1.split("..")[1])
                dr1_coords = (dr1_start, dr1_end)

                dr2 = dr2.replace("attR:", "").split("(")[0]
                dr2_start = int(dr2.split("..")[0])
                dr2_end = int(dr2.split("..")[1])
                dr2_coords = (dr2_start, dr2_end)

                icf_dr[mge_id] = (dr1_coords, dr2_coords)

            description = (
                "mobile_element_type="
                + ice_type
                + ";relaxase_type="
                + relaxase_type
                + ";mpf_systems="
                + mpf_systems
            )
            value = (contig, description, mge_coord)
            mge_data[mge_id] = value

    return (mge_data, icf_dr)
