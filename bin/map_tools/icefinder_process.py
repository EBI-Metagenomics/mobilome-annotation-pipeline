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


def icf_dr_control(icf_dr, mge_data):
    # Correcting wrong coordinates
    icf_dr_corr = {}
    for mge in icf_dr:
        start_1 = int(icf_dr[mge][0][0])
        end_1 = int(icf_dr[mge][0][1])
        len_1 = end_1 - start_1

        start_2 = int(icf_dr[mge][1][0])
        end_2 = int(icf_dr[mge][1][1])
        len_2 = end_2 - start_2
        if start_2 > end_2:
            mge_start = mge_data[mge][2][0]
            mge_end = mge_data[mge][2][1]
            new_start_2 = mge_end - len_1
            print("Corrected  coordinates:", start_2, "   ->   ", new_start_2)
            dr_1 = (start_1, end_1)
            new_dr_2 = (new_start_2, mge_end)
            icf_dr[mge] = (dr_1, new_dr_2)

    return icf_dr


def icf_parser(icf_dr_file, icf_results):
    ### Saving the ICEfinder sequences
    mge_counter = 0
    mge_data = {}
    icf_dr = {}

    ### Parsing ICEfinder result
    with open(icf_results, "r") as input_table:
        for line in input_table:
            (
                icefinder_Job_id,
                strain,
                genome_len,
                icefinder_output_name,
                description,
                coordinate,
                length,
                oriT,
                gc,
                genome_GC,
                delta_GC,
                arg,
                vf,
            ) = line.rstrip().split("\t")
            contig = icefinder_Job_id
            description = (
                description.replace(" ", "_")
                .replace("Putative_", "")
                .replace("_AICE", "AICE")
                .replace(":", "")
            )

            description = "mobile_element_type=" + description

            if "conjugative_region" not in description:
                mge_counter += 1
                mge_id = "icf_" + str(mge_counter)
                start = int(coordinate.split("..")[0])
                end = int(coordinate.split("..")[1])
                coord = (start, end)
                value = (contig, description, coord)
                mge_data[mge_id] = value

                seq_id = contig + "|" + coordinate

                composite_id = contig + "|" + icefinder_output_name
                if composite_id in icf_dr:
                    icf_dr[mge_id] = icf_dr.pop(composite_id)

    icf_dr = icf_dr_control(icf_dr, mge_data)

    return (mge_data, icf_dr)
