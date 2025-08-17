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


def mge_data_parser(mge_data):
    integrated_dic = {}
    types_list = [
        "IME",
        "ICE",
        "AICE",
        "integron",
        "TIR",
    ]
    for mge in mge_data:
        contig, description, coord = mge_data[mge]
        if "prophage" in description:
            if contig in integrated_dic:
                integrated_dic[contig].append(coord)
            else:
                integrated_dic[contig] = [coord]
        for mge_type in types_list:
            if mge_type in description:
                if contig in integrated_dic:
                    integrated_dic[contig].append(coord)
                else:
                    integrated_dic[contig] = [coord]

    return integrated_dic


def calculate_intersection(range1, range2):
    """
    Calculate the intersection between two ranges.
    
    :param range1: First range object
    :type range1: range
    :param range2: Second range object  
    :type range2: range
    :return: Length of intersection between the two ranges
    :rtype: int
    """
    return len(set(range1) & set(range2))


def outliers_parser(
    comp_bed, mge_data, rnas_coord, rna_cov_threshold=1.0, co_cov_threshold=0.75
):
    mge_counter = 0
    co_repeats = {}

    if len(comp_bed) == 0:
        return (mge_data, co_repeats)

    # Saving integrated mges in the mge_data
    integrated_dic = mge_data_parser(mge_data)

    # Parsing the compositional outliers (CO) file
    # There is one line per flanking repeat. That means two lines per CO
    # contig_1	62056	62067	contig_1_001:62056-67629|direct_1|GGCGCGCGGCC|gc:3.67|kmer:2.63|combined:3.67|conf:1.000
    # contig_1	67618	67629	contig_1_001:62056-67629|direct_2|GGCGCGCGGCC|gc:3.67|kmer:2.63|combined:3.67|conf:1.000

    with open(comp_bed, "r") as input_table:
        lines = (
            line.rstrip().split("\t")
            for line in input_table
            if not line.startswith("#")
        )
        for line1, line2 in zip(lines, lines):
            # Collecting repeats info
            contig, repeat1_start, repeat1_end, co_info1 = line1
            coord1 = (int(repeat1_start), int(repeat1_end))
            co_id, repeat_type1, repeat_seq1, gc, kmer, score, conf = co_info1.split(
                "|"
            )
            metainfo1 = [
                repeat_type1.replace("direct_1", "DR").replace("inverted_1", "IR"),
                repeat_seq1,
                coord1,
            ]

            contig, repeat2_start, repeat2_end, co_info2 = line2
            coord2 = (int(repeat2_start), int(repeat2_end))
            co_id, repeat_type2, repeat_seq2, gc, kmer, score, conf = co_info2.split(
                "|"
            )
            metainfo2 = [
                repeat_type2.replace("direct_2", "DR").replace("inverted_2", "IR"),
                repeat_seq2,
                coord2,
            ]

            # Collecting CO location info
            co_contig = co_id.split("_")
            co_contig.pop(-1)
            co_contig = "_".join(co_contig)
            co_start = int(co_id.split(":")[1].split("-")[0])
            co_end = int(co_id.split(":")[1].split("-")[1])
            co_coord = (co_start, co_end)
            co_len = co_end - co_start
            co_range = range(co_start, co_end + 1)
            co_attributes = ";".join(
                [
                    "mobile_element_type=compositional_outlier",
                    gc.replace(":", "="),
                    kmer.replace(":", "="),
                    conf.replace(":", "="),
                ]
            )


            new_co_id = (
                co_contig
                + "|compositional_outlier:"
                + str(co_start)
                + ":"
                + str(co_end)
            )

            # Removing CO overlapping rnas
            # RNA flag turn to 1 when one RNA is fully overlapping a CO
            rna_flag = 0
            if co_contig in rnas_coord:
                for rna_coord_pair in rnas_coord[co_contig]:
                    rna_start = rna_coord_pair[0]
                    rna_end = rna_coord_pair[1]
                    rna_len = rna_end - rna_start
                    rna_range = range(rna_start, rna_end + 1)
                    intersection = calculate_intersection(rna_range, co_range)
                    if intersection > 0:
                        rna_cov = float(intersection) / float(rna_len)
                        if rna_cov == rna_cov_threshold:
                            rna_flag += 1
                            print(co_id + " removed due to RNAs genes")

            # If not multiple RNAs in CO then evaluate overlapping with integrated MGEs
            # MGE flag turn to 1 when an MGE overlaps >= 75% of CO length
            if rna_flag <= 2:
                mge_flag = False
                if co_contig in integrated_dic:
                    for mge_coord_pair in integrated_dic[co_contig]:
                        mge_start = mge_coord_pair[0]
                        mge_end = mge_coord_pair[1]
                        mge_len = mge_end - mge_start
                        mge_range = range(mge_start, mge_end + 1)
                        intersection = calculate_intersection(mge_range, co_range)
                        if intersection > 0:
                            co_cov = float(intersection) / float(co_len)
                            if co_cov >= co_cov_threshold:
                                mge_flag = True
                                print(
                                    co_id
                                    + " removed due to overlapping with MGE with coverage "
                                    + str(co_cov)
                                )

                # Add the CO to mge_data dict if it is not redundant with MGE annotations and not contain RNAs
                if not mge_flag:
                    mge_counter += 1
                    mge_id = "co_" + str(mge_counter)
                    composite_val = (co_contig, co_attributes, co_coord)
                    mge_data[mge_id] = composite_val
                    # Save the corresponding flanking repeats
                    co_repeats[mge_id] = (metainfo1, metainfo2)

    return (mge_data, co_repeats)

