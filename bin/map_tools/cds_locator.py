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


def quality_filter(
    mge_data, mge_proteins, contigs_elements, output_file
):
    # Removing predictions of len<500 and with no CDS
    no_cds = []
    len_500 = []
    for element in mge_data:
        element_len = int(mge_data[element][2][1]) - int(mge_data[element][2][0])
        if element_len < 500:
            len_500.append(element)
        elif len(mge_proteins[element]) == 0:
            no_cds.append(element)

    with open(output_file, "a") as to_discard:
        for element in no_cds:
            to_discard.write(element + "\t" + mge_data[element][1] + "\tno_cds\n")
            del mge_data[element]
            for contig in contigs_elements:
                if element in contigs_elements[contig]:
                    contigs_elements[contig].remove(element)

        for element in len_500:
            to_discard.write(element + "\t" + mge_data[element][1] + "\tmge<500bp\n")
            del mge_data[element]
            for contig in contigs_elements:
                if element in contigs_elements[contig]:
                    contigs_elements[contig].remove(element)

    return (mge_data, contigs_elements)


def _parse_cds_loc(cds_loc: str, contig_prots: dict, prots_coord: dict):
    if os.stat(cds_loc).st_size == 0:
        return

    with open(cds_loc, "r") as input_table:
        for line in input_table:
            l_line = line.rstrip().split("\t")
            # Annotation lines have exactly 9 columns
            if len(l_line) == 9:
                contig = l_line[0]
                prot_source = l_line[1]
                start = int(l_line[3])
                end = int(l_line[4])
                attrib = l_line[8]
                prot_id = attrib.split(";")[0].replace("ID=", "")
                value = (start, end)
                prots_coord[prot_id] = value
                if contig in contig_prots:
                    contig_prots[contig].append(prot_id)
                else:
                    contig_prots[contig] = [prot_id]


def location_parser(cds_loc, mge_data, output_file="discarded_mge.txt"):
    ### Saving the CDS' coordinates and contigs location
    contig_prots = {}
    prots_coord = {}
    COV_THRESHOLD = 0.9

    _parse_cds_loc(cds_loc, contig_prots, prots_coord)

    ### Finding the CDS encoded in the mobilome
    mge_proteins = {}
    protein_mge = {}
    contigs_elements = {}
    for element in mge_data:
        mge_proteins[element] = []
        contig = mge_data[element][0]
        mge_start = mge_data[element][2][0]
        mge_end = mge_data[element][2][1]
        mge_range = range(mge_start, mge_end + 1)
        mge_len = mge_end - mge_start

        if contig in contigs_elements:
            contigs_elements[contig].append(element)
        else:
            contigs_elements[contig] = [element]

        if contig in contig_prots:
            for protein in contig_prots[contig]:
                prot_start = prots_coord[protein][0]
                prot_end = prots_coord[protein][1]
                prot_range = range(prot_start, prot_end + 1)
                prot_len = prot_end - prot_start
                intersection = len(list(set(mge_range) & set(prot_range)))
                if intersection > 0:
                    mge_cov = float(intersection) / float(mge_len)
                    prot_cov = float(intersection) / float(prot_len)
                    if any([mge_cov > COV_THRESHOLD, prot_cov > COV_THRESHOLD]):
                        mge_proteins[element].append(protein)
                        protein_mge[protein] = element

    ## Quality control of small and empty MGEs
    (mge_data, contigs_elements) = quality_filter(
        mge_data, mge_proteins, contigs_elements, output_file
    )

    return (mge_data, contigs_elements, protein_mge)
