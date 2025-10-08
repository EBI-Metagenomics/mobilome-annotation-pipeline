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

import gzip

def line_builder(id_to_print, source, seq_type, start, end, feat_id, description):
    attributes = feat_id + ";" + description
    to_print = "\t".join(
        [id_to_print, source, seq_type, str(start), str(end), ".", ".", ".", attributes]
    )
    return to_print


def id_collector(contig, start, end, mge_type):
    feature_id = "ID=" + contig + "|" + mge_type + "-" + str(start) + ":" + str(end)
    return feature_id


def gff_writer(
    names_equiv,
    contigs_elements,
    itr_sites,
    icf_dr,
    attC_site,
    co_repeats,
    mge_data,
    output_gff,
):
    ### Writing the mobilome outputs
    source_tools = {
        "icf": "ICEfinder",
        "inf": "IntegronFinder",
        "iss": "ISEScan",
        "vir1": "geNomad",
        "vir2": "VIRify",
        "plas": "geNomad",
        "phpl": "geNomad_VIRify",
        "co": "MAP",
    }

    with gzip.open(output_gff + ".gz", "wt") as to_gff:
        # Writing header
        to_gff.write('##gff-version 3' + '\n')

        for contig in contigs_elements:
            id_to_print = names_equiv[contig]
            for element in contigs_elements[contig]:
                source = source_tools[element.split("_")[0]]
                e_start = mge_data[element][2][0]
                e_end = mge_data[element][2][1]
                e_desc = mge_data[element][1]
                if "iss" in element:
                    seq_type = "insertion_sequence"
                    if element in itr_sites:
                        f_seq_type = "inverted_repeat_element"
                        f_start = itr_sites[element][0][0]
                        f_end = itr_sites[element][0][1]
                        feature_id = id_collector(
                            id_to_print, f_start, f_end, f_seq_type
                        )
                        gff_line = line_builder(
                            id_to_print,
                            source,
                            f_seq_type,
                            f_start,
                            f_end,
                            feature_id,
                            "flanking_site=IR_1",
                        )
                        to_gff.write(gff_line + "\n")

                        f_start = itr_sites[element][1][0]
                        f_end = itr_sites[element][1][1]
                        feature_id = id_collector(
                            id_to_print, f_start, f_end, f_seq_type
                        )
                        gff_line = line_builder(
                            id_to_print,
                            source,
                            f_seq_type,
                            f_start,
                            f_end,
                            feature_id,
                            "flanking_site=IR_2",
                        )   
                        to_gff.write(gff_line + "\n")

                elif "inf" in element:
                    seq_type = "integron"
                    if element in attC_site:
                        f_seq_type = "attC_site"
                        for attc_element in attC_site[element]:
                            f_start = attc_element[0]
                            f_end = attc_element[1]
                            feature_id = id_collector(
                            id_to_print, f_start, f_end, f_seq_type
                            )
                            gff_line = line_builder(
                                id_to_print,
                                source,
                                f_seq_type,
                                f_start,
                                f_end,
                                feature_id,
                                "recombination_site=attC",
                            )
                            to_gff.write(gff_line + "\n")

                elif "icf" in element:
                    if "IME" in mge_data[element][1]:
                        seq_type = "integron"
                    elif "ICE" in mge_data[element][1]:
                        seq_type = "conjugative_integron"
                    if element in icf_dr:
                        f_seq_type = "direct_repeat_element"
                        f_start = icf_dr[element][0][0]
                        f_end = icf_dr[element][0][1]
                        feature_id = id_collector(
                            id_to_print, f_start, f_end, f_seq_type
                        )   
                        gff_line = line_builder(
                            id_to_print,
                            source,
                            f_seq_type,
                            f_start,
                            f_end,
                            feature_id,
                            "flanking_site=DR_1",
                        )
                        to_gff.write(gff_line + "\n")

                        f_start = icf_dr[element][1][0]
                        f_end = icf_dr[element][1][1]
                        feature_id = id_collector(
                            id_to_print, f_start, f_end, f_seq_type
                        )
                        gff_line = line_builder(
                            id_to_print,
                            source,
                            f_seq_type,
                            f_start,
                            f_end,
                            feature_id,
                            "flanking_site=DR_2",
                        )
                        to_gff.write(gff_line + "\n")

                elif "vir" in element:
                    # Extracting the data type from the attributes line in VIRify and geNomad predictions:
                    # mobile_element_type=phage_linear;checkv_provirus=No;checkv_quality=Low-quality;checkv_miuvig_quality=Genome-fragment;checkv_kmer_freq=1.0;checkv_viral_genes=10;taxonomy=Prymnesiovirus%3BPhycodnaviridae%3BAlgavirales

                    vir_attributes = mge_data[element][1].split(";")
                    for vir_att in vir_attributes:
                        vir_key, vir_value = vir_att.split("=")
                        if vir_key == "mobile_element_type":
                            mobile_element_type = vir_value

                    if "prophage" in mobile_element_type:
                        seq_type = "prophage"
                    else:
                        seq_type = "viral_sequence"

                elif "plas" in element:                    
                    seq_type = "plasmid"

                elif "phpl" in element:
                    seq_type = "phage_plasmid"

                elif "co" in element:
                    seq_type = "compositional_outlier"
                    if element in co_repeats:
                        # Parsing CO flanking repreats info:
                        # metainfo1 = [ repeat_type1.replace('direct_1','DR').replace('inverted_1','IR'), repeat_seq1, coord1 ]
                        # metainfo2 = [ repeat_type2.replace('direct_2','DR').replace('inverted_2','IR'), repeat_seq2, coord2 ]
                        # co_repeats[mge_id] = (metainfo1, metainfo2)
                        element_counter = 0
                        for repeat_element in co_repeats[element]:
                            element_counter += 1
                            f_seq_prefix, repeat_seq, coords_tuple = list(
                                repeat_element
                            )
                            flanking_site = (
                                f_seq_prefix + "_" + str(element_counter)
                            )
                            f_seq_type = (
                                f_seq_prefix.replace(
                                    "DR", "direct"
                                ).replace("IR", "inverted")
                                + "_repeat_element"
                            )
                            extra_info = ";".join(
                                [
                                    "repeat_sequence=" + repeat_seq,
                                    "flanking_site=" + flanking_site,
                                ]
                            )

                            f_start = coords_tuple[0]
                            f_end = coords_tuple[1]
                            feature_id = id_collector(
                                id_to_print, f_start, f_end, f_seq_type
                            )
                            gff_line = line_builder(
                                id_to_print,
                                source,
                                f_seq_type,
                                f_start,
                                f_end,
                                feature_id,
                                extra_info,
                            )
                            to_gff.write(gff_line + "\n")

                feature_id = id_collector(
                    id_to_print, e_start, e_end, seq_type
                )
                gff_line = line_builder(
                    id_to_print,
                    source,
                    seq_type,
                    e_start,
                    e_end,
                    feature_id,
                    e_desc,
                )
                to_gff.write(gff_line + "\n")

