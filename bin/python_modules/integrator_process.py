#!/usr/bin/env python

import sys
import os.path

##### This module integrates the mobilome predictions
##### Alejandra Escobar, EMBL-EBI
##### October 18, 2023


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
    cds_gff,
    names_equiv,
    contigs_elements,
    itr_sites,
    icf_dr,
    attC_site,
    mge_data,
    proteins_mge,
    mog_annot,
    crispr_annot,
    virify_prots,
):

    ### Writing the mobilome outputs
    source_tools = {
        "icf": "ICEfinder",
        "inf": "IntegronFinder",
        "iss": "ISEScan",
        "pal": "PaliDIS",
        "vir1": "geNomad",
        "vir2": "VIRify",
        "plas": "geNomad",
        "phpl": "geNomad_VIRify",
    }

    output_gff = "mobilome_prokka.gff"
    used_contigs = []
    used_mges = {}
    with open(cds_gff, "r") as input_table, open(output_gff, "w") as to_gff:
        for line in input_table:
            l_line = line.rstrip().split("\t")
            if len(l_line) == 9:
                contig = l_line[0]
                id_to_print = names_equiv[contig]

                # Writing the MGEs found in the current contig
                if not contig in used_contigs:
                    used_contigs.append(contig)

                    if contig in crispr_annot:
                        for prediction in crispr_annot[contig]:
                            pred_start = prediction[1][0]
                            pred_end = prediction[1][1]
                            pred_attributes = prediction[0]
                            to_print = "\t".join(
                                [
                                    id_to_print,
                                    "crisprcasfinder",
                                    "repeat_region",
                                    str(pred_start),
                                    str(pred_end),
                                    ".",
                                    ".",
                                    ".",
                                    pred_attributes,
                                ]
                            )

                    if contig in contigs_elements:
                        for element in contigs_elements[contig]:
                            source = source_tools[element.split("_")[0]]
                            e_start = mge_data[element][2][0]
                            e_end = mge_data[element][2][1]
                            e_desc = mge_data[element][1]
                            if any(["iss" in element, "pal" in element]):
                                seq_type = "insertion_sequence"
                                if element in itr_sites:
                                    f_seq_type = "terminal_inverted_repeat_element"
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
                                        "flanking_site=TIR_1",
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
                                        "flanking_site=TIR_2",
                                    )
                                    to_gff.write(gff_line + "\n")

                            elif "inf" in element:
                                seq_type = "integron"
                                if element in attC_site:
                                    f_seq_type = "attC_site"
                                    f_start = attC_site[element][0]
                                    f_end = attC_site[element][1]
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
                                if "prophage" in mge_data[element][1]:
                                    seq_type = "prophage"
                                elif "viral_sequence" in mge_data[element][1]:
                                    seq_type = "viral_sequence"

                            elif "plas" in element:
                                seq_type = "plasmid"

                            elif "phpl" in element:
                                seq_type = "phage_plasmid"

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

                # Writing the rest of the annotations
                attrib = l_line[8]
                prot_id = attrib.split(";")[0].replace("ID=", "")

                if prot_id in proteins_mge:
                    attrib = attrib + ";location=mobilome"
                else:
                    attrib = attrib + ";location=chromosome"

                if prot_id in mog_annot:
                    function = mog_annot[prot_id].replace(" ", "_")
                    attrib = attrib + ";mobileOG=" + function

                start = int(l_line[3])
                end = int(l_line[4])
                prot_location = (contig, start, end)
                if prot_location in virify_prots:
                    function = virify_prots[prot_location]
                    attrib = attrib + ";" + function

                pkka_contig = l_line.pop(0)
                contig = names_equiv[pkka_contig]
                l_line.pop(-1)
                gff_line = [contig] + l_line + [attrib]
                to_gff.write("\t".join(gff_line) + "\n")

            elif line.startswith("##sequence-region"):
                tag, seqid, start, end = line.rstrip().split()
                id_to_print = names_equiv[seqid]
                to_gff.write(
                    " ".join(
                        [
                            tag,
                            id_to_print,
                            start,
                            end,
                        ]
                    )
                    + "\n"
                )

            elif line.startswith(">"):
                seqid = line.rstrip().replace(">", "")
                id_to_print = names_equiv[seqid]
                to_gff.write(">" + id_to_print + "\n")

            else:
                to_gff.write(line)
