#!/usr/bin/env python

import sys
import os.path

##### This module parse the VIRify output files
##### Alejandra Escobar, EMBL-EBI
##### October 9, 2023


def mge_data_parser(mge_data):
    plasmids_list = []
    prophages_dic, prophages_ids, viral_dic = {}, {}, {}
    for mge in mge_data:
        contig, description, coord = mge_data[mge]
        prefix = mge.split("_")[0]
        if prefix == "vir1":
            if "viral_sequence" in description.split(";")[0]:
                viral_dic[contig] = mge
            elif "prophage" in description.split(";")[0]:
                composite_key = (contig, coord)
                prophages_ids[composite_key] = mge
                if contig in prophages_dic:
                    prophages_dic[contig].append(coord)
                else:
                    prophages_dic[contig] = [coord]
        if prefix == "plas":
            plasmids_list.append(contig)

    return (viral_dic, prophages_dic, prophages_ids, plasmids_list)


def virify_reader(virify_gff, inv_names_equiv, mge_data):
    virify_predictions, virify_prots = {}, {}
    mge_counter = 0
    with open(virify_gff, "r") as input_table:
        for line in input_table:
            line_l = line.rstrip().split("\t")
            # Annotation lines have exactly 9 columns
            if len(line_l) == 9:
                (
                    contig,
                    seq_source,
                    seq_type,
                    start,
                    end,
                    score,
                    strand,
                    phase,
                    attr,
                ) = line.rstrip().split("\t")
                contig = inv_names_equiv[contig]
                coord = (int(start), int(end))

                if seq_source == "VIRify":
                    mge_counter += 1
                    mge_id = "vir2_" + str(mge_counter)

                    (
                        vir_mge_id,
                        gbkey,
                        mobile_element_type,
                        checkv_provirus,
                        checkv_quality,
                        miuvig_quality,
                        taxonomy,
                    ) = attr.split(";")

                    composite_val = (contig, attr, coord)
                    virify_predictions[mge_id] = composite_val

                # Saving protein predictions
                elif seq_source == "Prodigal":
                    (gene_id, gene_gbkey, viphog, viphog_taxonomy) = attr.split(";")
                    prot_viphog = viphog + ";" + viphog_taxonomy
                    prot_location = (contig, int(start), int(end))
                    virify_prots[prot_location] = prot_viphog

    ## Parsing the mge_data
    (viral_dic, prophages_dic, prophages_ids, plasmids_list) = mge_data_parser(mge_data)

    ## Removing redundancy on viral genome fragments
    to_discard = []
    phage_plasmids = []
    for phage in virify_predictions:
        v_contig, v_description, v_coord = virify_predictions[phage]

        # Finding and catching phage-plasmids
        if v_contig in plasmids_list:
            phage_plasmids.append(v_contig)

        # Finding redundancy on viral genome fragments
        if "viral_sequence" in v_description.split(";")[0]:
            if v_contig in viral_dic:
                genomad_id = viral_dic[v_contig]
                to_discard.append(genomad_id)

        # Finding redundancy on prophages
        elif "prophage" in v_description.split(";")[0]:
            v_start = v_coord[0]
            v_end = v_coord[1]
            v_len = v_end - v_start
            v_range = range(v_start, v_end + 1)
            if v_contig in prophages_dic:
                for g_coord_pair in prophages_dic[v_contig]:
                    g_start = g_coord_pair[0]
                    g_end = g_coord_pair[1]
                    g_id = prophages_ids[(v_contig, (g_start, g_end))]
                    g_len = g_end - g_start
                    g_range = range(g_start, g_end + 1)
                    intersection = len(list(set(v_range) & set(g_range)))
                    if intersection > 0:
                        v_cov = float(intersection) / float(v_len)
                        g_cov = float(intersection) / float(g_len)
                        if any([v_cov > 0.25, g_cov > 0.25]):
                            to_discard.append(g_id)

    ## Discarding genomad redundance from mge_data dictionary
    to_discard = list(set(to_discard))
    print("Number of geNomad predictions discarded: " + str(len(to_discard)))
    for phage_id in to_discard:
        if phage_id in mge_data:
            print(phage_id, mge_data[phage_id])
            del mge_data[phage_id]

    ## Adding missing viral predictions to mge_data
    print(
        "Number of VIRify predictions to be added: "
        + str(len(list(virify_predictions.keys())))
    )
    for phage in virify_predictions:
        description = virify_predictions[phage][1].split(";")
        description.pop(0)
        description.pop(0)
        description = ";".join(description)
        new_value = (
            virify_predictions[phage][0],
            description,
            virify_predictions[phage][2],
        )
        mge_data[phage] = new_value
        print(phage, new_value)

    ## Labelling phage-plasmids
    # Storing the plasmid IDs to access values
    plas_phage_ids = []
    for mge in mge_data:
        contig, description, coord = mge_data[mge]
        prefix = mge.split("_")[0]
        if prefix == "plas":
            if contig in phage_plasmids:
                plas_phage_ids.append(mge)

    # Replacing description
    mge_counter = 0
    for pp in plas_phage_ids:
        mge_counter += 1
        old_contig, old_description, old_coord = mge_data[pp]
        new_desc = "mobile_element_type=phage_plasmid"
        new_val = (contig, new_desc, old_coord)
        new_mge_id = "phpl_" + str(mge_counter)
        del mge_data[pp]
        mge_data[new_mge_id] = new_val

    return (mge_data, virify_prots)
