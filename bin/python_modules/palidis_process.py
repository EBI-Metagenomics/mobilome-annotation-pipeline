#!/usr/bin/env python

import sys
import os.path

##### This module parse palidis files and check redundancy with ISEscan
##### Alejandra Escobar, EMBL-EBI
##### October 4, 2023


def mge_data_parser(mge_data):
    iss_dir = {}
    for mge in mge_data:
        contig, description, coord = mge_data[mge]
        prefix = mge.split("_")[0]
        if prefix == "iss":
            if contig in iss_dir:
                iss_dir[contig].append(coord)
            else:
                iss_dir[contig] = [coord]

    return iss_dir


def palids_parser(pal_tsv, inv_names_equiv, mge_data, itr_sites):
    ### Parsing Palidis output
    mge_counter = 0
    palidis_predictions = {}
    with open(pal_tsv, "r") as input_table:
        next(input_table)
        for line in input_table:
            mge_counter += 1
            mge_id = "pal_" + str(mge_counter)
            (
                is_name,
                sample_id,
                contig,
                itr1_start_position,
                itr1_end_position,
                itr2_start_position,
                itr2_end_position,
                pal_description,
            ) = line.rstrip().split("\t")
            description = "mobile_element_type=IS_with_TIR"
            coord = (int(itr1_start_position), int(itr2_end_position))
            contig = inv_names_equiv[contig]
            value = (contig, description, coord)
            palidis_predictions[mge_id] = value
            ir_1 = (itr1_start_position, itr1_end_position)
            ir_2 = (itr2_start_position, itr2_end_position)
            itr_sites[mge_id] = (ir_1, ir_2)

    ## Parsing the mge_data
    iss_dir = mge_data_parser(mge_data)

    # Finding overlapping predictions
    to_discard = []
    for pal_id in palidis_predictions:
        pal_contig, pal_desc, pal_coord = palidis_predictions[pal_id]
        p_start = pal_coord[0]
        p_end = pal_coord[1]
        p_len = p_end - p_start
        p_range = range(p_start, p_end + 1)

        if pal_contig in iss_dir:
            for i_coord_pair in iss_dir[pal_contig]:
                i_start = i_coord_pair[0]
                i_end = i_coord_pair[1]
                i_len = i_end - i_start
                i_range = range(i_start, i_end + 1)
                intersection = len(list(set(p_range) & set(i_range)))
                if intersection > 0:
                    p_cov = float(intersection) / float(p_len)
                    i_cov = float(intersection) / float(i_len)
                    if any([p_cov > 0.1, i_cov > 0.1]):
                        to_discard.append(pal_id)

    ## Discarding redundancy from virify dictionary
    to_discard = list(set(to_discard))
    print("Number of PaliDis predictions discarded: " + str(len(to_discard)))
    for pal_id in to_discard:
        del palidis_predictions[pal_id]
    print(
        "Number of PaliDis predictions remained: "
        + str(len(list(palidis_predictions.keys())))
    )

    ## Adding missing viral predictions to mge_data
    for pal_mge in palidis_predictions:
        mge_data[pal_mge] = palidis_predictions[pal_mge]

    return (mge_data, itr_sites)
