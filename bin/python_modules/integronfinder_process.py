#!/usr/bin/env python

import os.path
from Bio import SeqIO

##### This module parse the IntegronFinder outputs
##### Alejandra Escobar, EMBL-EBI
##### October 4, 2023


def integron_parser(mge_data, integron_results, inf_gbks):
    ### Parsing IntegronFinder summary file
    gbk_files_list = []
    if os.stat(integron_results).st_size > 0:
        with open(integron_results, "r") as input_table:
            next(input_table)
            next(input_table)
            for line in input_table:
                id_replicon, calin, complete, in0, topology, size = line.rstrip().split("\t")

                # Saving gbk file names of complete integrons
                if int(complete) > 0:
                    gbk_file = id_replicon + ".gbk"
                    gbk_files_list.append(gbk_file)

    mge_counter = 0
    attC_site = {}
    for gbk_file in gbk_files_list:
        for gb_record in SeqIO.parse(gbk_file, "genbank"):

            # Get indexes for integron records and attC
            indexes_integron, indexes_attc = [], []
            for index, feature in enumerate(gb_record.features):
                if feature.type == "integron":
                    indexes_integron.append(index)
                elif feature.type == "attC":
                    indexes_attc.append(index)

            # process all integrons
            for index_integron in indexes_integron:
                integron_feature = gb_record.features[index_integron]
                if integron_feature.qualifiers["integron_type"][0] == "complete":
                    mge_counter += 1
                    mge_id = "inf_" + str(mge_counter)
                    mge_start = int(integron_feature.location.start)
                    if not mge_start:
                        mge_start = 1
                    mge_end = int(integron_feature.location.end)
                    mge_coord = (mge_start, mge_end)
                    description = "mobile_element_type=complete_integron"
                    id_replicon = gbk_file.replace('.gbk','')
                    value = (id_replicon, description, mge_coord)
                    mge_data[mge_id] = value

                    # check attC
                    for index_attc in indexes_attc:
                        attc_feature = gb_record.features[index_attc]
                        attc_start = str(attc_feature.location.start)
                        if not attc_start:
                            attc_start = 1
                        attc_end = str(attc_feature.location.end)
                        attc_coord = (attc_start, attc_end)

                        attC_site.setdefault(mge_id, [])
                        attC_site[mge_id].append(attc_coord)

    return (mge_data, attC_site)
