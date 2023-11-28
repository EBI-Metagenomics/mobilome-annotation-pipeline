#!/usr/bin/env python

import os.path
from Bio import SeqIO

##### This module parse the IntegronFinder outputs
##### Alejandra Escobar, EMBL-EBI
##### October 4, 2023


def integron_parser(mge_data, integron_results, inf_gbks):
    ### Parsing IntegronFinder output
    mge_counter = 0
    attC_site = {}
    if os.stat(integron_results).st_size > 0:
        with open(integron_results, "r") as input_table:
            next(input_table)
            next(input_table)
            for line in input_table:
                id_replicon, calin, complete, in0, topology, size = line.rstrip().split(
                    "\t"
                )
                if int(complete) > 0:
                    description = "mobile_element_type=complete_integron"
                    gbk_file = id_replicon + ".gbk"
                    if gbk_file in inf_gbks:
                        flag = 0
                        for gb_record in SeqIO.parse(gbk_file, "genbank"):
                            for feature in gb_record.features:
                                if feature.type == "integron":
                                    if (
                                        feature.qualifiers["integron_type"][0]
                                        == "complete"
                                    ):
                                        flag = 1
                                        mge_counter += 1
                                        mge_id = "inf_" + str(mge_counter)
                                        start = int(feature.location.start)
                                        end = int(feature.location.end)
                                        coord = (start, end)
                                        value = (id_replicon, description, coord)
                                        mge_data[mge_id] = value
                                if feature.type == "attC":
                                    if flag == 1:
                                        start = str(feature.location.start)
                                        end = str(feature.location.end)
                                        coord = (start, end)
                                        attC_site[mge_id] = coord
                                        flag = 0
    return (mge_data, attC_site)
