#!/usr/bin/env python

import sys
import os.path

##### This module write the report of overlapping long integrons
##### Alejandra Escobar, EMBL-EBI
##### October 12, 2023


def overlap_report(mge_data, names_equiv):
    # Grouping elements by contig
    contig_inmge = {}
    long_mges = [
        "icf",
        "inf",
        "vir1",
        "vir2",
        "plas",
        "phpl",
    ]

    for prediction in mge_data:
        prefix = prediction.split("_")[0]
        if prefix in long_mges:
            contig = mge_data[prediction][0]
            if contig in contig_inmge:
                contig_inmge[contig].append(prediction)
            else:
                contig_inmge[contig] = [prediction]

    output_nested = "overlapping_integrons.txt"
    with open(output_nested, "w") as to_nested:
        to_nested.write(
            "\t".join(
                [
                    "contig",
                    "coord_1",
                    "coord_2",
                    "coverage_1",
                    "coverage_2",
                ]
            )
            + "\n"
        )
        for contig in contig_inmge:
            if len(contig_inmge[contig]) > 1:
                tools_list = []
                for e1 in contig_inmge[contig]:
                    e1_coord = mge_data[e1][2]
                    e1_len = e1_coord[1] - e1_coord[0]
                    for e2 in contig_inmge[contig]:
                        if e1 != e2:
                            e2_coord = mge_data[e2][2]
                            e2_len = e2_coord[1] - e2_coord[0]
                            intersection = len(
                                list(
                                    set(range(e1_coord[0], e1_coord[1] + 1))
                                    & set(range(e2_coord[0], e2_coord[1] + 1))
                                )
                            )
                            if intersection > 1000:
                                e1_cov = float(intersection) / float(e1_len)
                                e2_cov = float(intersection) / float(e2_len)
                                to_nested.write(
                                    "\t".join(
                                        [
                                            names_equiv[contig],
                                            str(e1_coord[0]) + "-" + str(e1_coord[1]),
                                            str(e2_coord[0]) + "-" + str(e2_coord[1]),
                                            str(e1_cov),
                                            str(e2_cov),
                                        ]
                                    )
                                    + "\n"
                                )
