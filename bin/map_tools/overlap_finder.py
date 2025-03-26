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


def overlap_report(mge_data, names_equiv, output_file = "overlapping_integrons.txt"):
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

    with open(output_file, "w") as to_nested:
        to_nested.write(
            "\t".join(
                [
                    "contig",
                    "coord_1",
                    "coord_2",
                    "coverage_1",
                    "coverage_2",
                    "type_1",
                    "type_2",
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
                                            e1,
                                            e2,
                                        ]
                                    )
                                    + "\n"
                                )
