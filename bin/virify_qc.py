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

import argparse
import gzip
import os.path


def open_file(filename):
    """
    Open a file, handling both compressed (.gz) and uncompressed files.
    Returns a file handle that can be used for reading.
    """
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')  # 'rt' for text mode
    else:
        return open(filename, 'r')


def virify_parser(virify_gff, output_prefix):
    qc_passed = []
    all_proteins = {}
    with (
        open_file(virify_gff) as input_table,
        open(f"{output_prefix}_virify_hq.gff", "w") as output_gff,
    ):
        output_gff.write("##gff-version 3\n")
        for line in input_table:
            line = line.rstrip()
            line_l = line.split("\t")
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

                # Looking for HQ predictions
                features_list = attr.split(";")
                feature_id = features_list[0].split("=")[1]
                if seq_source == "VIRify":
                    virify_quality = features_list[1].split("=")[1]
                    checkv_quality = features_list[5].split("=")[1]
                    checkv_kmer_freq = float(features_list[7].split("=")[1])
                    checkv_viral_genes = int(features_list[8].split("=")[1])

                    # Keeping virify high-quality predictions
                    if virify_quality == "HC":
                        output_gff.write(line + "\n")
                        qc_passed.append(feature_id)
                    else:
                        if any(
                            [
                                checkv_quality == "High-quality",
                                checkv_quality == "Complete",
                            ]
                        ):
                            output_gff.write(line + "\n")
                            qc_passed.append(feature_id)
                        elif all(
                            [
                                checkv_viral_genes > 0,
                                checkv_quality != "Not-determined",
                                checkv_kmer_freq <= 1.0,
                            ]
                        ):
                            output_gff.write(line + "\n")
                            qc_passed.append(feature_id)

                # These are proteins. Saving and parsing later
                else:
                    id_spliced = feature_id.split("_")
                    id_spliced.pop(-1)
                    parent_feature = "_".join(id_spliced)

                    if not "prophage" in parent_feature:
                        parent_feature = parent_feature + "|viral_sequence"

                    if "|phage-circular" in parent_feature:
                        parent_feature = parent_feature.replace(
                            "phage-circular", "viral_sequence"
                        )

                    if parent_feature in all_proteins:
                        all_proteins[parent_feature].append(line)
                    else:
                        all_proteins[parent_feature] = [line]

        # Printing proteins belonging to HQ virify predictions
        for prediction_id in all_proteins:
            if prediction_id in qc_passed:
                for protein in all_proteins[prediction_id]:
                    output_gff.write(protein + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="This script parse the virify output for quality control."
    )
    parser.add_argument(
        "--virify_gff",
        type=str,
        help="VIRify v2.0 output in GFF format (08-final/gff/)",
    )
    parser.add_argument("--prefix", type=str, help="Output files prefix", required=True)

    args = parser.parse_args()

    virify_parser(args.virify_gff, args.prefix)


if __name__ == "__main__":
    main()
