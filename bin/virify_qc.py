#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2024 EMBL - European Bioinformatics Institute
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
import os.path


# TODO: this is not used anymore - move relevant pieces to the other method
def quality_parser(checkv_out):
    checkv_pass = []
    for summ_file in checkv_out:
        if os.path.exists(summ_file):
            head, tail = os.path.split(summ_file)
            current_qual = tail.split("_")[0]
            with open(summ_file, "r") as input_table:
                next(input_table)
                for line in input_table:
                    l_line = line.rstrip().split("\t")
                    phage_id = l_line[0].replace("prophage-0:", "prophage-1:")
                    if "|prophage" not in phage_id:
                        phage_id = phage_id + "|viral_sequence"
                    viral_genes = int(l_line[5])
                    checkv_quality = l_line[7]
                    kmer_freq = float(l_line[12])
                    if current_qual == "high":
                        checkv_pass.append(phage_id)
                    else:
                        if all([viral_genes > 0, checkv_quality != "Not-determined"]):
                            if any(
                                [
                                    checkv_quality == "Low-quality",
                                    checkv_quality == "Medium",
                                ]
                            ):
                                if kmer_freq <= 1.0:
                                    checkv_pass.append(phage_id)
                            else:
                                checkv_pass.append(phage_id)
        else:
            print("No checkV files found\n")

    return checkv_pass


def virify_parser(virify_gff, hq_list, output_prefix):
    hq_prophages = {}
    hq_viral_contigs = []
    for phage_id in hq_list:
        if "|prophage" in phage_id:
            prophage_contig = phage_id.split("|")[0]
            coordinates = phage_id.split("|")[1].split("-")[1]
            start = int(coordinates.split(":")[0])
            end = int(coordinates.split(":")[1])
            coord_tuple = (start, end)
            if prophage_contig in hq_prophages:
                hq_prophages[prophage_contig].append(coord_tuple)
            else:
                hq_prophages[prophage_contig] = [coord_tuple]
        else:
            hq_viral_contigs.append(phage_id.split("|")[0])

    with open(virify_gff, "r") as input_table, open(f"{output_prefix}_virify_hq.gff", "w") as output_gff:
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

                # Saving mobilome predictions
                if seq_source == "VIRify":
                    feature_id = attr.split(";")[0].replace("ID=", "")
                    if feature_id in hq_list:
                        output_gff.write(line + "\n")
                else:
                    if contig in hq_viral_contigs:
                        output_gff.write(line + "\n")
                    elif contig in hq_prophages:
                        prot_start = int(start)
                        prot_end = int(end)
                        prot_range = range(prot_start, prot_end + 1)
                        prot_len = prot_end - prot_start
                        for phage_loc in hq_prophages[contig]:
                            phage_start = phage_loc[0]
                            phage_end = phage_loc[1]
                            phage_range = range(phage_start, phage_end + 1)
                            intersection = len(list(set(phage_range) & set(prot_range)))
                            if intersection > 0:
                                prot_cov = float(intersection) / float(prot_len)
                                if prot_cov > 0.75:
                                    output_gff.write(line + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="This script parse the virify output for quality control."
    )
    parser.add_argument(
        "--virify_gff",
        type=str,
        help="VIRify v2.0 output in GFF format (08-final/gff/)",
    )
    parser.add_argument(
        "--prefix",
        type=str,
        help="Output files prefix",
        required=True
    )

    args = parser.parse_args()
    
    # TODO: this needs to be adjusted as the checkv files are not provided anymore

    virify_parser(args.virify_gff, args.prefix)


if __name__ == "__main__":
    main()
