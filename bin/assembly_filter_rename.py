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

from Bio import SeqIO
import argparse


def rename(input_file, prefix):
    output_1kb = prefix + "_1kb_contigs.fasta"
    output_5kb = prefix + "_5kb_contigs.fasta"
    output_100kb = prefix + "_100kb_contigs.fasta"
    output_map = prefix + "_contigID.map"
    counter = 0

    with open(output_1kb, "w") as to_1kb, open(output_5kb, "w") as to_5kb, open(output_map, "w") as to_map, open(output_100kb, "w") as to_100kb:
        for record in SeqIO.parse(input_file, "fasta"):
            counter += 1
            new_id = ">contig_" + str(counter)
            my_chain = str(record.seq).upper()
            to_map.write(new_id + "\t" + str(record.id) + "\n")
            if len(my_chain) > 1000:
                to_1kb.write(new_id + "\n")
                to_1kb.write(my_chain + "\n")
            if len(my_chain) > 5000:
                to_5kb.write(new_id + "\n")
                to_5kb.write(my_chain + "\n")
            if len(my_chain) >= 100000:
                to_100kb.write(new_id + "\n")
                to_100kb.write(my_chain + "\n")


def main():
    parser = argparse.ArgumentParser(
        description="This script filter by contig size an assembly and generates two files to feed the MoMofy pipeline: contigs > 1kb for insertion sequences prediction and contigs > 5 kb for integron prediction"
    )
    parser.add_argument(
        "--assembly",
        type=str,
        help="Input fasta file",
        required=True,
    )
    parser.add_argument(
        "--output",
        type=str,
        help="Output prefix",
        required=True,
    )
    args = parser.parse_args()

    rename(args.assembly, args.output)


if __name__ == "__main__":
    main()
