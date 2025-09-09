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

#
# This script contains algorithms and methods derived from or inspired by ICEfinder2
# Original ICEfinder2 work licensed under CC BY-NC-SA 4.0
# (http://creativecommons.org/licenses/by-nc-sa/4.0/)
# Modified and adapted for this pipeline by EMBL-EBI
#

import argparse
import gzip
import os
import sys
from collections import defaultdict

from Bio import SeqIO
from Bio.SearchIO import parse


def scanf(hmmlist):
    """
    Scan HMM list to determine if a contig contains essential ICE components.

    This function analyzes a list of HMM hits to determine
    whether a contig contains the essential components of an Integrative and
    Conjugative Element (ICE), including mobilization proteins, Type IV coupling
    proteins, integrases, and Type IV secretion systems.

    :param hmmlist: List of HMM hit names
    :type hmmlist: list
    :return: True if contig contains essential ICE components, False otherwise
    :rtype: bool

    .. note::
       The function requires at least: 1 MOB protein, 1 T4CP protein,
       1 integrase, and 5 T4SS proteins for a positive ICE classification.
    """
    ICEcount = []
    for line in hmmlist:
        if "MOB" in line:
            ICEcount.append("MOB")
        elif "t4cp" in line or "tcpA" in line:
            ICEcount.append("t4cp")
        elif "FA" in line:
            ICEcount.append("T4SS")
        elif line in [
            "Phage_integrase",
            "UPF0236",
            "Recombinase",
            "rve",
            "TIGR02224",
            "TIGR02249",
            "TIGR02225",
            "PB001819",
        ]:
            ICEcount.append("Int")
        else:
            ICEcount.append("T4SS")
    if (
        ICEcount.count("MOB")
        and ICEcount.count("t4cp")
        and ICEcount.count("Int")
        and ICEcount.count("T4SS") >= 5
    ):
        return True
    else:
        return False


def hmm_parser(hmm_out, evalue_threshold=0.00001):
    """
    Parse HMM search output to identify candidate ICE-containing contigs.

    This function processes HMM search results and groups hits by contig,
    filtering by E-value threshold. It then applies ICE component scanning
    to identify contigs that contain essential ICE components.

    :param hmm_out: Path to HMM search output file
    :type hmm_out: str
    :param evalue_threshold: E-value threshold for filtering HMM hits
    :type evalue_threshold: float
    :return: List of contig identifiers that contain ICE candidates
    :rtype: list
    .. note::
       The default E-value threshold of 0.00001 provides stringent filtering
       for high-quality HMM matches.
    """
    # Validate input file exists
    if not os.path.exists(hmm_out):
        raise FileNotFoundError(f"HMM output file not found: {hmm_out}")

    icedict = defaultdict(list)
    processed_ids = set()

    with open(hmm_out, "r") as outfile:
        for line_num, line in enumerate(outfile, 1):
            # Skip comment lines
            if line.startswith("#") or not line.strip():
                continue

            # Parse line with error handling
            fields = line.strip().split()
            if len(fields) < 5:
                raise ValueError(
                    f"Invalid format at line {line_num}: insufficient fields"
                )

            query_name, _, target_id, _, evalue_str = fields[:5]

            # Skip already processed targets
            if target_id in processed_ids:
                continue

            processed_ids.add(target_id)

            # Extract contig identifier (first two underscore-separated parts)
            id_parts = target_id.split("_")
            if len(id_parts) < 2:
                continue  # Skip malformed IDs

            contig_key = "_".join(id_parts[:2])

            # Filter by E-value threshold
            try:
                evalue = float(evalue_str)
                if evalue < evalue_threshold:
                    icedict[contig_key].append(query_name)
            except ValueError:
                continue  # Skip lines with invalid E-values

    # Apply ICE component scanning and return candidates
    return [
        contig_key for contig_key, components in icedict.items() if scanf(components)
    ]


def fasta_parser(assembly_file, candidates_list, output_prefix):
    """
    Extract candidate contig sequences from the assembly FASTA file.

    :param assembly_file: Path to assembly FASTA file
    :type assembly_file: str
    :param candidates_list: List of candidate contig identifiers
    :type candidates_list: list
    :param output_prefix: Output file prefix for candidate sequences
    :type output_prefix: str
    """
    candidates_set = set(candidates_list)

    # Determine if input file is gzip-compressed
    try:
        if assembly_file.endswith(".gz"):
            file_handle = gzip.open(assembly_file, "rt")
        else:
            file_handle = open(assembly_file, "r")

        with open(output_prefix + "_candidates.fasta", "w") as fasta_out:
            for record in SeqIO.parse(file_handle, "fasta"):
                seq_id = str(record.id)
                if seq_id in candidates_set:
                    fasta_out.write(f">{seq_id}\n")
                    fasta_out.write(f"{str(record.seq).upper()}\n")
    finally:
        file_handle.close()


def proteins_parser(proteins_fasta, candidates_list, output_prefix):
    """
    Extract protein sequences from candidate contigs.

    :param proteins_fasta: Path to protein FASTA file
    :type proteins_fasta: str
    :param candidates_list: List of candidate contig identifiers
    :type candidates_list: list
    :param output_prefix: Output file prefix for candidate protein sequences
    :type output_prefix: str
    """
    # Determine if input file is gzip-compressed
    try:
        if proteins_fasta.endswith(".gz"):
            file_handle = gzip.open(proteins_fasta, "rt")
        else:
            file_handle = open(proteins_fasta, "r")

        with open(output_prefix + "_candidates.faa", "w") as faa_out:
            for record in SeqIO.parse(file_handle, "fasta"):
                prot_id = str(record.id)
                seq = str(record.seq).upper().replace("*", "")
                id_split = prot_id.split("_")
                id_split.pop(-1)
                contig_id = "_".join(id_split)
                if contig_id in candidates_list:
                    faa_out.write(">" + prot_id + "\n")
                    faa_out.write(seq + "\n")
    finally:
        file_handle.close()


def main():
    parser = argparse.ArgumentParser(
        description="Process hmmscan output vs ICEs models for prescaning of contigs"
    )
    parser.add_argument(
        "--hmm_out", required=True, help="Path to all_systems.tsv file from macsyfinder"
    )
    parser.add_argument(
        "--proteins", required=True, help="Path to aminoacids FASTA file"
    )
    parser.add_argument("--assembly", required=True, help="Path to assembly FASTA file")
    parser.add_argument("--output", required=True, help="Output prefix name")
    parser.add_argument(
        "--evalue_threshold",
        type=float,
        default=0.00001,
        help="E-value threshold for filtering HMM hits (default: 0.00001)",
    )
    args = parser.parse_args()

    # Validate input files
    if not os.path.exists(args.hmm_out):
        print(f"Error: hmmer output file not found: {args.hmm_out}")
        sys.exit(1)

    if not os.path.exists(args.proteins):
        print(f"Error: Proteins file not found: {args.proteins}")
        sys.exit(1)

    if not os.path.exists(args.assembly):
        print(f"Error: Assembly file not found: {args.assembly}")
        sys.exit(1)

    # Running functions
    candidates_list = hmm_parser(args.hmm_out, args.evalue_threshold)

    if len(candidates_list) > 0:
        fasta_parser(args.assembly, candidates_list, args.output)
        proteins_parser(args.proteins, candidates_list, args.output)


if __name__ == "__main__":
    main()
