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
import os.path
from map_tools import cds_locator
from map_tools import genomad_parser
from map_tools import icefinder_process
from map_tools import integrator_process
from map_tools import integronfinder_process
from map_tools import isescan_process
from map_tools import mapping_names
from map_tools import mobileog_process
from map_tools import outliers_process
from map_tools import overlap_finder
from map_tools import prokka_process
from map_tools import virify_process


def main():
    parser = argparse.ArgumentParser(
        description="This script integrates the results for the Mobilome Annotation Pipeline"
    )
    parser.add_argument(
        "--pkka_gff",
        type=str,
        help="Prokka output GFF file",
        required=True,
    )
    parser.add_argument(
        "--map",
        type=str,
        help="Rename contigs file: contigID.map",
        required=True,
    )
    parser.add_argument(
        "--iss_tsv",
        type=str,
        help="ISEScan predictions table",
    )
    parser.add_argument(
        "--inf_tsv",
        type=str,
        help="IntegronFinder predictions table",
    )
    parser.add_argument(
        "--inf_gbks",
        nargs="*",
        help="Space separated list of IntegonFinder gbk files per contig",
    )
    parser.add_argument(
        "--icf_tsv",
        type=str,
        nargs='?',
        help="ICEfinder2-lite prediction file",
    )
    parser.add_argument(
        "--geno_out",
        type=str,
        help="geNomad virus summary table (virus_summary.tsv)",
    )
    parser.add_argument(
        "--geno_plas",
        type=str,
        help="geNomad plasmids summary table (plasmid_summary.tsv)",
    )
    parser.add_argument(
        "--comp_bed",
        type=str,
        nargs='?',
        help="Compositional outliers prediction in bed format",
    )
    parser.add_argument(
        "--virify_out",
        type=str,
        help="HQ virify results",
    )
    parser.add_argument("--prefix", type=str, help="The output prefix", required=True)
    args = parser.parse_args()

    ### Calling functions
    mge_data = {}

    ## Mapping contig names
    (names_equiv, inv_names_equiv) = mapping_names.names_map(args.map)

    ## Parsing results of mobilome prediction tools
    # Parsing ICEfinder results
    if args.icf_tsv and os.path.exists(args.icf_tsv):
        (mge_data, icf_dr) = icefinder_process.icf_parser(args.icf_tsv)
    else:
        icf_dr = {}

    # Parsing IntegronFinder results
    (mge_data, attC_site) = integronfinder_process.integron_parser(
        mge_data,
        args.inf_tsv,
        args.inf_gbks,
    )

    # Parsing ISEScan results
    (mge_data, itr_sites) = isescan_process.isescan_parser(
        mge_data,
        args.iss_tsv,
    )

    # Parsing geNomad results
    (mge_data) = genomad_parser.genomad_viral(args.geno_out, mge_data)
    (mge_data) = genomad_parser.plasmids_parser(args.geno_plas, mge_data)

    # Parsing VIRIfy results and solving redundancy with geNomad
    if args.virify_out and os.path.exists(args.virify_out):
        (mge_data, virify_prots) = virify_process.virify_reader(
            args.virify_out, inv_names_equiv, mge_data
        )
    else:
        virify_prots = {}

    # Parsing prokka rrnas and genes location
    contig_prots, prots_coord, rnas_coord = prokka_process.prokka_parser(args.pkka_gff)

    # Parsing compositional outliers and removing redundancy with other MGEs
    if args.comp_bed and os.path.exists(args.comp_bed):
        mge_data, co_repeats = outliers_process.outliers_parser(
            args.comp_bed, mge_data, rnas_coord
        )
    else:
        co_repeats = {}

    ## Overlapping report for long MGEs
    # Collecting integrons, virus and plasmids predicted per contig
    overlap_finder.overlap_report(
        mge_data, names_equiv, output_file=f"{args.prefix}_overlap_report.txt"
    )

    ## Tagging genes on the mobilome and MGEs
    mge_data, contigs_elements, proteins_mge = cds_locator.location_parser(
        contig_prots,
        prots_coord,
        mge_data,
        output_file=f"{args.prefix}_discarded_mge.txt",
    )

    ## Storing extra annotation results
    # Parsing mobileOG results
    # mog_annot = mobileog_process.mobileog_parser(args.mog_tsv)

    # Adding the mobilome annotation to the GFF file
    integrator_process.gff_writer(
        args.pkka_gff,
        names_equiv,
        contigs_elements,
        itr_sites,
        icf_dr,
        attC_site,
        co_repeats,
        mge_data,
        proteins_mge,
        virify_prots,
        f"{args.prefix}_mobilome_prokka.gff",
    )


if __name__ == "__main__":
    main()
