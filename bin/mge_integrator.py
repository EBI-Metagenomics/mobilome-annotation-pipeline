#!/usr/bin/env python


##### This script integrates the results for the Mobilome Annotation Pipeline
##### Alejandra Escobar, EMBL-EBI
##### January 11, 2023


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
        "--pal_tsv",
        type=str,
        help="PaliDIS predictions table",
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
        help="ICEfinder prediction files (concatenated)",
    )
    parser.add_argument(
        "--icf_lim",
        type=str,
        help="ICEfinder DR coordinates",
    )
    parser.add_argument(
        "--mog_tsv",
        type=str,
        help="Diamond output versus MobileOG-DB format 6",
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
        "--virify_out",
        type=str,
        help="HQ virify results",
    )
    parser.add_argument(
        "--crispr_out",
        type=str,
        help="CRISPRCasFinder results (tsv file)",
    )
    args = parser.parse_args()

    ### Calling functions
    ## Mapping contig names
    (names_equiv, inv_names_equiv) = mapping_names.names_map(args.map)

    ## Parsing results of mobilome prediction tools
    # Parsing ICEfinder results
    (mge_data, icf_dr,) = icefinder_process.icf_parser(
        args.icf_lim,
        args.icf_tsv,
    )

    # Parsing IntegronFinder results
    (mge_data, attC_site,) = integronfinder_process.integron_parser(
        mge_data,
        args.inf_tsv,
        args.inf_gbks,
    )

    # Parsing ISEScan results
    (mge_data, itr_sites,) = isescan_process.isescan_parser(
        mge_data,
        args.iss_tsv,
    )

    # Parsing PaliDIS results and removing redundancy with ISEScan
    if os.path.exists(args.pal_tsv):
        (mge_data, itr_sites) = palidis_process.palids_parser(
            args.pal_tsv, inv_names_equiv, mge_data, itr_sites
        )

    # Parsing geNomad results
    (mge_data) = genomad_parser.genomad_viral(args.geno_out, mge_data)
    (mge_data) = genomad_parser.plasmids_parser(args.geno_plas, mge_data)

    # Parsing VIRIfy results and solving redundancy with geNomad
    if os.path.exists(args.virify_out):
        (mge_data, virify_prots) = virify_process.virify_reader(
            args.virify_out, inv_names_equiv, mge_data
        )
    else:
        virify_prots = {}

    ## Overlapping report for long MGEs
    # Collecting integrons, virus and plasmids predicted per contig
    overlap_finder.overlap_report(mge_data, names_equiv)

    ## Catching proteins on the mobilome and MGEs QC
    (mge_data, contigs_elements, proteins_mge) = cds_locator.location_parser(
        args.pkka_gff, mge_data
    )

    ## Storing extra annotation results
    # Parsing mobileOG results
    mog_annot = mobileog_process.mobileog_parser(args.mog_tsv)

    # Parsing CRISPRCasFinder results
    if os.path.exists(args.crispr_out):
        crispr_annot = crispr_process.crispr_parser(args.crispr_out, args.pkka_gff)
    else:
        crispr_annot = {}

    # Adding the mobilome annotation to the GFF file
    integrator_process.gff_writer(
        args.pkka_gff,
        names_equiv,
        contigs_elements,
        itr_sites,
        icf_dr,
        attC_site,
        mge_data,
        proteins_mge,
        mog_annot,
        crispr_annot,
        virify_prots,
    )


if __name__ == "__main__":
    import argparse
    import os.path
    import sys

    sys.path.append("./python_modules")

    from python_modules import cds_locator
    from python_modules import crispr_process
    from python_modules import genomad_parser
    from python_modules import icefinder_process
    from python_modules import integrator_process
    from python_modules import integronfinder_process
    from python_modules import isescan_process
    from python_modules import mapping_names
    from python_modules import mobileog_process
    from python_modules import overlap_finder
    from python_modules import palidis_process
    from python_modules import virify_process

    main()
