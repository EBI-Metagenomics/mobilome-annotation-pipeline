#!/usr/bin/env python

import argparse
import os.path
import csv

##### This script integrates the results of amrfinderplus with the mobilome
##### Alejandra Escobar, EMBL-EBI
##### May 31, 2023


def user_gff_parser(user_gff):
    user_genes = {}
    if os.stat(user_gff).st_size == 0:
        return user_genes

    with open(user_gff, "r") as input_file:
        for line in input_file:
            l_line = line.rstrip().split("\t")
            # Annotation lines have exactly 9 columns
            if len(l_line) == 9:
                contig = l_line[0]
                start = l_line[3]
                end = l_line[4]
                strand = l_line[6]
                attrib = l_line[8]
                prot_id = attrib.split(";")[0].replace("ID=", "")
                composite_key = (contig, start, end, strand)
                user_genes[composite_key] = prot_id
    return user_genes


def names_parser(contigs_map):
    contig_names = {}
    with open(contigs_map, "r") as input_file:
        for line in input_file:
            pkka_name, fasta_name = line.rstrip().split("\t")
            contig_names[pkka_name.replace(">", "")] = fasta_name
    return contig_names


def arg_parser(amr_out, contig_names):
    """
    Parse the AMRFinderPlus CSV output and map contig IDs to contig names.

    This function reads a CSV file generated by AMRFinderPlus and processes
    its contents to extract AMR (antimicrobial resistance) predictions. It
    returns a dictionary where keys are protein IDs and values are the
    corresponding contig names.

    Args:
        amr_out (str): The file path to the AMRFinderPlus CSV output.
        contig_names (dict): A dictionary mapping contig IDs to contig names.

    Returns:
        dict: A dictionary containing AMR predictions with protein IDs as keys
        and contig names as values. If the file does not exist or is empty,
        an empty dictionary is returned.
    """

    amr_data = {}
    if not os.path.exists(amr_out) and os.stat(amr_out).st_size == 0:
        return arm_data

    with open(amr_out, "r") as input_csv:
        csv_reader = csv.DictReader(input_csv, delimiter="\t")
        for row in csv_reader:
            protein_id = row["Protein identifier"]
            contig_id = contig_names[row["Contig id"]]
            start = int(row["Start"])
            end = int(row["Stop"])
            strand = row["Strand"]
            gene_name = row["Gene symbol"]
            sequence_name = row["Sequence name"]
            aggregated_class = "|".join(
                [
                    row["Element type"],
                    row["Element subtype"],
                    row["Class"],
                    row["Subclass"],
                ]
            )
            coverage = row["% Coverage of reference sequence"]
            identity = row["% Identity to reference sequence"]
            descriptors = "\t".join(
                [coverage, identity, aggregated_class, gene_name, sequence_name]
            )
            composite_key = (protein_id, contig_id, start, end, strand)
            amr_data[composite_key] = descriptors
    return amr_data


def mob_parser(mobilome):
    MOB_THRES = 0.75
    ### Saving the proteins in the mobilome
    mob_coords, mob_types, prots_loc = {}, {}, {}
    mges_list = [
        "conjugative_integron",
        "insertion_sequence",
        "integron",
        "phage_plasmid",
        "plasmid",
        "prophage",
        "viral_sequence",
    ]
    with open(mobilome, "r") as input_gff:
        for line in input_gff:
            l_line = line.rstrip().split("\t")
            ## Annotation lines have exactly 9 columns
            if len(l_line) == 9:
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

                # Saving genes tagged as mobilome
                if seq_type == "CDS":
                    att_fields = attr.split(";")
                    protein_id = att_fields[0].replace("ID=", "")
                    for attribute in att_fields:
                        if attribute == "location=mobilome":
                            prots_loc[protein_id] = (contig, start, end)

                # Saving the MGEs coordinates
                elif seq_type in mges_list:
                    coords_tuple = (int(start), int(end))
                    if contig in mob_coords:
                        mob_coords[contig].append(coords_tuple)
                    else:
                        mob_coords[contig] = [coords_tuple]
                    composite_key = (contig, int(start), int(end))
                    mob_types[composite_key] = seq_type

    # Adding the MGE type to each protein
    mob_prots = {}
    for protein_id in prots_loc:
        mges_loc = []
        prot_contig = prots_loc[protein_id][0]
        prot_start = int(prots_loc[protein_id][1])
        prot_end = int(prots_loc[protein_id][2])
        prot_range = range(prot_start, prot_end + 1)
        prot_len = prot_end - prot_start
        if prot_contig in mob_coords:
            for mobile_element in mob_coords[prot_contig]:
                mge_start = int(mobile_element[0])
                mge_end = int(mobile_element[1])
                mge_range = range(mge_start, mge_end + 1)
                current_mge_type = mob_types[(prot_contig, mge_start, mge_end)]
                intersection = len(list(set(mge_range) & set(prot_range)))
                if intersection > 0:
                    prot_cov = float(intersection) / float(prot_len)
                    if prot_cov > MOB_THRES:
                        mges_loc.append(current_mge_type)

        if len(mges_loc) > 0:
            mob_prots[protein_id] = ";".join(mges_loc)
        else:
            print(
                f"No MGE in the protein coordinates for protein {protein_id} in contig {prot_contig}"
            )

    return mob_prots, mob_coords, mob_types


def location_parser(amr_data, mob_prots, mob_coords, mob_types, user_genes):
    AMR_GENE_THRES = 0.75
    with open("amr_location.tsv", "w", newline="") as to_output:
        writer = csv.writer(
            to_output, delimiter="\t", quoting=csv.QUOTE_MINIMAL, lineterminator="\n"
        )

        # Header
        headers_with_user_gene = [
            "user_gene_id",
            "prokka_gene_id",
            "contig_id",
            "gene_start",
            "gene_end",
            "ref_coverage",
            "ref_identity",
            "class_summary",
            "ref_gene_name",
            "ref_sequence_name",
            "location",
        ]
        headers_without_user_gene = [
            "gene_id",
            "contig_id",
            "gene_start",
            "gene_end",
            "ref_coverage",
            "ref_identity",
            "class_summary",
            "ref_gene_name",
            "ref_sequence_name",
            "location",
        ]

        if len(user_genes) > 0:
            writer.writerow(headers_with_user_gene)
        else:
            writer.writerow(headers_without_user_gene)

        # Data rows
        for composite_key in amr_data:
            protein_id, contig_id, amr_start, amr_end, strand = composite_key
            description = amr_data[composite_key].split("\t")
            amr_location_key = (contig_id, str(amr_start), str(amr_end), strand)
            user_prot_id = user_genes.get(amr_location_key, "NA")

            if protein_id == "NA":  # Special processing for NA protein IDs
                mges_loc = []
                if contig_id in mob_coords:
                    amr_range = range(amr_start, amr_end + 1)
                    amr_len = amr_end - amr_start + 1
                    for mge_start, mge_end in mob_coords[contig_id]:
                        mge_range = range(mge_start, mge_end + 1)
                        intersection_len = len(list(set(mge_range) & set(amr_range)))
                        if intersection_len > 0:
                            amr_cov = intersection_len / amr_len
                            if amr_cov > AMR_GENE_THRES:
                                mges_loc.append(
                                    mob_types[(contig_id, mge_start, mge_end)]
                                )
                location = (
                    "mobilome:" + ";".join(mges_loc) if mges_loc else "chromosome"
                )

                if len(user_genes) > 0:
                    row = (
                        [
                            user_prot_id,
                            protein_id,
                            contig_id,
                            str(amr_start),
                            str(amr_end),
                        ]
                        + description
                        + [location]
                    )
                else:
                    row = (
                        [protein_id, contig_id, str(amr_start), str(amr_end)]
                        + description
                        + [location]
                    )
                writer.writerow(row)

            else:  # Regular processing for valid protein IDs
                location = (
                    "mobilome:" + mob_prots[protein_id]
                    if protein_id in mob_prots
                    else "chromosome"
                )
                if len(user_genes) > 0:
                    row = (
                        [
                            user_prot_id,
                            protein_id,
                            contig_id,
                            str(amr_start),
                            str(amr_end),
                        ]
                        + description
                        + [location]
                    )
                else:
                    row = (
                        [protein_id, contig_id, str(amr_start), str(amr_end)]
                        + description
                        + [location]
                    )
                writer.writerow(row)


def main():
    parser = argparse.ArgumentParser(
        description="This script integrates the results of amrfinderplus with the mobilome. Please provide the relevant input files"
    )
    parser.add_argument(
        "--mobilome",
        type=str,
        help="Clean version of the output of the mobilome annotation pipeline in GFF format (mobilome_prokka.gff)",
    )
    parser.add_argument(
        "--amr_out",
        type=str,
        help="ARG prediction output (amrfinderplus.tsv)",
        required=True,
    )
    parser.add_argument(
        "--contigs_map",
        type=str,
        help="Mapping file with contig names (contigID.map)",
        required=True,
    )
    parser.add_argument(
        "--user_gff",
        type=str,
        help="User gff file",
        required=False,
    )
    args = parser.parse_args()

    ### Calling functions

    user_genes = {}
    if args.user_gff:
        user_genes = user_gff_parser(args.user_gff)

    contig_names = names_parser(args.contigs_map)
    amr_data = arg_parser(args.amr_out, contig_names)
    (mob_prots, mob_coords, mob_types) = mob_parser(args.mobilome)
    location_parser(amr_data, mob_prots, mob_coords, mob_types, user_genes)


if __name__ == "__main__":
    main()
