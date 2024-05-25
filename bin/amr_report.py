#!/usr/bin/env python

import argparse
import os.path

##### This script integrates the results of amrfinderplus with the mobilome
##### Alejandra Escobar, EMBL-EBI
##### May 31, 2023


def names_parser( contigs_map ):
    contig_names = {}
    with open(contigs_map, "r") as input_file:
        for line in input_file:
            pkka_name, fasta_name = line.rstrip().split("\t")
            contig_names[pkka_name.replace(">", "")] = fasta_name
    return contig_names


def arg_parser(amr_out, contig_names):
    ### Saving amr predictions
    amr_data = {}
    if os.path.exists(amr_out):
        if os.stat(amr_out).st_size > 0:
            with open(amr_out, "r") as input_csv:
                next(input_csv)
                for line in input_csv:
                    line_l = line.rstrip().split("\t")
                    protein_id = line_l[0]
                    contig_id = contig_names[line_l[1]]
                    start = int(line_l[2])
                    end = int(line_l[3])
                    strand = line_l[4]
                    gene_name = line_l[5]
                    sequence_name = line_l[6]
                    aggregated_class = '|'.join([
                        line_l[8],
                        line_l[9],
                        line_l[10],
                        line_l[11],
                    ])
                    coverage = line_l[15]
                    identity = line_l[16]
                    descriptors = '\t'.join([
                        coverage,
                        identity,
                        aggregated_class,
                        gene_name,
                        sequence_name
                    ])
                    composite_key = ( protein_id, contig_id, start, end, strand )
                    amr_data[composite_key] = descriptors
    return amr_data


def mob_parser(mobilome):
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
                            prots_loc[protein_id] = ( contig, start, end )

                # Saving the MGEs coordinates
                elif seq_type in mges_list: 
                    coords_tuple = ( int(start), int(end) )
                    if contig in mob_coords:
                        mob_coords[contig].append(coords_tuple)
                    else:
                        mob_coords[contig] = [coords_tuple]
                    composite_key = ( contig, int(start), int(end) )
                    mob_types[composite_key] = seq_type

    # Adding the MGE type to each protein
    mob_prots = {}
    for protein_id in prots_loc:
        mges_loc = []
        prot_contig = prots_loc[protein_id][0]
        prot_start = prots_loc[protein_id][1]
        prot_end = prots_loc[protein_id][2]
        print(prot_contig,prot_start,prot_end)
        prot_range = range(prot_start, prot_end + 1)
        prot_len = prot_end - prot_start
        if prot_contig in mob_coords:
            for mobile_element in mob_coords[prot_contig]:
                mge_start = mobile_element[0]
                mge_end = mobile_element[1]
                mge_range = range(mge_start, mge_end + 1)
                current_mge_type = mob_types[(prot_contig, mge_start, mge_end)]
                intersection = len(list(set(mge_range) & set(prot_range)))
                if intersection > 0:
                    prot_cov = float(intersection) / float(prot_len)
                    if prot_cov > 0.75:
                        mges_loc.append(current_mge_type)

        if len(mges_loc)>0:
            mob_prots[protein_id] = ';'.join(mges_loc)
        else:
            print("No MGE in the protein coordinates for protein "+protein_id+" in contig "+prot_contig)
    
    return ( mob_prots, mob_coords, mob_types )


def location_parser(amr_data, mob_prots, mob_coords, mob_types):
    with open("amr_location.txt", "w") as to_output:
        to_output.write(
            "\t".join([
                "gene_id",
                "contig_id", 
                "ref_coverage",
                "ref_identity",
                "class_summary",
                "ref_gene_name",
                "ref_sequence_name",
                "location"
            ]) + "\n"
        )
        
        for composite_key in amr_data:
            protein_id = composite_key[0]
            contig_id = composite_key[1]
            amr_start = composite_key[2]
            amr_end = composite_key[3]
            strand = composite_key[4]
            description = amr_data[composite_key]

            # When protein_id is NA, we need to check if the prediction coordinates are in the mobilome boundaries
            if protein_id == 'NA':
                mges_loc = []
                if contig_id in mob_coords:
                    amr_range = range( amr_start, amr_end + 1 )
                    amr_len = amr_end - amr_start
                    for mobile_element in mob_coords[contig_id]:
                        mge_start = mobile_element[0]
                        mge_end = mobile_element[1]
                        mge_range = range(mge_start, mge_end + 1)
                        current_mge_type = mob_types[(contig_id, mge_start, mge_end)]
                        intersection = len(list(set(mge_range) & set(amr_range)))
                        if intersection > 0:
                            amr_cov = float(intersection) / float(amr_len)
                            if amr_cov > 0.75:
                                mges_loc.append(current_mge_type)
                if len(mges_loc) > 0:
                    location = 'mobilome:'+';'.join(mges_loc)
                else:
                    location = 'chromosome'

                to_output.write("\t".join([
                    protein_id, 
                    contig_id,
                    description,
                    location,
                ])+ "\n")
            
            else:
                if protein_id in mob_prots:
                    location = 'mobilome:'+mob_prots[protein_id]
                    to_output.write("\t".join([
                        protein_id, 
                        contig_id,
                        description,
                        location,
                        ])+ "\n")
                else:
                    location = 'chromosome'
                    to_output.write("\t".join([
                        protein_id, 
                        contig_id,
                        description,
                        location,
                    ])+ "\n")
            

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
    args = parser.parse_args()

    ### Calling functions
    contig_names = names_parser( args.contigs_map )
    amr_data = arg_parser( args.amr_out, contig_names )
    ( mob_prots, mob_coords, mob_types ) = mob_parser( args.mobilome )
    location_parser( amr_data, mob_prots, mob_coords, mob_types )


if __name__ == "__main__":
    main()
