#!/usr/bin/env python

import re
import argparse
import sys
import os.path


##### This script integrates MOMOfy, proMGE, and alien_hunter predictions to report agreement score per gene
##### Alejandra Escobar, EMBL-EBI
##### Apr 24, 2023

def gff_parser(current_line):
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
    ) = current_line.split("\t")
    seq_id = attr.split(';')[0].replace('ID=','')
    data_list = [contig, int(start), int(end), seq_type, seq_id]
    return( data_list )


def momo_parser(momofy):
    ### Parsing MoMofy gff file
    momofy_dict = {}
    genes_dict = {}
    genes_coord = {}
    mges = [
        "insertion_sequence",
        "integron",
        "conjugative_transposon",
        "plasmid",
        "viral_sequence",
        "prophage",
    ]
    flanking = [
        "terminal_inverted_repeat_element",
        "attC_site",
        "direct_repeat",
    ]
    dummy_list = []
    with open(momofy, 'r') as input_file:
        for line in input_file:
            line = line.rstrip()
            # Annotation lines have exactly 9 columns
            if len(line.split("\t")) == 9:
                data_list = gff_parser( line )
                contig = data_list.pop(0)
                seq_id = data_list.pop(-1)
                data_list = tuple(data_list)
                if data_list[2] in mges:
                    if contig in momofy_dict:
                        momofy_dict[contig].append(data_list)
                    else:
                        momofy_dict[contig]=[data_list]
                elif data_list[2] in flanking:
                    dummy_list.append(contig)
                else:
                    genes_coord[seq_id] = data_list
                    if contig in genes_dict:
                        genes_dict[contig].append(seq_id)
                    else:
                        genes_dict[contig]=[seq_id]
    return(genes_dict, momofy_dict, genes_coord)


def promge_parser(promge, meta):
    mge_type=[
        "is_tn",
        "phage",
        "phage_like",
        "ce",
        "integron",
        "mi",
        "cellular",
    ]
    mge_desc={}
    with open(meta, 'r') as input_meta:
        next(input_meta)
        for line in input_meta:
            is_tn,phage,phage_like,ce,integron,mi,cellular,contig,start,end,size,n_genes,mgeR = line.rstrip().split("\t")
            new_id = contig+':'+start+'-'+end
            description = []
            for index in range(7):
                if int(line.rstrip().split("\t")[index])>0:
                    description.append(mge_type[index])
            description='|'.join(description)
            mge_desc[new_id]=description

    promge_dict = {}
    with open(promge, 'r') as input_gff:
        for line in input_gff:
            line = line.rstrip()
            # Annotation lines have exactly 9 columns
            if len(line.split("\t")) == 9:
                data_list = gff_parser( line )
                contig = data_list.pop(0)
                seq_id = data_list.pop(-1)
                new_id = re.sub(".*contig_", "contig_", seq_id)

                if mge_desc[new_id] != 'cellular':
                    if contig in promge_dict:
                        promge_dict[contig].append(data_list)
                    else:
                        promge_dict[contig]=[data_list]
    return(promge_dict)


def alien_parser(alien):
    alien_dict = {}
    with open(alien, 'r') as input_txt:
        next(input_txt)
        for line in input_txt:
            line_l = line.rstrip().split('\t')
            contig = line_l[2]
            start = int(line_l[5].split('..')[0])
            end = int(line_l[5].split('..')[1])
            coords = (start, end, 'CO')
            score = float(line_l[8])

            if score >= 20:
                if contig in alien_dict:
                    alien_dict[contig].append(coords)
                else:
                    alien_dict[contig]=[coords]
    return(alien_dict)


def blast_parser(mobog):
    mobog_list = []
    with open(mobog, 'r') as input_txt:
        for line in input_txt:
            line_l = line.rstrip().split('\t')
            gene_id = line_l[1].split(' ')[0]
            mobog_list.append(gene_id)
    mobog_list = list(set(mobog_list))
    return(mobog_list)

def mapper(gene_start, gene_end, mge_start, mge_end):
    match = '0'
    gene_len = gene_end - gene_start
    gene_range = range(gene_start , gene_end+1)
    mge_range = range(mge_start, mge_end+1)
    intersection = len(list(set(gene_range) & set(mge_range)))
    if intersection > 0:
        gene_cov = float(intersection) / float(gene_len)
        if gene_cov >= 0.75:
            match = '1'
    return(match)


def integrator(genes_dict, genes_coord, momofy_dict, promge_dict, alien_dict, mobog_dict):
    agreement = {
        0: 'Chromosome',
        1: 'Low', 
        2: 'Medium_low',
        3: 'Medium_high',
        4: 'High',
    }
    with open('per_gene_mobilome.tsv', 'w') as to_pergene:
        to_pergene.write(
            "\t".join(
                [
                    'contig',
                    'gene_id',
                    'start',
                    'end',
                    'type',
                    'mobileOG',
                    'momofy',
                    'promge',
                    'alien_hunter',
                    'score',
                    'agreement',
                ]
            ) + '\n'
        )

        for contig in genes_dict:
            for gene in genes_dict[contig]:
                to_print = [
                    contig,
                    gene,
                    str(genes_coord[gene][0]),
                    str(genes_coord[gene][1]),
                    genes_coord[gene][2]
                ]

                # Checking for matches in mobileOG-DB
                score = 0
                if gene in mobog_dict:
                    to_print.append('1')
                    score+=1
                else:
                    to_print.append('0')

                # Checking coverage in momofy, promge, and alien_hunter predictions
                gene_start = genes_coord[gene][0]
                gene_end = genes_coord[gene][1]

                # MoMofy
                if contig in momofy_dict:
                    for element in momofy_dict[contig]:
                        mge_start = element[0]
                        mge_end = element[1]
                        match = mapper(gene_start, gene_end, mge_start, mge_end)
                        if match == '1':
                            break
                    if match == '1':
                        score+=1
                    to_print.append(match)
                else:
                    to_print.append('0')

                # proMGE
                if contig in promge_dict:
                    for element in promge_dict[contig]:
                        mge_start = element[0]
                        mge_end = element[1]
                        match = mapper(gene_start, gene_end, mge_start, mge_end)
                        if match == '1':
                            break
                    if match == '1':
                        score+=1
                    to_print.append(match)
                else:
                    to_print.append('0')
                    
                # alie_hunter
                if contig in alien_dict:
                    for element in alien_dict[contig]:
                        mge_start = element[0]
                        mge_end = element[1]
                        match = mapper(gene_start, gene_end, mge_start, mge_end)
                        if match == '1':
                            break
                    if match == '1':
                        score+=1
                    to_print.append(match)
                else:
                    to_print.append('0')

                to_print.append(str(score))
                to_print.append(agreement[score])
                to_pergene.write("\t".join(to_print)+'\n')

def main():
    parser = argparse.ArgumentParser(
        description="This script integrates MOMOfy, proMGE, and alien_hunter predictions to report agreement score per gene"
    )
    parser.add_argument(
        "--momofy",
        type=str,
        help="MoMofy+Virify mobilome annotation file (mobilome_predictions.gff)",
        required=True,
    )
    parser.add_argument(
        "--proMGE",
        type=str,
        help="ProMGE annotation file (gff)",
        required=True,
    )
    parser.add_argument(
        "--meta",
        type=str,
        help="proMGE metadata annotation file (txt)",
        required=True,
    )
    parser.add_argument(
        "--alien",
        type=str,
        help="corrected alien_hunter output (tsv)",
        required=True,
    )
    parser.add_argument(
        "--mobileOG",
        type=str,
        help="Diamond output vs mobileOG-DB (tsv)",
        required=True,
    )


    args = parser.parse_args()

    ### Setting up variables
    momofy = args.momofy
    promge = args.proMGE
    meta = args.meta
    alien = args.alien
    mobog = args.mobileOG

    (genes_dict, momofy_dict, genes_coord) = momo_parser(momofy)
    promge_dict = promge_parser(promge, meta)
    alien_dict = alien_parser(alien)
    mobog_list = blast_parser(mobog)
    integrator(genes_dict, genes_coord, momofy_dict, promge_dict, alien_dict, mobog_list)

if __name__ == "__main__":
    main()

