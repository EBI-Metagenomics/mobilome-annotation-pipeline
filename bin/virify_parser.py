#!/usr/bin/env python

from Bio import SeqIO
import argparse
import sys
import os.path
import glob

##### This script parse the VIRify v2.0 outputs to be integrated to MoMofy results
##### Alejandra Escobar, EMBL-EBI
##### Aug 14, 2023


def names_map(ori_fasta, ren_fasta):
    ### Creating an ID conversion dict
    originals = []
    if os.stat(ori_fasta).st_size > 0:
        for record in SeqIO.parse(ori_fasta, "fasta"):
            originals.append(str(record.id))

    ren_fasta = ren_fasta[0]
    renamed = []
    if os.stat(ren_fasta).st_size > 0:
        for record in SeqIO.parse(ren_fasta, "fasta"):
            renamed.append(str(record.id))
    contig_names = dict(zip(renamed, originals))

    return(contig_names)


def pprmeta_parser(pprm, contig_names):
    ### Saving ppr-meta predictions
    pprm = pprm[0]
    with open('plasmids.list', "w") as list_out:
        if os.path.exists(pprm):
            with open(pprm, "r") as input_csv:
                next(input_csv)
                for line in input_csv:
                    (
                        header,
                        length,
                        phage_score,
                        chromosome_score,
                        plasmid_score,
                        possible_source,
                    ) = line.rstrip().split(",")
                    if possible_source == "plasmid":
                        contig = contig_names[header]
                        list_out.write(contig+'\t'+plasmid_score+'\n')



def phage_seq_save(viri_fa):
    with open ('all_virify.fasta','w') as fasta_out:
        for fasta_file in viri_fa:
            if os.path.exists(fasta_file):
                head, tail = os.path.split(fasta_file)
                quality = tail.split('_')[0]
                for record in SeqIO.parse(fasta_file, "fasta"):
                    phage_ID = str(record.description)
                    if 'prophage' in phage_ID:
                        phage_ID = phage_ID.replace(' ','|').replace('prophage-0:','prophage-1:')
                    else:
                        phage_ID = str(record.id)
                    phage_ID = phage_ID + ' ' + quality
                    fasta_out.write('>'+phage_ID+'\n')
                    fasta_out.write(str(record.seq)+'\n')



def main():
    parser = argparse.ArgumentParser(
        description="This script parse the VIRify v2.0 outputs to be integrated to MoMofy results. Please provide the full path of the GFF output file of VIRify"
    )
    parser.add_argument(
        "--virify_gff", 
        type=str, 
        help="Virify output in gff format"
    )
    parser.add_argument(
        "--original_assem",
        type=str, 
        help="Original assembly fasta file"
    )
    args = parser.parse_args()


    ### Setting up the paths of all the input files
    viri = args.virify_gff
    ori_fasta = args.original_assem

    path_prefix = args.virify_gff.split('/')[:-3]
    path_prefix = '/'.join(path_prefix)

    viri_fa = glob.glob(path_prefix + '/*_original.fasta')
    pprm = glob.glob(path_prefix + '/01-viruses/pprmeta/*_pprmeta.csv')
    ren_fasta = glob.glob(path_prefix + '/*_renamed.fasta' )

    ### Calling functions
    contig_names = names_map(ori_fasta, ren_fasta)
    pprmeta_parser(pprm, contig_names)
    phage_seq_save(viri_fa)


if __name__ == "__main__":
    main()

