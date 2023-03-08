#!/usr/bin/env python
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

from Bio import SeqIO
import argparse
import sys
import os.path
import glob


##### This script transfers the functional annotation from the annotation pipeline V5 to the MoMofy output
##### Alejandra Escobar, EMBL-EBI
##### Feb 13, 2023

parser = argparse.ArgumentParser(
        description='This script transfers the functional annotation from the annotation pipeline V5 to the MoMofy output. Please provide the annotation file on gff format, the aminoacid sequences in fasta, and the MoMofy file output on gff format')
parser.add_argument('--func_gff', type=str, help='V5 pipeline functional annotation file (gff)')
parser.add_argument('--aa', type=str, help='The protein sequences in fasta format')
parser.add_argument('--momo', type=str, help='momofy_predictions (gff)')

args = parser.parse_args()


### Setting up variables
annot=args.func_gff
aa_seq=args.aa
momo_results=args.momo


### Saving the aminoacid sequences
protein_seq={}
if os.path.isfile(aa_seq):
    for record in SeqIO.parse(aa_seq, "fasta"):
        my_chain=str(record.seq).upper()
        my_id=str(record.id)
        protein_seq[my_id]=my_chain

### Saving the functional annotation
protein_func={}
if os.path.isfile(annot):
    with open(annot,'r') as input_file:
        for line in input_file:
            if not line.startswith('#'):
                contig,seq_source,seq_type,start,end,score,strand,phase,attr=line.rstrip().split('\t')
                comp_key=(contig,start,end)

                functions=attr.split(';')

                protein_func[comp_key]={
                        'eggnog':'NA',
                        'cog':'NA',
                        'kegg':'NA'
                }
                for element in functions:
                    database=element.split('=')[0]
                    if database=='eggnog':
                        description=element.split('=')[1].replace('%20',' ')
                        if description==' ':
                            protein_func[comp_key]['eggnog']='NA'
                        else:
                            protein_func[comp_key]['eggnog']=element.split('=')[1]
                    elif database=='cog':
                        protein_func[comp_key]['cog']=element.split('=')[1]
                    elif database=='kegg':
                        protein_func[comp_key]['kegg']=element.split('=')[1]
                    
simple_func={}
for prot_coor in protein_func.keys():
    to_print=[]
    for database in protein_func[prot_coor].keys():
        to_print.append(database+'='+protein_func[prot_coor][database])
    to_print=';'.join(to_print)
    simple_func[prot_coor]=to_print

### Parsing the momofy output
with open('MoMofy_results/v5_momofy_predictions.gff','w') as new_annot, open('MoMofy_results/momofy_proteins.fasta','w') as to_fasta:
    if os.path.isfile(momo_results):
        with open(momo_results,'r') as input_table:
            for line in input_table:
                if not line.startswith('#'):
                    contig,seq_source,seq_type,start,end,score,strand,phase,attr=line.rstrip().split('\t')
    
                    if 'product=mobileOG_' in attr:
                        attr=attr.replace('product=mobileOG_','mobileOG=mobileOG_').replace(' ','_')
                    else:
                        attr=attr.replace('product=hypothetical_protein','mobileOG=NA')

                    if seq_type=='CDS':
                        comp_key=(contig,start,end)
                        if comp_key in simple_func.keys():
                            function=simple_func[comp_key]
                        else:
                            function='eggnog=NA;cog=NA;kegg=NA'

                        new_attr=attr+';'+function
                        updated_line=[contig,seq_source,seq_type,start,end,score,strand,phase,new_attr]
                        updated_line='\t'.join(updated_line)
                        new_annot.write(updated_line+'\n')

                        protein_id=new_attr.split(';')[0].replace('ID=','')
                        if protein_id in protein_seq.keys():
                            to_fasta.write('>'+protein_id+'\n')
                            to_fasta.write(protein_seq[protein_id]+'\n')
                        else:
                            print('No sequence for '+protein_id)

                    else:
                        new_annot.write(line)
                else:
                    new_annot.write(line)


