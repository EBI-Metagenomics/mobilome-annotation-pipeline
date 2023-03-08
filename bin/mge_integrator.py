#!/usr/bin/env python
import warnings
warnings.filterwarnings("ignore", category=DeprecationWarning)
warnings.filterwarnings("ignore", category=FutureWarning)

from Bio import SeqIO
import argparse
import sys
import os.path
import glob


##### This script integrates and parse the output of ICEfinder, IntegronFinder, ISEScan, and PaliDIS for MoMofy
##### Alejandra Escobar, EMBL-EBI
##### January 11, 2023

parser = argparse.ArgumentParser(
        description='This script integrates and parse the output of ICEfinder, IntegronFinder, ISEScan, and PaliDIS for MoMofy. Please provide the relevant input files')
parser.add_argument('--cgc_fa', type=str, help='CDS prediction fasta file')
parser.add_argument('--map', type=str, help='Rename contigs file: contigID.map')
parser.add_argument('--iss_fa', type=str, help='ISEScan fasta file')
parser.add_argument('--iss_tsv', type=str, help='ISEScan predictions table')
parser.add_argument('--pal_fa', type=str, help='PaliDIS fasta file')
parser.add_argument('--pal_tsv', type=str, help='PaliDIS predictions table')
parser.add_argument('--inf_tsv', type=str, help='IntegronFinder predictions table')
parser.add_argument('--inf_gbks', nargs='*', help='Space separated list of IntegonFinder gbk files per contig')
parser.add_argument('--icf_tsv', type=str, help='ICEfinder prediction files (concatenated)')
parser.add_argument('--icf_fa', type=str, help='ICEfinder fasta files (concatenated)')
parser.add_argument('--mog_tsv', type=str, help='Diamond output versus MobileOG-DB format 6')
args = parser.parse_args()

### For debugging
#to_test=open('test.out','w')

### Setting up variables
cgc_seqs=args.cgc_fa
map_file=args.map
iss_seqs=args.iss_fa
iss_results=args.iss_tsv
pal_seqs=args.pal_fa
pal_results=args.pal_tsv
integron_results=args.inf_tsv
inf_gbks=args.inf_gbks
icf_results=args.icf_tsv
icf_seqs=args.icf_fa
mog_results=args.mog_tsv

mge_data={}
mge_nuc={}
mge_counter=0


### Saving contig names equivalence
names_equiv={}
if os.stat(map_file).st_size > 0:
    with open(map_file,'r') as input_map:
        for line in input_map:
            new_name,old_name=line.rstrip().split('\t')
            names_equiv[new_name.replace('>','')]=old_name
    inv_names_equiv = {v: k for k, v in names_equiv.items()}


### Saving the ICEfinder sequences
icf_nuc={}
if os.stat(icf_seqs).st_size > 0:
    for record in SeqIO.parse(icf_seqs, "fasta"):
        my_chain=str(record.seq).upper()
        my_desc=str(record.description)
        my_id=my_desc.split(' ')[0]+'|'+my_desc.split(' ')[5]
        icf_nuc[my_id]=my_chain


    ### Parsing ICEfinder summaries
    with open(icf_results,'r') as input_table:
        for line in input_table:
            mge_counter+=1
            mge_id='icf_'+str(mge_counter)
            ICEfinder_Job_id,Strain,Genome_len,ICEfinder_output_name,Description,Coordinate,Length,oriT,GC,Genome_GC,Delta_GC,ARG,VF=line.rstrip().split('\t')
            contig=ICEfinder_Job_id.replace('.prokka','')
            Description=Description.replace(' ','_').replace('Putative_','').replace('_AICE','AICE').replace(':','')
            start=int(Coordinate.split('..')[0])
            end=int(Coordinate.split('..')[1])
            coord=(start,end)
            value=(contig,Description,coord)
            mge_data[mge_id]=value

            seq_id=contig+'|'+Coordinate
            if seq_id in icf_nuc.keys():
                my_id='>'+mge_id+'|'+contig+'|'+str(start)+':'+str(end)+'|'+Description
                mge_nuc[my_id]=icf_nuc[seq_id]


### Parsing IntegronFinder output
if os.stat(integron_results).st_size > 0:
    with open(integron_results,'r') as input_table:
        next(input_table)
        next(input_table)
        for line in input_table:
            ID_replicon,CALIN,complete,In0,topology,size=line.rstrip().split('\t')
            if int(complete)>0:
                description='Complete_integron'
                gbk_file=ID_replicon+'.gbk'
                if gbk_file in inf_gbks:
                    for gb_record in SeqIO.parse(gbk_file, "genbank"):
                        for feature in gb_record.features:
                            if feature.type=='integron':
                                mge_counter+=1
                                mge_id='inf_'+str(mge_counter)
                                start=int(feature.location.start)
                                end=int(feature.location.end)
                                coord=(start,end)
                                value=(ID_replicon,description,coord)
                                mge_data[mge_id]=value
                                my_id='>'+mge_id+'|'+ID_replicon+'|'+str(start)+':'+str(end)+'|'+description
                                my_chain=str(feature.extract(gb_record.seq)).upper()
                                mge_nuc[my_id]=my_chain


### Parsing ISEScan outputs
raw_iss={}
if os.stat(iss_results).st_size > 0:
    for record in SeqIO.parse(iss_seqs, "fasta"):
        my_chain=str(record.seq).upper()
        raw_iss[str(record.id)]=my_chain

    with open(iss_results,'r') as input_table:
        next(input_table)
        for line in input_table:
            mge_counter+=1
            mge_id='iss_'+str(mge_counter)
            contig,family,cluster,isBegin,isEnd,isLen,ncopy4is,start1,end1,start2,end2,score,irId,irLen,nGaps,orfBegin,orfEnd,strand,orfLen,E_value,E_value4copy,iss_type,ov,tir=line.rstrip().split('\t')

            if iss_type == 'c':
                if tir == '-:-':
                    description='without_TIR'
                else:
                    description='with_TIR'
                description=family+'_'+description
                coord=(int(isBegin),int(isEnd))
                value=(contig,description,coord)
                mge_data[mge_id]=value
                fasta_id=contig+'_'+isBegin+'_'+isEnd+'_'+strand
                if fasta_id in raw_iss.keys():
                    new_fasta_id='>'+mge_id+'|'+contig+'|'+isBegin+':'+isEnd+'|'+description
                    mge_nuc[new_fasta_id]=raw_iss[fasta_id]
                else:
                    print('No fasta sequence for: '+fasta_id)


### Parsing Palidis output
raw_pal={}
if os.path.isfile(pal_seqs):
    if os.stat(pal_seqs).st_size > 0:
        for record in SeqIO.parse(pal_seqs, "fasta"):
            my_chain=str(record.seq).upper()
            raw_pal[str(record.id)]=my_chain

        with open(pal_results,'r') as input_table:
            next(input_table)
            for line in input_table:
                mge_counter+=1
                mge_id='pal_'+str(mge_counter)
                IS_name,sample_id,contig,itr1_start_position,itr1_end_position,itr2_start_position,itr2_end_position,description=line.rstrip().split('\t')
                description='IS_with_TIR'
                coord=(int(itr1_start_position),int(itr2_end_position))
                value=(contig,description,coord)
                mge_data[mge_id]=value
    
                if IS_name in raw_pal.keys():
                    new_fasta_id='>'+mge_id+'|'+contig+'|'+itr1_start_position+':'+itr2_end_position+'|'+description
                    mge_nuc[new_fasta_id]=raw_pal[IS_name]
                else:
                    print('No fasta file for: '+IS_name)


### Cleaning redundant predictions
# Grouping elements by contig
contig_ismge={}
contig_inmge={}
for prediction in mge_data.keys():
    if any([ 'iss' in prediction , 'pal' in prediction ]):
        contig=mge_data[prediction][0]
        if contig in contig_ismge.keys():
            contig_ismge[contig].append(prediction)
        else:
            contig_ismge[contig]=[prediction]
    elif any([ 'icf' in prediction , 'inf' in prediction ]):
        contig=mge_data[prediction][0]
        if contig in contig_inmge.keys():
            contig_inmge[contig].append(prediction)
        else:
            contig_inmge[contig]=[prediction]


# Parsing insertion sequences overlapps
to_remove=[]
for contig in contig_ismge.keys():
    if len(contig_ismge[contig])>1:
        tools_list=[]
        for element in contig_ismge[contig]:
            prefix=element.split('_')[0]
            tools_list.append(prefix)
        tools_list=list(set(tools_list))
        if len(tools_list)>1:
            iss_list=[x for x in contig_ismge[contig] if 'iss' in x]
            pal_list=[x for x in contig_ismge[contig] if 'pal' in x]
            #print(contig)
            for pal_element in pal_list:
                pal_coord=mge_data[pal_element][2]
                pal_len=pal_coord[1]-pal_coord[0]
                for iss_element in iss_list:
                    iss_coord=mge_data[iss_element][2]
                    iss_len=iss_coord[1]-iss_coord[0]       
                    intersection=len(list(set(range(pal_coord[0], pal_coord[1]+1)) & set(range(iss_coord[0], iss_coord[1]+1))))
                    if intersection>0:
                        iss_cov=float(intersection)/float(iss_len)
                        pal_cov=float(intersection)/float(pal_len)
                        if any([ iss_cov>0.1 , pal_cov>0.1 ]):
                            to_remove.append(iss_element)

to_remove=list(set(to_remove))

# Printing discarded redundant predictions to file
output_discard='discarded_mge.txt'
with open(output_discard, 'w') as to_discard:
    to_discard.write('#element_id\tdescription\tcause\n')
    for element in to_remove:
        to_discard.write(element+'\t'+mge_data[element][1]+'\toverlapping\n')

# Removing redundant IS
for element in to_remove:
    nuc_llave=['>'+element,mge_data[element][0],str(mge_data[element][2][0])+':'+str(mge_data[element][2][1]),mge_data[element][1]]
    nuc_llave='|'.join(nuc_llave)
    del mge_data[element]
    del mge_nuc[nuc_llave]

# Parsing integron overlapps. Reporting only, not to be removed
output_nested='nested_integrons.txt'
with open(output_nested, 'w') as to_nested:
    to_nested.write('#contig\ticefinder_element\tintfinder_element\ticefinder_cov\tintfinder_cov\n')
    for contig in contig_inmge.keys():
        if len(contig_inmge[contig])>1:
            tools_list=[]
            for element in contig_inmge[contig]:
                prefix=element.split('_')[0]
                tools_list.append(prefix)
            tools_list=list(set(tools_list))
            if len(tools_list)>1:
                inf_list=[x for x in contig_inmge[contig] if 'inf' in x]
                icf_list=[x for x in contig_inmge[contig] if 'icf' in x]
                #print(contig)
                for inf_element in inf_list:
                    inf_coord=mge_data[inf_element][2]
                    inf_len=inf_coord[1]-inf_coord[0]
                    for icf_element in icf_list:
                        icf_coord=mge_data[icf_element][2]
                        icf_len=icf_coord[1]-icf_coord[0]
                        intersection=len(list(set(range(inf_coord[0], inf_coord[1]+1)) & set(range(icf_coord[0], icf_coord[1]+1))))
                        if intersection>0:
                            inf_cov=float(intersection)/float(inf_len)
                            icf_cov=float(intersection)/float(icf_len)
                            if any([ inf_cov>0.75 , icf_cov>0.75 ]):
                                to_nested.write(contig+'\t'+icf_element+'\t'+inf_element+'\t'+str(icf_cov)+'\t'+str(inf_cov)+'\n')


### Saving the CDS' coordinates and contigs location
contig_prots={}
prots_coord={}
gff_simplecontig_realcontig={}
if os.stat(cgc_seqs).st_size > 0:
    for record in SeqIO.parse(cgc_seqs, "fasta"):
        seq_id=str(record.id).split('-')
        seq_id.pop(0)
        seq_id='-'.join(seq_id)
        contig=seq_id.split('_')[0].replace('-','_').replace('.','_')
        contigtogff=str(record.id).split('_')[0]
        gff_simplecontig_realcontig[contig]=contigtogff

        protein_id=str(record.description)
        if ' ' in protein_id:
            start=int(protein_id.split(' ')[2])
            end=int(protein_id.split(' ')[4])
            strand=int(protein_id.split(' ')[6])
            if strand=='-1':
                strand='-'
            else:
                strand='+'
        else:
            start=int(protein_id.split('_')[1])
            end=int(protein_id.split('_')[2])
            strand=protein_id.split('_')[3]
    
        if contig not in contig_prots.keys():
            contig_prots[contig]=[protein_id]
        else:
            contig_prots[contig].append(protein_id)

        prots_coord[protein_id]=(start,end,strand)


### Finding the CDS encoded in the mobilome
mge_proteins={}
for element in mge_data.keys():
    mge_proteins[element]=[]
    contig=names_equiv[mge_data[element][0]].replace('.','_')
    mge_start=mge_data[element][2][0]
    mge_end=mge_data[element][2][1]
    mge_range=range(mge_start,mge_end+1)
    mge_len=mge_end-mge_start

    if contig in contig_prots.keys():
        for protein in contig_prots[contig]:
            prot_start=prots_coord[protein][0]
            prot_end=prots_coord[protein][1]
            prot_range=range(prot_start,prot_end+1)
            prot_len=prot_end-prot_start

            intersection=len(list(set(mge_range) & set(prot_range)))
            if intersection>0:
                mge_cov=float(intersection)/float(mge_len)
                prot_cov=float(intersection)/float(prot_len)
                if any([ mge_cov>0.75 , prot_cov>0.75 ]):
                    #if all(item in mge_range for item in prot_range):
                    mge_proteins[element].append(protein)

### Parsing the mobileOG annotation file
mog_annot={}
if os.stat(mog_results).st_size > 0:
    with open(mog_results) as input_table:
        for line in input_table:
            line_l=line.rstrip().split('\t')
            target_id=line_l[0]
            query_id=line_l[1]
            mog_id,gene_name,best_hit_id,major,minor,db,evidence=target_id.split('|')
            function=mog_id+'|'+major+'|'+minor

            if query_id not in mog_annot.keys():
                mog_annot[query_id]=function.replace(',','/')
                #to_test.write(function+'\t'+query_id+'\n')


### Writing the mobilome outputs
prefixes={
    'icf':'ICEfinder',
    'inf':'IntegronFinder',
    'iss':'ISEScan',
    'pal':'PaliDIS'
}

mgecontigs=[]
no_cds=[]
len_500=[]
for element in mge_data.keys():
    element_len=int(mge_data[element][2][1])-int(mge_data[element][2][0])
    if element_len<500:
        len_500.append(element)
    elif len(mge_proteins[element])==0:
        no_cds.append(element)
    else:
        contig=names_equiv[mge_data[element][0]]
        contig_key=contig.replace('.','_')
        contig_cds=gff_simplecontig_realcontig[contig_key]
        if contig_cds not in mgecontigs:
            mgecontigs.append(contig_cds)

# Removing predictions of len<500 and with no CDS
with open(output_discard, 'a') as to_discard:
    for element in no_cds:
        to_discard.write(element+'\t'+mge_data[element][1]+'\tno_cds\n')
        del mge_data[element]

    for element in len_500:
        to_discard.write(element+'\t'+mge_data[element][1]+'\tmge<500bp\n')
        del mge_data[element]


used_proteins=[]
output_fna='momofy_predictions.fna'
output_tsv='momofy_predictions.gff'
with open(output_fna, 'w') as to_fasta, open(output_tsv, 'w') as to_tsv:
    to_tsv.write('##gff-version 3\n')
    for contig in mgecontigs: 
        contig_len=contig.split('-')[4]
        to_tsv.write('##sequence-region '+contig+' 1 '+contig_len+'\n')

    for element in mge_data.keys():
        seqid=gff_simplecontig_realcontig[names_equiv[mge_data[element][0]].replace('.','_')]
        source=prefixes[element.split('_')[0]]

        if any( [ 'iss' in element , 'pal' in element ] ): 
            seq_type='insertion_sequence'
        elif 'inf' in element:
            seq_type='integron'
        else:
            if 'IME' in mge_data[element][1]:
                seq_type='integron'       

            elif 'ICE' in mge_data[element][1]:
                seq_type='conjugative_transposon'
            else:
                seq_type='gene_group'

        start=str(mge_data[element][2][0])
        end=str(mge_data[element][2][1])
        score='.'
        strand='.'
        phase='.'
        attributes='ID='+element+';gbkey=mobile_element;mobile_element_type='+mge_data[element][1]

        tsv_line=[seqid,source,seq_type,start,end,score,strand,phase,attributes]
        tsv_line='\t'.join(tsv_line)
        to_tsv.write(tsv_line+'\n')

        seqid_saved=['>'+element,mge_data[element][0],str(mge_data[element][2][0])+':'+str(mge_data[element][2][1]),mge_data[element][1]]
        seqid_saved='|'.join(seqid_saved)
        my_sequence=mge_nuc[seqid_saved]
        nuc_id=['>'+element,seqid,str(mge_data[element][2][0])+'..'+str(mge_data[element][2][1]),mge_data[element][1]]
        nuc_id='|'.join(nuc_id)
        to_fasta.write(nuc_id+'\n')
        to_fasta.write(my_sequence+'\n')

        for protein in mge_proteins[element]:
            if protein in mog_annot.keys():
                function=mog_annot[protein].replace(' ','_')
                if protein not in used_proteins:
                    used_proteins.append(protein)
                    if ' ' in protein:
                        ID=protein.split(' ')[0]
                    else:
                        ID=protein
                    source='CGC_V5'
                    seq_type='CDS'
                    start=str(prots_coord[protein][0])
                    end=str(prots_coord[protein][1])
                    score='.'
                    strand=prots_coord[protein][2]
                    phase='0'
                    attributes='ID='+ID+';gbkey=CDS;product='+function

                    tsv_line=[seqid,source,seq_type,start,end,score,strand,phase,attributes]
                    tsv_line='\t'.join(tsv_line)
                    to_tsv.write(tsv_line+'\n')


