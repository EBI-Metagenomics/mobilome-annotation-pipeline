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
parser.add_argument('--assem', type=str, help='Original input assembly')
parser.add_argument('--cds_gff', type=str, help='GFF prediction fasta file')
parser.add_argument('--map', type=str, help='Rename contigs file: contigID.map')
parser.add_argument('--iss_fa', type=str, help='ISEScan fasta file')
parser.add_argument('--iss_tsv', type=str, help='ISEScan predictions table')
parser.add_argument('--pal_fa', type=str, help='PaliDIS fasta file')
parser.add_argument('--pal_tsv', type=str, help='PaliDIS predictions table')
parser.add_argument('--inf_tsv', type=str, help='IntegronFinder predictions table')
parser.add_argument('--inf_gbks', nargs='*', help='Space separated list of IntegonFinder gbk files per contig')
parser.add_argument('--icf_tsv', type=str, help='ICEfinder prediction files (concatenated)')
parser.add_argument('--icf_fa', type=str, help='ICEfinder fasta files (concatenated)')
parser.add_argument('--icf_lim', type=str, help='ICEfinder DR coordinates')
parser.add_argument('--mog_tsv', type=str, help='Diamond output versus MobileOG-DB format 6')
args = parser.parse_args()

### For debugging
to_test=open('test.out','w')

### Setting up variables
ori_contigs=args.assem
cds_loc=args.cds_gff
map_file=args.map
iss_seqs=args.iss_fa
iss_results=args.iss_tsv
pal_seqs=args.pal_fa
pal_results=args.pal_tsv
integron_results=args.inf_tsv
inf_gbks=args.inf_gbks
icf_results=args.icf_tsv
icf_seqs=args.icf_fa
icf_dr_file=args.icf_lim
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


### Saving contigs len
assembly_len={}
if os.stat(ori_contigs).st_size > 0:
    for record in SeqIO.parse(ori_contigs, "fasta"):
        chain_len=len(str(record.seq))
        seq_id=str(record.id)
        assembly_len[seq_id]=str(chain_len)

### Saving the ICEfinder sequences
icf_nuc={}
icf_dr={}
if os.stat(icf_seqs).st_size > 0:
    for record in SeqIO.parse(icf_seqs, "fasta"):
        my_chain=str(record.seq).upper()
        my_desc=str(record.description)
        my_id=my_desc.split(' ')[0]+'|'+my_desc.split(' ')[5]
        icf_nuc[my_id]=my_chain

    ### Saving ICEfinder direct repeat coordinates
    with open(icf_dr_file,'r') as input_table:
        for line in input_table:
            line_l=line.rstrip().split()
            if len(line_l)==3:
                new_id=line_l[0].replace('result/','').replace(':DR:','').replace('/','|')
                dr_1_s=line_l[1].split('..')[0]
                dr_1_e=line_l[1].split('..')[1]
                dr_1=(dr_1_s,dr_1_e)
                dr_2_s=line_l[2].split('..')[0]
                dr_2_e=line_l[2].split('..')[1]
                dr_2=(dr_2_s,dr_2_e)
                icf_dr[new_id]=(dr_1,dr_2)

    ### Parsing ICEfinder summaries
    with open(icf_results,'r') as input_table:
        for line in input_table:
            ICEfinder_Job_id,Strain,Genome_len,ICEfinder_output_name,Description,Coordinate,Length,oriT,GC,Genome_GC,Delta_GC,ARG,VF=line.rstrip().split('\t')
            contig=ICEfinder_Job_id
            Description=Description.replace(' ','_').replace('Putative_','').replace('_AICE','AICE').replace(':','')

            if not 'conjugative_region' in Description:
                mge_counter+=1
                mge_id='icf_'+str(mge_counter)
                start=int(Coordinate.split('..')[0])
                end=int(Coordinate.split('..')[1])
                coord=(start,end)
                value=(contig,Description,coord)
                mge_data[mge_id]=value

                seq_id=contig+'|'+Coordinate
                if seq_id in icf_nuc.keys():
                    my_id='>'+mge_id+'|'+contig+'|'+str(start)+':'+str(end)+'|'+Description
                    mge_nuc[my_id]=icf_nuc[seq_id]

                composite_id=contig+'|'+ICEfinder_output_name
                if composite_id in icf_dr.keys():
                    icf_dr[mge_id]=icf_dr.pop(composite_id)


### Parsing IntegronFinder output
attC_site={}
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
                    flag=0
                    for gb_record in SeqIO.parse(gbk_file, "genbank"):
                        for feature in gb_record.features:
                            if feature.type == 'integron':
                                if feature.qualifiers['integron_type'][0] == 'complete':
                                    flag=1
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
                            if feature.type == 'attC':
                                if flag == 1:
                                    start=str(feature.location.start)
                                    end=str(feature.location.end)
                                    coord=(start,end)
                                    attC_site[mge_id]=coord
                                    flag=0

### Parsing ISEScan outputs
raw_iss={}
itr_sites={}
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
                if int(irLen) == 0:
                    description='without_TIR'
                else:
                    description='with_TIR'
                    ir_1=(start1,end1)
                    ir_2=(start2,end2)
                    itr_sites[mge_id]=(ir_1,ir_2)

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
                contig=inv_names_equiv[contig]
                value=(contig,description,coord)
                mge_data[mge_id]=value

                ir_1=(itr1_start_position,itr1_end_position)
                ir_2=(itr2_start_position,itr2_end_position)
                itr_sites[mge_id]=(ir_1,ir_2)
    
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
if os.stat(cds_loc).st_size > 0:
    with open(cds_loc,'r') as input_table:
        for line in input_table:
            l_line=line.rstrip().split('\t')
            if len(l_line)==9:
                if l_line[2]=='CDS':
                    contig=l_line[0]
                    prot_source=l_line[1]
                    start=int(l_line[3])
                    end=int(l_line[4])
                    strand=l_line[6]
                    attrib=l_line[8]
                    prot_id=attrib.split(';')[0]
                    if prot_id.startswith('ID='):
                        prot_id=prot_id.replace('ID=','')
                        value=(start,end,strand)
                        prots_coord[prot_id]=value
                        if contig not in contig_prots.keys():
                            contig_prots[contig]=[prot_id]
                        else:
                            contig_prots[contig].append(prot_id)
                    else:
                        print('CDS ID is expected as the first token in the attributes field with the key "ID" in line:\n'+line)


### Finding the CDS encoded in the mobilome
mge_proteins={}
mob_proteome=[]
for element in mge_data.keys():
    mge_proteins[element]=[]
    contig=names_equiv[mge_data[element][0]]
    mge_start=mge_data[element][2][0]
    mge_end=mge_data[element][2][1]
    mge_range=range(mge_start,mge_end+1)
    mge_len=mge_end-mge_start

    flag=0
    if contig in contig_prots.keys():
        flag=1
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
                    mob_proteome.append(protein)

    if flag==0:
        contig=mge_data[element][0]
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
                        mob_proteome.append(protein)

### Parsing the mobileOG annotation file
mog_annot={}
if os.stat(mog_results).st_size > 0:
    with open(mog_results) as input_table:
        for line in input_table:
            line_l=line.rstrip().split('\t')
            target_id=line_l[0]
            query_id=line_l[1].split(' ')[0]
            mog_id,gene_name,best_hit_id,major,minor,db,evidence=target_id.split('|')
            function=mog_id+'|'+major+'|'+minor

            if query_id in mob_proteome:
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

no_cds=[]
len_500=[]
for element in mge_data.keys():
    element_len=int(mge_data[element][2][1])-int(mge_data[element][2][0])
    if element_len<500:
        len_500.append(element)
    elif len(mge_proteins[element])==0:
        no_cds.append(element)

# Removing predictions of len<500 and with no CDS
with open(output_discard, 'a') as to_discard:
    for element in no_cds:
        to_discard.write(element+'\t'+mge_data[element][1]+'\tno_cds\n')
        del mge_data[element]

    for element in len_500:
        to_discard.write(element+'\t'+mge_data[element][1]+'\tmge<500bp\n')
        del mge_data[element]

contigs_elements={}
proteins_mge={}
for element in mge_data.keys():
    contig=mge_data[element][0]
    if contig not in contigs_elements.keys():
        contigs_elements[contig]=[element]
    else:
        contigs_elements[contig].append(element)

    for protein in mge_proteins[element]:
        proteins_mge[protein]=element

def flanking_data(ID_to_print, source, seq_type, start, end, flank_id):
    to_print=[ID_to_print,source,seq_type,start,end,'.','.','.',flank_id]
    to_print='\t'.join(to_print)
    return(to_print)

used_contigs=[]
output_fna='momofy_predictions.fna'
output_gff='momofy_predictions.gff'
with open(cds_loc,'r') as input_table, open(output_fna, 'w') as to_fasta, open(output_gff, 'w') as to_gff:
    for line in input_table:
        l_line=line.rstrip().split('\t')
        if len(l_line)==9:
            seqid=l_line[0]
            # If prokka genes are used, then contig name won't require transformation to find MGEs
            if seqid in names_equiv.keys():
                contig=seqid
                ID_to_print=names_equiv[contig]
            # If contig ID not in names_equiv dir, means that we have user contig IDs and require contig name transformation to find MGEs
            else:
                contig=inv_names_equiv[seqid]
                ID_to_print=seqid

            if not seqid in used_contigs:
                used_contigs.append(seqid)
                if contig in contigs_elements.keys():
                    for element in contigs_elements[contig]:
                        source=prefixes[element.split('_')[0]]
                        if any( [ 'iss' in element , 'pal' in element ] ): 
                            seq_type='insertion_sequence'
                            if element in itr_sites.keys():
                                tir_seq_type='terminal_inverted_repeat_element'
                                tir_1_id='ID=TIR_1:'+element
                                gff_line=flanking_data(ID_to_print,source,tir_seq_type,itr_sites[element][0][0],itr_sites[element][0][1],tir_1_id)
                                to_gff.write(gff_line+'\n')
                                tir_2_id='ID=TIR_2:'+element
                                gff_line=flanking_data(ID_to_print,source,tir_seq_type,itr_sites[element][1][0],itr_sites[element][1][1],tir_2_id)
                                to_gff.write(gff_line+'\n')
                        elif 'inf' in element:
                            seq_type='integron'
                            attc_seq_type='attC_site'
                            attc_id='ID=attC:'+element
                            gff_line=flanking_data(ID_to_print,source,attc_seq_type,attC_site[element][0],attC_site[element][1],attc_id)
                            to_gff.write(gff_line+'\n')
                        else:
                            if 'IME' in mge_data[element][1]:
                                seq_type='integron'
                                if element in icf_dr.keys():
                                    dr_seq_type='direct_repeat'
                                    dr_1_id='ID=DR_1:'+element
                                    gff_line=flanking_data(ID_to_print,source,dr_seq_type,icf_dr[element][0][0],icf_dr[element][0][1],dr_1_id)
                                    to_gff.write(gff_line+'\n')
                                    dr_2_id='ID=DR_2:'+element
                                    gff_line=flanking_data(ID_to_print,source,dr_seq_type,icf_dr[element][1][0],icf_dr[element][1][1],dr_2_id)
                                    to_gff.write(gff_line+'\n')
                            elif 'ICE' in mge_data[element][1]:
                                seq_type='conjugative_transposon'
                                if element in icf_dr.keys():
                                    dr_seq_type='direct_repeat'
                                    dr_1_id='ID=DR_1:'+element
                                    gff_line=flanking_data(ID_to_print,source,dr_seq_type,icf_dr[element][0][0],icf_dr[element][0][1],dr_1_id)
                                    to_gff.write(gff_line+'\n')
                                    dr_2_id='ID=DR_2:'+element
                                    gff_line=flanking_data(ID_to_print,source,dr_seq_type,icf_dr[element][1][0],icf_dr[element][1][1],dr_2_id)
                                    to_gff.write(gff_line+'\n')
                        start=str(mge_data[element][2][0])
                        end=str(mge_data[element][2][1])
                        score='.'
                        strand='.'
                        phase='.'
                        attributes='ID='+element+';gbkey=mobile_element;mobile_element_type='+mge_data[element][1]

                        gff_line=[ID_to_print,source,seq_type,start,end,score,strand,phase,attributes]
                        gff_line='\t'.join(gff_line)
                        to_gff.write(gff_line+'\n')

                        seqid_saved=['>'+element,mge_data[element][0],str(mge_data[element][2][0])+':'+str(mge_data[element][2][1]),mge_data[element][1]]
                        seqid_saved='|'.join(seqid_saved)
                        my_sequence=mge_nuc[seqid_saved]
                        nuc_id=['>'+element,seqid,str(mge_data[element][2][0])+'..'+str(mge_data[element][2][1]),mge_data[element][1]]
                        nuc_id='|'.join(nuc_id)
                        to_fasta.write(nuc_id+'\n')
                        to_fasta.write(my_sequence+'\n')

            attrib=l_line[8]
            prot_id=attrib.split(';')[0].replace('ID=','')
            if prot_id in proteins_mge.keys():
                parent_mge=proteins_mge[prot_id]
                if prot_id in mog_annot.keys():
                    function=mog_annot[prot_id].replace(' ','_')
                    attrib=attrib+';mobileOG='+function+';from_mge='+parent_mge
                else:
                    attrib=attrib+';mobileOG=-'+';from_mge='+parent_mge
        
            l_line.pop(-1)
            l_line.pop(0)
            gff_line=[ID_to_print]+l_line+[attrib]
            gff_line='\t'.join(gff_line)
            to_gff.write(gff_line+'\n')

        elif line.startswith('##sequence-region'):
            tag,seqid,start,end=line.rstrip().split()
            # If prokka genes are used, then contig name require transformation to be printed on gff
            if seqid in names_equiv.keys():
                ID_to_print=names_equiv[seqid]
            # If contig ID not in names_equiv dir, means that we have user contig IDs and contig can be printed directly on gff
            else:
                ID_to_print=seqid
            gff_line=[tag,ID_to_print,start,end]
            gff_line=' '.join(gff_line)
            to_gff.write(gff_line+'\n')

        elif line.startswith('>'):
            seqid=line.rstrip().replace('>','')
            # If prokka genes are used, then contig name require transformation to be printed on gff
            if seqid in names_equiv.keys():
                ID_to_print=names_equiv[seqid]
            # If contig ID not in names_equiv dir, means that we have user contig IDs and contig can be printed directly on gff
            else:
                ID_to_print=seqid
            gff_line='>'+ID_to_print
            to_gff.write(gff_line+'\n')

        elif '##gff_version 3' in line:
            to_gff.write('##gff-version 3\n')

        else:
            to_gff.write(line)


'''
    to_tsv.write('##gff-version 3\n')
        contig_len=assembly_len[contig]
        to_tsv.write('##sequence-region '+contig+' 1 '+contig_len+'\n')

            attrib=l_line[8]
            prot_id=attrib.split(';')[0].replace('ID=','')

            if prot_id in mog_annot.keys():
                function=mog_annot[protein].replace(' ','_')
                attrib=attrib+';'+function
            if protein not in used_proteins:
                used_proteins.append(protein)

        for protein in mge_proteins[element]:
            if protein in mog_annot.keys():
                function=mog_annot[protein].replace(' ','_')
                if protein not in used_proteins:
                    used_proteins.append(protein)
                    ID=protein
                    source=prot_source
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

'''
