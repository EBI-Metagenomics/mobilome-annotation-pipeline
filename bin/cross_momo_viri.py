#!/usr/bin/env python

from Bio import SeqIO
import argparse
import sys
import os.path
import glob

##### This script integrates the outputs of MoMofy and VIRify to generate a single gff file
##### Alejandra Escobar, EMBL-EBI
##### March 14, 2023


def names_map(ori_fasta, ren_fasta):
    ### Creating an ID conversion dict
    originals = []
    if os.stat(ori_fasta).st_size > 0:
        for record in SeqIO.parse(ori_fasta, "fasta"):
            originals.append(str(record.id))

    renamed = []
    if os.stat(ren_fasta).st_size > 0:
        for record in SeqIO.parse(ren_fasta, "fasta"):
            renamed.append(str(record.id))
    contig_names = dict(zip(renamed, originals))
    return(contig_names)

def pprmeta_parser(pprm, contig_names):
    ### Saving ppr-meta predictions
    plasmids = {}
    pl_counter = 0
    if os.path.exists(pprm):
        if os.stat(pprm).st_size > 0:
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
                        pl_counter += 1
                        plas_name = "plasmid_" + str(pl_counter)
                        contig = contig_names[header]
                        plasmids[contig] = (length, plas_name)
    return(plasmids)

def checkv_parser(checkv):
    ### Saving checkV metadata
    checkv_pass = []
    viriqual = {}
    for summ_file in checkv:
        if os.path.exists(summ_file):
            head, tail = os.path.split(summ_file)
            current_qual = tail.split("_")[0]
            with open(summ_file, "r") as input_table:
                next(input_table)
                for line in input_table:
                    l_line = line.rstrip().split("\t")
                    phage_id = l_line[0].replace('prophage-0:','prophage-1:')
                    viral_genes = int(l_line[5])
                    checkv_quality = l_line[7]
                    kmer_freq = float(l_line[12])
                    viriqual[phage_id] = current_qual
                    if current_qual == "high":
                        checkv_pass.append(phage_id)
                    else:
                        if all([viral_genes > 0, checkv_quality != "Not-determined"]):
                            if any(
                                [
                                    checkv_quality == "Low-quality",
                                    checkv_quality == "Medium",
                                ]
                            ):
                                if kmer_freq <= 1.0:
                                    checkv_pass.append(phage_id)
                            else:
                                checkv_pass.append(phage_id)
        else:
            print('No checkV files found\n')
    return(checkv_pass, viriqual)

def phage_seq_save(viri_fa):
    phage_seqs = {}
    for fasta_file in viri_fa:
        if os.path.exists(fasta_file):
            for record in SeqIO.parse(fasta_file, "fasta"):
                phage_ID = str(record.description).replace(' ','|').replace('prophage-0:','prophage-1:')
                phage_seqs[phage_ID] = str(record.seq)
    return(phage_seqs)


def virify_parser(viri, checkv_pass, viriqual, phage_seqs):
    ### Saving virify predictions
    phages_metadata = {}
    contig_phages = {}
    protcoord_protid = {}
    prot_phage = {}
    viralprot_annot = {}
    viralname_coord = {}
    if os.path.exists(viri):
        with open(viri, "r") as input_gff, open('virify_HQ.fasta','w') as to_fasta:
            for line in input_gff:
                if not line.startswith("#"):
                    (
                        contig,
                        source,
                        seq_type,
                        start,
                        end,
                        score,
                        strand,
                        phase,
                        attributes,
                    ) = line.rstrip().split("\t")
                    if seq_type == "prophage":
                        current_phage = attributes.split(";")[0].replace("ID=", "")
                        viralname_coord[(contig, start, end)] = current_phage
                        seq_id = contig + "|" + seq_type + "-" + start + ":" + end
                        if seq_id in checkv_pass:
                            phages_metadata[seq_id] = (
                                contig,
                                source,
                                seq_type,
                                start,
                                end,
                                score,
                                strand,
                                phase,
                                attributes,
                            )
                            if contig in contig_phages:
                                contig_phages[contig].append(seq_id)
                            else:
                                contig_phages[contig] = [seq_id]  

                            if seq_id in phage_seqs:
                                to_fasta.write('>'+seq_id+'\n')
                                to_fasta.write(str(record.seq)+'\n')

                    elif seq_type == "viral_sequence":
                        seq_id = contig
                        current_phage = attributes.split(";")[0].replace("ID=", "")
                        viralname_coord[(contig, start, end)] = current_phage
                        if seq_id in checkv_pass:
                            new_attr = "virify_quality=" + viriqual[seq_id] + "-quality"
                            attributes = attributes + ";" + new_attr
                            phages_metadata[seq_id] = (
                                contig,
                                source,
                                seq_type,
                                start,
                                end,
                                score,
                                strand,
                                phase,
                                attributes,
                            )
                            contig_phages[contig] = [seq_id]

                            if seq_id in phage_seqs:
                                to_fasta.write('>'+seq_id+'|viral_sequence\n')
                                to_fasta.write(str(record.seq)+'\n')

                    elif seq_type == "CDS":
                        protein_id = attributes.split(";")[0].replace("ID=", "")
                        seq_belong = protein_id.split("_")
                        seq_belong.pop(-1)
                        seq_belong = "_".join(seq_belong)
                        if seq_belong in phages_metadata:
                            prot_phage[protein_id] = current_phage
                            viralprot_annot[protein_id] = (
                                contig,
                                source,
                                seq_type,
                                start,
                                end,
                                score,
                                strand,
                                phase,
                                attributes,
                            )
                            protcoord_protid[(contig, start, end)] = protein_id
    else:
        print('No VIRIfy files found\n')

    return(
        phages_metadata,
        contig_phages,
        protcoord_protid,
        prot_phage,
        viralprot_annot,
        viralname_coord,
    )

def momo_parser(momo):
    ### Saving momofy predictions to find obvious misannotations
    momo_pred = {}
    with open(momo, "r") as input_table:
        for line in input_table:
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
                if any([seq_type == "integron", seq_type == "conjugative_transposon"]):
                    mge_id = attr.split(";")[0].replace("ID=", "")
                    momo_pred[mge_id] = [
                        contig,
                        seq_source,
                        seq_type,
                        start,
                        end,
                        score,
                        strand,
                        phase,
                        attr,
                    ]
    return(momo_pred)

def overlaps(momo_pred, phages_metadata, viralname_coord):
    ## Finding overlapping predictions and printing a list of candidates to remove
    bad_annot = []
    with open("integron_phage_overlapps.txt", "w") as to_over:
        to_over.write(
            "#contig\tintegron_ID\tphage_ID\tintegron_cov\tphage_cov\tintegron_len\tphage_len\tto_discard\n"
        )
        for integron in momo_pred:
            m_contig = momo_pred[integron][0]
            m_seqtype = momo_pred[integron][2]
            m_start = int(momo_pred[integron][3])
            m_end = int(momo_pred[integron][4])
            m_atrib = momo_pred[integron][8]
            m_range = range(m_start, m_end + 1)
            m_len = m_end - m_start

            for phage in phages_metadata:
                p_contig = phages_metadata[phage][0]
                p_seqtype = phages_metadata[phage][2]
                p_start = int(phages_metadata[phage][3])
                p_end = int(phages_metadata[phage][4])
                p_atrib = phages_metadata[phage][8]
                p_range = range(p_start, p_end + 1)
                p_len = p_end - p_start
                phage_name = viralname_coord[(p_contig, str(p_start), str(p_end))]

                if m_contig == p_contig:
                    intersection = len(list(set(m_range) & set(p_range)))
                    if intersection > 0:
                        m_cov = float(intersection) / float(m_len)
                        p_cov = float(intersection) / float(p_len)
                        if any([m_cov > 0.9, p_cov > 0.9]):
                            to_discard = "none"
                            if all([m_cov > p_cov, "without_identified_DR" in m_atrib]):
                                to_discard = integron
                            elif all(
                                [p_cov > m_cov, "checkv_quality=Low-quality" in p_atrib]
                            ):
                                to_discard = phage_name
                            to_over.write(
                                "\t".join(
                                    [
                                        m_contig,
                                        integron,
                                        phage_name,
                                        str(m_cov),
                                        str(p_cov),
                                        str(m_len),
                                        str(p_len),
                                        to_discard,
                                    ]
                                )
                                + "\n"
                            )
                            bad_annot.append(to_discard)
    return(bad_annot)


def integra(
    momo, 
    bad_annot, 
    contig_phages, 
    plasmids, 
    protcoord_protid, 
    phages_metadata, 
    viralname_coord,
    viralprot_annot
    ):

    ### Parsing momofy output and adding virify and ppr-meta predictions
    ## Restoring protein ID on virify proteins
    output_gff = "mobilome_predictions.gff"
    pl_counter = 0
    used_contigs = []
    momo_labels = [
        "insertion_sequence",
        "terminal_inverted_repeat_element",
        "integron",
        "attC_site",
        "conjugative_transposon",
        "direct_repeat",
    ]

    with open(momo, "r") as input_table, open(output_gff, "w") as to_gff:
        for line in input_table:
            l_line = line.rstrip().split("\t")
            ## Annotation lines have exactly 9 columns
            if len(l_line) == 9:
                other = 0
                seqid,source,seq_type,start,end,score,strand,phase,attrib = line.rstrip().split("\t")

                # Printing mobilome predictions
                if seq_type in momo_labels:
                    other = 1
                    mge_id = attrib.split(";")[0].replace("ID=", "")
                    if ":" in mge_id:
                        mge_id = mge_id.split(":")[1]
                    if not mge_id in bad_annot:
                        to_gff.write(line)

                if all([seqid in contig_phages, not seqid in used_contigs]):
                    used_contigs.append(seqid)
                    for phage in contig_phages[seqid]:
                        if not phage in bad_annot:
                            to_gff.write(
                                "\t".join(
                                    phages_metadata[phage]
                                )
                                + "\n"
                            )

                if all([seqid in plasmids, not seqid in used_contigs]):
                    used_contigs.append(seqid)
                    plas_id = plasmids[seqid][1]
                    plas_attrib = (
                        "ID="
                        + plas_id
                        + ";gbkey=mobile_element;mobile_element_type=plasmid"
                    )
                    to_gff.write(
                        "\t".join(
                            [ 
                                seqid,
                                "PPR-meta",
                                "plasmid",
                                "1",
                                plasmids[seqid][0],
                                ".",
                                ".",
                                ".",
                                plas_attrib,
                            ]
                        )
                        + "\n"
                    )
                        
                # Printing the CDS
                if seq_type == "CDS":
                    other = 1
                    new_atrr = attrib
                    attr_l = attrib.split(";")
                    prot_id = attr_l[0].replace("ID=", "")

                    if "from_mge" in new_atrr:
                        attr_l = []
                        for element in new_atrr.split(";"):
                            attr_id, attr_val = element.split("=")
                            if attr_val not in bad_annot:
                                attr_l.append(element)
                        new_atrr = ";".join(attr_l)

                    if seqid in contig_phages:
                        flag = 0
                        gene_start = int(start)
                        gene_end = int(end)
                        gene_len = gene_end - gene_start
                        gene_range = range(gene_start, gene_end + 1)
                        for phage in contig_phages[seqid]:
                            phage_start = int(phages_metadata[phage][3])
                            phage_end = int(phages_metadata[phage][4])
                            composite_key = (seqid, str(phage_start), str(phage_end))
                            virus_name = viralname_coord[composite_key]
                            if not virus_name in bad_annot:
                                phage_range = range(phage_start, phage_end + 1)
                                intersection = len(list(set(phage_range) & set(gene_range)))
                                if intersection > 0:
                                    gene_cov = float(intersection) / float(gene_len)
                                    if gene_cov >= 0.75:
                                        flag = 1
                                        phag_atrr = []
                                        for element in attr_l:
                                            if element.split("=")[0] == "from_mge":
                                                if seqid in plasmids:
                                                    added_phag = (
                                                        element
                                                        + ","
                                                        + plasmids[seqid][1]
                                                        + ","
                                                        + virus_name
                                                    )
                                                else:
                                                    added_phag = element + "," + virus_name
                                                phag_atrr.append(added_phag)
                                            else:
                                                phag_atrr.append(element)
                                        new_atrr = ";".join(phag_atrr)
                                        if not "from_mge" in new_atrr:
                                            new_atrr = new_atrr + ";from_mge=" + virus_name
    
                                        # Checking if exact protein is in virify annotation list
                                        match_coord = (
                                            seqid,
                                            str(gene_start),
                                            str(gene_end),
                                        )
                                        extra_att = []
                                        if match_coord in protcoord_protid:
                                            for element in viralprot_annot[
                                                protcoord_protid[match_coord]
                                            ][8].split(";"):
                                                if "viphog" in element:
                                                    extra_att.append(element)
                                        if len(extra_att) > 0:
                                            extra_att = ";".join(extra_att)
                                            new_atrr = new_atrr + ";" + extra_att
                        if flag == 0:
                            if seqid in plasmids:
                                extra_att = ";from_mge=" + plasmids[seqid][1]
                                new_atrr = attrib + extra_att
                            else:
                                new_atrr = attrib
                    if seqid in plasmids:
                        plas_atrr = []
                        for element in attr_l:
                            if element.split("=")[0] == "from_mge":
                                added_plas = element + "," + plasmids[seqid][1]
                                plas_atrr.append(added_plas)
                            else:
                                plas_atrr.append(element)
                        new_atrr = ";".join(plas_atrr)
                        if not "from_mge" in new_atrr:
                            new_atrr = new_atrr + ";from_mge=" + plasmids[seqid][1]
                    
                    to_gff.write(
                        "\t".join(
                            [
                                seqid,
                                source,
                                seq_type,
                                start,
                                end,
                                score,
                                strand,
                                phase,
                                new_atrr,
                            ]
                        )
                        + "\n"
                    )
                    
                # Printing anything with other label: ncRNA, rRNA, tRNA, CRISPR, etc...
                if other == 0:
                    to_gff.write(line)
            else:
                to_gff.write(line)



def main():
    parser = argparse.ArgumentParser(
        description="This script integrates the outputs of MoMofy and VIRify to generate a single gff file. It also use the PPR-meta output generated by VIRify to label plasmids and the MobileOG-DB annotation generated by MoMofy. Please provide the relevant input files"
    )
    parser.add_argument(
        "--virify_gff", 
        type=str, 
        help="Virify output in gff format"
    )
    parser.add_argument(
        "--virify_fasta",
        nargs="*",
        help="Virify outputs in fasta format"
    )
    parser.add_argument(
        "--checkv_summ", 
        nargs="*", 
        help="checkv original summary files"
    )
    parser.add_argument(
        "--momofy_gff", 
        type=str, 
        help="MoMofy output in gff format", 
        required=True
    )
    parser.add_argument(
        "--pprmeta",
        type=str,
        help="Find the file in VIRify output 01-viruses/pprmeta/pprmeta.csv",
    )
    parser.add_argument(
        "--original_assem",
        type=str, 
        help="Original assembly fasta file"
    )
    parser.add_argument(
        "--renamed_virify",
        type=str,
        help="Renamed assembly fasta file"
    )
    args = parser.parse_args()

    ### Setting up variables
    viri = args.virify_gff
    viri_fa = args.virify_fasta
    checkv = args.checkv_summ
    momo = args.momofy_gff
    pprm = args.pprmeta
    ori_fasta = args.original_assem
    ren_fasta = args.renamed_virify


    ### Calling functions
    contig_names = names_map(ori_fasta, ren_fasta)

    plasmids = pprmeta_parser(pprm, contig_names)

    (checkv_pass, viriqual) = checkv_parser(checkv)

    phage_seqs = phage_seq_save(viri_fa)

    (
        phages_metadata,
        contig_phages,
        protcoord_protid,
        prot_phage,
        viralprot_annot,
        viralname_coord,
    ) = virify_parser(
        viri, 
        checkv_pass, 
        viriqual,
        phage_seqs,
        )

    momo_pred = momo_parser(momo)

    bad_annot = overlaps(
        momo_pred, 
        phages_metadata, 
        viralname_coord,
        )

    integra(
        momo,
        bad_annot, 
        contig_phages, 
        plasmids, 
        protcoord_protid, 
        phages_metadata, 
        viralname_coord,
        viralprot_annot,
        )    


if __name__ == "__main__":
    main()

