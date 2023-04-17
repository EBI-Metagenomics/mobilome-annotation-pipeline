#!/usr/bin/env python

from Bio import SeqIO
import argparse
import sys
import os.path
import glob


##### This script integrates and parse the output of ICEfinder, IntegronFinder, ISEScan, and PaliDIS for MoMofy
##### Alejandra Escobar, EMBL-EBI
##### January 11, 2023


def names_map(map_file):
    ### Saving contig names equivalence
    names_equiv = {}
    with open(map_file, "r") as input_map:
        for line in input_map:
            new_name, old_name = line.rstrip().split("\t")
            names_equiv[new_name.replace(">", "")] = old_name
    inv_names_equiv = {v: k for k, v in names_equiv.items()}
    return(names_equiv, inv_names_equiv)


def icf_parser(icf_seqs, icf_dr_file, icf_results):
    ### Saving the ICEfinder sequences
    mge_counter = 0
    mge_data = {}
    mge_nuc = {}
    icf_nuc = {}
    icf_dr = {}
    if os.stat(icf_seqs).st_size > 0:
        for record in SeqIO.parse(icf_seqs, "fasta"):
            my_chain = str(record.seq).upper()
            my_desc = str(record.description)
            my_id = my_desc.split(" ")[0] + "|" + my_desc.split(" ")[5]
            icf_nuc[my_id] = my_chain

        ### Saving ICEfinder direct repeat coordinates
        with open(icf_dr_file, "r") as input_table:
            for line in input_table:
                line_l = line.rstrip().split()
                if len(line_l) == 3:
                    new_id = (
                        line_l[0]
                        .replace("result/", "")
                        .replace(":DR:", "")
                        .replace("/", "|")
                    )
                    dr_1_s = line_l[1].split("..")[0]
                    dr_1_e = line_l[1].split("..")[1]
                    dr_1 = (dr_1_s, dr_1_e)
                    dr_2_s = line_l[2].split("..")[0]
                    dr_2_e = line_l[2].split("..")[1]
                    dr_2 = (dr_2_s, dr_2_e)
                    icf_dr[new_id] = (dr_1, dr_2)

        ### Parsing ICEfinder summaries
        with open(icf_results, "r") as input_table:
            for line in input_table:
                (
                    icefinder_Job_id,
                    strain,
                    genome_len,
                    icefinder_output_name,
                    description,
                    coordinate,
                    length,
                    oriT,
                    gc,
                    genome_GC,
                    delta_GC,
                    arg,
                    vf,
                ) = line.rstrip().split("\t")
                contig = icefinder_Job_id
                description = (
                    description.replace(" ", "_")
                    .replace("Putative_", "")
                    .replace("_AICE", "AICE")
                    .replace(":", "")
                )

                if not "conjugative_region" in description:
                    mge_counter += 1
                    mge_id = "icf_" + str(mge_counter)
                    start = int(coordinate.split("..")[0])
                    end = int(coordinate.split("..")[1])
                    coord = (start, end)
                    value = (contig, description, coord)
                    mge_data[mge_id] = value

                    seq_id = contig + "|" + coordinate
                    if seq_id in icf_nuc:
                        my_id = (
                            ">"
                            + mge_id
                            + "|"
                            + contig
                            + "|"
                            + str(start)
                            + ":"
                            + str(end)
                            + "|"
                            + description
                        )
                        mge_nuc[my_id] = icf_nuc[seq_id]

                    composite_id = contig + "|" + icefinder_output_name
                    if composite_id in icf_dr:
                        icf_dr[mge_id] = icf_dr.pop(composite_id)
    return(mge_data, mge_nuc, icf_dr)


def integron_parser(mge_data, mge_nuc, integron_results, inf_gbks):
    ### Parsing IntegronFinder output
    mge_counter = 0
    attC_site = {}
    if os.stat(integron_results).st_size > 0:
        with open(integron_results, "r") as input_table:
            next(input_table)
            next(input_table)
            for line in input_table:
                id_replicon, calin, complete, in0, topology, size = line.rstrip().split("\t")
                if int(complete) > 0:
                    description = "Complete_integron"
                    gbk_file = id_replicon + ".gbk"
                    if gbk_file in inf_gbks:
                        flag = 0
                        for gb_record in SeqIO.parse(gbk_file, "genbank"):
                            for feature in gb_record.features:
                                if feature.type == "integron":
                                    if feature.qualifiers["integron_type"][0] == "complete":
                                        flag = 1
                                        mge_counter += 1
                                        mge_id = "inf_" + str(mge_counter)
                                        start = int(feature.location.start)
                                        end = int(feature.location.end)
                                        coord = (start, end)
                                        value = (id_replicon, description, coord)
                                        mge_data[mge_id] = value
                                        my_id = (
                                            ">"
                                            + mge_id
                                            + "|"
                                            + id_replicon
                                            + "|"
                                            + str(start)
                                            + ":"
                                            + str(end)
                                            + "|"
                                            + description
                                        )
                                        my_chain = str(feature.extract(gb_record.seq)).upper()
                                        mge_nuc[my_id] = my_chain
                                if feature.type == "attC":
                                    if flag == 1:
                                        start = str(feature.location.start)
                                        end = str(feature.location.end)
                                        coord = (start, end)
                                        attC_site[mge_id] = coord
                                        flag = 0
    return(mge_data, mge_nuc, attC_site)

def isescan_parser(mge_data, mge_nuc, iss_seqs, iss_results):
    ### Parsing ISEScan outputs
    mge_counter = 0
    raw_iss = {}
    itr_sites = {}
    if os.stat(iss_results).st_size > 0:
        for record in SeqIO.parse(iss_seqs, "fasta"):
            my_chain = str(record.seq).upper()
            raw_iss[str(record.id)] = my_chain

        with open(iss_results, "r") as input_table:
            next(input_table)
            for line in input_table:
                mge_counter += 1
                mge_id = "iss_" + str(mge_counter)
                (
                    contig,
                    family,
                    cluster,
                    isBegin,
                    isEnd,
                    isLen,
                    ncopy4is,
                    start1,
                    end1,
                    start2,
                    end2,
                    score,
                    irId,
                    irLen,
                    nGaps,
                    orfBegin,
                    orfEnd,
                    strand,
                    orfLen,
                    e_value,
                    e_value4copy,
                    iss_type,
                    ov,
                    tir,
                ) = line.rstrip().split("\t")

                ## Keeping complete predictions only
                if iss_type == "c":
                    if int(irLen) == 0:
                        description = "without_TIR"
                    else:
                        description = "with_TIR"
                        ir_1 = (start1, end1)
                        ir_2 = (start2, end2)
                        itr_sites[mge_id] = (ir_1, ir_2)

                    description = family + "_" + description
                    coord = (int(isBegin), int(isEnd))
                    value = (contig, description, coord)
                    mge_data[mge_id] = value
                    fasta_id = contig + "_" + isBegin + "_" + isEnd + "_" + strand
                    if fasta_id in raw_iss:
                        new_fasta_id = (
                            ">"
                            + mge_id
                            + "|"
                            + contig
                            + "|"
                            + isBegin
                            + ":"
                            + isEnd
                            + "|"
                            + description
                        )
                        mge_nuc[new_fasta_id] = raw_iss[fasta_id]
                    else:
                        print("No fasta sequence for: " + fasta_id)
    return(mge_data, mge_nuc, itr_sites)


def palids_parser(mge_data, mge_nuc, itr_sites, pal_seqs, pal_results, inv_names_equiv):
    ### Parsing Palidis output
    mge_counter = 0
    raw_pal = {}
    if os.path.isfile(pal_seqs):
        if os.stat(pal_seqs).st_size > 0:
            for record in SeqIO.parse(pal_seqs, "fasta"):
                my_chain = str(record.seq).upper()
                raw_pal[str(record.id)] = my_chain

            with open(pal_results, "r") as input_table:
                next(input_table)
                for line in input_table:
                    mge_counter += 1
                    mge_id = "pal_" + str(mge_counter)
                    (
                        is_name,
                        sample_id,
                        contig,
                        itr1_start_position,
                        itr1_end_position,
                        itr2_start_position,
                        itr2_end_position,
                        description,
                    ) = line.rstrip().split("\t")
                    description = "IS_with_TIR"
                    coord = (int(itr1_start_position), int(itr2_end_position))
                    contig = inv_names_equiv[contig]
                    value = (contig, description, coord)
                    mge_data[mge_id] = value

                    ir_1 = (itr1_start_position, itr1_end_position)
                    ir_2 = (itr2_start_position, itr2_end_position)
                    itr_sites[mge_id] = (ir_1, ir_2)
    
                    if is_name in raw_pal:
                        new_fasta_id = (
                            ">"
                            + mge_id
                            + "|"
                            + contig
                            + "|"
                            + itr1_start_position
                            + ":"
                            + itr2_end_position
                            + "|"
                            + description
                        )
                        mge_nuc[new_fasta_id] = raw_pal[is_name]
                    else:
                        print("No fasta file for: " + is_name)
    return(mge_data, mge_nuc, itr_sites)

def is_grouping(mge_data):
    ### Cleaning redundant predictions
    # Grouping elements by contig
    contig_ismge = {}
    contig_inmge = {}
    for prediction in mge_data:
        if any(["iss" in prediction, "pal" in prediction]):
            contig = mge_data[prediction][0]
            if contig in contig_ismge:
                contig_ismge[contig].append(prediction)
            else:
                contig_ismge[contig] = [prediction]
        elif any(["icf" in prediction, "inf" in prediction]):
            contig = mge_data[prediction][0]
            if contig in contig_inmge:
                contig_inmge[contig].append(prediction)
            else:
                contig_inmge[contig] = [prediction]
    return(contig_ismge, contig_inmge)

def is_overlap(contig_ismge, mge_data, mge_nuc):
    # Parsing insertion sequences overlapps
    to_remove = []
    for contig in contig_ismge:
        if len(contig_ismge[contig]) > 1:
            tools_list = []
            for element in contig_ismge[contig]:
                prefix = element.split("_")[0]
                tools_list.append(prefix)
            tools_list = list(set(tools_list))
            if len(tools_list) > 1:
                iss_list = [x for x in contig_ismge[contig] if "iss" in x]
                pal_list = [x for x in contig_ismge[contig] if "pal" in x]
                for pal_element in pal_list:
                    pal_coord = mge_data[pal_element][2]
                    pal_len = pal_coord[1] - pal_coord[0]
                    for iss_element in iss_list:
                        iss_coord = mge_data[iss_element][2]
                        iss_len = iss_coord[1] - iss_coord[0]
                        intersection = len(
                            list(
                                set(range(pal_coord[0], pal_coord[1] + 1))
                                & set(range(iss_coord[0], iss_coord[1] + 1))
                            )
                        )
                        if intersection > 0:
                            iss_cov = float(intersection) / float(iss_len)
                            pal_cov = float(intersection) / float(pal_len)
                            if any([iss_cov > 0.1, pal_cov > 0.1]):
                                to_remove.append(iss_element)

    to_remove = list(set(to_remove))

    # Printing discarded redundant predictions to file
    output_discard = "discarded_mge.txt"
    with open(output_discard, "w") as to_discard:
        to_discard.write(
            "\t".join(
                [
                    '#element_id',
                    'description',
                    'cause',
                ]
            )
            + "\n"
        )
        for element in to_remove:
            to_discard.write(
                "\t".join(
                    [
                        element,
                        mge_data[element][1],
                        'overlapping',
                    ]
                )
                + "\n"
            )
                    
    # Removing redundant IS
    for element in to_remove:
        nuc_llave = [
            ">" + element,
            mge_data[element][0],
            str(mge_data[element][2][0]) + ":" + str(mge_data[element][2][1]),
            mge_data[element][1],
        ]
        nuc_llave = "|".join(nuc_llave)
        del mge_data[element]
        del mge_nuc[nuc_llave]

    return(mge_data, mge_nuc)


def int_overlap(contig_inmge):
    # Parsing integron overlapps. Reporting only, not to be removed
    output_nested = "nested_integrons.txt"
    with open(output_nested, "w") as to_nested:
        to_nested.write(
            "\t".join(
                [
                    '#contig',
                    'icefinder_element',
                    'intfinder_element',
                    'icefinder_cov',
                    'intfinder_cov',
                ]
            )
            + '\n'
        )
        for contig in contig_inmge:
            if len(contig_inmge[contig]) > 1:
                tools_list = []
                for element in contig_inmge[contig]:
                    prefix = element.split("_")[0]
                    tools_list.append(prefix)
                tools_list = list(set(tools_list))
                if len(tools_list) > 1:
                    inf_list = [x for x in contig_inmge[contig] if "inf" in x]
                    icf_list = [x for x in contig_inmge[contig] if "icf" in x]
                    for inf_element in inf_list:
                        inf_coord = mge_data[inf_element][2]
                        inf_len = inf_coord[1] - inf_coord[0]
                        for icf_element in icf_list:
                            icf_coord = mge_data[icf_element][2]
                            icf_len = icf_coord[1] - icf_coord[0]
                            intersection = len(
                                list(
                                    set(range(inf_coord[0], inf_coord[1] + 1))
                                    & set(range(icf_coord[0], icf_coord[1] + 1))
                                )
                            )
                            if intersection > 0:
                                inf_cov = float(intersection) / float(inf_len)
                                icf_cov = float(intersection) / float(icf_len)
                                if any([inf_cov > 0.75, icf_cov > 0.75]):
                                    to_nested.write(
                                        "\t".join(
                                            [
                                                contig,
                                               icf_element,
                                               inf_element,
                                               str(icf_cov),
                                               str(inf_cov),
                                            ]
                                        )
                                        + "\n"
                                    )

def location_parser(cds_loc):
    ### Saving the CDS' coordinates and contigs location
    contig_prots = {}
    prots_coord = {}
    if os.stat(cds_loc).st_size > 0:
        with open(cds_loc, "r") as input_table:
            for line in input_table:
                l_line = line.rstrip().split("\t")
                if len(l_line) == 9:
                    if l_line[2] == "CDS":
                        contig = l_line[0]
                        prot_source = l_line[1]
                        start = int(l_line[3])
                        end = int(l_line[4])
                        strand = l_line[6]
                        attrib = l_line[8]
                        prot_id = attrib.split(";")[0]
                        if prot_id.startswith("ID="):
                            prot_id = prot_id.replace("ID=", "")
                            value = (start, end, strand)
                            prots_coord[prot_id] = value
                            if contig in contig_prots:
                                contig_prots[contig].append(prot_id)
                            else:
                                contig_prots[contig] = [prot_id]
                        else:
                            sys.exit('CDS ID is expected as the first token in the attributes field')
                            
    return(contig_prots, prots_coord)

def mob_cds(mge_data, names_equiv, contig_prots, user_gff, prots_coord):
    ### Finding the CDS encoded in the mobilome
    mob_proteome = []
    mge_proteins = {}
    proteins_mge = {}
    contigs_elements = {}
    for element in mge_data:
        mge_proteins[element] = []
        mge_start = mge_data[element][2][0]
        mge_end = mge_data[element][2][1]
        mge_range = range(mge_start, mge_end + 1)
        mge_len = mge_end - mge_start

        if mge_data[element][0] in contigs_elements:
            contigs_elements[mge_data[element][0]].append(element)
        else:
            contigs_elements[mge_data[element][0]] = [element]

        if user_gff == "T":
            # We are using user contig ID. It is the same as original assembly. We need to transform to find MGE
            contig = names_equiv[mge_data[element][0]]
        else:
            # We are using prokka contig ID. It is the same as MGEs
            contig = mge_data[element][0]

        if contig in contig_prots:
            for protein in contig_prots[contig]:
                prot_start = prots_coord[protein][0]
                prot_end = prots_coord[protein][1]
                prot_range = range(prot_start, prot_end + 1)
                prot_len = prot_end - prot_start
                intersection = len(list(set(mge_range) & set(prot_range)))
                if intersection > 0:
                    mge_cov = float(intersection) / float(mge_len)
                    prot_cov = float(intersection) / float(prot_len)
                    if any([mge_cov > 0.75, prot_cov > 0.75]):
                        mge_proteins[element].append(protein)
                        mob_proteome.append(protein)
                        proteins_mge[protein] = element

    return(mge_proteins, mob_proteome, contigs_elements, proteins_mge)


def mobileog_parser(mog_results, mob_proteome):
    ### Parsing the mobileOG annotation file
    mog_annot = {}
    if os.stat(mog_results).st_size > 0:
        with open(mog_results) as input_table:
            for line in input_table:
                line_l = line.rstrip().split("\t")
                target_id = line_l[0]
                query_id = line_l[1].split(" ")[0]
                (
                    mog_id,
                    gene_name,
                    best_hit_id,
                    major,
                    minor,
                    db,
                    evidence,
                ) = target_id.split("|")
                function = mog_id + "|" + major + "|" + minor

                if query_id in mob_proteome:
                    if query_id not in mog_annot:
                        mog_annot[query_id] = function.replace(",", "/")
    return(mog_annot)


def quality_filter(mge_data, mge_proteins, contigs_elements):
    # Removing predictions of len<500 and with no CDS
    no_cds = []
    len_500 = []
    for element in mge_data:
        element_len = int(mge_data[element][2][1]) - int(mge_data[element][2][0])
        if element_len < 500:
            len_500.append(element)
        elif len(mge_proteins[element]) == 0:
            no_cds.append(element)

    with open("discarded_mge.txt", "a") as to_discard:
        for element in no_cds:
            to_discard.write(element + "\t" + mge_data[element][1] + "\tno_cds\n")
            del mge_data[element]
            for contig in contigs_elements:
                if element in contigs_elements[contig]:
                    contigs_elements[contig].remove(element)

        for element in len_500:
            to_discard.write(element + "\t" + mge_data[element][1] + "\tmge<500bp\n")
            del mge_data[element]
            for contig in contigs_elements:
                if element in contigs_elements[contig]:
                    contigs_elements[contig].remove(element)


    return(mge_data, contigs_elements)


def flanking_data(id_to_print, source, seq_type, start, end, flank_id):
    to_print = "\t".join(
        [
            id_to_print, 
            source, 
            seq_type, 
            start, 
            end, 
            ".", 
            ".", 
            ".", 
            flank_id
        ]
    )
    return to_print


def writing_gff(
    cds_loc, 
    names_equiv, 
    inv_names_equiv, 
    contigs_elements, 
    itr_sites, 
    icf_dr,
    attC_site,
    mge_data, 
    mge_nuc, 
    proteins_mge, 
    mog_annot,
    user_gff,
    ):

    ### Writing the mobilome outputs
    prefixes = {
        "icf": "ICEfinder",
        "inf": "IntegronFinder",
        "iss": "ISEScan",
        "pal": "PaliDIS",
    }
    used_contigs = []
    output_fna = "momofy_predictions.fna"
    output_gff = "momofy_predictions.gff"
    with open(cds_loc, "r") as input_table, \
        open(output_fna, "w") as to_fasta, \
        open(output_gff, "w") as to_gff:
        for line in input_table:
            l_line = line.rstrip().split("\t")
            if len(l_line) == 9:
                seqid = l_line[0]
                if user_gff == "T":
                    # We are using user contig ID. It is the same as original assembly. We need to transform to find MGE
                    if seqid in inv_names_equiv:
                        contig = inv_names_equiv[seqid]
                        id_to_print = seqid
                    else:
                        next(input_table)
                else:
                    # We are using prokka contig ID. It is the same as MGEs, but we need to transform to print gff
                    contig = seqid
                    id_to_print = names_equiv[contig]

                if not seqid in used_contigs:
                    used_contigs.append(seqid)
                    if contig in contigs_elements:
                        for element in contigs_elements[contig]:
                            source = prefixes[element.split("_")[0]]
                            if any(["iss" in element, "pal" in element]):
                                seq_type = "insertion_sequence"
                                if element in itr_sites:
                                    tir_seq_type = "terminal_inverted_repeat_element"
                                    tir_1_id = "ID=TIR_1:" + element
                                    gff_line = flanking_data(
                                        id_to_print,
                                        source,
                                        tir_seq_type,
                                        itr_sites[element][0][0],
                                        itr_sites[element][0][1],
                                        tir_1_id,
                                    )
                                    to_gff.write(gff_line + "\n")
                                    tir_2_id = "ID=TIR_2:" + element
                                    gff_line = flanking_data(
                                        id_to_print,
                                        source,
                                        tir_seq_type,
                                        itr_sites[element][1][0],
                                        itr_sites[element][1][1],
                                        tir_2_id,
                                    )
                                    to_gff.write(gff_line + "\n")
                            elif "inf" in element:
                                seq_type = "integron"
                                attc_seq_type = "attC_site"
                                attc_id = "ID=attC:" + element
                                gff_line = flanking_data(
                                    id_to_print,
                                    source,
                                    attc_seq_type,
                                    attC_site[element][0],
                                    attC_site[element][1],
                                    attc_id,
                                )
                                to_gff.write(gff_line + "\n")
                            else:
                                if "IME" in mge_data[element][1]:
                                    seq_type = "integron"
                                    if element in icf_dr:
                                        dr_seq_type = "direct_repeat"
                                        dr_1_id = "ID=DR_1:" + element
                                        gff_line = flanking_data(
                                            id_to_print,
                                            source,
                                            dr_seq_type,
                                            icf_dr[element][0][0],
                                            icf_dr[element][0][1],
                                            dr_1_id,
                                        )
                                        to_gff.write(gff_line + "\n")
                                        dr_2_id = "ID=DR_2:" + element
                                        gff_line = flanking_data(
                                            id_to_print,
                                            source,
                                            dr_seq_type,
                                            icf_dr[element][1][0],
                                            icf_dr[element][1][1],
                                            dr_2_id,
                                        )
                                        to_gff.write(gff_line + "\n")
                                elif "ICE" in mge_data[element][1]:
                                    seq_type = "conjugative_transposon"
                                    if element in icf_dr:
                                        dr_seq_type = "direct_repeat"
                                        dr_1_id = "ID=DR_1:" + element
                                        gff_line = flanking_data(
                                            id_to_print,
                                            source,
                                            dr_seq_type,
                                            icf_dr[element][0][0],
                                            icf_dr[element][0][1],
                                            dr_1_id,
                                        )
                                        to_gff.write(gff_line + "\n")
                                        dr_2_id = "ID=DR_2:" + element
                                        gff_line = flanking_data(
                                            id_to_print,
                                            source,
                                            dr_seq_type,
                                            icf_dr[element][1][0],
                                            icf_dr[element][1][1],
                                            dr_2_id,
                                        )
                                        to_gff.write(gff_line + "\n")
                            start = str(mge_data[element][2][0])
                            end = str(mge_data[element][2][1])
                            score = "."
                            strand = "."
                            phase = "."
                            attributes = (
                                "ID="
                                + element
                                + ";gbkey=mobile_element;mobile_element_type="
                                + mge_data[element][1]
                            )

                            to_gff.write(
                                "\t".join(
                                    [
                                        id_to_print,
                                        source,
                                        seq_type,
                                        start,
                                        end,
                                        score,
                                        strand,
                                        phase,
                                        attributes,
                                    ]
                                )
                                + "\n"
                            )

                            seqid_saved = [
                                ">" + element,
                                mge_data[element][0],
                                str(mge_data[element][2][0])
                                + ":"
                                + str(mge_data[element][2][1]),
                                mge_data[element][1],
                            ]
                            seqid_saved = "|".join(seqid_saved)
                            my_sequence = mge_nuc[seqid_saved]
                            nuc_id = [
                                ">" + element,
                                id_to_print,
                                str(mge_data[element][2][0])
                                + ".."
                                + str(mge_data[element][2][1]),
                                mge_data[element][1],
                            ]
                            nuc_id = "|".join(nuc_id)
                            to_fasta.write(nuc_id + "\n")
                            to_fasta.write(my_sequence + "\n")

                attrib = l_line[8]
                prot_id = attrib.split(";")[0].replace("ID=", "")
                if prot_id in proteins_mge:
                    parent_mge = proteins_mge[prot_id]
                    if prot_id in mog_annot:
                        function = mog_annot[prot_id].replace(" ", "_")
                        attrib = (
                            attrib 
                            + ";mobileOG=" 
                            + function 
                            + ";from_mge=" 
                            + parent_mge
                        )
                    else:
                        attrib = ( 
                            attrib 
                            + ";mobileOG=-" 
                            + ";from_mge=" 
                            + parent_mge
                        )

                l_line.pop(-1)
                l_line.pop(0)
                gff_line = [id_to_print] + l_line + [attrib]
                to_gff.write(
                    "\t".join(
                        gff_line
                    )
                    + "\n"
                )

            elif line.startswith("##sequence-region"):
                tag, seqid, start, end = line.rstrip().split()
                if user_gff == "T":
                    # We are using user contig ID. It is the same as original assembly
                    id_to_print = seqid
                else:
                    # We are using prokka contig ID. We need to transform to print gff
                    id_to_print = names_equiv[seqid]

                to_gff.write(
                    " ".join(
                        [
                            tag,
                            id_to_print,
                            start,
                            end,
                        ]
                    )
                    + "\n"
                )

            elif line.startswith(">"):
                seqid = line.rstrip().replace(">", "")
                if user_gff == "T":
                    # We are using user contig ID. It is the same as original assembly
                    id_to_print = seqid
                else:
                    # We are using prokka contig ID. We need to transform to print gff
                    id_to_print = names_equiv[seqid]
                to_gff.write(
                    ">" 
                    + id_to_print 
                    + "\n"
                )

            elif "##gff_version 3" in line:
                to_gff.write("##gff-version 3\n")

            else:
                to_gff.write(line)



def main():
    parser = argparse.ArgumentParser(
        description="This script integrates and parse the output of ICEfinder, IntegronFinder, ISEScan, and PaliDIS for MoMofy. Please provide the relevant input files"
    )
    parser.add_argument(
        "--user",
        type=str,
        help="Flag indicating if gff files are provided by the user",
        required=True,
    )
    parser.add_argument(
        "--cds_gff",
        type=str, 
        help="GFF prediction fasta file", 
        required=True,
    )
    parser.add_argument(
        "--map", 
        type=str, 
        help="Rename contigs file: contigID.map", 
        required=True,
    )
    parser.add_argument(
        "--iss_fa", 
        type=str, 
        help="ISEScan fasta file",
    )
    parser.add_argument(
        "--iss_tsv", 
        type=str, 
        help="ISEScan predictions table",
    )
    parser.add_argument(
        "--pal_fa", 
        type=str, 
        help="PaliDIS fasta file",
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
        "--icf_fa", 
        type=str, 
        help="ICEfinder fasta files (concatenated)",
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
    args = parser.parse_args()


    ### Setting up variables
    user_gff = args.user
    cds_loc = args.cds_gff
    map_file = args.map
    iss_seqs = args.iss_fa
    iss_results = args.iss_tsv
    pal_seqs = args.pal_fa
    pal_results = args.pal_tsv
    integron_results = args.inf_tsv
    inf_gbks = args.inf_gbks
    icf_results = args.icf_tsv
    icf_seqs = args.icf_fa
    icf_dr_file = args.icf_lim
    mog_results = args.mog_tsv

    ## Calling functions
    # Mapping contig names
    (names_equiv, inv_names_equiv) = names_map(map_file)

    # Parsing ICEfinder results
    (
        mge_data,
        mge_nuc,
        icf_dr,
    ) = icf_parser(
        icf_seqs, 
        icf_dr_file, 
        icf_results,
    )

    # Parsing IntegronFinder results
    (
        mge_data, 
        mge_nuc, 
        attC_site,
    ) = integron_parser(
        mge_data, 
        mge_nuc, 
        integron_results, 
        inf_gbks,
    )

    # Parsing ISEScan results
    (
        mge_data, 
        mge_nuc, 
        itr_sites,
    ) = isescan_parser(
        mge_data, 
        mge_nuc, 
        iss_seqs, 
        iss_results,
    )

    # Parsing PaliDIS results
    (
        mge_data, 
        mge_nuc, 
        itr_sites
    )= palids_parser(
        mge_data, 
        mge_nuc, 
        itr_sites, 
        pal_seqs, 
        pal_results,
        inv_names_equiv,
    )

    # Collecting insertion sequences predicted per contig
    (contig_ismge, contig_inmge) = is_grouping(mge_data)

    # Removing overlapping ISs
    (mge_data, mge_nuc) = is_overlap(
        contig_ismge, 
        mge_data, 
        mge_nuc,
    )

    # Reporting overlapping integrons
    int_overlap(contig_inmge)

    # Saving protein coordinates
    (contig_prots, prots_coord) = location_parser(cds_loc)

    # Saving mobilome proteins
    (
        mge_proteins, 
        mob_proteome, 
        contigs_elements, 
        proteins_mge,
    ) = mob_cds(
        mge_data, 
        names_equiv, 
        contig_prots,
        user_gff,
        prots_coord,
    )

    # Parsing mobileOG results
    mog_annot = mobileog_parser(mog_results, mob_proteome)

    # Filtering out low-quality results
    (mge_data , contigs_elements ) = quality_filter(mge_data, mge_proteins, contigs_elements)

    # Generating the final output
    writing_gff(
        cds_loc, 
        names_equiv, 
        inv_names_equiv, 
        contigs_elements, 
        itr_sites, 
        icf_dr,
        attC_site,
        mge_data, 
        mge_nuc, 
        proteins_mge, 
        mog_annot,
        user_gff,
    )


if __name__ == "__main__":
    main()


