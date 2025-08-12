#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2025 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import re
import argparse
from Bio import SeqIO


def parse_blast_uniprot(uniprot_annot):
    """
    Parse UniProt BLAST annotation file to create a query-to-product mapping dictionary.
    
    This function reads a tab-separated BLAST annotation file and extracts the mapping
    between query IDs and their corresponding UniProt product descriptions. The file
    is expected to follow the BLAST settings and quality filtering that that Prokka applies.
    
    :param uniprot_annot: Path to the UniProt annotation file
    :type uniprot_annot: str
    :return: Dictionary mapping query IDs to product descriptions
    :rtype: dict
    
    .. note::
       Expects BLAST results with quality filters matching prokka settings:
       "-evalue 1E-9 -qcov_hsp_perc 80 -num_descriptions 1 -num_alignments 1 -seg no"
    """
    uniprot_annot_dict = {}
    with open(uniprot_annot, "r") as input_file:
        for line in input_file:
            l_line = line.rstrip().split("\t")
            # Match qulity filters set up during blast process according with prokka settings:
            # "-evalue 1E-9 -qcov_hsp_perc 80 -num_descriptions 1 -num_alignments 1 -seg no"
            query_id = l_line[0]
            product = l_line[2]
            uniprot_annot_dict[query_id] = product
    return uniprot_annot_dict


def parse_merged_gff(gff_file, uniprot_annot_dict):
    """
    Parse the merged GFF file and extract the feature information with UniProt annotations.
    
    This function processes a GFF3 file containing genomic features and extracts
    information about tRNA/tmRNA features and coding sequences (CDS). It enriches CDS
    features with UniProt product annotations and organizes data for downstream analysis.
    
    :param gff_file: Path to the merged GFF3 file
    :type gff_file: str
    :param uniprot_annot_dict: Dictionary mapping query IDs to UniProt product descriptions
    :type uniprot_annot_dict: dict
    :return: Tuple containing multiple dictionaries for names mapping, tRNA data, 
             position data, total counts, locus mapping, and protein-contig associations
    :rtype: tuple
    
    .. note::
       The function processes both tRNA/tmRNA and CDS features, creating separate
       data structures for each feature type based on the GFF3 format.
    """
    names_map, trnadict, posdict, totalnum_dict, locusdict, prots_contigs = {}, {}, {}, {}, {}, {}
    valid_rnas = ["tRNA", "tmRNA"]
    header = ''
    with open(gff_file, "r") as input_gff:
        for line in input_gff:
            l_line = line.rstrip().split("\t")
            # Annotation lines have exactly 9 columns
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

                if contig in totalnum_dict:
                    totalnum_dict[contig] += 1
                else:
                    totalnum_dict[contig] = 1

                if seq_type in valid_rnas:
                    for attribute in attr.split(";"):
                        key, value = attribute.split("=")
                        if key == "ID":
                            ids = value
                            locusdict[ids] = ids
                        if key == "product":
                            product = value
                    pos = [int(start), int(end), strand, product]
                    trnadict[ids] = pos

                elif seq_type == "CDS":
                    for attribute in attr.split(";"):
                        key, value = attribute.split("=")
                        if key == "ID":
                            ids = value
                            locusdict[ids] = ids
                            header = ids.split("_")
                            header.pop(-1)
                            header = "_".join(header)
                        if key == "ori_id":
                            blastp_id = value

                            if blastp_id in uniprot_annot_dict:
                                product = uniprot_annot_dict[blastp_id]
                            else:
                                product = ""

                    names_map[blastp_id] = ids
                    pos = [int(start), int(end), strand, product]

                posdict[ids] = pos
                prots_contigs[ids] = contig

    return names_map, trnadict, posdict, header, totalnum_dict, locusdict, prots_contigs


def get_DR(dr_file):
    dr_dict = {}
    with open(dr_file, "r") as input_tsv:
        for line in input_tsv.readlines():
            contig, dr1_start, dr1_end, dr2_start, dr2_end = line.strip().split()
            dr = "|".join([dr1_start, dr1_end, dr2_start, dr2_end])
            if contig in dr_dict:
                dr_dict[contig].append(dr)
            else:
                dr_dict[contig] = [dr]
    return dr_dict


def ICE_filter(macsyfinder_out):
    filtered_lines = []
    with open(macsyfinder_out, "r") as ICEin:
        ICEfdict = {"ICE": []}
        IMEfdict = {}
        IMEgendict = {}
        AICEfdict = {}
        fICE = []
        fAICE = []
        for line in ICEin.readlines():
            if "Chromosome" in line:
                lines = line.strip().split("\t")
                if lines[7] != "1":
                    continue
                else:
                    filtered_lines.append(line)
                    IDtag = lines[3]
                    ICEtag = lines[5]
                    if "IME" in ICEtag:
                        if ICEtag not in IMEfdict:
                            IMEfdict[ICEtag] = [IDtag]
                            IMEgendict[ICEtag] = [lines[2]]
                        else:
                            IMEfdict[ICEtag] += [IDtag]
                            IMEgendict[ICEtag] += [lines[2]]
                    elif "AICE" in ICEtag:
                        if ICEtag not in AICEfdict:
                            AICEfdict[ICEtag] = [IDtag]
                        else:
                            AICEfdict[ICEtag] += [IDtag]
                        if ICEtag not in fAICE:
                            fAICE.append(ICEtag)
                    else:
                        ICEfdict["ICE"] += [IDtag]
                        if ICEtag not in fICE:
                            fICE.append(ICEtag)

    Indict = {
        "Phage_integrase": "Integrase",
        "UPF0236": "Integrase",
        "Recombinase": "Integrase",
        "rve": "Integrase",
        "TIGR02224": "Integrase",
        "TIGR02249": "Integrase",
        "TIGR02225": "Integrase",
        "PB001819": "Integrase",
    }

    IMEgenlist = []
    for k, v in IMEgendict.items():
        genecount = []
        for line in v:
            if "Relaxase_" in line or "T4SS_MOB" in line:
                genecount.append("MOB")
            elif line in Indict:
                genecount.append("Int")
        if len(set(genecount)) == 2:
            IMEgenlist.append(k)

    fIME = []
    for k, v in IMEfdict.items():
        for k1, v1 in ICEfdict.items():
            if not set(v).issubset(set(v1)):
                if k not in fIME and k in IMEgenlist:
                    fIME.append(k)

    concat_array = fICE + fIME + fAICE

    return concat_array, filtered_lines


def get_feat(feat):
    featuredict = {
        "Phage_integrase": "Integrase",
        "UPF0236": "Integrase",
        "Recombinase": "Integrase",
        "rve": "Integrase",
        "TIGR02224": "Integrase",
        "TIGR02249": "Integrase",
        "TIGR02225": "Integrase",
        "PB001819": "Integrase",
        "RepSAv2": "Rep",
        "DUF3631": "Rep",
        "Prim-Pol": "Rep",
        "FtsK_SpoIIIE": "Tra",
    }

    if feat in featuredict:
        return featuredict[feat] + "@" + feat
    elif "T4SS_MOB" in feat:
        tag = feat.split("_")[1]
        return "Relaxase@" + tag
    elif "Relaxase_" in feat:
        tag = feat.split("_")[1:]
        return "Relaxase@" + "_".join(tag)
    elif "t4cp" in feat:
        tag = feat.split("_")[1]
        return "T4CP@" + tag
    elif "tcpA" in feat:
        tag = feat.split("_")[1]
        return "T4CP@" + tag
    elif "FATA_" in feat or "FA_" in feat:
        tag = feat.split("_")[1]
        return "T4SS@" + tag
    else:
        return "T4SS@" + feat.replace("T4SS_", "")


def getnum(gene_id):
    id_num = int(gene_id.split("_")[-1].lstrip("0"))
    return id_num


def zill(header, num):
    restored_id = header + "_" + str(num).zfill(5)
    return restored_id


def find_max_distance(numbers):
    max_distance = -1
    max_distance_index = -1

    for i in range(len(numbers) - 1):
        distance = abs(numbers[i] - numbers[i + 1])
        if distance > max_distance:
            max_distance = distance
            max_distance_index = i

    if max_distance_index == -1:
        return None

    return numbers[max_distance_index], numbers[max_distance_index + 1]


def pos_tag(pos, posdict, ICE, final, totalnum, dirtag):
    tICE = ICE
    tfinal = final
    for k, v in posdict.items():
        vstart, vend = int(v[0]), int(v[1])
        if int(pos) <= vend:
            if dirtag == "s":
                tICE = getnum(k)
                tfinal = max(1, tICE - 5)
            else:
                if vstart > int(pos):
                    tICE = getnum(k) - 1
                else:
                    tICE = getnum(k)
                tfinal = min(totalnum, tICE + 5)
            break
    return tICE, tfinal


def merge_tRNA(ice_id, ICEdict, DR_dict, listgff, prots_contigs):

    [trnadict, posdict, header, total_dict, locusdict] = listgff

    contig = ''
    for key in ICEdict:
        if key in prots_contigs:
            contig = prots_contigs[key]
    totalnum = total_dict[contig]
    
    fICE = getnum(next(iter(ICEdict)))
    eICE = getnum(list(ICEdict.keys())[-1])
    nfICEnum = max(1, fICE - 5)
    neICEnum = min(totalnum, eICE + 5)

    ICEtagnum = [nfICEnum, neICEnum]
    trnalist = []
    for key, value in trnadict.items():
        if nfICEnum <= getnum(key) <= neICEnum:
            ICEtagnum.append(getnum(key))
            trnalist.append(value)

    DRlist = DR_dict[contig]

    ICEtagnum.sort()
    finalstart, finalend = find_max_distance(ICEtagnum)

    myDR1 = posdict[zill(header, fICE)][0]
    myDR2 = ""
    myDR3 = ""
    myDR4 = posdict[zill(header, eICE)][1]

    if trnalist:
        if finalend == neICEnum:
            fICE = finalstart
            finalstart = max(1, finalstart - 5)
            myDR1 = posdict[zill(header, fICE)][0]
            for line in DRlist:
                DRs = line.split("|")
                if int(DRs[3]) - int(DRs[0]) > 500000:
                    continue
                if int(DRs[3]) - int(DRs[0]) < 5000:
                    continue
                if (
                    int(posdict[zill(header, fICE)][0])
                    < int(DRs[0])
                    < int(posdict[zill(header, fICE)][1])
                ):
                    checktrna = 0
                    for key, value in trnadict.items():
                        if int(DRs[0]) <= int(value[0]) <= int(DRs[3]) and int(
                            DRs[0]
                        ) <= int(value[1]) <= int(DRs[3]):
                            checktrna += 1
                    if checktrna >= 2:
                        break

                    eICE, finalend = pos_tag(
                        DRs[3], posdict, eICE, finalend, totalnum, "e"
                    )
                    myDR1 = DRs[0]
                    myDR2 = DRs[1]
                    myDR3 = DRs[2]
                    myDR4 = DRs[3]
                    break

        elif finalstart == nfICEnum:
            eICE = finalend
            finalend = min(totalnum, finalend + 5)
            myDR4 = posdict[zill(header, eICE)][1]
            for line in DRlist:
                DRs = line.split("|")
                if int(DRs[3]) - int(DRs[0]) > 500000:
                    continue
                if int(DRs[3]) - int(DRs[0]) < 5000:
                    continue
                if (
                    int(posdict[zill(header, eICE)][0])
                    < int(DRs[3])
                    < int(posdict[zill(header, eICE)][1])
                ):
                    checktrna = 0
                    for key, value in trnadict.items():
                        if int(DRs[0]) <= int(value[0]) <= int(DRs[3]) and int(
                            DRs[0]
                        ) <= int(value[1]) <= int(DRs[3]):
                            checktrna += 1
                    if checktrna >= 2:
                        break

                    fICE, finalstart = pos_tag(
                        DRs[0], posdict, fICE, finalstart, totalnum, "s"
                    )
                    myDR1 = DRs[0]
                    myDR2 = DRs[1]
                    myDR3 = DRs[2]
                    myDR4 = DRs[3]
                    break

    return (
        contig,
        myDR1,
        myDR2,
        myDR3,
        myDR4,
        fICE,
        eICE,
        finalstart,
        finalend,
        posdict,
        header,
        trnalist,
        locusdict,
    )


def get_ICE(
    macsyfinder_out,
    dr_dict,
    trnadict,
    posdict,
    header,
    totalnum_dict,
    locusdict,
    names_map,
    prots_contigs,
):
    # Based in get_ICE, ICE_filter, get_feat, and merge_tRNA functions in single.py code
    # Filtering ICE predictions macsyfinder output file all_systmes.tsv
    ftag, filtered_lines = ICE_filter(macsyfinder_out)

    # Parsing again the macsyfinder info
    genes_icedict = {}
    infodict = {}
    for line in filtered_lines:
        lines = line.strip().split("\t")
        if lines[5] not in ftag:
            continue
        else:
            gbname = names_map[lines[1]]
            tags = get_feat(lines[2])

            if "T4SS" in lines[4]:
                mpf = lines[4].split("/")[-1].split("_")[1]
            else:
                mpf = ""
            if "Relaxase@" in tags:
                mob = tags.split("@")[1]
            else:
                mob = ""
            ICEtag = "ICE" + lines[5].split("_")[-1]

            if "IME" in lines[5]:
                ICEtag = "IME" + lines[5].split("_")[-1]
            elif "AICE" in lines[5]:
                ICEtag = "AICE" + lines[5].split("_")[-1]
            else:
                ICEtag = "ICE" + lines[5].split("_")[-1]

            genes_icedict.setdefault(ICEtag, {})[gbname] = tags

            if ICEtag not in infodict:
                infodict[ICEtag] = {"mob": [], "mpf": []}
            if mob not in infodict[ICEtag]["mob"]:
                if mob:
                    infodict[ICEtag]["mob"].append(mob)
            if mpf not in infodict[ICEtag]["mpf"]:
                infodict[ICEtag]["mpf"].append(mpf)

    drs_ice_dict = {}
    listgff = [trnadict, posdict, header, totalnum_dict, locusdict]
    for key, value in genes_icedict.items():
        (
            contig,
            myDR1,
            myDR2,
            myDR3,
            myDR4,
            fICE,
            eICE,
            finalstart,
            finalend,
            posdict,
            header,
            trnalist,
            locusdict,
        ) = merge_tRNA(key, value, dr_dict, listgff, prots_contigs)

        drs_ice_dict[key] = [
            contig,
            myDR1,
            myDR2,
            myDR3,
            myDR4,
            fICE,
            eICE,
            finalstart,
            finalend,
            trnalist,
            locusdict,
        ]

    return drs_ice_dict, genes_icedict, posdict, header, infodict


def get_map(drs_ice_dict, genes_icedict, posdict, header, infodict, assembly, prefix):

    ice_output_file = prefix + "_ices.tsv"
    genes_output_file = prefix + "_ice_genes.tsv"

    with open(genes_output_file, "w") as genes_output:
        genes_header = "\t".join(
            ["contig", "ice_id", "gene_id", "gene_coords", "gene_annot", "ice_feature"]
        )
        genes_output.write(genes_header + "\n")
        # Finding gene context per ICE
        ices_summary = []
        for key, value in drs_ice_dict.items():
            genelist = []
            [
                contig,
                myDR1,
                myDR2,
                myDR3,
                myDR4,
                fICE,
                eICE,
                finalstart,
                finalend,
                trnalist,
                locusdict,
            ] = value
            regi = contig + "_" + key
            start = finalstart
            while start < fICE:
                gene = zill(header, start)
                s, e, strand, pro = posdict[gene]
                pos = str(s) + ".." + str(e) + " [" + strand + "]"
                feature = "Flank"
                product = pro
                start += 1
                content = {
                    "gene": locusdict[gene],
                    "pos": pos,
                    "prod": product,
                    "featu": feature,
                }
                genelist.append(content)

            mov = fICE
            while mov <= eICE:
                gene = zill(header, mov)
                s, e, strand, pro = posdict[gene]
                pos = str(s) + ".." + str(e) + " [" + strand + "]"

                if gene in genes_icedict[key]:
                    [feature, pro11] = genes_icedict[key][gene].split("@")
                else:
                    feature, pro11 = "", ""

                if pro11:
                    if pro == "hypothetical protein":
                        product = pro11
                    else:
                        product = pro + ", " + pro11
                else:
                    product = pro

                mov += 1
                content = {
                    "gene": locusdict[gene],
                    "pos": pos,
                    "prod": product,
                    "featu": feature,
                }
                genelist.append(content)

            while mov <= finalend:
                gene = zill(header, mov)
                s, e, strand, pro = posdict[gene]
                pos = str(s) + ".." + str(e) + "[" + strand + "]"
                feature = "Flank"
                product = pro

                mov += 1
                content = {
                    "gene": locusdict[gene],
                    "pos": pos,
                    "prod": product,
                    "featu": feature,
                }
                genelist.append(content)

            # Print gene information into output file
            for gene in genelist:
                raw_product = gene["prod"]
                raw_product = re.sub(r"^, ", "", raw_product)
                to_print = "\t".join(
                    [
                        contig,
                        key,
                        gene["gene"],
                        gene["pos"],
                        raw_product,
                        gene["featu"],
                    ]
                )
                genes_output.write(to_print + "\n")

            # Integrating ICEs meta info
            sgene = zill(header, fICE)
            egene = zill(header, eICE)
            s1, e1, strand1, pro1 = posdict[sgene]
            s2, e2, strand2, pro2 = posdict[egene]
            if myDR1 == "0":
                myDR1 = "1"

            gcc = get_gc(get_sequence(assembly, contig, int(myDR1) - 1, int(myDR4)))

            if myDR2:
                DR1 = get_sequence(assembly, contig, int(myDR1), int(myDR2))
                DR2 = get_sequence(assembly, contig, int(myDR3), int(myDR4))
                DRw = (
                    "attL:"
                    + myDR1
                    + ".."
                    + myDR2
                    + "("
                    + DR1
                    + "),"
                    + "attR:"
                    + myDR3
                    + ".."
                    + myDR4
                    + "("
                    + DR2
                    + ")"
                )
            else:
                DRw = "-"
            if trnalist:
                trnaout = (
                    trnalist[0][3]
                    + "("
                    + str(trnalist[0][0])
                    + ".."
                    + str(trnalist[0][1])
                    + ")["
                    + trnalist[0][2]
                    + "]"
                )
            else:
                trnaout = "-"

            if "IME" in regi:
                typeIE = "IME"
            elif "AICE" in regi:
                typeIE = "AICE"
            else:
                typeIE = "T4SS-type ICE"

            ICEinfo = {
                "contig": contig,
                "ice_id": regi,
                "type": typeIE,
                "location": str(myDR1) + ".." + str(myDR4),
                "length": str(int(myDR4) - int(myDR1) + 1),
                "gc": gcc,
                "drs": DRw,
                "relaxase": ",".join(infodict[key]["mob"]),
                "mating_sys": ",".join(infodict[key]["mpf"]),
                "tRNA": trnaout,
            }
            ices_summary.append(ICEinfo)

    with open(ice_output_file, "w") as ice_output:
        summ_header = "\t".join(
            [
                "contig",
                "ice_id",
                "ice_type",
                "ice_location",
                "ice_length",
                "gc_content",
                "direct_repeats",
                "relaxase_type",
                "mating_pair_formation_systems",
                "close_to_RNA",
            ]
        )
        ice_output.write(summ_header + "\n")
        for ice in ices_summary:
            to_print = "\t".join(
                [
                    ice["contig"],
                    ice["ice_id"],
                    ice["type"],
                    ice["location"],
                    ice["length"],
                    ice["gc"],
                    ice["drs"],
                    ice["relaxase"],
                    ice["mating_sys"],
                    ice["tRNA"],
                ]
            )
            ice_output.write(to_print + "\n")


def get_sequence(fasta_file, contig, start, end):
    for record in SeqIO.parse(fasta_file, "fasta"):
        if str(record.id) == contig:
            sequence = record.seq[start:end]
    return str(sequence)


def get_gc(seq):
    seq = seq.upper()
    gc_count = seq.count("G") + seq.count("C")
    if len(seq) > 0:
        gcs = str("{:.2f}".format(gc_count / len(seq)))
    else:
        gcs = "0.00"
    return gcs


def main():
    parser = argparse.ArgumentParser(
        description="This script process the following input files: the (meta)genomic assembly, macsyfinder results, vmatch direct repeats predictions, prodigal GFF file integrated with aragorn annotation of tRNAs, and proteins annotation vs uniprotKB as generated by prokka. The outputs are two tsv files, a summary of ICEs predictions and the info per gene"
    )
    parser.add_argument(
        "--assembly",
        type=str,
        help="Assembly file in fasta format",
        required=True,
    )
    parser.add_argument(
        "--gff_file",
        type=str,
        help="Prodigal gff file with aragorn trnas interleaved (merged_gff)",
        required=True,
    )
    parser.add_argument(
        "--macsyfinder_out",
        type=str,
        help="Result of macsyfinder tool in tsv format (all_systems.tsv)",
        required=True,
    )
    parser.add_argument(
        "--uniprot_annot",
        type=str,
        help="Result of blastp vs uniprotkb-sp in prokka style",
        required=True,
    )
    parser.add_argument(
        "--drs_tsv",
        type=str,
        help="Result of Vmatch tool in tsv format",
        required=True,
    )
    parser.add_argument(
        "--prefix",
        type=str,
        help="Output files prefix",
        required=True,
    )
    args = parser.parse_args()

    # We have to generate the the following data structure according with single.py script
    # listgff = [trnadict,posdict,header,totalnum,locusdict]
    # We are recycling as much as possible code to avoid introducing bugs

    # Parsing uniprotkb blast result
    print("Parsing uniprotkb blast result...")
    uniprot_annot_dict = parse_blast_uniprot(args.uniprot_annot)

    # Parsing gff annotation to save the protein names correspondance with per-protein annotation results: macsyfinder and uniprotkb
    print("Parsing merged gff file...")
    names_map, trnadict, posdict, header, totalnum_dict, locusdict, prots_contigs = parse_merged_gff(
        args.gff_file, uniprot_annot_dict
    )

    # Parsing direct repeats prediction:
    print("Parsing direct repeats prediction file...")
    dr_dict = get_DR(args.drs_tsv)

    # Parsing macsyfinder results and processing ICEs
    print("Parsing macsyfinder results...")
    drs_ice_dict, genes_icedict, posdict, header, infodict = get_ICE(
        args.macsyfinder_out,
        dr_dict,
        trnadict,
        posdict,
        header,
        totalnum_dict,
        locusdict,
        names_map,
        prots_contigs,
    )

    # Integrating per gene info and ICE summaries
    print("Integrating and writing outputs...")
    get_map(
        drs_ice_dict,
        genes_icedict,
        posdict,
        header,
        infodict,
        args.assembly,
        args.prefix,
    )

    print("Processing complete!")


if __name__ == "__main__":
    main()
