#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2025-2026 EMBL - European Bioinformatics Institute
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


CHECKV_FIELDS = [
    "checkv_kmer_freq",
    "checkv_miuvig_quality",
    "checkv_provirus",
    "checkv_quality",
    "checkv_viral_genes"
]

def mge_data_parser(mge_data):
```
    Extracting genomad predictions in three different structures:
        viral_dic -> Viral genomes in a single contig
        prophages_dic -> Prophages
        plasmids_list -> Plasmid contigs
```
    plasmids_list = []
    prophages_dic, prophages_ids, viral_dic = {}, {}, {}
    for mge in mge_data:
        contig, description, coord = mge_data[mge]
        prefix = mge.split("_")[0]
        if prefix == "vir1":
            if "viral_sequence" in description:
                viral_dic[contig] = mge
            elif "prophage" in description:
                composite_key = (contig, coord)
                prophages_ids[composite_key] = mge
                if contig in prophages_dic:
                    prophages_dic[contig].append(coord)
                else:
                    prophages_dic[contig] = [coord]
        if prefix == "plas":
            plasmids_list.append(contig)

    return (viral_dic, prophages_dic, prophages_ids, plasmids_list)


def virify_reader(virify_gff, inv_names_equiv, mge_data):
    ```
    Parsing virify predictions. We are storing the virify_prots
    and virify_with_viphogs to use them to replace genomad overlapping predictions.
    ```
    virify_predictions, virify_prots = {}, {}
    mge_counter = 0
    virify_with_viphogs = set()
    contig_to_mge_id = {}
    names_equiv = {v: k for k, v in inv_names_equiv.items()}

    with open(virify_gff, "r") as input_table:
        for line in input_table:
            line_l = line.rstrip().split("\t")
            # Annotation lines have exactly 9 columns
            if len(line_l) == 9:
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
                attr = attr.replace("gbkey=mobile_element;", "")
                contig = inv_names_equiv[contig]
                coord = (int(start), int(end))

                if seq_source == "VIRify":
                    mge_counter += 1
                    mge_id = "vir2_" + str(mge_counter)
                    composite_val = (contig, attr, coord)
                    virify_predictions[mge_id] = composite_val
                    contig_to_mge_id[contig] = mge_id

                # Saving protein predictions having viphog annotation
                # ID=<cds_id>;virify_quality=HC;gbkey=CDS;viphog=ViPhOG18043;viphog_taxonomy=Andromedavirus
                elif seq_type == "CDS" and 'viphog' in attr:
                    parent_mge_id = contig_to_mge_id.get(contig)
                    if parent_mge_id:
                        virify_with_viphogs.add(parent_mge_id)
                    (
                        gene_id,
                        virify_quality,
                        gbkey,
                        viphog,
                        viphog_taxonomy,
                    ) = attr.split(";")
                    prot_viphog = viphog + ";" + viphog_taxonomy
                    prot_location = (contig, int(start), int(end))
                    virify_prots[prot_location] = prot_viphog

    ## Parsing the mge_data to retrieve genomad predictions
    (viral_dic, prophages_dic, prophages_ids, plasmids_list) = mge_data_parser(mge_data)

    ## Removing redundancy on viral genomes to keep only one entry as plasmid_phage
    ## This applies only for viral genomes, not for prophages
    to_discard = []
    virify_plasmids = {}

    for phage in virify_predictions:
        v_contig, v_description, v_coord = virify_predictions[phage]

        # Catching phage-plasmids in viral predictions. No prophages
        if (
            v_description.split(";")[0].split("|")[1] == "viral_sequence"
            and v_contig in plasmids_list
        ):
            virify_plasmids[v_contig] = virify_predictions[phage]
            to_discard.append(phage)

    # Removing viral-phages from virify list to avoid double entry in the final GFF
    for phage in to_discard:
        del virify_predictions[phage]

    print("Number of plasmid-phages detected: " + str(len(to_discard)))

    # Assessing redundancy between genomad and virify
    # We are replacing genomad with virify predictions just when overlappiing virify have viphog matches
    # and we are adding genomad taxonomy to virify overlapping predictions
    to_discard_genomad = []
    to_discard_virify = []
    for phage in virify_predictions:
        v_contig, v_description, v_coord = virify_predictions[phage]

        # Finding redundancy on viral genomes (whole contig)
        if v_description.split(";")[0].split("|")[1] == "viral_sequence":
            if v_contig in viral_dic:
                # This is the structure of the mge_data dict:
                # mge_data[mge_id]  = (contig, description, coord)
                # This is the content of the description variable:
                # description = ( "mobile_element_type=viral_sequence;" + "taxonomy=" + taxonomy)

                if phage in virify_with_viphogs:
                    # Keeping virify prediction with genomad taxonomy
                    genomad_id = viral_dic[v_contig]
                    genomad_taxonomy_attr = mge_data[genomad_id][1].split(';')[1].replace('taxonomy=', 'genomad_taxonomy=')
                    v_extended_description = v_description + ";" + genomad_taxonomy_attr
                    virify_predictions[phage] = (v_contig, v_extended_description, v_coord)
                    to_discard_genomad.append(genomad_id)
                else:
                    # Keeping genomad prediction
                    to_discard_virify.append(phage)

            if v_contig in prophages_dic:
                # This means that virify predicted a viral genome and genomad predicted prophage(s)
                if phage in virify_with_viphogs:
                    # Keeping virify viral genome and discarding all the genomad prophages
                    # We are not keeping genomad taxonomy here as matches won't be identical 
                    for g_coord_pair in prophages_dic[v_contig]:
                        g_start = g_coord_pair[0]
                        g_end = g_coord_pair[1]
                        g_id = prophages_ids[(v_contig, (g_start, g_end))]
                        to_discard_genomad.append(g_id)
                else:
                    # Discarding virify viral genome due to lack of viphogs and overlapping with genomad phages
                    to_discard_virify.append(phage)

        # Finding redundancy on prophages predicted by genomad and virify
        elif v_description.split(";")[0].split("|")[1].split("-")[0] == "prophage":
            v_start = v_coord[0]
            v_end = v_coord[1]
            v_len = v_end - v_start
            v_range = range(v_start, v_end + 1)

            if v_contig in prophages_dic:
                for g_coord_pair in prophages_dic[v_contig]:
                    g_start = g_coord_pair[0]
                    g_end = g_coord_pair[1]
                    g_id = prophages_ids[(v_contig, (g_start, g_end))]
                    g_len = g_end - g_start
                    g_range = range(g_start, g_end + 1)
                    intersection = len(list(set(v_range) & set(g_range)))

                    if intersection > 0:
                        v_cov = float(intersection) / float(v_len)
                        g_cov = float(intersection) / float(g_len)

                        if any([v_cov > 0.25, g_cov > 0.25]):
                            if phage in virify_with_viphogs:
                                # We keep virify prophage and discard genomad's
                                to_discard_genomad.append(g_id)
                                # We append genomad taxonomy if predictions are identical
                                if g_cov == 1 and v_cov == 1:
                                    genomad_taxonomy_attr = mge_data[g_id][1].split(';')[1].replace('taxonomy=', 'genomad_taxonomy=')
                                    v_extended_description = v_description + ";" + genomad_taxonomy_attr
                                    virify_predictions[phage] = (v_contig, v_extended_description, v_coord)
                            else:
                                # Overlapping virify predictions without viphogs are discarded
                                to_discard_virify.append(phage)

            if v_contig in viral_dic:
                # This means that virify predicted a prophage in a contig that genomad predicted as a viral genome
                # In this case we are discarding virify as the prophage will be fully contained in the viral genome
                to_discard_virify.append(phage)

    ## Removing genomad predictions saved on to_discard_genomad list from mge_data dictionary
    to_discard_genomad = list(set(to_discard_genomad))
    print("Number of geNomad predictions discarded: " + str(len(to_discard_genomad)))
    for phage_id in to_discard_genomad:
        if phage_id in mge_data:
            del mge_data[phage_id]

    ## Removing virify predictions saved in the to_discard_virify list from virify_predictions dict
    to_discard_virify = list(set(to_discard_virify))
    print("Number of VIRify predictions discarded: " + str(len(to_discard_virify)))
    for phage_id in to_discard_virify:
        if phage_id in virify_predictions:
            del virify_predictions[phage_id]

    ## Adding Virify predictions to mge_data
    print(
        "Number of VIRify predictions to be added: "
        + str(len(list(virify_predictions.keys())))
    )
    for phage in virify_predictions:
        description = virify_predictions[phage][1].split(";")
        # Removing the prediction ID. This will be regenerated later on. Indentical predictions to genomad will have genomad_taxonomy
        # ['ID=contig_id|viral_sequence', 'virify_quality=HC', 'mobile_element_type=phage_linear', 'checkv_provirus=No', 'checkv_quality=Medium-quality', 'checkv_miuvig_quality=Genome-fragment', 'checkv_kmer_freq=1.0', 'checkv_viral_genes=11', 'virify_taxonomy=Bruynoghevirus']
        description.pop(0)
        description = ";".join(description)
        new_value = (
            virify_predictions[phage][0],
            description,
            virify_predictions[phage][2],
        )
        mge_data[phage] = new_value

    ## Labelling phage-plasmids
    # Storing the plasmid IDs to access values
    plas_phage_ids = []
    for mge in mge_data:
        contig, description, coord = mge_data[mge]
        prefix = mge.split("_")[0]
        if prefix == "plas":
            if contig in virify_plasmids:
                plas_phage_ids.append(mge)

    # Replacing descriptions in phage-plasmids
    mge_counter = 0
    useful_info = CHECKV_FIELDS + [
        "virify_taxonomy",
        "virify_quality",
    ]
    for pp in plas_phage_ids:
        mge_counter += 1
        old_contig, old_description, old_coord = mge_data[pp]
        new_met = "mobile_element_type=phage_plasmid"
        viral_info = virify_plasmids[old_contig][1].split(";")
        viral_desc = []
        for info in viral_info:
            key = info.split("=")[0]
            if key in useful_info:
                viral_desc.append(info)

        viral_desc = ";".join(viral_desc)
        new_desc = new_met + ";" + viral_desc
        new_val = (old_contig, new_desc, old_coord)
        new_mge_id = "phpl_" + str(mge_counter)
        del mge_data[pp]
        mge_data[new_mge_id] = new_val

    # for data in mge_data:
    #    print(data,mge_data[data])

    return (mge_data, virify_prots)
