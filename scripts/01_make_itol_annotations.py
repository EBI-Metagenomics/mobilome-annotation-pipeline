#!/usr/bin/env python3
"""
Figure 2A - generate iTOL annotation datasets for the GTDB/IQ-TREE tree of the 30
plasmid-carrying genomes. Drag each output .txt onto the tree in iTOL.

Inputs : supplementary_1.xlsx  (sheet 'genomes': genome_id, genome_type, plasmids_num,
                                 vf_in_mobilome, arg_in_mobilome, lineage)
Outputs: itol_taxonomy_order.txt, itol_assembly_type.txt, itol_plasmid_count.txt,
         itol_arg_mobilome.txt, itol_vf_mobilome.txt, itol_selected_hosts.txt, itol_labels.txt

Usage:  python 01_make_itol_annotations.py --xlsx supplementary_1.xlsx --out itol_annotations [--suffix .fna]
Note:   tree leaf labels are <genome_id><suffix>; GTDB-Tk names leaves '<genome_id>.fna'.
"""
import argparse, os, openpyxl

ORDER_COLORS={'Burkholderiales':'#1f78b4','Enterobacterales':'#33a02c','Rhizobiales':'#a6cee3',
 'Pseudomonadales':'#b2df8a','Mycobacteriales':'#ff7f00','Streptomycetales':'#6a3d9a',
 'Streptosporangiales':'#b15928','Paenibacillales':'#fb9a99'}
SELECTED={'MGYG000509863':'Pantoea agglomerans','MGYG000499014':'Rahnella variigena'}

def rank(l,p):
    for tok in str(l).split(';'):
        if tok.startswith(p+'__'): return tok[3:] or '(unclassified)'
    return '(none)'

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--xlsx', default='supplementary_1.xlsx')
    ap.add_argument('--out',  default='itol_annotations')
    ap.add_argument('--suffix', default='.fna', help='appended to genome_id to match tree leaf labels')
    a=ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    S=a.suffix
    g=list(openpyxl.load_workbook(a.xlsx, data_only=True)['genomes'].iter_rows(values_only=True))[1:]
    def lid(x): return f"{x}{S}"

    orders=[o for o in ORDER_COLORS if any(rank(r[6],'o')==o for r in g)]
    with open(f'{a.out}/itol_taxonomy_order.txt','w') as f:
        f.write("DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\tTaxonomic order\nCOLOR\t#333333\n")
        f.write("COLOR_BRANCHES\t1\nSTRIP_WIDTH\t30\nMARGIN\t4\nSHOW_INTERNAL\t0\n")
        f.write("LEGEND_TITLE\tTaxonomic order\nLEGEND_SHAPES\t"+"\t".join("1" for _ in orders)+"\n")
        f.write("LEGEND_COLORS\t"+"\t".join(ORDER_COLORS[o] for o in orders)+"\n")
        f.write("LEGEND_LABELS\t"+"\t".join(orders)+"\nDATA\n")
        for r in g: o=rank(r[6],'o'); f.write(f"{lid(r[0])}\t{ORDER_COLORS.get(o,'#999999')}\t{o}\n")

    with open(f'{a.out}/itol_assembly_type.txt','w') as f:
        f.write("DATASET_COLORSTRIP\nSEPARATOR TAB\nDATASET_LABEL\tAssembly type\nCOLOR\t#444444\n")
        f.write("STRIP_WIDTH\t20\nMARGIN\t4\nSHOW_INTERNAL\t0\nLEGEND_TITLE\tAssembly type\n")
        f.write("LEGEND_SHAPES\t1\t1\nLEGEND_COLORS\t#444444\t#bdbdbd\nLEGEND_LABELS\tMAG\tIsolate\nDATA\n")
        for r in g: f.write(f"{lid(r[0])}\t{'#444444' if r[1]=='MAG' else '#bdbdbd'}\t{r[1]}\n")

    def bar(fn,label,color,idx):
        with open(f'{a.out}/{fn}','w') as f:
            f.write(f"DATASET_SIMPLEBAR\nSEPARATOR TAB\nDATASET_LABEL\t{label}\nCOLOR\t{color}\n")
            f.write("WIDTH\t120\nMARGIN\t6\nHEIGHT_FACTOR\t0.9\nBORDER_WIDTH\t0.5\nDATA\n")
            for r in g: f.write(f"{lid(r[0])}\t{int(r[idx])}\n")
    bar('itol_plasmid_count.txt','Plasmids per genome','#1f78b4',3)
    bar('itol_arg_mobilome.txt','ARG genes in mobilome','#e31a1c',5)
    bar('itol_vf_mobilome.txt','VF genes in mobilome','#6a3d9a',4)

    with open(f'{a.out}/itol_selected_hosts.txt','w') as f:
        f.write("DATASET_SYMBOL\nSEPARATOR TAB\nDATASET_LABEL\tHosts of selected plasmids\nCOLOR\t#e31a1c\n")
        f.write("MAXIMUM_SIZE\t18\nLEGEND_TITLE\tSelected for deep-dive\nLEGEND_SHAPES\t3\n")
        f.write("LEGEND_COLORS\t#e31a1c\nLEGEND_LABELS\tSelected plasmid host\nDATA\n")
        for gid,sp in SELECTED.items(): f.write(f"{lid(gid)}\t3\t16\t#e31a1c\t1\t1\t{sp}\n")

    with open(f'{a.out}/itol_labels.txt','w') as f:
        f.write("LABELS\nSEPARATOR TAB\nDATA\n")
        for r in g:
            sp=rank(r[6],'s'); genus=rank(r[6],'g')
            name=sp if sp not in ('','(none)','(unclassified)') else genus
            f.write(f"{lid(r[0])}\t{r[0]} | {name}\n")
    print("wrote iTOL annotation files to", a.out)

if __name__=='__main__': main()
