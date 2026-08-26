#!/usr/bin/env python3
"""
Build a CGView (Proksee) JSON from a MAP mobilome GFF for one plasmid.
Import the resulting .json in Proksee via  Map > Import > CGView JSON.

Usage:  python build_cgview_json.py MGYG000509863_2.gff [more.gff ...]
Output: <name>_cgview.json  next to each input.
"""
import sys, os, json, re

# palette (matches the iTOL tree: ARG red, VF purple, plasmids blue)
COL = {
    "ARG":                  "#e31a1c",
    "VF":                   "#6a3d9a",
    "Other CDS":            "#bdbdbd",
    "BGC":                  "#ff7f00",
    "Compositional outlier":"#1b9e77",
    "Insertion sequence":   "#8c6d31",
}
DECOR = {"ARG":"arrow","VF":"arrow","Other CDS":"arrow",
         "BGC":"arc","Compositional outlier":"arc","Insertion sequence":"arc","Repeat":"arc"}

def attrs(col9):
    d={}
    for kv in col9.strip().split(';'):
        if '=' in kv:
            k,v=kv.split('=',1); d[k]=v
    return d

def cds_legend(tags):
    t=set(tags.split(',')) if tags else set()
    if 'arg' in t: return "ARG"
    if 'vf'  in t: return "VF"
    return "Other CDS"

SRC = {  # feature type -> (track source group)
    "CDS":"genes",
    "bgc":"BGC",
    "compositional_outlier":"MGE",
    "insertion_sequence":"MGE",
}
TYPE_LEGEND = {
    "bgc":"BGC",
    "compositional_outlier":"Compositional outlier",
    "insertion_sequence":"Insertion sequence",
}

def build(gff):
    name=os.path.basename(gff).replace('.gff','')
    length=None; feats=[]; used=set()
    with open(gff) as fh:
        for line in fh:
            if line.startswith('##sequence-region'):
                parts=line.split()
                if len(parts)>=4: length=int(parts[3])
                continue
            if line.startswith('#') or not line.strip(): continue
            c=line.rstrip('\n').split('\t')
            if len(c)<9 or not c[2]: continue
            typ=c[2]
            if typ=="plasmid":
                if length is None: length=int(c[4]); continue
                continue
            a=attrs(c[8]); start=int(c[3]); stop=int(c[4])
            strand=1 if c[6]=='+' else (-1 if c[6]=='-' else 1)
            if typ=="CDS":
                leg=cds_legend(a.get('pathofact2',''))
                fid=a.get('ID','CDS')
                feats.append({"name":"","type":"CDS","start":start,"stop":stop,
                              "strand":strand,"source":"genes","legend":leg})
                used.add(leg)
            elif typ in SRC:
                leg=TYPE_LEGEND[typ]
                feats.append({"name":typ.replace('_',' '),"type":typ,"start":start,"stop":stop,
                              "strand":strand,"source":SRC[typ],"legend":leg})
                used.add(leg)
    if length is None:
        length=max(f["stop"] for f in feats)

    legend_order=["ARG","VF","Other CDS","BGC","Compositional outlier","Insertion sequence"]
    legend_items=[{"name":n,"swatchColor":COL[n],"decoration":DECOR[n]}
                  for n in legend_order if n in used]

    tracks=[]
    if any(f["source"]=="BGC" for f in feats):
        tracks.append({"name":"BGCs","position":"outside","separateFeaturesBy":"none",
                       "thicknessRatio":1,"dataType":"feature","dataMethod":"source","dataKeys":"BGC"})
    tracks.append({"name":"Cargo CDS","position":"both","separateFeaturesBy":"strand",
                   "thicknessRatio":2,"dataType":"feature","dataMethod":"source","dataKeys":"genes"})
    if any(f["source"]=="MGE" for f in feats):
        tracks.append({"name":"Mobile elements","position":"inside","separateFeaturesBy":"none",
                       "thicknessRatio":1,"dataType":"feature","dataMethod":"source","dataKeys":"MGE"})

    cgv={"cgview":{
        "version":"1.7.0","name":name,"geneticCode":11,
        "settings":{"format":"circular","showShading":True},
        "backbone":{"color":"#bbbbbb","thickness":5},
        "ruler":{"font":"sans-serif, plain, 10"},
        "sequence":{"length":length},
        "legend":{"position":"top-right","defaultFont":"sans-serif, plain, 12","items":legend_items},
        "captions":[{"name":f"{name}  ({length:,} bp)","position":"bottom-center",
                     "font":"sans-serif, bold, 16","fontColor":"#333333"}],
        "tracks":tracks,
        "features":feats,
    }}
    out=os.path.join(os.path.dirname(os.path.abspath(gff)),name+"_cgview.json")
    with open(out,"w") as fh: json.dump(cgv,fh,indent=2)
    from collections import Counter
    cc=Counter(f["legend"] for f in feats)
    print(f"{name}: length={length:,} bp, {len(feats)} features -> {out}")
    print("   ", dict(cc))
    return out

for g in sys.argv[1:]:
    build(g)
