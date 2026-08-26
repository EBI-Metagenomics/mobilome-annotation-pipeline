#!/usr/bin/env python3
"""
Figure 2B data - parse MGnify Branchwater results for the 43 plasmid contigs and build
the biome-distribution table (distinct BioProjects per biome, containment >= 0.90).

Inputs : --bwdir  directory with one '<plasmid>_bw.csv' per plasmid (Branchwater export)
         --xlsx   supplementary_1.xlsx (cargo counts + genome types)
Outputs (in --out): P2_plasmid_biome_matrix.csv          (raw sample counts per biome)
                    P2_plasmid_biome_matrix_bioproject.csv (distinct BioProjects per biome)
                    supplementary_2.xlsx                  (formatted table + README)

Usage: python 02_branchwater_biome_matrix.py --bwdir results --xlsx supplementary_1.xlsx --out .
"""
import argparse, csv, glob, os
from collections import Counter, defaultdict
import openpyxl

BIOMES=['soil','plant/rhizosphere','food/silage','host-associated','wastewater','compost/digester',
        'freshwater','sediment','marine','peat/wetland','lichen/moss','other','unclassified']
KW=[('soil','soil'),('permafrost','soil'),('rhizosphere','plant/rhizosphere'),('root','plant/rhizosphere'),
 ('leaf','plant/rhizosphere'),('phyllosphere','plant/rhizosphere'),('seed','plant/rhizosphere'),
 ('plant','plant/rhizosphere'),('paddy','plant/rhizosphere'),('endophyte','plant/rhizosphere'),
 ('silage','food/silage'),('food','food/silage'),('cheese','food/silage'),('milk','food/silage'),
 ('wine','food/silage'),('fermented','food/silage'),('compost','compost/digester'),('digester','compost/digester'),
 ('biogas','compost/digester'),('bioreactor','compost/digester'),('wastewater','wastewater'),('sludge','wastewater'),
 ('sewage','wastewater'),('gut','host-associated'),('fec','host-associated'),('human','host-associated'),
 ('homo sapiens','host-associated'),('rumen','host-associated'),('bovine','host-associated'),('gallus','host-associated'),
 ('poultry','host-associated'),('clinical','host-associated'),('hospital','host-associated'),('bat','host-associated'),
 ('freshwater','freshwater'),('lake','freshwater'),('river','freshwater'),('groundwater','freshwater'),
 ('drinking','freshwater'),('pond','freshwater'),('sediment','sediment'),('marine','marine'),('seawater','marine'),
 ('seagrass','marine'),('coral','marine'),('estuary','marine'),('mangrove','marine'),('peat','peat/wetland'),
 ('bog','peat/wetland'),('wetland','peat/wetland'),('lichen','lichen/moss'),('moss','lichen/moss')]
def biome(o):
    o=(o or '').lower()
    if o in ('','uncalculated','na','not provided','metagenome'): return 'unclassified'
    for k,v in KW:
        if k in o: return v
    return 'other'
def load(f):
    l=open(f).read().splitlines(); hi=[i for i,x in enumerate(l) if x.startswith('acc,')][0]
    return list(csv.DictReader(l[hi:]))
def unit(r):
    bp=(r['bioproject'] or '').strip(); return bp if bp else 'acc:'+r['acc']

def main():
    ap=argparse.ArgumentParser()
    ap.add_argument('--bwdir', required=True); ap.add_argument('--xlsx', default='supplementary_1.xlsx')
    ap.add_argument('--out', default='.'); ap.add_argument('--min-containment', type=float, default=0.90)
    a=ap.parse_args()
    wb=openpyxl.load_workbook(a.xlsx, data_only=True)
    cargo={r[0]:{'len':int(r[1] or 0),'ARG':int(r[2] or 0),'VF':int(r[3] or 0),'BGC':int(r[4] or 0),'tot':int(r[5] or 0)}
           for r in list(wb['plasmids'].iter_rows(values_only=True))[1:]}
    gtype={r[0]:r[1] for r in list(wb['genomes'].iter_rows(values_only=True))[1:]}
    SEL={'MGYG000509863_2','MGYG000499014_6'}
    rows=[]
    for f in sorted(glob.glob(os.path.join(a.bwdir,'*_bw.csv'))):
        pid=os.path.basename(f).replace('_bw.csv','')
        hits=[r for r in load(f) if float(r['containment'])>=a.min_containment]
        samp=Counter(); bpset=defaultdict(set); allbp=set()
        for r in hits:
            b=biome(r['organism']); samp[b]+=1; u=unit(r); bpset[b].add(u); allbp.add(u)
        bp={b:len(bpset[b]) for b in BIOMES}
        g=pid.rsplit('_',1)[0]
        nsb=len([b for b in BIOMES if b not in('soil','unclassified') and bp[b]>0])
        rows.append(dict(pid=pid,g=g,type=gtype.get(g,''),cargo=cargo.get(pid,{}),
                         samp=samp,bp=bp,tot_s=len(hits),tot_bp=len(allbp),nsb=nsb))
    rows.sort(key=lambda x:(x['nsb'], sum(v for b,v in x['bp'].items() if b not in('soil','unclassified')), x['tot_bp']), reverse=True)

    def write_csv(path,key):
        with open(path,'w',newline='') as fh:
            w=csv.writer(fh)
            w.writerow(['plasmid','selected','host_genome','host_type','ARG','VF','BGC','total_ARG_VF',
                        'total_near_complete','nonsoil_biomes']+BIOMES)
            for r in rows:
                cg=r['cargo']; tot=r['tot_bp'] if key=='bp' else r['tot_s']
                w.writerow([r['pid'],'yes' if r['pid'] in SEL else '',r['g'],r['type'],
                            cg.get('ARG',''),cg.get('VF',''),cg.get('BGC',''),cg.get('tot',''),
                            tot,r['nsb']]+[r[key][b] for b in BIOMES])
    write_csv(os.path.join(a.out,'P2_plasmid_biome_matrix.csv'),'samp')
    write_csv(os.path.join(a.out,'P2_plasmid_biome_matrix_bioproject.csv'),'bp')

    out=openpyxl.Workbook(); ws=out.active; ws.title='branchwater_biome_distribution'
    ws.append(['plasmid_contig','host_genome','host_type','plasmid_len_bp','ARG_genes','VF_genes','BGC_genes',
               'total_ARG_VF','recovered_beyond_source','n_nonsoil_biomes','n_samples_ge0.90','n_bioprojects_ge0.90']+BIOMES)
    for r in rows:
        cg=r['cargo']
        ws.append([r['pid'],r['g'],r['type'],cg.get('len',''),cg.get('ARG',''),cg.get('VF',''),cg.get('BGC',''),
                   cg.get('tot',''),'yes' if r['tot_s']>0 else 'no',r['nsb'],r['tot_s'],r['tot_bp']]+[r['bp'][b] for b in BIOMES])
    rd=out.create_sheet('README')
    for line in ["Supplementary Table 2. Biome distribution of the 43 cargo-rich plasmids (MGnify Branchwater).","",
      "Branchwater v0.4.0 (database v2024-11-28); near-complete matches only (containment >= 0.90).",
      "Biome columns = distinct BioProjects per biome (repeated sampling within a study collapsed), matching Figure 2B.",
      "n_samples_ge0.90 = raw matching runs; n_bioprojects_ge0.90 = distinct BioProjects across all biomes.",
      "Cargo counts are from Supplementary Table 1.",
      "recovered_beyond_source = found near-complete in >=1 metagenome besides the source genome."]:
        rd.append([line])
    out.save(os.path.join(a.out,'supplementary_2.xlsx'))
    print("wrote P2 matrices and supplementary_2.xlsx to", a.out, "(%d plasmids)"%len(rows))

if __name__=='__main__': main()
