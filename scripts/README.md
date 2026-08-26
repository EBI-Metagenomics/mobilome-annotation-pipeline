# Figure-generation scripts (MAP paper, use-case figures)

Reproduces the analysis figures. Requires Python 3 with `openpyxl`, `pandas`,
`numpy`, `matplotlib`. Input data: `supplementary_1.xlsx` (genomes/plasmids +
cargo), the per-plasmid MGnify Branchwater exports (`*_bw.csv`), and the GTDB-Tk
`bac120` tree.

## Figure 2A - phylogeny + annotations
The tree was built from the 30 plasmid-carrying genomes with GTDB-Tk v2.7.1
(`identify` then `align`, default parameters, bac120 marker set), then IQ-TREE v2.2.0.3:

    gtdbtk identify --genome_dir genomes --out_dir gtdb_id_results --extension gz --cpus 8
    gtdbtk align --identify_dir gtdb_id_results/identify --out_dir gtdb_align_results --cpus 8
    iqtree2 -s gtdb_align_results/align/gtdbtk.bac120.user_msa.fasta.gz -m LG+R10 -pre gtdbtk_bac120

Decorations are generated with:

    python 01_make_itol_annotations.py --xlsx supplementary_1.xlsx --out itol_annotations

Load the tree in iTOL (https://itol.embl.de) and drag each `itol_annotations/*.txt`
file onto it (order strip, MAG/isolate strip, plasmid/ARG/VF bars, selected-host stars).

## Figure 2B - biome distribution heatmap
Build the tables (raw + BioProject-collapsed) and Supplementary Table 2:

    python 02_branchwater_biome_matrix.py --bwdir <dir with *_bw.csv> --xlsx supplementary_1.xlsx --out .

Then render the heatmap (reads P2_plasmid_biome_matrix_bioproject.csv from its own dir):

    python 03_plot_biome_heatmap.py            # -> Figure1B_biome_heatmap.pdf

## Figure 2C - plasmid genomic context (Proksee/CGView)
From the MAP mobilome GFF of each selected plasmid:

    python 04_build_cgview_json.py MGYG000509863_2.gff MGYG000499014_6.gff

Import each `*_cgview.json` in Proksee (https://proksee.ca): Map > Import > CGView JSON.

## Notes
- Cargo bars in Figure 2A are per-genome mobilome totals; the two selected plasmids'
  per-plasmid counts in the text/legend differ by design.
- Branchwater biome counts are collapsed to distinct BioProjects to avoid sampling-effort bias.
