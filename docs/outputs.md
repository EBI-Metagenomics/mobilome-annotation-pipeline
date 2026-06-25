# MAP Outputs

Results are written to `--outdir` (default: `results/`).

## Directory structure

```
sample/
├── sample_combined_report.tsv
├── sample_discarded_mge.txt
├── sample_mobilome.fasta.gz
├── sample_overlap_report.txt
├── gff/
│   ├── sample_mobilome.gff.gz
│   │   # When user proteins+GFF are provided, the three derived files are named
│   │   # sample_user_mobilome_{clean,extra,full}.gff.gz (baseline: the user GFF).
│   │   # Otherwise they are named sample_mobilome_{clean,extra,full}.gff.gz
│   │   # (baseline: the Prodigal/tRNA genes GFF). Each has matching .csi and .gzi indexes.
│   ├── sample_user_mobilome_clean.gff.gz      # mobilome + CDSs that fall inside an MGE
│   ├── sample_user_mobilome_clean.gff.gz.csi
│   ├── sample_user_mobilome_clean.gff.gz.gzi
│   ├── sample_user_mobilome_extra.gff.gz      # mobilome + functionally annotated passengers
│   ├── sample_user_mobilome_extra.gff.gz.csi
│   ├── sample_user_mobilome_extra.gff.gz.gzi
│   ├── sample_user_mobilome_full.gff.gz       # mobilome + every CDS from the genes GFF
│   ├── sample_user_mobilome_full.gff.gz.csi
│   └── sample_user_mobilome_full.gff.gz.gzi
├── prediction/
│   ├── amr_genes/
│   │   ├── integrated_sample.gff
│   │   ├── amrfinderplus/
│   │   │   └── sample.tsv
│   │   ├── deeparg/
│   │   │   └── sample.mapping.ARG
│   │   └── rgi/
│   │       └── sample.txt
│   ├── bgcs/
│   │   ├── sample_bgcs.gff
│   │   ├── sample_bgcs.json
│   │   ├── antismash/
│   │   │   └── sample_antismash.gff
│   │   ├── gecco/
│   │   │   └── sample.gff
│   │   └── sanntis/
│   │       └── sample_sanntis.gff.gz
│   ├── compositional_outliers/
│   │   └── sample_100kb_contigs.1.bed
│   ├── genomad/
│   │   ├── sample_5kb_contigs_plasmid_summary.tsv
│   │   └── sample_5kb_contigs_virus_summary.tsv
│   ├── icefinder2lite/
│   │   ├── sample_ice_genes.tsv
│   │   └── sample_ices.tsv
│   ├── integronfinder/
│   │   ├── sample_100kb_contigs.summary
│   │   └── contig_1.gbk
│   ├── interproscan/
│   │   └── sample.tsv.gz
│   ├── isescan/
│   │   └── sample_1kb_contigs.fasta.tsv
│   ├── virify_filter/
│   │   └── sample_virify_hq.gff
│   └── virulence/
│       └── sample_pathofact2.gff
└── preprocessing/
    ├── sample_1kb_contigs.fasta
    ├── sample_5kb_contigs.fasta
    ├── sample_100kb_contigs.fasta
    └── sample_contigID.map
```

## Derived mobilome GFFs (clean / extra / full)

The `gff/` directory holds the mobilome GFF (`sample_mobilome.gff.gz`) plus three derived
GFFs produced by `gff_mapping.py`. Each merges the mobilome features with a baseline genes
GFF: the user-provided proteins GFF when one is supplied (outputs named
`sample_user_mobilome_*`), otherwise the pipeline's Prodigal/tRNA genes GFF (outputs named
`sample_mobilome_*`). All three carry the **same set of mobilome features**; they differ in
which genes from the baseline GFF they additionally include:

- **`full`** — the superset: every mobilome feature plus **every** feature from the baseline
  genes GFF. It preserves the baseline GFF's header. Where available, each CDS gains its
  VIRify ViPhOG attributes (`viphog` / `viphog_taxonomy`) and a `pathofact2=<summary_string>`
  attribute (from the [combined report](#combined-report)).
- **`clean`** — mobilome features plus only the "passenger" CDSs that fall **inside** a mobile
  element (>75% of the CDS length overlapping an MGE on the same contig). Each passenger CDS
  additionally carries an `mge_location=` attribute, plus ViPhOG and `pathofact2=` attributes
  when present.
- **`extra`** — a subset of `clean`: the mobilome features plus only those **passenger** CDSs
  that also carry a **functional annotation** — a VIRify ViPhOG hit (`viphog` /
  `viphog_taxonomy`) and/or a `pathofact2=` summary. Rows use the same format as in `clean`
  (including `mge_location=`). A passenger CDS with no functional annotation appears in
  `clean` but not in `extra`; a functionally-annotated CDS that is not a passenger appears in
  `full` but not in `extra`.

`clean` and `extra` carry a minimal `##gff-version 3` header; `full` preserves the baseline
GFF's full header. All three are bgzip-compressed with matching `.csi` and `.gzi` indexes.

## Discarded predictions

`sample_discarded_mge.txt` lists predictions excluded during integration and the reason:

1. `mge < 500bp` — discarded by length
2. `no_cds` — no coding sequences within the prediction
3. `tRNAs_in_window` — tRNA genes present inside a compositional outlier
4. `CO_overlap_with_MGE` — compositional outlier overlapping another MGE

`sample_overlap_report.txt` lists long-MGEs with overlapping coordinates. No predictions are discarded for this reason; it is informational only.

## Feature identifiers and types

Mobilome feature IDs follow this pattern: `contig_id|mge_type-start:end`.

GFF feature types and their Sequence Ontology mappings:

| Type | SO ID | Element | Tool |
|---|---|---|---|
| `insertion_sequence` | [SO:0000973](http://www.sequenceontology.org/browser/current_svn/term/SO:0000973) | Insertion sequence | ISEScan |
| `inverted_repeat_element` | [SO:0000481](http://www.sequenceontology.org/browser/current_svn/term/SO:0000481) | Inverted repeat flanking IS or COD | ISEScan, MAP |
| `integron` | [SO:0000365](http://www.sequenceontology.org/browser/current_svn/term/SO:0000365) | Integrative mobilisable element | IntegronFinder, ICEfinder |
| `attC_site` | [SO:0000950](http://www.sequenceontology.org/browser/current_svn/term/SO:0000950) | Integron integration site | IntegronFinder |
| `conjugative_integron` | [SO:0000371](http://www.sequenceontology.org/browser/current_svn/term/SO:0000371) | Integrative conjugative element | ICEfinder |
| `direct_repeat` | [SO:0000314](http://www.sequenceontology.org/browser/current_svn/term/SO:0000314) | Flanking region of mobilisable element | ICEfinder, MAP |
| `prophage` | [SO:0001006](http://www.sequenceontology.org/browser/current_svn/term/SO:0001006) | Temperate phage | geNomad, VIRify |
| `viral_sequence` | [SO:0001041](http://www.sequenceontology.org/browser/current_svn/term/SO:0001041) | Viral genome fragment | geNomad, VIRify |
| `plasmid` | [SO:0000155](http://www.sequenceontology.org/browser/current_svn/term/SO:0000155) | Plasmid | geNomad |
| `compositional_outlier` | — | Non-autonomous element on contigs ≥100 kb | MAP |

## Combined report

`sample_combined_report.tsv` is produced by the PathFact2 integrator step, which merges PathoFact2, AMR, mobilome, biosynthetic gene cluster (BGC), and InterProScan (IPS) annotations.

**Seed proteins** — rows in the report — come from PathoFact2 and/or AMR GFF files. If both inputs are absent or empty, no report is generated. BGC annotations are only resolved for seed proteins; BGC-only or MGE-only proteins are not included as rows. MGE assignment requires ≥90% CDS overlap with a mobilome feature on the same contig. IPS input is optional; only SignalP entries are retained. Missing values are reported as `-`.

Each protein's `summary_string` is also carried into the derived `*_mobilome_{clean,extra,full}.gff.gz` files as a `pathofact2=<summary_string>` attribute on the matching CDS.

| Column | Description |
|---|---|
| `protein_id` | Protein identifier from Prodigal or the user-provided GFF |
| `contig_id` | Contig the protein is located on |
| `summary_string` | Condensed annotation tags for the gene, comma-joined in fixed order `vf,arg,mge,bgc`: `vf` (virulence — VFDB hit or PATHOFACT2 toxin/VF), `arg` (AMR), `mge` (within a mobile element), `bgc` (within a BGC). Every row carries at least `vf` or `arg` |
| `vfdb_hit` | Best VFDB hit accession (DIAMOND blastp), `-` if no hit |
| `vfdb_blastp_eval` | E-value of the VFDB hit |
| `pathofact2_tox_prob` | PATHOFACT2 toxin probability (0–1); `-` if PATHOFACT2 was not run |
| `pathofact2_vf_prob` | PATHOFACT2 virulence factor probability (0–1); `-` if PATHOFACT2 was not run |
| `cdd_annotation` | CDD domain annotation from local CDsearch or InterProScan |
| `amr_drug_class` | Drug class of the AMR hit (e.g. `beta-lactam`, `tetracycline`) |
| `amr_tool` | Tool reporting the AMR hit (`amrfinderplus`, `deeparg`, or `rgi`) |
| `amr_tool_ident` | Percent identity to the AMR reference sequence |
| `mge_type` | MGE type when the protein has ≥90% overlap with a mobilome feature (e.g. `prophage`, `insertion_sequence`); `-` otherwise |
| `bgc_type` | BGC class (e.g. `T3PKS`, `Saccharide`, `RiPP-like`); `-` if no BGC |
| `bgc_tools` | Tool(s) that called the BGC (`antismash`, `sanntis`, `gecco`; comma-separated if multiple) |
| `signalp` | SignalP annotation from InterProScan if an IPS TSV was provided or IPS was run internally with `--interpro_licensed_software true`; `-` otherwise |
