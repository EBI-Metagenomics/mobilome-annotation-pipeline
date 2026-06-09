# MAP Outputs

Results are written to `--outdir` (default: `results/`).

## Directory structure

```
{sample}/
├── {sample}_combined_report.tsv
├── {sample}_discarded_mge.txt
├── {sample}_mobilome.fasta.gz
├── {sample}_overlap_report.txt
├── gff/
│   ├── {sample}_mobilome.gff.gz
│   ├── {sample}_user_mobilome_clean.gff.gz      # mobilome + matching CDSs
│   ├── {sample}_user_mobilome_clean.gff.gz.csi
│   ├── {sample}_user_mobilome_clean.gff.gz.gzi
│   ├── {sample}_user_mobilome_extra.gff.gz      # mobilome + VIRify ViPhOG-annotated genes
│   ├── {sample}_user_mobilome_extra.gff.gz.csi
│   ├── {sample}_user_mobilome_extra.gff.gz.gzi
│   ├── {sample}_user_mobilome_full.gff.gz       # mobilome + all features from user GFF
│   ├── {sample}_user_mobilome_full.gff.gz.csi
│   └── {sample}_user_mobilome_full.gff.gz.gzi
├── prediction/
│   ├── amr_genes/
│   ├── bgcs/
│   │   ├── {sample}_bgcs.gff
│   │   └── {sample}_bgcs.json
│   ├── compositional_outliers/
│   │   └── {sample}_100kb_contigs.1.bed
│   ├── genomad/
│   │   ├── {sample}_5kb_contigs_plasmid_summary.tsv
│   │   └── {sample}_5kb_contigs_virus_summary.tsv
│   ├── integronfinder/
│   │   ├── {sample}_100kb_contigs.summary
│   │   └── contig_1.gbk
│   ├── isescan/
│   │   └── {sample}_1kb_contigs.fasta.tsv
│   └── virulence/
│       └── {sample}_pathofact2.gff
└── preprocessing/
    ├── {sample}_1kb_contigs.fasta
    ├── {sample}_5kb_contigs.fasta
    ├── {sample}_100kb_contigs.fasta
    └── {sample}_contigID.map
```

## Discarded predictions

`{sample}_discarded_mge.txt` lists predictions excluded during integration and the reason:

1. `mge < 500bp` — discarded by length
2. `no_cds` — no coding sequences within the prediction
3. `tRNAs_in_window` — tRNA genes present inside a compositional outlier
4. `CO_overlap_with_MGE` — compositional outlier overlapping another MGE

`{sample}_overlap_report.txt` lists long-MGEs with overlapping coordinates. No predictions are discarded for this reason; it is informational only.

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

When functional annotation is enabled, `{sample}_combined_report.tsv` is produced. Each row is one protein. Only proteins with a virulence/toxin annotation from PATHOFACT2 or an AMR annotation are included. BGC context (`bgc_type`, `bgc_tools`) and mobilome context (`mge_type`) are added as additional columns when relevant, but BGC-only or MGE-only proteins are not included as rows.

| Column | Description |
|---|---|
| `protein_id` | Protein identifier from Prodigal or the user-provided GFF |
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
| `signalP` | SignalP annotation from InterProScan if an IPS TSV was provided or IPS was run internally with `--interpro_licensed_software true`; `-` otherwise |
