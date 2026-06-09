[![Run nf-tests for modules](https://github.com/EBI-Metagenomics/mobilome-annotation-pipeline/actions/workflows/test_modules.yml/badge.svg)](https://github.com/EBI-Metagenomics/mobilome-annotation-pipeline/actions/workflows/test_modules.yml)

# Mobilome Annotation Pipeline (MAP)

<p align="center" width="100%">
   <img src="media/mges_schema.png" width="100%"/>
</p>

Bacteria can acquire genetic material through horizontal gene transfer, allowing them to rapidly adapt to changing environmental conditions. These mobile genetic elements can be classified into three main categories: plasmids, phages, and integrative elements. Plasmids are mostly extrachromosomal; phages can be found extrachromosomal or as temperate phages (prophages); whereas integrons are stably inserted in the chromosome. Autonomous elements are those integrative elements capable of excising themselves from the chromosome and reintegrating elsewhere. They can use a transposase (like insertion sequences and transposons) or an integrase/excisionase (like ICEs and IMEs).

The Mobilome Annotation Pipeline integrates the output of different tools designed for the prediction of plasmids, phages, insertion sequences, integrative mobile genetic elements (ICEs, IMEs), integrons, and non-autonomous mobile genetic elements in prokaryotic genomes and metagenomes. The primary output is a compressed GFF3 file of the mobilome annotation. Since v5, MAP also runs functional annotation subworkflows for antimicrobial resistance (AMR), virulence factors, and biosynthetic gene clusters (BGCs), producing a combined per-protein report.

## Contents

- [Workflow](#wf)
- [Install and dependencies](#install)
- [Usage](#usage)
- [Functional annotation](#functional)
- [Outputs](#out)
- [Tests](#test)
- [License and Attribution](#license)
- [Citation](#cite)

<a name="wf"></a>

## Workflow

<p align="center" width="100%">
   <img src="media/map_pathofact2.png" width="90%"/>
</p>

The pipeline has five main stages:

**1. Preprocessing** — Contigs are filtered by length and renamed to short IDs (`RENAME`). CDS are annotated with Prodigal; tRNAs with ARAGORN.

**2. MGE prediction** (parallel) — geNomad (plasmids/phages), ICEfinder2-lite (ICEs/IMEs), IntegronFinder (integrons), ISEScan (insertion sequences), and a compositional outlier detection subworkflow (contigs ≥100 kb).

**3. Integration** — All predictions are merged into a single `{sample}_mobilome.gff.gz`. Predictions <500 bp or with no CDS are discarded.

**4. Functional annotation** — Three independent subworkflows annotate virulence factors, antimicrobial resistance genes (ARG) and biosynthetic gene clusters (BGC):
- **PATHOFACT2**: toxin and virulence factor prediction via machine learning and DIAMOND vs VFDB.
- **AMR_ANNOTATION**: antimicrobial resistance gene detection with AMRFinderPlus, DeepARG, and RGI (CARD).
- **BGC_ANNOTATION**: biosynthetic gene cluster prediction and overlps merging with SanntiS, GECCO, and antiSMASH.

**5. Pathofact2-style report** — Virulence factor and ARG report in the context of BGCs and MGEs.


<a name="install"></a>

## Install and downloading dependencies

The only prerequisites are [Nextflow](https://www.nextflow.io/) and a container tool such as [Docker](https://www.docker.com/) or [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html).

If this the first time running nextflow please refer to [this page](https://www.nextflow.io/index.html#GetStarted)

The first time you run the pipeline you will need to set up the required databases. MAP includes a built-in download subworkflow that handles this automatically.

### Downloading all databases

Run the pipeline with `--download_dbs` pointing to a target directory:

```bash
nextflow run EBI-Metagenomics/mobilome-annotation-pipeline \
    --download_dbs /path/to/dbs \
    -profile singularity
```

This downloads and sets up all databases in parallel:

| Database | Tool | Purpose |
|---|---|---|
| geNomad v1.9 | geNomad | Plasmid/phage prediction |
| ICEfinder2-lite | ICEfinder2 | ICE/IME prediction |
| PATHOFACT2 models | PATHOFACT2 | Toxin/virulence ML models |
| VFDB (`VFDB_setB_pro.dmnd`) | DIAMOND + PATHOFACT2 | Virulence factor search |
| CDD | local-cd-search | Domain annotation (used when no IPS provided) |
| AMRFinderPlus DB | AMRFinderPlus | AMR gene detection |
| DeepARG DB | DeepARG | AMR gene detection |
| CARD | RGI | AMR gene detection |
| antiSMASH DB | antiSMASH | BGC prediction |

> **InterProScan** (~100 GB) is not included due to its size. Download it manually following the [InterProScan documentation](https://interproscan-docs.readthedocs.io/en/v5/HowToDownload.html) if you need it. It is only required for SanntiS BGC prediction; the pipeline can run without it.

On completion, the pipeline prints a ready-to-paste config block with the exact paths for your system. Save it to a file (e.g. `my_paths.config`) and pass it on every run with `-c my_paths.config`:

```nextflow
params {
    // Mobilome
    genomad_db                   = "/path/to/dbs/genomad_db_v1.9"
    icefinder_macsyfinder_models = "/path/to/dbs/icf2_dbs/macsydata"
    icefinder_hmm_models         = "/path/to/dbs/icf2_dbs/icehmm/icescan"
    icefinder_prokka_uniprot_db  = "/path/to/dbs/icf2_dbs/icefinder_prokka_uniprot"

    // PATHOFACT2
    pathofact_models             = "/path/to/dbs/Models.tar.gz"
    virulencefactors_db          = "/path/to/dbs/VFDB_setB_pro.dmnd"
    ncbi_cdd                     = "/path/to/dbs/database"

    // AMR
    amrfinderplus_db             = "/path/to/dbs/amrfinderdb"
    deeparg_db                   = "/path/to/dbs/db"
    rgi_db                       = "/path/to/dbs/card_dir"

    // BGC
    antismash_db                 = "/path/to/dbs/antismash_db"
    // SanntiS requires InterProScan output (provided via samplesheet or run internally)
}
```

<a name="usage"></a>

## Usage

Prepare a samplesheet with your input data:

```csv
sample,assembly,proteins_gff,proteins_faa,virify_gff,interproscan_tsv
sample1,/PATH/assembly.fasta,,,,
sample2,/PATH/assembly.fasta,/PATH/proteins.gff,/PATH/proteins.faa,,
sample3,/PATH/assembly.fasta,/PATH/proteins.gff,/PATH/proteins.faa,,/PATH/ips.tsv
```

Only `sample` and `assembly` are required. Optional columns:

| Column | Description |
|---|---|
| `proteins_gff` | Pre-computed CDS annotation GFF (Prodigal or equivalent). If absent, MAP runs Prodigal. |
| `proteins_faa` | Protein FASTA matching `proteins_gff`. |
| `virify_gff` | VIRify ≥3.0.0 output GFF. Prophage predictions are incorporated into the mobilome. |
| `interproscan_tsv` | InterProScan TSV. Used by PATHOFACT2 (SignalP entries) instead of local CDsearch. |

Basic run:

```bash
nextflow run ebi-metagenomics/mobilome-annotation-pipeline \
    --input samplesheet.csv \
    -c my_paths.config \
    -profile singularity
```

<a name="functional"></a>

## Functional annotation

### Running mobilome-only (v4 behaviour)

The functional annotation subworkflows are enabled by default but can be disabled individually. To run MAP in mobilome-only mode (equivalent to v4 behaviour), skip all functional annotation:

```bash
nextflow run /PATH/mobilome-annotation-pipeline/main.nf \
    --input samplesheet.csv \
    --skip_virulence true \
    --skip_amrfinderplus true \
    --skip_deeparg true \
    --skip_rgi true \
    --skip_sanntis true \
    --skip_gecco true \
    --skip_antismash true
```

### Functional annotation flags

| Parameter | Default | Description |
|---|---|---|
| `--skip_virulence` | `false` | Skip PATHOFACT2 toxin/virulence annotation and VFDB search |
| `--skip_amrfinderplus` | `false` | Skip AMRFinderPlus |
| `--skip_deeparg` | `false` | Skip DeepARG |
| `--skip_rgi` | `false` | Skip RGI (CARD) |
| `--skip_sanntis` | `false` | Skip SanntiS BGC prediction |
| `--skip_gecco` | `false` | Skip GECCO BGC prediction |
| `--skip_antismash` | `false` | Skip antiSMASH BGC prediction |

For reusing pre-computed annotation outputs from the MGnify Genomes Catalogues pipeline, see [docs/annotation_manifest.md](docs/annotation_manifest.md).

<a name="out"></a>

## Outputs

Results are written to `--outdir` (default: `results/`).

### Mobilome outputs

```
sample/
├── sample_combined_report.tsv
├── sample_discarded_mge.txt
├── sample_mobilome.fasta.gz
├── sample_overlap_report.txt
├── gff/
│   ├── sample_mobilome.gff.gz
│   ├── sample_user_mobilome_clean.gff.gz      # mobilome + matching CDSs
│   ├── sample_user_mobilome_clean.gff.gz.csi
│   ├── sample_user_mobilome_clean.gff.gz.gzi
│   ├── sample_user_mobilome_extra.gff.gz      # mobilome + VIRify ViPhOG-annotated genes
│   ├── sample_user_mobilome_extra.gff.gz.csi
│   ├── sample_user_mobilome_extra.gff.gz.gzi
│   ├── sample_user_mobilome_full.gff.gz       # mobilome + all features from user GFF
│   ├── sample_user_mobilome_full.gff.gz.csi
│   └── sample_user_mobilome_full.gff.gz.gzi
├── prediction/
│   ├── amr_genes/
│   ├── bgcs/
│   │   ├── sample_bgcs.gff
│   │   └── sample_bgcs.json
│   ├── compositional_outliers/
│   │   └── sample_100kb_contigs.1.bed
│   ├── genomad/
│   │   ├── sample_5kb_contigs_plasmid_summary.tsv
│   │   └── sample_5kb_contigs_virus_summary.tsv
│   ├── integronfinder/
│   │   ├── sample_100kb_contigs.summary
│   │   └── contig_1.gbk
│   ├── isescan/
│   │   └── sample_1kb_contigs.fasta.tsv
│   └── virulence/
│       └── sample_pathofact2.gff
└── preprocessing/
    ├── sample_1kb_contigs.fasta
    ├── sample_5kb_contigs.fasta
    ├── sample_100kb_contigs.fasta
    └── sample_contigID.map
```

The `sample_discarded_mge.txt` file lists predictions excluded during integration and the reason:

1. `mge < 500bp` — discarded by length
2. `no_cds` — no coding sequences within the prediction
3. `tRNAs_in_window` — tRNA genes present inside a compositional outlier
4. `CO_overlap_with_MGE` — compositional outlier overlapping another MGE

The `_overlap_report.txt` lists long-MGEs with overlapping coordinates (no predictions are discarded for this reason).

Mobilome feature IDs follow this pattern: `contig_id|mge_type-start:end`.

The GFF feature types and their Sequence Ontology mappings:

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

### Combined report

When functional annotation is enabled, `sample_combined_report.tsv` is produced. Each row is one protein. Only proteins with a virulence/toxin annotation from PATHOFACT2 or an AMR annotation are included. BGC context (`bgc_type`, `bgc_tools`) and mobilome context (`mge_type`) are added as additional columns for those proteins when relevant, but BGC-only or MGE-only proteins are not included as rows.

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
| `signalP` | SignalP annotation from InterProScan if an IPS TSV was provided; `-` otherwise |

<a name="test"></a>

## Tests

Nextflow integration tests use [nf-test](https://github.com/askimed/nf-test). The full suite takes ~10 minutes.

```bash
cd mobilome-annotation-pipeline/
nf-test test --profile test,singularity
```

Individual test tags: `positive`, `negative`, `pathofact`.

Python unit tests for `bin/` scripts use [pytest](https://pytest.org/):

```bash
task setup-venv   # once
task test
```

## Development Tasks

This project uses [Task](https://taskfile.dev/) and [uv](https://docs.astral.sh/uv/) for development workflow management.

```bash
task setup-venv              # Bootstrap Python virtual environment
task test                    # Run all pytest tests
task test-verbose            # Run pytest with verbose output
task test-coverage           # Run pytest with HTML coverage report
task test-specific -- <pat>  # Run a specific test file or pattern
task clean                   # Remove virtual environment
```

<a name="license"></a>

## License and Attribution

### ICEfinder2 Attribution

This pipeline includes scripts derived from or inspired by ICEfinder2 algorithms and methods. The following scripts contain code adapted from ICEfinder2:

- `bin/ice_boundary_refinement.py` — ICE boundary refinement and direct repeat processing
- `bin/map_tools/icefinder_process.py` — ICE result processing and data formatting
- `bin/prescan_to_fasta.py` — ICE prescanning and candidate detection methods

**ICEfinder2 License**: These algorithms are used under the Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International License (CC BY-NC-SA 4.0).

- Original work: ICEfinder2 (http://creativecommons.org/licenses/by-nc-sa/4.0/)
- Modifications: Licensed under Apache 2.0 by EMBL-EBI

<a name="cite"></a>

## Citation

The Mobilome Annotation Pipeline integrates the following tools and databases (alphabetical):

**MGE prediction**
- geNomad v1.11.1 — [Camargo et al., Nature Biotechnology, 2023](https://doi.org/10.1038/s41587-023-01953-y)
- ICEfinder v2.0 — [Wang et al., Nucleic Acids Res, 2024](https://academic.oup.com/nar/article/52/D1/D732/7327075)
- IntegronFinder2 v2.0.6 — [Néron et al., Microorganisms, 2022](https://doi.org/10.3390/microorganisms10040700)
- ISEScan v1.7.3 — [Xie et al., Bioinformatics, 2017](https://doi.org/10.1093/bioinformatics/btx433)
- Prodigal v2.6.3 — [Hyatt et al., Bioinformatics, 2010](https://doi.org/10.1186/1471-2105-11-119)
- VIRify v3.0.0 — [Rangel-Pineros et al., PLoS Comput Biol, 2023](https://doi.org/10.1371/journal.pcbi.1011422)

**Functional annotation**
- AMRFinderPlus — [Feldgarden et al., Scientific Reports, 2021](https://doi.org/10.1038/s41598-021-91456-0)
- antiSMASH — [Blin et al., Nucleic Acids Res, 2023](https://doi.org/10.1093/nar/gkad344)
- DeepARG — [Arango-Argoty et al., Microbiome, 2018](https://doi.org/10.1186/s40168-018-0401-z)
- GECCO — [Carroll et al., bioRxiv, 2021](https://doi.org/10.1101/2021.05.03.442509)
- PATHOFACT2 — [Delgado et al., GigaScience, 2026](https://doi.org/10.1093/gigascience/giag062)
- RGI / CARD — [Alcock et al., Nucleic Acids Res, 2023](https://doi.org/10.1093/nar/gkac920)
- SanntiS — [Sanchez et al., bioRxiv, 2023](https://doi.org/10.1101/2023.08.14.552309)
- VFDB — [Liu et al., Nucleic Acids Res, 2022](https://doi.org/10.1093/nar/gkab1107)
