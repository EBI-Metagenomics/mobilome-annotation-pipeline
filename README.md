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

The pipeline has four main stages:

**1. Preprocessing** — Contigs are filtered by length and renamed to short IDs (`RENAME`). CDS are annotated with Prodigal; tRNAs with ARAGORN.

**2. MGE prediction** (parallel) — geNomad (plasmids/phages), ICEfinder2-lite (ICEs/IMEs), IntegronFinder (integrons), ISEScan (insertion sequences), and a compositional outlier detection subworkflow (contigs ≥100 kb).

**3. Integration** — All predictions are merged into a single `{sample}_mobilome.gff.gz`. Predictions <500 bp or with no CDS are discarded.

**4. Functional annotation** (runs after mobilome prediction) — Three independent subworkflows annotate virulence factors, ARG and BGC:
- **PATHOFACT2**: toxin and virulence factor prediction via machine learning and DIAMOND vs VFDB.
- **AMR_ANNOTATION**: antimicrobial resistance gene detection with AMRFinderPlus, DeepARG, and RGI (CARD).
- **BGC_ANNOTATION**: biosynthetic gene cluster prediction and overlps merging with SanntiS, GECCO, and antiSMASH.

**5. Pathofact2-style report** — Virulence factor and ARG report in the context of BGCs and MGEs.


<a name="install"></a>

## Install and downloading dependencies

The only prerequisites are [Nextflow](https://www.nextflow.io/) and a container tool such as [Docker](https://www.docker.com/) or [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html).

If this the first time running nextflow please refer to [this page](https://www.nextflow.io/index.html#GetStarted)

The first time you run the pipeline you will need to set up the following databases:

### Databases to run mobilome prediction
1. Download and extract the geNomad database:
```bash
wget https://zenodo.org/records/14886553/files/genomad_db_v1.9.tar.gz
tar -xvf genomad_db_v1.9.tar.gz
```

2. Download and extract the ICEfinder2-lite databases:
```bash
wget https://ftp.ebi.ac.uk/pub/databases/metagenomics/pipelines/tool-dbs/icefinder2lite/icf2_dbs.tar.gz
tar -xvf icf2_dbs.tar.gz
```

### Required databases for Pathofact2-style functional annotation
Some databases set up require using the tool iteself. You can do it using docker
```bash
docker run --rm -v $(pwd):/data -w /data quay.io/biocontainers/TOOL:VERSION
```
Or singularity:
```bash
singularity run https://depot.galaxyproject.org/singularity/TOOL:VERSION
```

1. Download and extract the pathofact2 models
```bash
wget https://zenodo.org/records/18223764/files/Models.tar.gz?download=1 -O Models.tar.gz
tar -xvf Models.tar.gz
```
2. Download and extract VFDB, and format using diamond with docker or singularity
```bash
wget https://www.mgc.ac.cn/VFs/Down/VFDB_setB_pro.fas.gz
gzip -d VFDB_setB_pro.fas.gz
diamond:2.1.16--h13889ed_0 diamond makedb --in VFDB_setB_pro.fas -d VFDB_setB_pro
```
3. Download CDD database using local-cd-search tool with docker or singularity. Optional; used when no IPS is provided
```bash
mkdir cdd_database
local-cd-search:0.3.0--pyhdfd78af_0 local-cd-search download cdd_database/ cdd
```

4. Download and decompress AMRfinderPlus database with docker or singularity
```bash
ncbi-amrfinderplus:4.2.7--hf69ffd2_0 amrfinder_update -d amrfinderdb
tar czvf amrfinderdb.tar.gz -C amrfinderdb/\$(readlink amrfinderdb/latest) ./
```

5. Download deeparg database with docker or singularity
```bash
# if docker add: -v $(which bash):/usr/local/lib/python2.7/site-packages/Theano-0.8.2-py2.7.egg-info/PKG-INFO
# If singularity add: -B $(which bash):/usr/local/lib/python2.7/site-packages/Theano-0.8.2-py2.7.egg-info/PKG-INFO
mkdir -p theano
export THEANO_FLAGS="base_compiledir=\$PWD/theano"
deeparg:1.0.4--pyhdfd78af_0 deeparg download_data -o db/
```

6. Download CARD database for RGI
```bash
mkdir CARD_db
wget https://card.mcmaster.ca/latest/data
tar -xvf data ./card.json
mv card.json CARD_db
```

7. Download antismash database with docker or singularity
```bash
antismash:8.0.1--pyhdfd78af_0 download-antismash-databases --database-dir antismash_db
```

8. Download [InterProScan](https://interproscan-docs.readthedocs.io/en/v5/HowToDownload.html) databases
```bash
mkdir interproscan && cd interproscan
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.76-107.0/interproscan-5.76-107.0-64-bit.tar.gz
wget https://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/5.76-107.0/interproscan-5.76-107.0-64-bit.tar.gz.md5
md5sum -c interproscan-5.76-107.0-64-bit.tar.gz.md5
```

Once downloaded, we recommend creating a config file with all paths and passing it with `-c my_paths.config`:

```nextflow
params {
    // Mobilome
    genomad_db                   = "/PATH/genomad_db_v1.9"
    icefinder_macsyfinder_models = "/PATH/icf2_dbs/macsydata/"
    icefinder_hmm_models         = "/PATH/icf2_dbs/icehmm/icescan.hmm"
    icefinder_prokka_uniprot_db  = "/PATH/icf2_dbs/icefinder_prokka_uniprot/"

    // PATHOFACT2
    pathofact_models   = "/PATH/pathofact2_models"
    virulecefactors_db = "/PATH/VFDB_setB_pro.dmnd"
    ncbi_cdd           = "/PATH/cdd_database"   // optional; used when no IPS is provided

    // AMR
    amrfinderplus_db   = "/PATH/amrfinderplus_db"
    deeparg_db         = "/PATH/deeparg_db"
    rgi_db             = "/PATH/CARD_db"

    // BGC
    antismash_db       = "/PATH/antismash_db"
    // SanntiS require InterProScan output (provided via samplesheet or run internally)
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
| `proteins_gff` | Pre-computed CDS annotation GFF (Prodigal or equivalent). If absent, MAP runs Prodigal internally. |
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
    -c my_paths.config \
    -profile singularity \
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

For reusing pre-computed annotation outputs from the EBI Genomes Catalogues pipeline, see [docs/annotation_manifest.md](docs/annotation_manifest.md).

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
