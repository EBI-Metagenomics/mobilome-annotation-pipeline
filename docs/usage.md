# MAP Usage Guide

## Contents

- [Install and dependencies](#install)
- [Downloading databases](#databases)
- [Samplesheet format](#samplesheet)
- [Running the pipeline](#running)
- [Functional annotation options](#functional)
- [Annotation manifest (reuse mode)](#manifest)

<a name="install"></a>

## Install and dependencies

The only prerequisites are [Nextflow](https://www.nextflow.io/) >=24.04.0 and a container tool such as [Docker](https://www.docker.com/) or [Singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html).

If this is the first time running Nextflow, refer to [this page](https://www.nextflow.io/index.html#GetStarted).

<a name="databases"></a>

## Downloading databases

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

<a name="samplesheet"></a>

## Samplesheet format

Prepare a CSV with your input data:

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
| `interproscan_tsv` | Pre-computed InterProScan TSV. See below. |

### InterProScan input

MAP uses InterProScan (IPS) for two purposes:

1. **SanntiS BGC prediction** — SanntiS requires IPS output to identify BGCs. If no IPS is available, MAP runs IPS internally.
2. **SignalP annotation in the combined report** — SignalP entries from the IPS TSV are passed directly into the PATHOFACT2 combined report as the `signalP` column. They are not used for MGE prediction.

When MAP runs IPS internally, it uses the following applications:

| Application | Purpose |
|---|---|
| CDD | Conserved domain detection |
| TIGRFAM | TIGR protein family HMMs |
| GENE3D | CATH structural annotation |
| PRINTS | Protein fingerprint motifs |
| PROSITEPATTERNS | PROSITE pattern matches |
| PFAM | Pfam domain annotation |

[SignalP](https://services.healthtech.dtu.dk/services/SignalP-6.0/) is added automatically when `--interpro_licensed_software true` is set (SignalP is a licensed database not distributed with IPS by default).

If you provide a pre-computed `interproscan_tsv` in the samplesheet, MAP skips the internal IPS run entirely and uses your file for both SanntiS and the SignalP column.

<a name="running"></a>

## Running the pipeline

```bash
nextflow run EBI-Metagenomics/mobilome-annotation-pipeline \
    --input samplesheet.csv \
    -c my_paths.config \
    -profile singularity
```

<a name="functional"></a>

## Functional annotation options

### Mobilome-only mode (v4 behaviour)

The functional annotation subworkflows are enabled by default but can be disabled individually. To run MAP in mobilome-only mode, skip all functional annotation:

```bash
nextflow run EBI-Metagenomics/mobilome-annotation-pipeline \
    --input samplesheet.csv \
    --skip_virulence \
    --skip_amrfinderplus \
    --skip_deeparg \
    --skip_rgi \
    --skip_sanntis \
    --skip_gecco \
    --skip_antismash 
```

### Skip flags

| Parameter | Default | Description |
|---|---|---|
| `--skip_virulence` | `false` | Skip PATHOFACT2 toxin/virulence annotation and VFDB search |
| `--skip_amrfinderplus` | `false` | Skip AMRFinderPlus |
| `--skip_deeparg` | `false` | Skip DeepARG |
| `--skip_rgi` | `false` | Skip RGI (CARD) |
| `--skip_sanntis` | `false` | Skip SanntiS BGC prediction |
| `--skip_gecco` | `false` | Skip GECCO BGC prediction |
| `--skip_antismash` | `false` | Skip antiSMASH BGC prediction |

<a name="manifest"></a>

## Annotation manifest (reuse mode)

When running MAP on assemblies that were previously processed by the following MGnify pipelines: [Genomes Catalogues pipeline](https://github.com/EBI-Metagenomics/genomes-catalogue-pipeline), [Assembly Analysis Pipeline](https://github.com/EBI-Metagenomics/assembly-analysis-pipeline), and [Mettanotator](https://github.com/EBI-Metagenomics/mettannotator); pre-computed annotation outputs (IPS, AMRFinderPlus, antiSMASH, GECCO, SanntiS) can be reused directly to avoid redundant computation.

See [annotation_manifest.md](annotation_manifest.md) for the manifest format and column reference.
