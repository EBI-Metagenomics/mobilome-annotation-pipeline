# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Pipeline Does

The **Mobilome Annotation Pipeline (MAP)** is a Nextflow DSL2 pipeline that integrates multiple tools to predict and annotate mobile genetic elements (MGEs) — plasmids, phages, insertion sequences, integrons, ICEs/IMEs — in prokaryotic genomes and metagenomes. The primary output is a compressed GFF3 file (`mobilome.gff.gz`) and a FASTA of mobilome sequences.

The pipeline also runs functional annotation subworkflows: **PATHOFACT2** (toxins/virulence), **AMR_ANNOTATION** (AMRFinderPlus, DeepARG, RGI), and **BGC_ANNOTATION** (SanntiS, GECCO, antiSMASH).

## Commands

### Python unit tests (bin/ scripts)

```bash
task setup-venv              # Create .venv with uv (run once)
task test                    # Run all pytest tests
task test-verbose            # Run pytest -v
task test-coverage           # Run pytest with HTML coverage report
task test-specific -- tests/unit/test_ice_boundary_refinement.py  # Run one file
```

Python tests live in `tests/unit/` and target scripts in `bin/`. The `bin/` directory is on `PYTHONPATH` (set in `pyproject.toml`).

### Nextflow integration tests

```bash
nf-test test --profile test,singularity          # Full nf-test suite (~9 min)
nf-test test tests/default.nf.test               # Single test file
```

Test databases are checked in under `tests/reference_databases/`. The `NFT_WORKDIR` env var controls where nf-test writes temp files.

### Running the pipeline

```bash
nextflow run main.nf --input samplesheet.csv -c my_paths.config -profile singularity
```

Required database params (pass via `-c my_paths.config` or flags):
- `--genomad_db` — geNomad DB directory
- `--icefinder_hmm_models` — ICEfinder HMM file (path without extension, used as glob `${path}.*`)
- `--icefinder_macsyfinder_models` — MaCSyFinder models directory
- `--icefinder_prokka_uniprot_db` — Prokka UniProt BLAST DB directory
- `--pathofact_models`, `--virulecefactors_db`, `--ncbi_cdd` — PATHOFACT2 databases
- `--amrfinderplus_db`, `--deeparg_db`, `--rgi_db` — AMR databases
- `--antismash_db`, `--ips_db` — BGC databases

Samplesheet columns: `sample,assembly,user_proteins_gff,user_proteins_fasta,virify_gff,ips_tsv` (only `sample` and `assembly` are mandatory).

## Architecture

### Pipeline stages (in `workflows/mobilomeannotation.nf`)

1. **Preprocessing** — `RENAME` filters contigs by length (1kb/5kb/100kb cutoffs) and renames them to short IDs, writing a mapping file. `PRODIGAL` annotates CDS; `ARAGORN` finds tRNAs; `TRNAS_INTEGRATOR` merges them into a single GFF/FAA.

2. **MGE Prediction** (all run in parallel):
   - `GENOMAD` — plasmid/phage/virus prediction on 5kb+ contigs
   - `ICEFINDER2_LITE` subworkflow — ICE/IME prediction via HMMSCAN prescanning → MACSYFINDER + BLASTP + VMATCH → `REFINE_BOUNDARIES`
   - `INTEGRONFINDER` — integron detection on 100kb+ contigs
   - `ISESCAN` — insertion sequence detection on 1kb+ contigs
   - `COMPOSITIONAL_OUTLIER_DETECTION` subworkflow — FASTA split into 20-sequence chunks, outlier scoring, BED merge; runs on 100kb+ contigs only

3. **Integration** — `INTEGRATOR` runs `mge_integrator.py`, consolidating all predictions into `{prefix}_mobilome.gff.gz`. Predictions <500bp or with no CDS are discarded.

4. **Postprocessing**:
   - `FASTA_WRITER` — extracts mobilome sequences
   - `GT_GFF3VALIDATOR` — validates GFF3 (skippable via `--gff_validation false`)
   - `GFF_MAPPING_COMPRESSION_AND_INDEXING` subworkflow — when user proteins are provided, produces `mobilome_clean/extra/full.gff.gz` with `.gzi` and `.csi` indexes

5. **Functional annotation** (runs after mobilome prediction):
   - `PATHOFACT2` subworkflow — toxin/virulence annotation via HMMER models + DIAMOND vs VFDB + local CDsearch (or IPS if provided)
   - `AMR_ANNOTATION` subworkflow — AMRFinderPlus + DeepARG + RGI; each tool skippable
   - `BGC_ANNOTATION` subworkflow — SanntiS (needs IPS or runs InterProScan) + GECCO + antiSMASH; each tool skippable
   - `COMBINEREPORTER` — merges mobilome GFF with functional annotation GFFs into final report

### Key Python scripts (`bin/`)

| Script | Purpose |
|--------|---------|
| `mge_integrator.py` | Core integration logic; merges all tool outputs into a single mobilome GFF |
| `assembly_filter_rename.py` | Filters contigs by length, renames to short IDs |
| `ice_boundary_refinement.py` | Boundary refinement for ICEs using direct repeats (ICEfinder2-derived) |
| `prescan_to_fasta.py` | Filters HMM hits to extract ICE candidate FASTAs |
| `gff_mapping.py` | Maps user CDSs onto mobilome GFF to produce clean/extra/full GFFs |
| `trnas_integrator.py` | Merges Prodigal GFF + ARAGORN tRNA table |
| `pathofact2_report.py` | Generates PATHOFACT2-style combined report |
| `bin/map_tools/` | Shared Python library: parsers for each tool (geNomad, ICEfinder, ISEScan, IntegronFinder, VIRify) plus overlap logic |

### Module organisation

- `modules/local/` — pipeline-specific Nextflow modules
- `modules/nf-core/` — standard nf-core modules (prodigal, BLAST, tabix, etc.)
- `modules/ebi-metagenomics/` — EBI-specific modules (PATHOFACT2 steps, AMR integrator, BGC mapper)
- `subworkflows/local/` — local multi-step workflows
- `subworkflows/ebi-metagenomics/` — EBI subworkflows (`pathofact2`, `amr_annotation`, `bgc_annotation`)

### Optional inputs and Nextflow channel handling

Nextflow does not support optional `path` inputs natively. The pipeline uses the community pattern of passing `[]` (empty list) for missing optional files, combined with `.filter{}` and `.join(..., remainder: true)` to route only samples that have a given input into the relevant processes. See the `INTEGRATOR` input construction in `workflows/mobilomeannotation.nf` for the canonical example.

### Contig length filtering

Contigs are filtered into three sets during `RENAME`:
- `1kb` — used by ISEScan, Prodigal, ARAGORN
- `5kb` — used by geNomad, ICEfinder2, BGC annotation
- `100kb` — used by IntegronFinder, compositional outlier detection

### Resource configuration

`conf/base.config` sets per-label resource defaults. `conf/modules.config` overrides per-process publish paths and tool arguments. The `local` profile limits resources to 8 CPUs / 12 GB.
