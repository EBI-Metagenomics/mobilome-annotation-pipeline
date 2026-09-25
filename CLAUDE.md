# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## What This Pipeline Does

The **Mobilome Annotation Pipeline (MAP)** is a Nextflow DSL2 pipeline (nf-core template, `org: ebi-metagenomics`, not an nf-core pipeline) that integrates multiple tools to predict and annotate mobile genetic elements (MGEs) — plasmids, phages, insertion sequences, integrons, ICEs/IMEs — in prokaryotic genomes and metagenomes. The primary output is a compressed GFF3 file (`{sample}_mobilome.gff.gz`) plus a FASTA of mobilome sequences and a combined per-protein report.

Since v5 the pipeline also runs functional annotation subworkflows: **PATHOFACT2** (toxins/virulence), **AMR_ANNOTATION** (AMRFinderPlus, DeepARG, RGI), and **BGC_ANNOTATION** (SanntiS, GECCO, antiSMASH).

## Commands

### Python unit tests (bin/ scripts)

```bash
task setup-venv              # Create .venv with uv (run once)
task test                    # Run all pytest tests
task test-verbose            # Run pytest -v
task test-coverage           # Run pytest with HTML coverage report
task test-specific -- tests/unit/test_ice_boundary_refinement.py  # Run one file
task clean                   # Remove the venv
```

Python tests live in `tests/unit/` and target scripts in `bin/`. The `bin/` directory is on `PYTHONPATH` (set in `pyproject.toml`). CI (`.github/workflows/pytest.yml`) just runs `uv run pytest`.

### Nextflow integration tests

```bash
task fetch-icefinder-test-db                     # ONCE: download ICEfinder2-lite test DBs (61 MB, EBI FTP)
nf-test test                                     # Full suite (~10 min), defaults to profile test,docker
nf-test test --profile test,singularity          # Same with Singularity
nf-test test tests/default.nf.test               # Single test file
nf-test test --tag positive                      # Tags: positive, negative, schema, pathofact
nf-test test --update-snapshot                   # Regenerate tests/default.nf.test.snap
```

`nf-test.config` sets `testsDir "."`, loads `tests/nextflow.config`, and ignores `modules/**/tests/*`. `NFT_WORKDIR` controls where nf-test writes temp files. Most test databases are checked in under `tests/reference_databases/`, **except** the ICEfinder2-lite DBs, which must be fetched with the task above (CI does this explicitly).

Gotcha: a `withName` block in `tests/nextflow.config` or `conf/test.config` replaces the `publishDir` selector from `conf/modules.config` rather than merging with it — use `withLabel` for test-only resource bumps.

### Linting / formatting

Prettier (`.prettierrc.yml`, `.prettierignore`) and `.editorconf` govern formatting. `nf-core lint` is configured via `.nf-core.yml`, which lists the template files intentionally absent from this repo — update that `lint.files_exist` list rather than adding placeholder files.

### Running the pipeline

```bash
# One-off: fetch every database (except InterProScan) and print a ready-to-paste params block
nextflow run main.nf --download_dbs /path/to/dbs -profile singularity

# Normal run
nextflow run main.nf --input samplesheet.csv -c my_paths.config -profile singularity
```

`main.nf` branches on `params.download_dbs`: set → `DOWNLOAD_DATABASES` subworkflow only; unset → the `MOBILOMEANNOTATION` workflow.

Required database params (pass via `-c my_paths.config` or flags):
- `--genomad_db`, `--checkv_db` — geNomad and CheckV DB directories
- `--icefinder_hmm_models` — ICEfinder HMM file (path without extension, used as glob `${path}.*`)
- `--icefinder_macsyfinder_models` — MaCSyFinder models directory
- `--icefinder_prokka_uniprot_db` — Prokka UniProt BLAST DB directory
- `--pathofact_models`, `--virulencefactors_db`, `--ncbi_cdd` — PATHOFACT2 databases
- `--amrfinderplus_db`, `--deeparg_db`, `--rgi_db` — AMR databases
- `--antismash_db`, `--interproscan_db` — BGC databases (InterProScan must be installed manually, ~100 GB; only needed for SanntiS)

Other notable params: `--gff_validation` (default true), `--publish_all` (default true; gates most `publishDir` blocks in `conf/modules.config`), `--annotation_manifest`, and the skip flags `--skip_virulence`, `--skip_amrfinderplus`, `--skip_deeparg`, `--skip_rgi`, `--skip_sanntis`, `--skip_gecco`, `--skip_antismash`.

### Inputs

Samplesheet (`--input`, validated by `assets/schema_input.json`):
`sample,assembly,proteins_gff,proteins_faa,virify_gff,interproscan_tsv` — only `sample` and `assembly` are mandatory; `proteins_gff` and `proteins_faa` are mutually `dependentRequired`.

Annotation manifest (`--annotation_manifest`, validated by `assets/schema_manifest.json`): optional CSV of pre-computed annotations from the genomes-catalogue-pipeline — `sample,amrfinder_tsv,antismash_gff,gecco_gff,sanntis_gff`. `PARSE_ANNOTATION_MANIFEST` errors out on sample IDs not present in the samplesheet, and each populated column makes MAP reuse that file instead of running the tool. InterProScan input lives only in the samplesheet (`interproscan_tsv`), never in the manifest.

## Architecture

### Pipeline stages (in `workflows/mobilomeannotation.nf`)

1. **Preprocessing** — `RENAME` filters contigs by length (1kb/5kb/100kb cutoffs) and renames them to short IDs, writing a mapping file. `PRODIGAL` annotates CDS; `ARAGORN` finds tRNAs; `TRNAS_INTEGRATOR` merges them into a single GFF/FAA. User-supplied proteins GFF/FAA replace the Prodigal output when given.

2. **MGE Prediction** (all run in parallel):
   - `GENOMAD` — plasmid/phage/virus prediction on 5kb+ contigs; `CHECKV_ENDTOEND` then scores the geNomad virus FASTA and its quality summary feeds `INTEGRATOR`
   - `ICEFINDER2_LITE` subworkflow — ICE/IME prediction via HMMSCAN prescanning → MACSYFINDER + BLASTP + VMATCH → `REFINE_BOUNDARIES`
   - `INTEGRONFINDER` — integron detection on 100kb+ contigs
   - `ISESCAN` — insertion sequence detection on 1kb+ contigs
   - `COMPOSITIONAL_OUTLIER_DETECTION` subworkflow — FASTA split into 20-sequence chunks, outlier scoring, BED merge; runs on 100kb+ contigs only
   - `VIRIFY_QC` — ingests an optional user-provided VIRify GFF

3. **Integration** — `INTEGRATOR` runs `mge_integrator.py`, consolidating all predictions into `{prefix}_mobilome.gff.gz`. Predictions <500bp or with no CDS are discarded.

4. **Postprocessing**:
   - `FASTA_WRITER` — extracts mobilome sequences
   - `GT_GFF3VALIDATOR` — validates GFF3 (skippable via `--gff_validation false`)
   - `GFF_MAPPING_COMPRESSION_AND_INDEXING` subworkflow — when user proteins are provided, produces `mobilome_clean/extra/full.gff.gz` with `.gzi` and `.csi` indexes

5. **Functional annotation** (runs after mobilome prediction):
   - `PATHOFACT2` subworkflow — toxin/virulence annotation via HMMER models + DIAMOND vs VFDB + local CDsearch (or IPS if provided)
   - `AMR_ANNOTATION` subworkflow — AMRFinderPlus + DeepARG + RGI; each tool skippable
   - `BGC_ANNOTATION` subworkflow — SanntiS (needs IPS or runs InterProScan) + GECCO + antiSMASH; each tool skippable
   - `COMBINEREPORTER` — merges mobilome GFF with functional annotation GFFs into `{sample}_combined_report.tsv`

### Key Python scripts (`bin/`)

| Script | Purpose |
|--------|---------|
| `mge_integrator.py` | Core integration logic; merges all tool outputs into a single mobilome GFF |
| `assembly_filter_rename.py` | Filters contigs by length, renames to short IDs |
| `ice_boundary_refinement.py` | Boundary refinement for ICEs using direct repeats (ICEfinder2-derived) |
| `prescan_to_fasta.py` | Filters HMM hits to extract ICE candidate FASTAs |
| `fast_composition_analyzer.py` | Compositional outlier scoring for the CO subworkflow |
| `gff_mapping.py` | Maps user CDSs onto mobilome GFF to produce clean/extra/full GFFs |
| `trnas_integrator.py` | Merges Prodigal GFF + ARAGORN tRNA table |
| `pathofact2_report.py` | Generates the PATHOFACT2-style combined report |
| `amr_report.py` / `merge_results.py` | AMR result reformatting and report merging |
| `bin/map_tools/` | Shared Python library: parsers for each tool (geNomad, ICEfinder, ISEScan, IntegronFinder, VIRify, mobileOG, outliers) plus overlap/CDS-locator logic |

Scripts under `bin/` are on the container PATH for every process. Additionally `nextflow.enable.moduleBinaries = true` (set in `nextflow.config`) puts module-scoped scripts in `modules/*/resources/usr/bin` on the PATH for their own module only — used by the PATHOFACT2 modules.

`scripts/` holds one-off analysis/figure scripts for the MAP paper, not pipeline code.

### Module organisation

- `modules/local/` — pipeline-specific Nextflow modules
- `modules/nf-core/` — standard nf-core modules (prodigal, checkv, BLAST, tabix, etc.), tracked in `modules.json`
- `modules/ebi-metagenomics/` — EBI-specific modules (PATHOFACT2 steps, antiSMASH, InterProScan, RGI, SanntiS, AMR integrator, toolkit)
- `subworkflows/local/` — local multi-step workflows
- `subworkflows/ebi-metagenomics/` — EBI subworkflows (`pathofact2`, `amr_annotation`, `bgc_annotation`)

### Optional inputs and Nextflow channel handling

Nextflow does not support optional `path` inputs natively. The pipeline uses the community pattern of passing `[]` (empty list) for missing optional files, combined with `.filter{}` and `.join(..., remainder: true)` to route only samples that have a given input into the relevant processes. See the `INTEGRATOR` input construction in `workflows/mobilomeannotation.nf` for the canonical example.

Several `.map` closures downstream of those joins take a single `row ->` parameter and index positionally instead of destructuring. That is deliberate: a chain of `remainder: true` joins yields variable-arity tuples, so destructured parameter lists would break. Don't "clean them up".

When splitting one channel into several per-tool channels, use `multiMap` (see `PARSE_ANNOTATION_MANIFEST`) — multiple `.map{}` calls on the same queue channel round-robin items between consumers instead of broadcasting.

### Contig length filtering

Contigs are filtered into three sets during `RENAME`:
- `1kb` — used by ISEScan, Prodigal, ARAGORN
- `5kb` — used by geNomad, ICEfinder2, BGC annotation
- `100kb` — used by IntegronFinder, compositional outlier detection

### Resource and container configuration

`conf/base.config` sets per-label resource defaults. `conf/modules.config` overrides per-process publish paths and tool arguments. `conf/download_dbs.config` and `conf/test.config` cover the DB-download and test profiles. The `local` profile limits resources to 8 CPUs / 12 GB. All container engines default to the `quay.io` registry (set in `nextflow.config`) because bare `biocontainers/*` names otherwise resolve to Docker Hub and 401.

## License note

`bin/ice_boundary_refinement.py`, `bin/map_tools/icefinder_process.py`, and `bin/prescan_to_fasta.py` contain code adapted from ICEfinder2 under CC BY-NC-SA 4.0; modifications are Apache 2.0. Keep the attribution intact when editing these.
