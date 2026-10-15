# ebi-metagenomics/mobilome-annotation-pipeline: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### Added

- Two columns to the combined report (`*_combined_report.tsv`): `contig_id` (the contig
  each protein is on) and `summary_string`, a condensed per-gene tag string with comma-joined
  tokens `vf,arg,mge,bgc` (virulence, AMR, mobile element, BGC).
- The `summary_string` is now also carried into the derived GFFs as a
  `pathofact2=<summary_string>` attribute on each matching CDS.
- The derived `clean`/`extra`/`full` GFFs are now produced even when no user proteins are
  provided, built from the Prodigal/tRNA genes GFF and named
  `*_mobilome_{clean,extra,full}.gff.gz` (without the `_user` infix). With user proteins the
  outputs keep their `*_user_mobilome_*` names.

### Changed

- `--download_dbs`: Zenodo downloads (geNomad DB, PathoFact2 models) now go through a new local
  `ZENODO_DOWNLOAD` module using the mgnify-pipelines-toolkit `download_zenodo` script. The
  ICEfinder2 download is split into its own `DB_DOWNLOAD_ICEFINDER2` module (formerly part of
  `DB_DOWNLOAD_MOBILOME_DBS`).
- Database sources are now params in `nextflow.config` instead of being hardcoded in modules:
  `genomad_db_version` (`1.9`), `genomad_zenodo_id` (`14886553`), `checkv_db_version` (`1.5`),
  `icefinder2_db_url`, and `pathofact_zenodo_id` (`18223764`). **Breaking:**
  `pathofact_zenodo_id` replaces `zenodo_id`.
- When `pathofact_models` is not given, MAP now downloads the PathoFact2 models with the same
  `ZENODO_DOWNLOAD` module as `--download_dbs`, instead of the vendored
  `PATHOFACT2_DOWNLOADDATA`.

### Fixed

- `--download_dbs`: VFDB is no longer fetched with `curl`, which is absent from the DIAMOND
  container. Nextflow now stages the VFDB FASTA and the nf-core `DIAMOND_MAKEDB` module builds
  the database, replacing the local `DB_DOWNLOAD_VFDB` module (#105).
- `--download_dbs`: the extracted geNomad directory is now renamed from `genomad_db/` to
  `genomad_db_v1.9/`, fixing the "Missing output file(s)" error (#106).
- `deeparg_db_version`, `deeparg_model` and `deeparg_tool_version` now default to `2`, `LS`
  and `1.0.4` instead of `null`, which crashed `AMR_ANNOTATION` (#107).
- The config block printed by `--download_dbs` and the usage docs now give
  `icefinder_hmm_models` as `.../icehmm/icescan.hmm`, matching the files in the database
  bundle (#108).
- `--download_dbs`: the ready-to-paste config block is now actually printed on completion
  (the `onComplete` handler previously failed with a `NullPointerException`).
- `--download_dbs`: the CheckV database is now published at `checkV/checkv-db-v*`, matching
  the printed `checkv_db` path, instead of `checkV/checkv_db/checkv-db-v*`.
- `*_user_mobilome_clean.gff` now starts with a `##gff-version 3` header and contains only
  genes covered by a mobile element (previously it had no header and included every gene).
- `*_user_mobilome_extra.gff` header is reduced to just `##gff-version 3`; the
  `##sequence-region` lines are kept only in the `full` output.

## v5.0.0 - [2026-06-15]

Major release adding a full functional annotation layer (toxins/virulence, AMR, BGC) on
top of mobilome prediction. This release contains **breaking changes** relative to
`v4.2.3` — see the upgrade guide below.

### ⚠️ Breaking changes

- **Input samplesheet schema reworked.** The header changed from
  `sample,assembly,user_proteins_gff,virify_gff` to
  `sample,assembly,proteins_gff,proteins_faa,virify_gff,interproscan_tsv`:
  - `user_proteins_gff` was **renamed** to `proteins_gff`.
  - A new `proteins_faa` column was added, and user proteins now require **both** the GFF
    and the FASTA — the schema enforces `proteins_gff` ⇄ `proteins_faa` as a co-required
    pair. Previously a proteins GFF could be supplied on its own.
  - Protein and VIRify GFF inputs must now use the `.gff`/`.gff.gz` extension. The `.gff3`
    extension is **no longer accepted**.
  - A new optional `interproscan_tsv` column was added for precomputed InterProScan results.
- **Functional annotation runs by default.** The pipeline now executes the PATHOFACT2,
  AMR_ANNOTATION, and BGC_ANNOTATION subworkflows in addition to mobilome prediction. The
  corresponding `--skip_*` flags all default to `false`, so a minimal `v4.2.3`-style
  invocation (geNomad + ICEfinder databases only) will attempt these new stages and **fail
  without their databases**. To keep the previous mobilome-only behavior, set the relevant
  `--skip_*` flags, or provide the new databases (see below).
- **New databases required for the default functional annotation stages.** Unless the
  matching `--skip_*` flag is set (or `--annotation_manifest` / `--download_dbs` is used),
  the following must now be provided: `--pathofact_models`, `--virulencefactors_db`,
  `--ncbi_cdd` (PATHOFACT2); `--amrfinderplus_db`, `--deeparg_db`, `--rgi_db` (AMR);
  `--antismash_db`, `--interproscan_database` (BGC).
- **Contig length filtering is now inclusive.** The 1 kb / 5 kb / 100 kb cutoffs changed
  from strict greater-than (`>`) to greater-than-or-equal (`>=`). Contigs exactly at a
  threshold are now included, which can change outputs compared to `v4.2.3`.

### Added

- Functional annotation subworkflows: **PATHOFACT2** (toxins/virulence), **AMR_ANNOTATION**
  (AMRFinderPlus, DeepARG, RGI), and **BGC_ANNOTATION** (SanntiS, GECCO, antiSMASH), each
  individually skippable via `--skip_*` flags.
- `--annotation_manifest`: reuse precomputed per-tool outputs from the MGnify
  genomes-catalogue-pipeline (AMRFinderPlus, antiSMASH, GECCO, SanntiS), bypassing the
  corresponding internal tools.
- `--download_dbs`: helper subworkflow to fetch the required databases.
- Pathofact2-style combined report consolidating the mobilome GFF with the functional
  annotation outputs.

### Upgrade guide (from v4.2.3)

1. In your samplesheet, rename the `user_proteins_gff` column to `proteins_gff` and add the
   now-required paired `proteins_faa` column. Update the header to
   `sample,assembly,proteins_gff,proteins_faa,virify_gff,interproscan_tsv`.
2. Convert any `.gff3` protein or VIRify inputs to `.gff` (`.gff.gz` is also accepted).
3. To reproduce the previous mobilome-only behavior, add the skip flags:
   `--skip_virulence --skip_amrfinderplus --skip_deeparg --skip_rgi --skip_sanntis
   --skip_gecco --skip_antismash`. Otherwise, provide the new functional annotation
   databases, or fetch them with `--download_dbs`.
