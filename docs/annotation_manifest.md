# Annotation Manifest

The `--annotation_manifest` option allows functional annotation outputs produced by the [EBI Genomes Catalogues pipeline](https://github.com/EBI-Metagenomics/genomes-catalogue-pipeline) to be reused directly, bypassing the corresponding internal tools. This avoids redundant compute when processing genomes that have already been annotated in a catalogue context.

## When to use it

Use `--annotation_manifest` when:

- You are running MAP on assemblies that are part of an EBI Genomes Catalogue run and you already have InterProScan, AMRFinderPlus, antiSMASH, GECCO, or SanntiS outputs for those assemblies.
- You want to skip re-running one or more annotation tools to reduce runtime or avoid database requirements.

When a manifest is provided, the pipeline skips the corresponding internal runs (e.g. does not invoke InterProScan or AMRFinderPlus internally) and instead routes the manifest files directly into the downstream integration and reporting steps. Functional annotation tools not covered by the manifest still run normally.

## CSV format

```
sample,ips_tsv,amrfinder_tsv,antismash_gff,gecco_gff,sanntis_gff
```

- `sample` — required; must match the `sample` column in the input samplesheet exactly.
- All other columns are optional. Leave the field empty (`,`) to skip that tool's manifest for a given sample.

Supported file extensions: `.tsv` / `.tsv.gz` for tabular files, `.gff` / `.gff.gz` for GFF files.

### Example

```csv
sample,ips_tsv,amrfinder_tsv,antismash_gff,gecco_gff,sanntis_gff
MGYG000528118,/data/catalogues/MGYG000528118_InterProScan.tsv,/data/catalogues/MGYG000528118_amrfinderplus.tsv,/data/catalogues/MGYG000528118_antismash.gff,/data/catalogues/MGYG000528118_gecco.gff,/data/catalogues/MGYG000528118_sanntis.gff
MGYG000012345,/data/catalogues/MGYG000012345_InterProScan.tsv,,,, 
```

In the second row, only the IPS TSV is provided; all BGC tools will run internally for that sample.

## Column details

| Column | Source tool | Used by |
|---|---|---|
| `ips_tsv` | InterProScan | PATHOFACT2 (CDS search step), COMBINEREPORTER |
| `amrfinder_tsv` | AMRFinderPlus (catalogue layout) | AMR_ANNOTATION subworkflow (REFORMAT step), COMBINEREPORTER |
| `antismash_gff` | antiSMASH | BGC_ANNOTATION subworkflow, COMBINEREPORTER |
| `gecco_gff` | GECCO | BGC_ANNOTATION subworkflow, COMBINEREPORTER |
| `sanntis_gff` | SanntiS | BGC_ANNOTATION subworkflow, COMBINEREPORTER |

### AMRFinderPlus column layout

The AMRFinderPlus output produced by the genomes-catalogue-pipeline uses a column layout that differs from the MAP-internal run. MAP automatically normalises the catalogue-format TSV via an internal reformat step before integration, so no manual preprocessing is required.

## Interaction with skip flags

Manifest columns and skip flags are independent:

- If a column is populated in the manifest, the corresponding internal tool is bypassed regardless of skip flags.
- If a column is empty and the skip flag is set (e.g. `--skip_amrfinderplus`), the tool is skipped and no results for that tool appear in the report.
- If a column is empty and the skip flag is not set, the tool runs internally as usual.

## Invocation

```bash
nextflow run ebi-metagenomics/mobilome-annotation-pipeline \
    --input samplesheet.csv \
    --annotation_manifest manifest.csv \
    -c my_paths.config \
    -profile singularity
```

The manifest is validated against `assets/schema_manifest.json` at startup. Sample IDs that appear in the manifest but not in the samplesheet are ignored with a warning.
