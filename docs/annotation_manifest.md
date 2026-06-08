# Annotation Manifest

The `--annotation_manifest` option allows functional annotation outputs produced by the [MGnify genomes catalogue pipeline](https://github.com/EBI-Metagenomics/genomes-catalogue-pipeline) to be reused directly, bypassing the corresponding internal tools. This avoids redundant compute when processing genomes that have already been annotated in a catalogue context.

## When to use it

Use `--annotation_manifest` when:

- You are running MAP on assemblies that are part of an MGnify genomes catalogue pipeline run and you already have InterProScan, AMRFinderPlus, antiSMASH, GECCO, or SanntiS outputs for those assemblies.
- You want to skip re-running one or more annotation tools to reduce runtime or avoid database requirements.

When a manifest is provided, the pipeline skips the corresponding internal runs (e.g. does not invoke InterProScan or AMRFinderPlus) and instead routes the manifest files directly into the downstream integration and reporting steps. Functional annotation tools not covered by the manifest still run normally.

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

The manifest and skip flags operate on different tool sets:

- **Tools covered by the manifest** (InterProScan, AMRFinderPlus, antiSMASH, GECCO, SanntiS): when `--annotation_manifest` is provided, these tools are always bypassed — skip flags have no effect. If a column is populated, the manifest file is used directly; if a column is empty, no results are produced for that tool.
- **Tools not covered by the manifest** (DeepARG, RGI): skip flags (`--skip_deeparg`, `--skip_rgi`) always work normally, regardless of whether a manifest is provided.

### Examples

**Skip a BGC tool — no manifest**

Reduce runtime by skipping antiSMASH when you have no precomputed results:

```bash
nextflow run ebi-metagenomics/mobilome-annotation-pipeline \
    --input samplesheet.csv \
    --skip_antismash \
    -profile singularity
```

**Partial manifest + skip an unrelated tool**

You have precomputed antiSMASH and GECCO outputs, and want to skip RGI:

```bash
nextflow run ebi-metagenomics/mobilome-annotation-pipeline \
    --input samplesheet.csv \
    --annotation_manifest manifest.csv \
    --skip_rgi \
    -profile singularity
```

The `antismash_gff` and `gecco_gff` columns are used directly; SanntiS, AMRFinderPlus, and IPS run internally. `--skip_rgi` takes effect because RGI has no manifest column.

**Skip flag alongside a populated manifest column — no effect**

`--skip_antismash` is redundant when `antismash_gff` is populated in the manifest; the tool is already bypassed even if there's no entry in the annotation mannifest:

```bash
# --skip_antismash has no effect here — the manifest already bypasses antiSMASH
nextflow run ebi-metagenomics/mobilome-annotation-pipeline \
    --input samplesheet.csv \
    --annotation_manifest manifest.csv \
    --skip_antismash \
    -profile singularity
```

## Invocation

```bash
nextflow run ebi-metagenomics/mobilome-annotation-pipeline \
    --input samplesheet.csv \
    --annotation_manifest manifest.csv \
    -c my_paths.config \
    -profile singularity
```

The manifest is validated against `assets/schema_manifest.json` at startup. Sample IDs that appear in the manifest but not in the samplesheet cause the pipeline to exit with an error.
