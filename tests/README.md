# Test dataset

The nf-test suite (`tests/default.nf.test`) runs several independent samples covering positive, negative, schema-validation, and functional-annotation (PathoFact2) scenarios. Each scenario is driven by a samplesheet CSV in this directory.

## Samplesheets

| Samplesheet | Sample | Inputs supplied | Backing nf-test | Purpose |
|---|---|---|---|---|
| [`pos_test_samplesheet.csv`](https://raw.githubusercontent.com/EBI-Metagenomics/mobilome-annotation-pipeline/refs/heads/dev/tests/pos_test_samplesheet.csv) | `pos_test` | assembly `pos_test_v2.fasta` (8 contigs) + `virify_gff` | `tag "positive"` | Mobilome-only run exercising IS, integron, ICE, plasmid, prophage (VIRify), and compositional outlier detection |
| [`neg_test_samplesheet.csv`](https://raw.githubusercontent.com/EBI-Metagenomics/mobilome-annotation-pipeline/refs/heads/dev/tests/neg_test_samplesheet.csv) | `neg_test` | assembly `neg_test.fasta.gz` (2 contigs) | `tag "negative"` | No mobile genetic elements; expects an empty mobilome GFF output |
| [`pathofact_samplesheet.csv`](https://raw.githubusercontent.com/EBI-Metagenomics/mobilome-annotation-pipeline/refs/heads/dev/tests/pathofact_samplesheet.csv) | `MGYG000528118_sub` | assembly `.fna` + `proteins_gff` + `proteins_faa` + `interproscan_tsv` | `tag "pathofact"` | Full functional-annotation run producing a PathoFact2 report; paired with `annotation_manifest.csv` (see below) and run with `skip_deeparg=true` |
| [`invalid_proteins_samplesheet.csv`](https://raw.githubusercontent.com/EBI-Metagenomics/mobilome-annotation-pipeline/refs/heads/dev/tests/invalid_proteins_samplesheet.csv) | `neg_test` | assembly `neg_test.fasta.gz` + `proteins_gff` but **no** `proteins_faa` | `tag "schema"` | Negative schema test: `proteins_gff` requires `proteins_faa` (`dependentRequired`), so the workflow is expected to fail validation |

### Positive test

The positive test assembly (`pos_test_v2.fasta`) is a compressed FASTA of eight contigs concatenated from different sources, between them covering every type of mobile genetic element the pipeline can detect. A VIRify v3.0.0 GFF (`virify_gff`) is supplied for the viral predictions; protein-coding genes are predicted by the pipeline's own Prodigal step rather than provided in the samplesheet.

### Negative test

The negative test assembly (`neg_test.fasta.gz`) is a compressed FASTA of two contigs with no mobile genetic elements, so the run is expected to produce an empty mobilome GFF.

## Samplesheet columns

All samplesheets share the header:

```
sample,assembly,proteins_gff,proteins_faa,virify_gff,interproscan_tsv
```

Only `sample` and `assembly` are mandatory; the remaining columns are optional per-sample inputs. All paths point to raw GitHub URLs pinned to the `dev` branch.

## annotation_manifest.csv

`tests/data/MGYG000528118_sub/annotation_manifest.csv` is a companion file (not a samplesheet) with the header `sample,amrfinder_tsv,antismash_gff,gecco_gff,sanntis_gff`. It provides pre-computed AMR and BGC annotations and is wired into the PathoFact2 test via the `annotation_manifest` parameter, letting that test run without invoking the AMR/BGC tools directly.

## Databases

The databases have been manually curated for testing purposes and aggressively trimmed down, with the exception of the ICEfinder2 databases for HMMER and MaCSyFinder.
