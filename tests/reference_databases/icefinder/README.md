# Test Databases for ICEFinder 2 Lite

ICEFinder 2 lite requires three reference databases (HMM models, MacSyFinder models,
and a Prokka-formatted UniProt BLAST DB). These are **not committed** to the repo.

The full DBs are downloaded on the fly from the EBI FTP (the same 61 MB `icf2_dbs.tar.gz`
the production `DB_DOWNLOAD_MOBILOME_DBS` module uses), and extracted into `icf2_dbs/`
here (git-ignored). Fetch them before running the nf-test suite:

```bash
task fetch-icefinder-test-db
```

This downloads and unpacks:

```
icf2_dbs/
├── icehmm/                    # icescan.hmm.* (hmmpressed) + ICEfinder.hmm.*
├── macsydata/ICEscan/         # full MacSyFinder model set (T4SS typeB..typeG, AICE, IME)
└── icefinder_prokka_uniprot/  # prokka_uniprot_sprot.fasta.* BLAST DB
```

`conf/test.config` points the `icefinder_*` params at these paths, and CI fetches the
same tarball before `nf-test test` (see `.github/workflows/full_pipeline_test.yml`).

## Why the full DB instead of a mini-DB?

The previous hand-curated mini-DB only modelled AICE-type systems, so it could never
detect the T4SS/MOBH/typeG ICE present in the positive-test assembly — `ices.tsv` was
always empty. Using the full model set (the same one the cluster runs use) gives a real
positive ICE detection end-to-end, and removes ~60 MB of binaries from the repo tree at
the cost of a 61 MB runtime fetch.
