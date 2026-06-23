# CheckV mini test database

`checkv_minimal_db/` is a custom subset of the official **CheckV database v1.5**
(Jan 2023), built specifically so the `pos_test_v2` positive nf-test runs a *real* CheckV
`end_to_end` and yields genuine viral matches.

## Why not the nf-core smoke DB?

The nf-core minimal CheckV DB (one HMM profile, two reference genomes) reports
`viral_genes=0` / `checkv_quality=Not-determined` for every contig. CheckV gates geNomad
viral predictions through `quality_decision()` (`bin/virify_qc.py`), so that DB would filter
out **all** viral predictions and make the "positive" test non-positive. This subset instead
contains exactly the reference proteins and HMM profiles that `pos_test_v2` matches.

## Contents

```
checkv_minimal_db/
├── genome_db/
│   ├── checkv_reps.dmnd     # DIAMOND db rebuilt from the subset .faa
│   ├── checkv_reps.faa      # the 12,216 reference proteins hit by pos_test_v2
│   ├── checkv_reps.tsv      # genome metadata for the 8,833 hit reference genomes
│   └── checkv_error.tsv     # completeness confidence lookup (kept whole, small)
└── hmm_db/
    ├── checkv_hmms/
    │   └── checkv_hmms.hmm  # the 71 per-gene best-hit HMM profiles
    ├── checkv_hmms.tsv      # HMM → category metadata for those 71 profiles
    └── genome_lengths.tsv   # HMM completeness stats for those 71 profiles
```

Total ~22 MB. The full DIAMOND/AAI reference set is kept (it drives the AAI-based
completeness that keeps the High-quality contig surviving `quality_decision()`), while the
HMM set is reduced to per-gene best hits — these alone determine each gene's viral/host
category, so `viral_genes` counts are preserved, and keeping the full hmmsearch hit set
would bloat the HMM file to ~164 MB for no behavioural gain.

## Rebuilding

`build_minidb.sh` regenerates this directory from the full CheckV DB plus the `tmp/` work
dir of a real-DB CheckV run on `tests/data/pos_test_v2.fasta`. It needs `diamond`,
`hmmfetch` and `awk`, all bundled in the CheckV biocontainer. Example invocation is in the
header comment of that script.
