# Tiny VFDB test database

Used by the `pathofact` nf-test (via `virulencefactors_db` in `conf/test.config`) so the test
does not download the full VFDB from www.mgc.ac.cn on every run.

It holds the 43 VFDB setB proteins that the test proteins
(`tests/data/MGYG000528118_sub/MGYG000528118_sub.faa`) hit when searched against the full
VFDB with the pipeline's `DIAMOND_BLASTP` settings. The mini DB reproduces exactly the same
query/subject hit pairs as the full DB.

## How it was built

Using DIAMOND 2.1.12 (`quay.io/biocontainers/diamond:2.1.12--hdb4b4cc_1`, the same version as
the `DIAMOND_BLASTP` module), from VFDB_setB_pro.fas.gz downloaded on 2026-10-14
(30,215 sequences).

```bash
# 1. Search the test proteins against the full VFDB with the pipeline's DIAMOND_BLASTP args
diamond makedb --in VFDB_setB_pro.fas.gz -d full
diamond blastp -q MGYG000528118_sub.faa -d full.dmnd -o hits.tsv \
    --outfmt 6 qseqid sseqid pident length qlen slen evalue bitscore \
    --max-target-seqs 25 --evalue 1e-5 --id 30 --query-cover 60 --subject-cover 50

# 2. Keep only the hit subjects (byte-wise: some VFDB headers are not valid UTF-8)
cut -f2 hits.tsv | sort -u > ids.txt
zcat VFDB_setB_pro.fas.gz | python3 -c '
import sys
ids = set(open("ids.txt").read().split())
keep = False
for line in sys.stdin.buffer:
    if line.startswith(b">"):
        keep = line[1:].split()[0].decode() in ids
    if keep:
        sys.stdout.buffer.write(line)
' > VFDB_setB_pro_mini.fas

# 3. Build the DIAMOND DB
diamond makedb --in VFDB_setB_pro_mini.fas -d VFDB_setB_pro_mini
```

If the test proteins or the `DIAMOND_BLASTP` args in `conf/modules.config` change, rebuild it.
