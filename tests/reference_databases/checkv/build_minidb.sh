#!/usr/bin/env bash
#
# Build a custom CheckV mini-database for the `pos_test_v2` positive test.
#
# CheckV's full v1.5 database is ~12 GB and cannot be committed. This script subsets it
# down to only the reference proteins and HMM profiles that the `pos_test_v2` assembly
# actually matches, so the positive nf-test exercises a *real* CheckV run (genuine viral
# matches surviving `quality_decision()`), not the all-zero nf-core smoke DB.
#
# The hit lists are derived from a real-database CheckV run on tests/data/pos_test_v2.fasta:
#   - genome_db: every protein target that appears in tmp/diamond.tsv (rebuilds .dmnd) so
#     AAI-based completeness/quality is reproduced exactly (this drives the High-quality
#     contig that must survive `quality_decision()`).
#   - hmm_db   : only the per-gene *best-hit* HMM profiles (tmp/gene_features.tsv col 9).
#     These fully determine each gene's viral/host category, so `viral_genes` counts are
#     preserved; keeping the full hmmsearch hit set instead would bloat the HMM file to
#     ~164 MB for no behavioural gain (best-hit only is ~12 MB).
#
# Requires `diamond`, `hmmfetch` and `awk` on PATH — all bundled in the CheckV biocontainer.
# Run it inside that container, e.g.:
#
#   docker run --rm --platform linux/amd64 \
#     -v /path/to/checkv-db-v1.5:/db:ro \
#     -v /path/to/checkv_results/pos_test_v2:/work:ro \
#     -v "$PWD/tests/reference_databases/checkv":/out \
#     -v /path/to/scratch:/scratch \
#     -e FULL_DB=/db -e WORKDIR=/work -e OUT=/out/checkv_minimal_db -e SCRATCH=/scratch \
#     quay.io/biocontainers/checkv:1.0.3--pyhdfd78af_0 \
#     bash /out/build_minidb.sh
#
set -euo pipefail

FULL_DB="${FULL_DB:?set FULL_DB to the full checkv-db-v1.5 directory}"
WORKDIR="${WORKDIR:?set WORKDIR to the CheckV run dir for pos_test_v2 (contains tmp/)}"
OUT="${OUT:?set OUT to the target checkv_minimal_db directory}"
SCRATCH="${SCRATCH:-/tmp}"

echo ">> Output: $OUT"
rm -rf "$OUT"
mkdir -p "$OUT/genome_db" "$OUT/hmm_db/checkv_hmms" "$SCRATCH"

############################################################################
# 1. genome_db — DIAMOND reference subset
############################################################################
echo ">> Collecting protein + genome hit lists from diamond.tsv"
cut -f2 "$WORKDIR/tmp/diamond.tsv" | sort -u > "$SCRATCH/protlist.txt"
sed 's/_[0-9]*$//' "$SCRATCH/protlist.txt" | sort -u > "$SCRATCH/genomelist.txt"
echo "   proteins: $(wc -l < "$SCRATCH/protlist.txt")  genomes: $(wc -l < "$SCRATCH/genomelist.txt")"

echo ">> Subsetting checkv_reps.faa"
awk 'NR==FNR{keep[$1]=1; next} /^>/{id=substr($1,2); p=(id in keep)} p' \
    "$SCRATCH/protlist.txt" "$FULL_DB/genome_db/checkv_reps.faa" > "$OUT/genome_db/checkv_reps.faa"
echo "   kept $(grep -c '^>' "$OUT/genome_db/checkv_reps.faa") proteins"

echo ">> Building DIAMOND db"
diamond makedb --in "$OUT/genome_db/checkv_reps.faa" -d "$OUT/genome_db/checkv_reps" --quiet

echo ">> Subsetting checkv_reps.tsv (genome metadata)"
awk -F'\t' 'NR==FNR{keep[$1]=1; next} FNR==1{print; next} ($1 in keep)' \
    "$SCRATCH/genomelist.txt" "$FULL_DB/genome_db/checkv_reps.tsv" > "$OUT/genome_db/checkv_reps.tsv"

echo ">> Copying checkv_error.tsv (small confidence lookup table, kept whole)"
cp "$FULL_DB/genome_db/checkv_error.tsv" "$OUT/genome_db/checkv_error.tsv"

############################################################################
# 2. hmm_db — HMM profile subset
############################################################################
echo ">> Collecting per-gene best-hit HMM names from gene_features.tsv (col 9)"
cut -f9 "$WORKDIR/tmp/gene_features.tsv" | tail -n +2 | grep -v '^NA$' | sort -u > "$SCRATCH/hmmlist.txt"
echo "   HMM profiles: $(wc -l < "$SCRATCH/hmmlist.txt")"

echo ">> Concatenating + indexing full HMM set (this is the slow/large step)"
cat "$FULL_DB"/hmm_db/checkv_hmms/*.hmm > "$SCRATCH/all.hmm"
hmmfetch --index "$SCRATCH/all.hmm" >/dev/null
echo ">> Fetching the matched profiles"
hmmfetch -f -o "$OUT/hmm_db/checkv_hmms/checkv_hmms.hmm" "$SCRATCH/all.hmm" "$SCRATCH/hmmlist.txt"

echo ">> Subsetting checkv_hmms.tsv (HMM->category, col3=hmm)"
awk -F'\t' 'NR==FNR{keep[$1]=1; next} FNR==1{print; next} ($3 in keep)' \
    "$SCRATCH/hmmlist.txt" "$FULL_DB/hmm_db/checkv_hmms.tsv" > "$OUT/hmm_db/checkv_hmms.tsv"

echo ">> Subsetting genome_lengths.tsv (HMM completeness stats, col1=hmm)"
awk -F'\t' 'NR==FNR{keep[$1]=1; next} FNR==1{print; next} ($1 in keep)' \
    "$SCRATCH/hmmlist.txt" "$FULL_DB/hmm_db/genome_lengths.tsv" > "$OUT/hmm_db/genome_lengths.tsv"

echo ">> Cleaning scratch index files"
rm -f "$SCRATCH/all.hmm" "$SCRATCH/all.hmm.ssi"

echo ">> Done. Tree:"
( cd "$OUT" && find . -type f -exec ls -lh {} \; )
