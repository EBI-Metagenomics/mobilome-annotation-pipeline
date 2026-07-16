"""Generate subsetted annotation test files for the annotation_manifest nf-test."""
import re

SRC = "/hps/nobackup/rdf/metagenomics/service-team/projects/infectome/map_pathofact_dev/test/data"
OUT = "/hps/nobackup/rdf/metagenomics/service-team/projects/infectome/map_pathofact_dev/mobilome-annotation-pipeline/tests/data/MGYG000528118_sub"

WIN_START = 1_000_000
WIN_END   = 1_100_000
OFFSET    = WIN_START - 1       # 999_999  (subtract to get 1-based coords in new contig)
OLD_CONTIG = "MGYG000528118_1"
NEW_CONTIG = "MGYG000528118_1_sub"

# ── 1. Collect protein IDs that fall entirely within the window ─────────────
in_window = set()
with open(f"{SRC}/MGYG000528118.gff") as fh:
    for line in fh:
        if line.startswith("#") or not line.strip():
            continue
        parts = line.split("\t")
        if parts[0] != OLD_CONTIG or parts[2] != "CDS":
            continue
        start, end = int(parts[3]), int(parts[4])
        if start >= WIN_START and end <= WIN_END:
            m = re.search(r"locus_tag=([^;]+)", parts[8])
            if m:
                in_window.add(m.group(1))

print(f"Proteins in window: {len(in_window)}")

# ── 2. Subset + offset user proteins GFF (all feature types) ───────────────
with open(f"{SRC}/MGYG000528118.gff") as fh, \
     open(f"{OUT}/MGYG000528118_sub.gff", "w") as out:
    for line in fh:
        if line.startswith("#"):
            out.write(line)
            continue
        if not line.strip():
            continue
        parts = line.split("\t")
        if parts[0] != OLD_CONTIG:
            continue
        start, end = int(parts[3]), int(parts[4])
        if start < WIN_START or end > WIN_END:
            continue
        parts[0] = NEW_CONTIG
        parts[3] = str(start - OFFSET)
        parts[4] = str(end   - OFFSET)
        # Update contig references inside the attributes column
        parts[8] = parts[8].replace(OLD_CONTIG, NEW_CONTIG)
        out.write("\t".join(parts))

print("GFF done")

# ── 3. Subset user proteins FAA ────────────────────────────────────────────
with open(f"{SRC}/MGYG000528118.faa") as fh, \
     open(f"{OUT}/MGYG000528118_sub.faa", "w") as out:
    write = False
    for line in fh:
        if line.startswith(">"):
            protein_id = line[1:].split()[0]
            write = protein_id in in_window
        if write:
            out.write(line)

print("FAA done")

# ── 4. Subset InterProScan TSV (protein ID in col 0, no header) ───────────
with open(f"{SRC}/MGYG000528118_InterProScan.tsv") as fh, \
     open(f"{OUT}/MGYG000528118_sub_InterProScan.tsv", "w") as out:
    for line in fh:
        if not line.strip():
            continue
        protein_id = line.split("\t")[0]
        if protein_id in in_window:
            out.write(line)

print("IPS TSV done")

# ── 5. AMRFinder TSV — header only (no hits in this region) ───────────────
with open(f"{SRC}/MGYG000528118_amrfinderplus.tsv") as fh, \
     open(f"{OUT}/MGYG000528118_sub_amrfinderplus.tsv", "w") as out:
    out.write(next(fh))   # write header line only

print("AMRFinder TSV done")

# ── 6. Subset + offset antiSMASH GFF (T3PKS region2 + its children only) ──
ANTISMASH_REGION = "MGYG000528118_1_region2"
with open(f"{SRC}/MGYG000528118_antismash.gff") as fh, \
     open(f"{OUT}/MGYG000528118_sub_antismash.gff", "w") as out:
    for line in fh:
        if line.startswith("#"):
            out.write(line)
            continue
        if not line.strip():
            continue
        parts = line.split("\t")
        if parts[0] != OLD_CONTIG:
            continue
        # Keep only region2 (T3PKS) and its gene children
        attrs = parts[8]
        is_region2 = f"ID={ANTISMASH_REGION}" in attrs
        is_child   = f"Parent={ANTISMASH_REGION}" in attrs
        if not (is_region2 or is_child):
            continue
        start, end = int(parts[3]), int(parts[4])
        parts[0] = NEW_CONTIG
        parts[3] = str(start - OFFSET)
        parts[4] = str(end   - OFFSET)
        parts[8] = attrs.replace(OLD_CONTIG, NEW_CONTIG)
        out.write("\t".join(parts))

print("antiSMASH GFF done")

# ── 7. GECCO GFF — header only (copy as-is) ────────────────────────────────
with open(f"{SRC}/MGYG000528118_gecco.gff") as fh, \
     open(f"{OUT}/MGYG000528118_sub_gecco.gff", "w") as out:
    out.write(fh.read())

print("GECCO GFF done")

# ── 8. SanntiS GFF — manufacture a CLUSTER matching the T3PKS window ───────
# T3PKS region in new coords: 1053466-OFFSET=53467  to  1094620-OFFSET=94621
sanntis_start = 1_053_466 - OFFSET
sanntis_end   = 1_094_620 - OFFSET
with open(f"{OUT}/MGYG000528118_sub_sanntis.gff", "w") as out:
    out.write("##gff-version 3\n")
    attrs = (
        f"ID={NEW_CONTIG}_sanntis_1;"
        "nearest_MiBIG=BGC0000792;"
        "nearest_MiBIG_class=Saccharide;"
        "nearest_MiBIG_diceDistance=0.353;"
        "score=0.850;"
        "partial=00"
    )
    out.write(
        f"{NEW_CONTIG}\tSanntiSv0.9.3.1\tCLUSTER\t"
        f"{sanntis_start}\t{sanntis_end}\t.\t.\t.\t{attrs}\n"
    )

print("SanntiS GFF done")
print("\nAll files generated.")
