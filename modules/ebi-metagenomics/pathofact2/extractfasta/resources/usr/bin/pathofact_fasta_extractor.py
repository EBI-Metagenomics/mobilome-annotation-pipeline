#!/usr/bin/env python

# Copyright 2024-2026 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

import argparse
import csv
import gzip
import logging
import os
from collections import namedtuple

from Bio import SearchIO, SeqIO

logger = logging.getLogger(__name__)

# A lightweight container for best-hit info (replaces the custom class)
BestBlastHit = namedtuple("BestBlastHit", ["sseqid", "evalue", "bitscore"])


def _open_text_maybe_gzip(path):
    """
    Return a text-mode file handle for plain or gzipped files.
    Caller is responsible for closing it (use a context manager).
    """
    if path.endswith(".gz"):
        return gzip.open(path, "rt", encoding="utf-8", errors="ignore")
    return open(path, encoding="utf-8", errors="ignore")


def _file_is_readable(path):
    """Return True if path exists and is non-empty."""
    if not path:
        return False
    if not os.path.exists(path):
        logger.warning("File not found: %s", path)
        return False
    if os.path.getsize(path) == 0:
        logger.info("Skipping empty file: %s", path)
        return False
    return True


def read_fasta_to_dict(fasta_path):
    """
    Read a FASTA file (optionally .gz) and return a dict of {sequence_id: SeqRecord}.
    """
    sequences = {}
    if not _file_is_readable(fasta_path):
        return sequences

    with _open_text_maybe_gzip(fasta_path) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            sequences[record.id] = record

    logger.info("Loaded %s sequences from FASTA: %s", len(sequences), fasta_path)
    return sequences


def parse_blastp_best_hits(blastp_out):
    """
    Parse DIAMOND blastp format 6 output with columns:
    qseqid sseqid pident length qlen slen evalue bitscore

    Select best hit per query:
      - prefer highest bitscore
      - then lowest evalue

    Store: {qseqid: BestBlastHit}
    """
    best = {}
    if not _file_is_readable(blastp_out):
        return best

    fields = [
        "qseqid",
        "sseqid",
        "pident",
        "length",
        "qlen",
        "slen",
        "evalue",
        "bitscore",
    ]

    with _open_text_maybe_gzip(blastp_out) as handle:
        for qresult in SearchIO.parse(handle, "blast-tab", fields=fields):
            qseqid = qresult.id
            current = best.get(qseqid)

            for hit in qresult:
                sseqid = (hit.id or "").split("(")[0]
                for hsp in hit.hsps:
                    candidate = BestBlastHit(
                        sseqid=sseqid,
                        evalue=hsp.evalue,
                        bitscore=hsp.bitscore,
                    )

                    if current is None:
                        best[qseqid] = candidate
                        current = candidate
                        continue

                    if candidate.bitscore > current.bitscore:
                        best[qseqid] = candidate
                        current = candidate
                    elif (
                        candidate.bitscore == current.bitscore
                        and candidate.evalue < current.evalue
                    ):
                        best[qseqid] = candidate
                        current = candidate

    logger.info(
        "Selected best BLASTP hits for %s queries from: %s", len(best), blastp_out
    )
    return best


def parse_pathofact2_predictions_tsv(tsv_path, threshold):
    """
    Parse PathoFact2 TSV output and filter by a probability threshold.

    Expected header:
      Sequence    Prediction    Probability

    Returns:
      {sequence_id: probability} for predictions with probability >= threshold

    Raises:
      ValueError: if the file structure is not as expected.
    """
    preds = {}
    if not _file_is_readable(tsv_path):
        return preds

    expected_header = ["Sequence", "Prediction", "Probability"]

    with _open_text_maybe_gzip(tsv_path) as handle:
        reader = csv.reader(handle, delimiter="\t")

        header = next(reader, None)
        if header is None:
            raise ValueError(f"Empty PathoFact2 TSV (no header): {tsv_path}")

        header = [h.strip() for h in header]
        if header[:3] != expected_header:
            raise ValueError(
                f"Unexpected header in {tsv_path}. "
                f"Expected first 3 columns {expected_header} but got {header[:3]}"
            )

        for line_no, row in enumerate(reader, start=2):
            if len(row) != 3:
                raise ValueError(
                    f"Malformed row in '{tsv_path}' at line {line_no}: {row}"
                )

            seq_id = row[0]
            try:
                probability = float(row[2])
            except ValueError as e:
                raise ValueError(
                    f"Invalid Probability value '{row[2]}' "
                    f"in '{tsv_path}' at line {line_no}"
                ) from e

            if probability >= threshold:
                preds[seq_id] = probability

    logger.info(
        "Loaded %s PathoFact2 predictions from: %s (threshold=%s)",
        len(preds),
        tsv_path,
        threshold,
    )
    return preds


def collect_detected_sequence_ids(blast_hits, tox_preds, vf_preds):
    """Union of all sequence IDs detected by any method."""
    ids = set()
    ids.update(blast_hits.keys())
    ids.update(tox_preds.keys())
    ids.update(vf_preds.keys())
    return ids


def write_detected_fasta(sequences, detected_ids, output_fasta):
    """
    Write a FASTA containing only sequences whose IDs are in detected_ids.

    Raises:
        ValueError: If any detected protein IDs are missing from the FASTA file.
    """
    records = []
    missing_ids = []

    for seq_id in sorted(detected_ids):
        rec = sequences.get(seq_id)
        if rec is None:
            missing_ids.append(seq_id)
            continue
        records.append(rec)

    # Check for missing sequences - this should never happen
    if missing_ids:
        # Show first 10 missing IDs, then summarize if more exist
        max_display = 10
        displayed_ids = missing_ids[:max_display]
        id_list = "\n".join(f"  - {seq_id}" for seq_id in displayed_ids)

        if len(missing_ids) > max_display:
            id_list += f"\n  ... and {len(missing_ids) - max_display} more"

        error_msg = (
            f"FATAL ERROR: {len(missing_ids)} predicted protein(s) not found in FASTA file.\n"
            f"This indicates a mismatch between prediction results and input sequences.\n"
            f"Missing protein IDs:\n{id_list}"
        )
        logger.error(error_msg)
        raise ValueError(error_msg)

    SeqIO.write(records, output_fasta, "fasta")
    logger.info("Wrote %s sequences to: %s", len(records), output_fasta)


def _fmt_prob(value, ndigits=6):
    """
    Format probability floats neatly for TSV output.
    Keeps a fixed number of decimals for column alignment.
    """
    return ("{0:." + str(ndigits) + "f}").format(float(value))


def write_support_table(blast_hits, tox_preds, vf_preds, output_tsv):
    """
    Write support TSV with header:
    sequence_id    detection_method    support_value_type    support_value    vfdb_hit

    detection_method:
      - blastp
      - pathofact2_tox
      - pathofact2_vf

    support_value_type:
      - blastp -> evalue
      - pathofact2_tox / pathofact2_vf -> probability
    """
    with open(output_tsv, "w", encoding="utf-8", newline="") as out:
        writer = csv.writer(out, delimiter="\t", lineterminator="\n")

        writer.writerow(
            [
                "sequence_id",
                "detection_method",
                "support_value_type",
                "support_value",
                "vfdb_hit",
            ]
        )

        for seq_id in sorted(blast_hits.keys()):
            hit = blast_hits[seq_id]
            writer.writerow([seq_id, "blastp", "evalue", hit.evalue, hit.sseqid])

        for seq_id in sorted(tox_preds.keys()):
            writer.writerow(
                [
                    seq_id,
                    "pathofact2_tox",
                    "probability",
                    _fmt_prob(tox_preds[seq_id]),
                    "",
                ]
            )

        for seq_id in sorted(vf_preds.keys()):
            writer.writerow(
                [
                    seq_id,
                    "pathofact2_vf",
                    "probability",
                    _fmt_prob(vf_preds[seq_id]),
                    "",
                ]
            )

    logger.info("Wrote support table to: %s", output_tsv)


def main():
    parser = argparse.ArgumentParser(
        description="Extract PathoFact2 candidate proteins from FASTA and report support from BLASTP and PathoFact2."
    )
    parser.add_argument(
        "-f", "--fasta", required=True, help="Protein sequences FASTA (optionally .gz)"
    )
    parser.add_argument(
        "-b", "--blastp_out", required=True, help="DIAMOND blastp tabular output"
    )
    parser.add_argument(
        "-t", "--pathofact2_tox", required=True, help="PathoFact2 toxins TSV"
    )
    parser.add_argument(
        "-v", "--pathofact2_vf", required=True, help="PathoFact2 VF TSV"
    )
    parser.add_argument("-o", "--output_prefix", required=True, help="Output prefix")

    parser.add_argument(
        "--threshold_vf",
        type=float,
        default=0.9,
        help="Minimum probability threshold for VF predictions (default: 0.9)",
    )
    parser.add_argument(
        "--threshold_tox",
        type=float,
        default=0.6,
        help="Minimum probability threshold for toxin predictions (default: 0.6)",
    )

    parser.add_argument(
        "--verbose", action="store_true", default=False, help="Enable DEBUG logging"
    )
    args = parser.parse_args()

    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[logging.StreamHandler()],
    )

    logger.info("Reading FASTA")
    seqs = read_fasta_to_dict(args.fasta)

    logger.info("Parsing BLASTP output (best hit per query)")
    blast_hits = parse_blastp_best_hits(args.blastp_out)

    logger.info(
        "Parsing PathoFact2 toxin predictions (threshold=%s)", args.threshold_tox
    )
    tox_preds = parse_pathofact2_predictions_tsv(
        args.pathofact2_tox, args.threshold_tox
    )

    logger.info("Parsing PathoFact2 VF predictions (threshold=%s)", args.threshold_vf)
    vf_preds = parse_pathofact2_predictions_tsv(args.pathofact2_vf, args.threshold_vf)

    detected_ids = collect_detected_sequence_ids(blast_hits, tox_preds, vf_preds)
    logger.info("Total detected sequence IDs (union of methods): %s", len(detected_ids))

    out_fasta = f"{args.output_prefix}_pathofact2.fasta"
    out_tsv = f"{args.output_prefix}_support.tsv"

    logger.info("Writing outputs")
    write_detected_fasta(seqs, detected_ids, out_fasta)
    write_support_table(blast_hits, tox_preds, vf_preds, out_tsv)

    logger.info("Done")


if __name__ == "__main__":
    main()
