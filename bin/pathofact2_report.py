#!/usr/bin/env python3

"""
Integrate PathoFact2, AMR, mobilome, BGC, and InterProScan annotations.

Output columns:
    protein_id
    vfdb_hit
    vfdb_blastp_eval
    pathofact2_tox_prob
    pathofact2_vf_prob
    amr_drug_class
    amr_tool
    amr_tool_ident
    mge_type
    bgc_type
    bgc_tools
    signalP

Rules
-----
- Seed proteins come from PathoFact2 and/or AMR GFF files.
- If both PathoFact2 and AMR inputs are absent/empty, the script exits
  without generating any output.
- BGC annotations are parsed only for proteins in the seed protein set.
- MGE assignment is based on >= 90% CDS overlap with a mobilome feature
  in the same contig.
- InterProScan is optional and only SignalP entries are retained.
- Missing values are reported as '-'.

The script accepts plain text or gzip-compressed inputs.
"""

from __future__ import annotations

import argparse
import csv
import fileinput
import logging
from collections import defaultdict
from pathlib import Path
from typing import DefaultDict


OUTPUT_HEADER = [
    "protein_id",
    "vfdb_hit",
    "vfdb_blastp_eval",
    "pathofact2_tox_prob",
    "pathofact2_vf_prob",
    "amr_drug_class",
    "amr_tool",
    "amr_tool_ident",
    "mge_type",
    "bgc_type",
    "bgc_tools",
    "signalP",
]

IGNORE_MOBILOME_FEATURES = {
    "inverted_repeat_element",
    "attC_site",
    "direct_repeat",
    "terminal_inverted_repeat_element",
}


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Parse PathoFact2, AMR, mobilome, BGC, and InterProScan files "
            "and produce a combined TSV report."
        )
    )
    parser.add_argument(
        "--pathofact2",
        required=False,
        default="",
        help="Optional PathoFact2 GFF file",
    )
    parser.add_argument(
        "--amr",
        required=False,
        default="",
        help="Optional AMR GFF file",
    )
    parser.add_argument(
        "--mobilome",
        required=False,
        default="",
        help="Optional mobilome GFF file",
    )
    parser.add_argument(
        "--bgc",
        required=False,
        default="",
        help="Optional BGC GFF file",
    )
    parser.add_argument(
        "--interproscan",
        required=False,
        default="",
        help="Optional InterProScan TSV file",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output TSV file",
    )
    parser.add_argument(
        "--log-level",
        default="INFO",
        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
        help="Logging level (default: INFO)",
    )
    return parser.parse_args()


def setup_logging(level: str) -> None:
    logging.basicConfig(
        level=getattr(logging, level),
        format="%(asctime)s - %(levelname)s - %(message)s",
    )


def path_is_missing_or_empty(path_str: str) -> bool:
    if path_str == "":
        return True

    path = Path(path_str)
    if not path.exists():
        return True

    if path.stat().st_size == 0:
        return True

    return False


def parse_attributes(attributes_field: str) -> dict[str, str]:
    attributes: dict[str, str] = {}

    for item in attributes_field.strip().split(";"):
        if not item:
            continue
        if "=" not in item:
            continue
        key, value = item.split("=", 1)
        attributes[key] = value

    return attributes


def split_csv_value(value: str) -> list[str]:
    if value in {"", "-"}:
        return []
    return [item.strip() for item in value.split(",") if item.strip()]


def merge_unique_preserving_order(existing: list[str], new_values: list[str]) -> list[str]:
    seen = set(existing)
    merged = list(existing)

    for value in new_values:
        if value not in seen:
            merged.append(value)
            seen.add(value)

    return merged


def join_or_dash(values: list[str]) -> str:
    return ",".join(values) if values else "-"


def parse_pathofact2_gff(
    gff_file: str,
) -> tuple[dict[str, tuple[str, int, int]], dict[str, tuple[str, str, str, str]]]:
    """
    Parse PathoFact2 GFF.

    Returns
    -------
    prots_coords
        protein_id -> (contig, start, end)
    pathofact_data
        protein_id -> (vfdb, blastp_eval, pathofact2_tox_prob, pathofact2_vf_prob)
    """
    logging.info("Parsing PathoFact2 GFF: %s", gff_file)

    prots_coords: dict[str, tuple[str, int, int]] = {}
    pathofact_data: dict[str, tuple[str, str, str, str]] = {}

    with fileinput.hook_compressed(gff_file, "r", encoding="utf-8", errors="ignore") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            columns = line.split("\t")
            if len(columns) != 9:
                logging.warning(
                    "Skipping malformed PathoFact2 GFF line %s with %s columns",
                    line_number,
                    len(columns),
                )
                continue

            contig, _, feature_type, start, end, _, _, _, attributes_field = columns
            if feature_type != "CDS":
                continue

            attributes = parse_attributes(attributes_field)
            protein_id = attributes.get("ID", "")
            if protein_id == "":
                logging.warning("Skipping PathoFact2 line %s without ID", line_number)
                continue

            start_i = int(start)
            end_i = int(end)
            prots_coords[protein_id] = (contig, start_i, end_i)

            vfdb = attributes.get("vfdb", "-")
            blastp_eval = attributes.get("blastp_eval", "-")
            tox_prob = attributes.get("pathofact2_tox_prob", "-")
            vf_prob = attributes.get("pathofact2_vf_prob", "-")
            pathofact_data[protein_id] = (vfdb, blastp_eval, tox_prob, vf_prob)

    logging.info(
        "Parsed PathoFact2: %s proteins with coordinates, %s proteins with PathoFact2 data",
        len(prots_coords),
        len(pathofact_data),
    )
    return prots_coords, pathofact_data


def parse_amr_gff(
    gff_file: str,
    prots_coords: dict[str, tuple[str, int, int]],
) -> dict[str, tuple[str, str, str]]:
    """Parse AMR GFF and extend prots_coords for missing proteins."""
    logging.info("Parsing AMR GFF: %s", gff_file)

    amr_data: dict[str, tuple[str, str, str]] = {}

    with fileinput.hook_compressed(gff_file, "r", encoding="utf-8", errors="ignore") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            columns = line.split("\t")
            if len(columns) != 9:
                logging.warning(
                    "Skipping malformed AMR GFF line %s with %s columns",
                    line_number,
                    len(columns),
                )
                continue

            contig, _, feature_type, start, end, _, _, _, attributes_field = columns
            if feature_type != "CDS":
                continue

            attributes = parse_attributes(attributes_field)
            protein_id = attributes.get("ID", "")
            if protein_id == "":
                logging.warning("Skipping AMR line %s without ID", line_number)
                continue

            if protein_id not in prots_coords:
                prots_coords[protein_id] = (contig, int(start), int(end))

            drug_class = attributes.get("drug_class", "-")
            amr_tool = attributes.get("amr_tool", "-")
            amr_tool_ident = attributes.get("amr_tool_ident", "-")
            amr_data[protein_id] = (drug_class, amr_tool, amr_tool_ident)

    logging.info(
        "Parsed AMR: %s proteins with AMR data, %s total proteins in seed set",
        len(amr_data),
        len(prots_coords),
    )
    return amr_data


def parse_mobilome_gff(gff_file: str) -> dict[str, list[tuple[int, int, str]]]:
    """Parse mobilome GFF into contig-level MGE intervals."""
    logging.info("Parsing mobilome GFF: %s", gff_file)

    mobilome_data: DefaultDict[str, list[tuple[int, int, str]]] = defaultdict(list)

    with fileinput.hook_compressed(gff_file, "r", encoding="utf-8", errors="ignore") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            columns = line.split("\t")
            if len(columns) != 9:
                logging.warning(
                    "Skipping malformed mobilome GFF line %s with %s columns",
                    line_number,
                    len(columns),
                )
                continue

            contig, _, feature_type, start, end, _, _, _, attributes_field = columns
            if feature_type in IGNORE_MOBILOME_FEATURES:
                continue

            attributes = parse_attributes(attributes_field)
            mobile_element_type = attributes.get("mobile_element_type", feature_type)
            mobilome_data[contig].append((int(start), int(end), mobile_element_type))

    logging.info(
        "Parsed mobilome: %s contigs with at least one retained MGE",
        len(mobilome_data),
    )
    return dict(mobilome_data)


def parse_bgc_gff(
    gff_file: str,
    prots_coords: dict[str, tuple[str, int, int]],
) -> dict[str, tuple[str, str]]:
    """
    Parse BGC GFF for CDS rows matching proteins in prots_coords.

    Returns
    -------
    bgc_data
        protein_id -> (bgc_tools, bgc_type)
    """
    logging.info("Parsing BGC GFF: %s", gff_file)

    bgc_tools_map: DefaultDict[str, list[str]] = defaultdict(list)
    bgc_type_map: DefaultDict[str, list[str]] = defaultdict(list)

    with fileinput.hook_compressed(gff_file, "r", encoding="utf-8", errors="ignore") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line or line.startswith("#"):
                continue

            columns = line.split("\t")
            if len(columns) != 9:
                logging.warning(
                    "Skipping malformed BGC GFF line %s with %s columns",
                    line_number,
                    len(columns),
                )
                continue

            _, _, feature_type, _, _, _, _, _, attributes_field = columns
            if feature_type != "CDS":
                continue

            attributes = parse_attributes(attributes_field)
            protein_id = attributes.get("ID", "")
            if protein_id == "":
                logging.warning("Skipping BGC line %s without ID", line_number)
                continue
            if protein_id not in prots_coords:
                continue

            bgc_tools_values = split_csv_value(attributes.get("bgc_tools", "-"))
            bgc_tools_map[protein_id] = merge_unique_preserving_order(
                bgc_tools_map[protein_id],
                bgc_tools_values,
            )

            bgc_type_values: list[str] = []
            bgc_type_values.extend(split_csv_value(attributes.get("nearest_MiBIG_class", "-")))
            bgc_type_values.extend(split_csv_value(attributes.get("antismash_product", "-")))
            bgc_type_values.extend(split_csv_value(attributes.get("gecco_bgc_type", "-")))

            bgc_type_map[protein_id] = merge_unique_preserving_order(
                bgc_type_map[protein_id],
                bgc_type_values,
            )

    bgc_data = {
        protein_id: (
            join_or_dash(bgc_tools_map[protein_id]),
            join_or_dash(bgc_type_map[protein_id]),
        )
        for protein_id in set(bgc_tools_map) | set(bgc_type_map)
    }

    logging.info("Parsed BGC data for %s proteins", len(bgc_data))
    return bgc_data


def parse_interproscan_signalp(tsv_file: str) -> dict[str, str]:
    """Parse optional InterProScan TSV and keep only SignalP annotations."""
    logging.info("Parsing InterProScan TSV: %s", tsv_file)

    signalp_map: DefaultDict[str, list[str]] = defaultdict(list)

    with fileinput.hook_compressed(tsv_file, "r", encoding="utf-8", errors="ignore") as handle:
        for line_number, raw_line in enumerate(handle, start=1):
            line = raw_line.rstrip("\n")
            if not line:
                continue

            columns = line.split("\t")
            if len(columns) < 5:
                logging.warning(
                    "Skipping malformed InterProScan line %s with %s columns",
                    line_number,
                    len(columns),
                )
                continue

            protein_id = columns[0]
            analysis = columns[3]
            signalp_annot = columns[4]

            if not analysis.startswith("SignalP"):
                continue
            if signalp_annot in {"", "-"}:
                continue

            signalp_map[protein_id] = merge_unique_preserving_order(
                signalp_map[protein_id],
                [signalp_annot],
            )

    signalp_data = {
        protein_id: join_or_dash(values)
        for protein_id, values in signalp_map.items()
    }
    logging.info("Parsed SignalP data for %s proteins", len(signalp_data))
    return signalp_data


def calculate_overlap_length(
    cds_start: int,
    cds_end: int,
    feature_start: int,
    feature_end: int,
) -> int:
    overlap_start = max(cds_start, feature_start)
    overlap_end = min(cds_end, feature_end)

    if overlap_end < overlap_start:
        return 0

    return overlap_end - overlap_start + 1


def assign_mge_types(
    prots_coords: dict[str, tuple[str, int, int]],
    mobilome_data: dict[str, list[tuple[int, int, str]]],
    min_fraction: float = 0.90,
) -> dict[str, str]:
    """
    Assign MGE type(s) to proteins based on interval overlap.

    A protein is considered part of an MGE when at least min_fraction of the
    protein length is covered by the MGE interval.
    """
    logging.info("Assigning MGE types to proteins using overlap threshold %.2f", min_fraction)

    protein_mge_map: dict[str, str] = {}

    for protein_id, (contig, cds_start, cds_end) in prots_coords.items():
        protein_length = cds_end - cds_start + 1
        assigned_mges: list[str] = []

        for mge_start, mge_end, mge_type in mobilome_data.get(contig, []):
            overlap_length = calculate_overlap_length(cds_start, cds_end, mge_start, mge_end)
            overlap_fraction = overlap_length / protein_length

            if overlap_fraction >= min_fraction:
                assigned_mges = merge_unique_preserving_order(assigned_mges, [mge_type])

        protein_mge_map[protein_id] = join_or_dash(assigned_mges)

    logging.info("Assigned MGE context to %s proteins", len(protein_mge_map))
    return protein_mge_map


def build_rows(
    prots_coords: dict[str, tuple[str, int, int]],
    pathofact_data: dict[str, tuple[str, str, str, str]],
    amr_data: dict[str, tuple[str, str, str]],
    protein_mge_map: dict[str, str],
    bgc_data: dict[str, tuple[str, str]],
    signalp_data: dict[str, str],
) -> list[dict[str, str]]:
    logging.info("Building output rows")

    rows: list[dict[str, str]] = []

    for protein_id in sorted(prots_coords):
        vfdb_hit, vfdb_blastp_eval, tox_prob, vf_prob = pathofact_data.get(
            protein_id,
            ("-", "-", "-", "-"),
        )
        amr_drug_class, amr_tool, amr_tool_ident = amr_data.get(
            protein_id,
            ("-", "-", "-"),
        )
        bgc_tools, bgc_type = bgc_data.get(
            protein_id,
            ("-", "-"),
        )
        signalp = signalp_data.get(protein_id, "-")
        mge_type = protein_mge_map.get(protein_id, "-")

        row = {
            "protein_id": protein_id,
            "vfdb_hit": vfdb_hit,
            "vfdb_blastp_eval": vfdb_blastp_eval,
            "pathofact2_tox_prob": tox_prob,
            "pathofact2_vf_prob": vf_prob,
            "amr_drug_class": amr_drug_class,
            "amr_tool": amr_tool,
            "amr_tool_ident": amr_tool_ident,
            "mge_type": mge_type,
            "bgc_type": bgc_type,
            "bgc_tools": bgc_tools,
            "signalP": signalp,
        }
        rows.append(row)

    logging.info("Built %s output rows", len(rows))
    return rows


def write_output(rows: list[dict[str, str]], output_file: str) -> None:
    logging.info("Writing output TSV: %s", output_file)

    output_path = Path(output_file)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=OUTPUT_HEADER, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)

    logging.info("Finished writing %s rows to %s", len(rows), output_file)


def main() -> None:
    args = parse_args()
    setup_logging(args.log_level)

    logging.info("Starting annotation integration")

    pathofact_missing = path_is_missing_or_empty(args.pathofact2)
    amr_missing = path_is_missing_or_empty(args.amr)
    mobilome_missing = path_is_missing_or_empty(args.mobilome)
    bgc_missing = path_is_missing_or_empty(args.bgc)
    interproscan_missing = path_is_missing_or_empty(args.interproscan)

    if pathofact_missing and amr_missing:
        logging.info(
            "PathoFact2 and AMR files are both absent or empty; no seed proteins found. "
            "Skipping downstream parsing and not generating an output file."
        )
        return

    prots_coords: dict[str, tuple[str, int, int]] = {}
    pathofact_data: dict[str, tuple[str, str, str, str]] = {}
    amr_data: dict[str, tuple[str, str, str]] = {}

    if not pathofact_missing:
        prots_coords, pathofact_data = parse_pathofact2_gff(args.pathofact2)
    else:
        logging.info("PathoFact2 GFF absent or empty; continuing with AMR proteins only")

    if not amr_missing:
        amr_data = parse_amr_gff(args.amr, prots_coords)
    else:
        logging.info("AMR GFF absent or empty; continuing with PathoFact2 proteins only")

    if len(prots_coords) == 0:
        logging.info(
            "No proteins were collected from PathoFact2 or AMR after parsing. "
            "Skipping downstream parsing and not generating an output file."
        )
        return

    if not mobilome_missing:
        mobilome_data = parse_mobilome_gff(args.mobilome)
    else:
        logging.info("Mobilome GFF absent or empty; mge_type will be filled with '-'")
        mobilome_data = {}

    if not bgc_missing:
        bgc_data = parse_bgc_gff(args.bgc, prots_coords)
    else:
        logging.info("BGC GFF absent or empty; bgc_type and bgc_tools will be filled with '-'")
        bgc_data = {}

    if not interproscan_missing:
        signalp_data = parse_interproscan_signalp(args.interproscan)
    else:
        logging.info("InterProScan TSV absent or empty; signalP will be filled with '-'")
        signalp_data = {}

    protein_mge_map = assign_mge_types(
        prots_coords,
        mobilome_data,
        min_fraction=0.90,
    )

    rows = build_rows(
        prots_coords=prots_coords,
        pathofact_data=pathofact_data,
        amr_data=amr_data,
        protein_mge_map=protein_mge_map,
        bgc_data=bgc_data,
        signalp_data=signalp_data,
    )

    write_output(rows, args.output)
    logging.info("Annotation integration completed successfully")


if __name__ == "__main__":
    main()
