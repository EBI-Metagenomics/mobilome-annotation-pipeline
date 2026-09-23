#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Copyright 2025 EMBL - European Bioinformatics Institute
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
import logging
import sys
import os.path
import gzip

from map_tools import mapping_names

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
COV_THRESHOLD = 0.75

# Attribute keys carried over from per-protein annotations (e.g. VIRify viphog hits)
EXTRA_ANNOT_KEYS = [
    "viphog",
    "viphog_taxonomy",
]

def normalize_attributes(attributes):
    """
    Drop trailing ';' from a GFF attributes column.

    Tools commonly end column 9 with a separator (``ID=x;``). Appending more
    attributes to such a line yields an empty field (``;;``), which is malformed
    GFF3, so the separator is removed before anything is appended.

    :param attributes: The attributes column (column 9) of a GFF line
    :return: The attributes column without trailing separators
    """
    return attributes.rstrip(";")

def open_file(filename, mode='r'):
    """
    Open a file, handling both compressed (.gz) and uncompressed files.
    
    Args:
        filename: Path to the file
        mode: File mode ('r' for read, 'w' for write, 'a' for append)
    
    Returns:
        File handle that can be used for reading/writing
    """
    if filename.endswith('.gz'):
        if 'w' in mode or 'a' in mode:
            return gzip.open(filename, mode + 't')  # Text mode for writing
        else:
            return gzip.open(filename, 'rt')  # Text mode for reading
    else:
        return open(filename, mode)

def is_file_empty(filepath):
    """
    Check if a file (compressed or uncompressed) has content.

    :param filepath: Path to the file to check
    :return: True if file is empty or doesn't exist, False if it has content
    """
    if not os.path.exists(filepath):
        logger.error(f"File not found: {filepath}")
        return True

    # Check file size first for uncompressed files
    if not filepath.endswith('.gz') and os.stat(filepath).st_size == 0:
        return True

    # For compressed files or to double-check content
    with open_file(filepath) as f:
        first_char = f.read(1)
        return len(first_char) == 0

def mobilome_parser(mobilome_clean):
    """Parse mobilome predictions from GFF file (handles compressed files)."""
    
    # Check if file exists
    if not os.path.exists(mobilome_clean):
        logger.error(f"Mobilome file not found: {mobilome_clean}")
        sys.exit(1)
    
    # Check if file has content (works with compressed files)
    if is_file_empty(mobilome_clean):
        logger.warning(f"Mobilome file is empty: {mobilome_clean}")
        return ({}, {}, {}, {})
    
    # Parsing the mobilome prediction
    proteins_annot, mobilome_annot, mges_dict, mob_types = {}, {}, {}, {}
    
    source_tools = [
        "ICEfinder",
        "IntegronFinder",
        "ISEScan",
        "geNomad",
        "VIRify",
        "geNomad",
        "geNomad_VIRify",
        "MAP",
    ]

    with open_file(mobilome_clean) as input_table:
        logger.info(f"Successfully opened mobilome file: {mobilome_clean}")
        
        for line in input_table:
            l_line = line.rstrip().split("\t")
            
            # Annotation lines have exactly 9 columns
            if len(l_line) == 9:
                # These rows are written verbatim into the outputs, so drop any
                # trailing ';' before storing them.
                l_line[8] = normalize_attributes(l_line[8])
                mge_line = "\t".join(l_line)

                contig = l_line[0]
                annot_source = l_line[1]
                seq_type = l_line[2]
                start = int(l_line[3])
                end = int(l_line[4])
                strand = l_line[6]
                coordinates = (start, end)

                if annot_source in source_tools:
                    composite_key = (contig, start, end)
                    mob_types[composite_key] = seq_type

                    if contig in mobilome_annot:
                        mobilome_annot[contig].append(mge_line)
                        mges_dict[contig].append(coordinates)
                    else:
                        mobilome_annot[contig] = [mge_line]
                        mges_dict[contig] = [coordinates]
                else:
                    str_composite_key = (contig, str(start), str(end), strand)
                    attrib = l_line[8]
                    extra_list = []
                    
                    for attr in attrib.split(";"):
                        if "=" in attr:  # Ensure attr has the expected format
                            att_key = attr.split("=")[0]
                            if att_key in EXTRA_ANNOT_KEYS:
                                extra_list.append(attr)
                    
                    if len(extra_list) > 0:
                        extra_val = ";".join(extra_list)
                        proteins_annot[str_composite_key] = extra_val
            else:
                logger.debug(
                    f"Skipping line: incorrect number of columns ({len(l_line)})"
                )
    
    # Log parsing statistics
    logger.info("Mobilome parsing completed:")
    logger.info(f"  - Unique contigs with mobilome annotations: {len(mobilome_annot)}")
    logger.info(f"  - Total protein annotations: {len(proteins_annot)}")
    
    return (proteins_annot, mobilome_annot, mges_dict, mob_types)

def parse_combined_report(report_file):
    """
    Parse the PathoFact2 combined report into a {protein_id: summary_string} map.

    Column positions are resolved from the header row so the parser is robust to
    column re-ordering. Returns an empty dict if the file is missing, empty, or does
    not contain the expected columns.
    """
    if is_file_empty(report_file):
        logger.warning(f"Combined report is empty or missing: {report_file}")
        return {}

    summary_map = {}
    with open_file(report_file) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        try:
            id_idx = header.index("protein_id")
            summary_idx = header.index("summary_string")
        except ValueError:
            logger.error(
                "Combined report missing 'protein_id' or 'summary_string' column; "
                "skipping pathofact2 annotation"
            )
            return {}

        for line in handle:
            cols = line.rstrip("\n").split("\t")
            if len(cols) <= max(id_idx, summary_idx):
                continue
            summary_map[cols[id_idx]] = cols[summary_idx]

    logger.info(f"Parsed combined report: {len(summary_map)} proteins with summary strings")
    return summary_map


def gff_updater(
    user_gff, output_prefix, proteins_annot, mobilome_annot, mges_dict, mob_types,
    summary_map=None, output_infix="_user_mobilome_", names_equiv=None,
):
    """Adding the mobilome predictions to the user file (handles compressed input/output).

    output_infix controls the output file naming: "_user_mobilome_" for user-provided
    genes, "_mobilome_" when the baseline is the Prodigal/tRNA genes GFF.

    names_equiv maps renamed contig ids -> original ids. When provided (the Prodigal/tRNA
    baseline still uses the internal renamed ids), each genes-GFF feature's contig is
    translated to its original name so it aligns with the already-renamed-back mobilome GFF.
    """

    summary_map = summary_map or {}
    names_equiv = names_equiv or {}

    extra_file = f"{output_prefix}{output_infix}extra.gff"
    full_file = f"{output_prefix}{output_infix}full.gff"
    clean_file = f"{output_prefix}{output_infix}clean.gff"

    # Check if input file exists
    if not os.path.exists(user_gff):
        logger.error(f"User GFF file not found: {user_gff}")
        sys.exit(1)

    # Check if input file has content
    if is_file_empty(user_gff):
        logger.warning(f"User GFF file is empty: {user_gff}")
        # Still create empty output files
        for output_file in (extra_file, full_file, clean_file):
            with open_file(output_file, 'w') as f:
                pass  # Create empty file
        logger.info(f"Created empty output files with prefix: {output_prefix}")
        return

    logger.info(f"Starting GFF update process with file: {user_gff}")

    used_contigs = set()
    flushed_contigs = set()
    processed_lines = 0
    annotation_lines = 0
    proteins_with_extra_annot = 0
    passenger_proteins = 0

    with open_file(user_gff) as input_table, \
         open_file(extra_file, "w") as output_extra, \
         open_file(full_file, "w") as output_full, \
         open_file(clean_file, "w") as output_clean:

        logger.info(f"Output files created with prefix: {output_prefix}")

        # clean and extra carry a minimal header; full preserves the user GFF's full header.
        output_clean.write("##gff-version 3\n")
        output_extra.write("##gff-version 3\n")

        # Rows of the contig currently being read, as (start, row, to_clean, to_extra).
        # Flushed sorted by start whenever the contig changes or the input ends, so contig
        # order follows the genes GFF while every entry of a contig -- genes and injected
        # mobilome alike -- comes out ascending by position.
        buffer = []
        current_contig = None

        def flush_contig(contig):
            """Write the buffered rows of one contig, ascending by start position."""
            # Stable, so a mobilome feature sharing a start with a CDS keeps the position
            # it was buffered in (first), which is the order a GFF3 parent belongs in.
            buffer.sort(key=lambda entry: entry[0])
            for _, row, to_clean, to_extra in buffer:
                output_full.write(row + "\n")
                if to_clean:
                    output_clean.write(row + "\n")
                    if to_extra:
                        output_extra.write(row + "\n")
            buffer.clear()
            if contig is not None:
                flushed_contigs.add(contig)

        for line in input_table:
            processed_lines += 1
            l_line = line.rstrip().split("\t")

            # Annotation lines have exactly 9 columns
            if len(l_line) == 9:
                annotation_lines += 1
                # Translate the contig from the internal renamed id to the original name so
                # it matches the mobilome GFF.
                if names_equiv:
                    l_line[0] = names_equiv.get(l_line[0], l_line[0])
                # Drop any trailing ';' so appending attributes below cannot yield ';;'.
                l_line[8] = normalize_attributes(l_line[8])
                contig = l_line[0]
                start = l_line[3]
                end = l_line[4]
                strand = l_line[6]
                composite_val = (contig, start, end, strand)

                if contig != current_contig:
                    flush_contig(current_contig)
                    current_contig = contig
                    if contig in flushed_contigs:
                        logger.warning(
                            f"Contig {contig} reappears after its block was written; its "
                            "rows are not contiguous in the genes GFF, which will break "
                            "tabix indexing of the outputs"
                        )
                    if contig not in used_contigs:
                        used_contigs.add(contig)
                        # Seed the buffer with this contig's mobilome entries so they sort
                        # in among the genes entries by position.
                        for mge in mobilome_annot.get(contig, []):
                            buffer.append((int(mge.split("\t")[3]), mge, True, True))

                # Append the PathoFact2 summary string (e.g. vf,mge,bgc) from the combined
                # report when this protein is present in it, as a `;pathofact2=...` attribute.
                protein_id = ""
                for attr in l_line[8].split(";"):
                    if attr.startswith("ID="):
                        protein_id = attr[3:]
                        break
                pf_suffix = (
                    f";pathofact2={summary_map[protein_id]}"
                    if protein_id in summary_map
                    else ""
                )

                has_viphog = composite_val in proteins_annot
                viphog_attr = proteins_annot[composite_val] if has_viphog else ""
                if has_viphog:
                    proteins_with_extra_annot += 1

                # Finding the mobilome proteins (passengers) in the user file. Done before the
                # full write so the mobile_element_type attribute is available for every output.
                u_prot_start = int(start)
                u_prot_end = int(end)
                u_prot_range = range(u_prot_start, u_prot_end + 1)
                u_prot_len = u_prot_end - u_prot_start
                passenger_flag = 0
                mge_loc = []

                if contig in mobilome_annot:
                    for coordinates in mges_dict[contig]:
                        mge_start = coordinates[0]
                        mge_end = coordinates[1]
                        mge_range = range(mge_start, mge_end + 1)
                        mge_label = mob_types[(contig, mge_start, mge_end)]
                        intersection = len(list(set(mge_range) & set(u_prot_range)))

                        if intersection > 0:
                            u_prot_cov = float(intersection) / float(u_prot_len)
                            if u_prot_cov > COV_THRESHOLD:
                                passenger_flag = 1
                                mge_loc.append(mge_label)

                # Shared attribute suffix: viphog (when available), mobile_element_type (for
                # passenger CDS), then pathofact2 (pf_suffix already begins with ";" or is "").
                extra_attrs = ""
                if has_viphog:
                    extra_attrs += ";" + viphog_attr
                if passenger_flag == 1:
                    extra_attrs += ";" + "mobile_element_type=" + ",".join(mge_loc)
                extra_attrs += pf_suffix

                # Rebuild column 9 rather than concatenating onto the raw line. The lstrip
                # guards a feature whose attributes were empty or a bare ';', which would
                # otherwise emit a leading separator; '.' is the GFF3 spelling for absent.
                l_line[8] = (l_line[8] + extra_attrs).lstrip(";") or "."
                out_line = "\t".join(l_line)

                # full keeps every feature, carrying the viphog, mobile_element_type and
                # pathofact2 attributes wherever they are available. clean keeps every
                # MGE-covered (passenger) CDS; extra keeps the subset of those passengers
                # that carry a functional annotation (viphog and/or pathofact2).
                if passenger_flag == 1:
                    passenger_proteins += 1
                buffer.append(
                    (
                        int(start),
                        out_line,
                        passenger_flag == 1,
                        bool(has_viphog or pf_suffix),
                    )
                )
            else:
                # Header/comment lines from the user GFF go to full only; clean and extra
                # use the minimal header written above. The FASTA block has to come after
                # every feature, so the pending contig is flushed before it starts.
                stripped = line.rstrip()
                if stripped.startswith("##FASTA"):
                    flush_contig(current_contig)
                    current_contig = None
                output_full.write(stripped + "\n")

        # Whatever contig the input ended on
        flush_contig(current_contig)
    
    # Log processing statistics
    logger.info("GFF update completed:")
    logger.info(f"  - Total lines processed: {processed_lines}")
    logger.info(f"  - Annotation lines (9 columns): {annotation_lines}")
    logger.info(f"  - Unique contigs processed: {len(used_contigs)}")
    logger.info(f"  - Proteins with extra annotations: {proteins_with_extra_annot}")
    logger.info(f"  - Passenger proteins identified: {passenger_proteins}")
    logger.info(
        f"  - Output files created: {output_prefix}{output_infix}[extra|full|clean].gff.gz"
    )
    logger.info(
        "Contig order follows the genes GFF; entries within each contig are sorted by start"
    )

def main():
    parser = argparse.ArgumentParser(
        description="This script adds extra annotations to the user GFF file. "
                   "Supports compressed (.gz) input files and generates uncompressed output files."
    )
    parser.add_argument(
        "--mobilome_gff",
        type=str,
        help="Mobilome prediction GFF file (can be compressed with .gz extension)",
        required=True,
    )
    parser.add_argument(
        "--user_gff",
        type=str,
        help="User GFF file (can be compressed with .gz extension)",
        required=False,
    )
    parser.add_argument(
        "--prefix",
        type=str,
        help="Output files prefix (outputs will be uncompressed .gff files)",
        required=True
    )
    parser.add_argument(
        "--combined_report",
        type=str,
        help="Optional PathoFact2 combined report TSV. When given, the summary_string of "
             "each protein is appended to its CDS as a `pathofact2=` attribute.",
        required=False,
    )
    parser.add_argument(
        "--user_proteins",
        action="store_true",
        help="The genes GFF comes from the user. Outputs are named "
             "`{prefix}_user_mobilome_*`; without this flag (Prodigal/tRNA baseline) they "
             "are named `{prefix}_mobilome_*`.",
    )
    parser.add_argument(
        "--contig_map",
        type=str,
        help="Optional contigID.map (renamed<TAB>original). When given, genes-GFF contig "
             "ids are translated back to their original names so they align with the "
             "mobilome GFF (used for the Prodigal/tRNA baseline).",
        required=False,
    )
    args = parser.parse_args()

    ## Calling functions
    # Storing the mobilome predictions
    (proteins_annot, mobilome_annot, mges_dict, mob_types) = mobilome_parser(
        args.mobilome_gff
    )

    # Optional per-protein summary strings from the combined report
    summary_map = parse_combined_report(args.combined_report) if args.combined_report else {}

    # Optional renamed -> original contig mapping (Prodigal/tRNA baseline only)
    names_equiv = mapping_names.names_map(args.contig_map)[0] if args.contig_map else {}

    # Adding the mobilome predictions to the genes GFF
    output_infix = "_user_mobilome_" if args.user_proteins else "_mobilome_"
    if args.user_gff:
        gff_updater(
            args.user_gff,
            args.prefix,
            proteins_annot,
            mobilome_annot,
            mges_dict,
            mob_types,
            summary_map,
            output_infix,
            names_equiv,
        )

if __name__ == "__main__":
    main()
