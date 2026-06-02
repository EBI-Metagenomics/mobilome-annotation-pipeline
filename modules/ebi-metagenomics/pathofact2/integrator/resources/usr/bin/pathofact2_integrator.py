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
import fileinput
import logging
import os

# Set up logger
logger = logging.getLogger(__name__)


def _pred_support_has_data_rows(path: str) -> bool:
    """
    Return True if pred_support TSV has at least one data row (beyond the header).

    Expected header:
      sequence_id, detection_method, support_value_type, support_value, vfdb_hit
    """
    with fileinput.hook_compressed(path, "r", encoding="utf-8", errors="ignore") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)

        # Empty file
        if header is None:
            return False

        # Header-only file => no data rows
        for row in reader:
            # Skip blank lines
            if not row:
                continue
            # Skip rows that are effectively empty (e.g. ["", "", ...])
            if all(cell.strip() == "" for cell in row):
                continue
            return True

    return False


def _has_gff_record_9cols(path: str, max_lines: int = 200) -> bool:
    """
    Return True if at least one non-comment, non-empty line with exactly 9 tab-separated
    columns exists within the first `max_lines` lines.
    """
    with fileinput.hook_compressed(path, "r", encoding="utf-8", errors="ignore") as fh:
        for _, line in zip(range(max_lines), fh):
            s = line.strip()
            if not s or s.startswith("#"):
                continue
            if len(s.split("\t")) == 9:
                return True
    return False


def validate_inputs(to_validate_inputs: dict[str, str | None]) -> dict[str, str | None]:
    """
    Validation policy:

      - annotation: must exist (content may be empty/header-only)
      - pred_support: must contain at least one data row (header-only is invalid)
      - gff: validated only if pred_support has data; then requires >=1 9-col record
    """
    non_valid_inputs: dict[str, str | None] = {}

    # ---- Basic existence checks ----
    for name, path in to_validate_inputs.items():
        if not path or not path.strip():
            non_valid_inputs[name] = "No input provided"
            continue

        if not os.path.exists(path):
            logger.warning("File not found for '%s' → %s", name, path)
            non_valid_inputs[name] = "File not found"
            continue

    # If any required path is missing, stop early
    if non_valid_inputs:
        return non_valid_inputs

    # ---- pred_support must have at least one predicted protein ----
    pred_support_path = to_validate_inputs["pred_support"]
    if pred_support_path is None:
        non_valid_inputs["pred_support"] = "No input provided"
        return non_valid_inputs

    if not _pred_support_has_data_rows(pred_support_path):
        logger.info("pred_support is empty or header-only → %s", pred_support_path)
        non_valid_inputs["pred_support"] = (
            "No predicted proteins (empty or header-only)"
        )
        return non_valid_inputs

    # ---- Validate GFF only if pred_support has predicted proteins ----
    gff_path = to_validate_inputs["gff"]
    if gff_path is None:
        non_valid_inputs["gff"] = "No input provided"
        return non_valid_inputs

    if not _has_gff_record_9cols(gff_path):
        logger.info("GFF has no valid 9-column feature records → %s", gff_path)
        non_valid_inputs["gff"] = "No valid GFF records (no 9-column features found)"

    # ---- annotation: existence already checked; no content validation ----
    return non_valid_inputs


def parse_pathofact_support(input_file: str) -> dict[str, str]:
    """
    Parse Pathofact2 prediction result in tsv format

    Args:
        input_file: Path to the tsv input file
                    This file has 5 columns:
                    sequence_id     detection_method        support_value_type      support_value   vfdb_hit
    Returns:
        pathofact_attributes defaultdict with protein id as key and list of annotations as values
    """
    pathofact_attributes: dict[str, str] = {}
    with fileinput.hook_compressed(input_file, "r", encoding="utf-8") as input_table:
        # Skipping the header line
        next(input_table)
        for line in input_table:
            line_l: list[str] = line.rstrip().split("\t")
            protein_id: str = line_l[0]
            method: str = line_l[1]
            support_value: str = line_l[3]

            if method == "blastp":
                vfdb_hit: str = line_l[4]
                annotation_line = ";".join(
                    [
                        ("vfdb=" + vfdb_hit),
                        ("blastp_eval=" + support_value),
                    ]
                )
            else:
                annotation_line = method + "_prob=" + support_value

            if protein_id in pathofact_attributes:
                stored_line = pathofact_attributes[protein_id]
                new_line = stored_line + ";" + annotation_line
                pathofact_attributes[protein_id] = new_line
            else:
                pathofact_attributes[protein_id] = annotation_line

    return pathofact_attributes


def _add_annotation_to_protein(
    protein_id: str,
    signature_acc: str,
    signature_desc: str,
    pathofact_attributes: dict[str, str],
    proteins_annotation: dict[str, str],
) -> None:
    """
    Add CDD annotation to a protein's annotation dictionary.

    This helper function updates the proteins_annotation dictionary by either:
    - Appending to existing CDD annotations (if protein already has some)
    - Creating new CDD annotation entry with pathofact attributes

    Args:
        protein_id: Protein identifier
        signature_acc: Signature accession (e.g., CDD accession)
        signature_desc: Signature description (e.g., CDD short name)
        pathofact_attributes: Dictionary with protein IDs and pathofact annotations
        proteins_annotation: Dictionary to update with new annotations (modified in place)
    """
    if protein_id in pathofact_attributes:
        current_line = signature_acc + ":" + signature_desc
        if protein_id in proteins_annotation:
            stored_line = proteins_annotation[protein_id]
            updated_line = stored_line + "," + current_line
            proteins_annotation[protein_id] = updated_line
        else:
            pathofact_values = pathofact_attributes[protein_id]
            proteins_annotation[protein_id] = pathofact_values + ";cdd=" + current_line


def parse_cdd(pathofact_attributes: dict[str, str], input_file: str) -> dict[str, str]:
    """
    Parse protein annotation file generated by local-cd-search

    Args:
        pathofact_attributes: Dictionary with protein IDs and pathofact annotations
        input_file: Path to the tsv input file
                    This file has 11 columns:
                    query hit_type pssm_id from to evalue bitscore accession short_name incomplete superfamily_pssm_id

    Returns:
        proteins_annotation defaultdict with protein id as key and cdd accession and shortname concatenated in a single string
    """
    proteins_annotation: dict[str, str] = {}
    with fileinput.hook_compressed(input_file, "r", encoding="utf-8") as input_table:
        # Skipping the header line
        next(input_table)
        for line in input_table:
            line_l: list[str] = line.rstrip().split("\t")
            protein_id: str = line_l[0]
            cdd_acc: str = line_l[7]
            cdd_short: str = line_l[8]

            _add_annotation_to_protein(
                protein_id,
                cdd_acc,
                cdd_short,
                pathofact_attributes,
                proteins_annotation,
            )

    return proteins_annotation


def parse_ips(pathofact_attributes: dict[str, str], input_file: str) -> dict[str, str]:
    """
     Parse the Interproscan protein annotation file provided by the user to extract CDD annotations

     Args:
         pathofact_attributes: Dictionary with protein IDs and pathofact annotations
         input_file: Path to the tsv input file
                     This file has 15 columns with no header:
    1. protein_id - Protein Accession (e.g. P51587)
    2. protein_md5 - Sequence MD5 digest (e.g. 14086411a2cdf1c4cba63020e1622579)
    3. protein_length - Sequence Length (e.g. 3418)
    4. analysis - Analysis (e.g. Pfam / PRINTS / Gene3D)
    5. signature_id - Signature Accession (e.g. PF09103 / G3DSA:2.40.50.140)
    6. signature_desc - Signature Description (e.g. BRCA2 repeat profile)
    7. start - Start location
    8. end - Stop location
    9. score - is the e-value (or score) of the match reported by member database method (e.g. 3.1E-52)
    10. status - is the status of the match (T: true)
    11. date - is the date of the run
    12. interpro_id - (InterPro annotations - accession (e.g. IPR002093) - optional column; only displayed if -iprlookup option is switched on)
    13. interpro_desc - (InterPro annotations - description (e.g. BRCA2 repeat) - optional column; only displayed if -iprlookup option is switched on)
    14. go_terms - (GO annotations (e.g. GO:0005515) - optional column; only displayed if –goterms option is switched on)
    15. pathway_terms - (Pathways annotations (e.g. REACT_71) - optional column; only displayed if –pathways option is switched on)

     Returns:
         proteins_annotation defaultdict with protein id as key and cdd accession and shortname concatenated in a single string
    """
    proteins_annotation: dict[str, str] = {}
    with fileinput.hook_compressed(input_file, "r", encoding="utf-8") as input_table:
        for line in input_table:
            line_l: list[str] = line.rstrip().split("\t")
            protein_id: str = line_l[0]
            analysis: str = line_l[3]
            signature_acc: str = line_l[4]
            signature_desc: str = line_l[5]

            if analysis == "CDD":
                _add_annotation_to_protein(
                    protein_id,
                    signature_acc,
                    signature_desc,
                    pathofact_attributes,
                    proteins_annotation,
                )

    return proteins_annotation


def parse_gff(
    cds_gff: str, output_file: str, proteins_annotation: dict[str, str]
) -> None:
    """
    Parse GFF file and integrate Pathofact2 annotations into output GFF.

    Args:
        cds_gff: Path to input GFF file containing CDS coordinates
        output_file: Path to output GFF file
        proteins_annotation: Dictionary mapping protein IDs to attribute lists
    """
    with (
        fileinput.hook_compressed(cds_gff, "r", encoding="utf-8") as input_table,
        open(output_file, "w") as output_gff,
    ):
        output_gff.write("##gff-version 3\n")
        for line in input_table:
            line = line.rstrip()
            line_l: list[str] = line.split("\t")
            # Annotation lines have exactly 9 columns
            if len(line_l) == 9:
                (
                    contig,
                    seq_source,
                    seq_type,
                    start,
                    end,
                    score,
                    strand,
                    phase,
                    attr,
                ) = line.split("\t")
                feature_id = attr.split(";")[0].replace("ID=", "")
                if feature_id in proteins_annotation:
                    new_attribute: str = (
                        "ID=" + feature_id + ";" + proteins_annotation[feature_id]
                    )
                    line_l.pop(-1)
                    line_l.append(new_attribute)
                    to_print: str = "\t".join(line_l)
                    output_gff.write(to_print + "\n")


def main() -> None:
    """
    Main function to orchestrate Pathofact2 results
    """
    parser: argparse.ArgumentParser = argparse.ArgumentParser(
        description="Integration of toxins and virulence factors predicted and annotated using Pathofact2 into a single gff3 file"
    )
    parser.add_argument(
        "-g",
        "--gff",
        dest="cds_gff",
        help="GFF file containing the coordinates of the CDSs used for Pathofact2 annotation",
        required=True,
        default=None,
    )
    parser.add_argument(
        "-a",
        "--annot",
        dest="proteins_annot",
        help="Tab-delimited file of functional annotation. It can be the interproscan result for all the proteins provided by the user or the CDD annotation for relevant proteins generated as part of the pathofact subworkflow",
        required=True,
        default=None,
    )
    parser.add_argument(
        "-s",
        "--support",
        dest="pred_support",
        help="Tab-delimited file corresponding to the Pathofact2 result containing the blastp hits evalues and the probability of the predictions",
        required=True,
        default=None,
    )
    parser.add_argument(
        "-t",
        "--type",
        dest="annot_type",
        choices=["cdd", "ips"],
        help="Annotation type to be parsed. Valid strings are 'cdd' for results generated in the subworkflow by local-cd-search module or 'ips' for user provided interproscan annotation",
        required=True,
        default=None,
    )
    parser.add_argument(
        "-o",
        "--output",
        dest="output",
        help="Name of the output file",
        required=False,
        default="integrated_result.gff",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        dest="verbose",
        help="Enable verbose logging (DEBUG level)",
        action="store_true",
        default=False,
    )
    args: argparse.Namespace = parser.parse_args()

    # Configure logging
    log_level = logging.DEBUG if args.verbose else logging.INFO
    logging.basicConfig(
        level=log_level,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
        handlers=[
            logging.StreamHandler(),  # Console output
        ],
    )

    to_validate_inputs: dict[str, str | None] = {
        "gff": args.cds_gff,
        "annotation": args.proteins_annot,
        "pred_support": args.pred_support,
    }

    non_valid_inputs: dict[str, str | None] = validate_inputs(to_validate_inputs)
    if len(non_valid_inputs) == 0:
        logger.info("Provided inputs are valid")
        logger.info(f"Parsing Pathofact2 predictions support file: {args.pred_support}")
        pathofact_attributes: dict[str, str] = parse_pathofact_support(
            args.pred_support
        )
        logger.debug(
            f"Pathofact2 predicted {len(pathofact_attributes)} VF and toxin proteins"
        )

        if args.annot_type == "cdd":
            logger.info(f"Parsing CDD proteins annotation file: {args.proteins_annot}")
            proteins_annotation: dict[str, str] = parse_cdd(
                pathofact_attributes, args.proteins_annot
            )

        elif args.annot_type == "ips":
            logger.info(f"Parsing IPS proteins annotation file: {args.proteins_annot}")
            proteins_annotation: dict[str, str] = parse_ips(
                pathofact_attributes, args.proteins_annot
            )

        logger.info("Parsing gff file and writing output file")
        parse_gff(args.cds_gff, args.output, proteins_annotation)
        logger.info(f"Output written to: {args.output}")

    else:
        logger.info(f"{len(non_valid_inputs)} invalid input files detected")
        for key, value in non_valid_inputs.items():
            logger.info(f"{key} file is invalid due to: {value}")
        logger.info("The output file is not generated")


if __name__ == "__main__":
    main()
