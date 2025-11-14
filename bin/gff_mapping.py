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

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)
COV_THRESHOLD = 0.75

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

def sort_gff_file(filepath):
    """
    Sort GFF file in-place by contig and start position.
    Headers (lines starting with #) are preserved at the top.

    :param filepath: Path to GFF file to sort
    """
    if not os.path.exists(filepath):
        return

    headers = []
    entries = []

    with open_file(filepath) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith('#'):
                headers.append(line)
            elif line:
                entries.append(line)

    # Sort by contig, then start position
    entries.sort(key=lambda x: (x.split('\t')[0], int(x.split('\t')[3])))

    with open_file(filepath, 'w') as f:
        for header in headers:
            f.write(header + '\n')
        for entry in entries:
            f.write(entry + '\n')

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
    
    extra_annot = [
        "viphog",
        "viphog_taxonomy",
    ]
    
    with open_file(mobilome_clean) as input_table:
        logger.info(f"Successfully opened mobilome file: {mobilome_clean}")
        
        for line in input_table:
            l_line = line.rstrip().split("\t")
            
            # Annotation lines have exactly 9 columns
            if len(l_line) == 9:
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
                        mobilome_annot[contig].append(line.rstrip())
                        mges_dict[contig].append(coordinates)
                    else:
                        mobilome_annot[contig] = [line.rstrip()]
                        mges_dict[contig] = [coordinates]
                else:
                    str_composite_key = (contig, str(start), str(end), strand)
                    attrib = l_line[8]
                    extra_list = []
                    
                    for attr in attrib.split(";"):
                        if "=" in attr:  # Ensure attr has the expected format
                            att_key = attr.split("=")[0]
                            if att_key in extra_annot:
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

def gff_updater(
    user_gff, output_prefix, proteins_annot, mobilome_annot, mges_dict, mob_types
):
    """Adding the mobilome predictions to the user file (handles compressed input/output)."""
    
    # Check if input file exists
    if not os.path.exists(user_gff):
        logger.error(f"User GFF file not found: {user_gff}")
        sys.exit(1)
    
    # Check if input file has content
    if is_file_empty(user_gff):
        logger.warning(f"User GFF file is empty: {user_gff}")
        # Still create empty output files
        output_files = [
            f"{output_prefix}_user_mobilome_extra.gff",
            f"{output_prefix}_user_mobilome_full.gff",
            f"{output_prefix}_user_mobilome_clean.gff"
        ]
        for output_file in output_files:
            with open_file(output_file, 'w') as f:
                pass  # Create empty file
        logger.info(f"Created empty output files with prefix: {output_prefix}")
        return
    
    logger.info(f"Starting GFF update process with file: {user_gff}")
    
    used_contigs = []
    processed_lines = 0
    annotation_lines = 0
    proteins_with_extra_annot = 0
    passenger_proteins = 0
    
    with open_file(user_gff) as input_table, \
         open_file(f"{output_prefix}_user_mobilome_extra.gff", "w") as output_extra, \
         open_file(f"{output_prefix}_user_mobilome_full.gff", "w") as output_full, \
         open_file(f"{output_prefix}_user_mobilome_clean.gff", "w") as output_clean:

        logger.info(f"Output files created with prefix: {output_prefix}")
        
        for line in input_table:
            processed_lines += 1
            l_line = line.rstrip().split("\t")
            
            # Annotation lines have exactly 9 columns
            if len(l_line) == 9:
                annotation_lines += 1
                contig = l_line[0]
                start = l_line[3]
                end = l_line[4]
                strand = l_line[6]
                composite_val = (contig, start, end, strand)
                
                if contig not in used_contigs:
                    used_contigs.append(contig)
                    
                    # Writing the mobilome entries in every output
                    if contig in mobilome_annot:
                        for mge in mobilome_annot[contig]:
                            output_clean.write(mge + "\n")
                            output_extra.write(mge + "\n")
                            output_full.write(mge + "\n")
                
                # Writing to extra and full outputs
                if composite_val in proteins_annot:
                    proteins_with_extra_annot += 1
                    extra_annot = proteins_annot[composite_val]
                    output_extra.write(line.rstrip() + ";" + extra_annot + "\n")
                    output_full.write(line.rstrip() + ";" + extra_annot + "\n")
                else:
                    output_full.write(line.rstrip() + "\n")
                
                # Finding mobilome proteins in the user file and writing to clean output
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
                
                if passenger_flag == 1:
                    passenger_proteins += 1
                    mge_loc = "mge_location=" + ",".join(mge_loc)
                    if composite_val in proteins_annot:
                        extra_annot = proteins_annot[composite_val]
                        output_clean.write(
                            line.rstrip() + ";" + extra_annot + ";" + mge_loc + "\n"
                        )
                    else:
                        output_clean.write(line.rstrip() + ";" + mge_loc + "\n")
                else:
                    output_clean.write(line.rstrip() + "\n")
            else:
                # Non-annotation lines (headers, comments, etc.)
                output_extra.write(line.rstrip() + "\n")
                output_full.write(line.rstrip() + "\n")
    
    # Log processing statistics
    logger.info("GFF update completed:")
    logger.info(f"  - Total lines processed: {processed_lines}")
    logger.info(f"  - Annotation lines (9 columns): {annotation_lines}")
    logger.info(f"  - Unique contigs processed: {len(used_contigs)}")
    logger.info(f"  - Proteins with extra annotations: {proteins_with_extra_annot}")
    logger.info(f"  - Passenger proteins identified: {passenger_proteins}")
    logger.info(
        f"  - Output files created: {output_prefix}_user_mobilome_[extra|full|clean].gff"
    )

    # Sort output files
    sort_gff_file(f"{output_prefix}_user_mobilome_extra.gff")
    sort_gff_file(f"{output_prefix}_user_mobilome_full.gff")
    sort_gff_file(f"{output_prefix}_user_mobilome_clean.gff")
    logger.info("Output files sorted by contig and position")

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
    args = parser.parse_args()

    ## Calling functions
    # Storing the mobilome predictions
    (proteins_annot, mobilome_annot, mges_dict, mob_types) = mobilome_parser(
        args.mobilome_gff
    )

    # Adding the mobilome predictions to the user file
    if args.user_gff:
        gff_updater(
            args.user_gff,
            args.prefix,
            proteins_annot,
            mobilome_annot,
            mges_dict,
            mob_types,
        )

if __name__ == "__main__":
    main()
