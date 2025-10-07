#!/usr/bin/env python3
"""
Python script to integrate Prodigal CDS predictions with Aragorn tRNA predictions
using standard bioinformatics overlap detection and resolution algorithms.

This implementation uses common practices for feature integration in genome annotation:
- Non-overlapping features are retained
- Overlapping CDS features are filtered when they conflict with RNA features
- Features are sorted by genomic coordinates for consistent output

This is compatible with Prokka
"""

import re
import gzip
import sys
import argparse
from Bio import SeqIO
from collections import defaultdict
from dataclasses import dataclass
from typing import Dict, List, Optional, Tuple


@dataclass
class Feature:
    """Represents a genomic feature (CDS or tRNA)"""

    contig_id: str
    feature_type: str  # 'CDS' or 'tRNA' or 'tmRNA'
    start: int
    end: int
    strand: str  # '+' or '-'
    source: str
    attributes: Dict[str, str]
    original_id: Optional[str] = None  # For tracking original prodigal IDs


def parse_aragorn_output(file_path: str) -> Dict[str, List[Feature]]:
    """
    Parse aragorn output following prokka rules.
    Returns dictionary with contig names as keys and tRNA features as values.
    """
    trnas_per_contig = defaultdict(list)
    current_contig = None

    with open(file_path, "r") as f:
        for line_num, line in enumerate(f, 1):
            line = line.strip()

            # Check for contig header (starts with >)
            if line.startswith(">"):
                current_contig = line[1:].split()[0]
                continue

            # Skip empty lines and gene count lines
            if not line or "genes found" in line:
                continue

            # Split line by whitespace (prokka rule)
            fields = line.split()

            # Apply prokka validation rules
            if len(fields) != 5:
                continue

            if not fields[0].isdigit():
                continue

            gene_num = fields[0]
            trna_type = fields[1]  # e.g., "tRNA-Glu"
            coordinates = fields[2]  # e.g., "[17805,17880]" or "c[405250,405325]"
            length = fields[3]
            anticodon = fields[4]  # e.g., "(ttc)"

            # Skip pseudo/wacky genes (prokka rule)
            if "?" in trna_type:
                print(f"Line {line_num}: Skipping pseudo/wacky gene: {trna_type}")
                continue

            # Parse coordinates using prokka's regex pattern
            coord_match = re.match(r"(c)?\[-?(\d+),(\d+)\]", coordinates)
            if not coord_match:
                print(f"Line {line_num}: Could not parse coordinates: {coordinates}")
                continue

            is_reverse = coord_match.group(1) == "c"
            start = int(coord_match.group(2))
            end = int(coord_match.group(3))

            # Apply prokka validation rules
            if start > end:
                print(
                    f"Line {line_num}: Skipping tRNA with start > end: {start} > {end}"
                )
                continue

            # Skip if too large (prokka rule: >500bp)
            if abs(end - start) > 500:
                print(
                    f"Line {line_num}: Skipping tRNA/tmRNA too big (>500bp): {abs(end-start)}bp"
                )
                continue

            # Determine strand and feature type
            strand = "-" if is_reverse else "+"

            # Handle tmRNA special case (prokka rule)
            if trna_type.startswith("tmRNA"):
                feature_type = "tmRNA"
                product = "tmRNA"
                gene_name = "ssrA"
            else:
                feature_type = "tRNA"
                product = f"{trna_type}"
                gene_name = None

            # Create feature following prokka structure
            attributes = {
                "product": product,
            }
            if gene_name:
                attributes["gene"] = gene_name

            feature = Feature(
                contig_id=current_contig,
                feature_type=feature_type,
                start=start,
                end=end,
                strand=strand,
                source="Aragorn",
                attributes=attributes,
            )

            if current_contig:
                trnas_per_contig[current_contig].append(feature)

    return dict(trnas_per_contig)


def parse_prodigal_output(file_path: str) -> Dict[str, List[Feature]]:
    """
    Parse prodigal output following prokka rules.
    Returns dictionary with contig names as keys and CDS features as values.
    """
    cds_per_contig = defaultdict(list)

    with gzip.open(file_path, "rt") as f:
        for line_num, line in enumerate(f, 1):
            l_line = line.rstrip().split("\t")
            # Annotation lines have exactly 9 columns
            if len(l_line) == 9:
                (
                    contig,
                    seq_source,
                    seq_type,
                    start,
                    end,
                    score,
                    strand_char,
                    phase,
                    attr,
                ) = line.rstrip().split("\t")

                original_id = (
                    contig + "_" + attr.split(";")[0].replace("ID=", "").split("_")[-1]
                )
                strand = "+" if strand_char == "+" else "-"

                attributes = {"ori_id": "" + original_id}

                feature = Feature(
                    contig_id=contig,
                    feature_type="CDS",
                    start=int(start),
                    end=int(end),
                    strand=strand,
                    source="Prodigal",
                    attributes=attributes,
                    original_id=f"{original_id}_{start}_{end}_{strand_char}",
                )

                cds_per_contig[contig].append(feature)

    return dict(cds_per_contig)


def features_overlap(feature1: Feature, feature2: Feature) -> bool:
    """
    Check if two features overlap (prokka rule).
    """
    if feature1.contig_id != feature2.contig_id:
        return False

    # Check for overlap
    return not (feature1.end < feature2.start or feature2.end < feature1.start)


def integrate_features(
    cds_dict: Dict[str, List[Feature]],
    trna_dict: Dict[str, List[Feature]],
    allow_cds_rna_overlap: bool = False,
) -> Tuple[Dict[str, List[Feature]], List[Dict]]:
    """
    Integrate CDS and tRNA features following prokka rules.
    Returns integrated features and list of excluded CDS.
    """
    integrated_features = defaultdict(list)
    excluded_cds = []

    # Get all contigs
    all_contigs = set(cds_dict.keys()) | set(trna_dict.keys())

    for contig_id in all_contigs:
        contig_features = []

        # Add tRNA features first (prokka processes these before CDS)
        if contig_id in trna_dict:
            contig_features.extend(trna_dict[contig_id])

        # Add CDS features, checking for overlaps with RNA features
        if contig_id in cds_dict:
            rna_features = trna_dict.get(contig_id, [])

            for cds in cds_dict[contig_id]:
                # Check for overlaps with RNA features (prokka rule)
                overlapping_rna = None
                for rna in rna_features:
                    if features_overlap(cds, rna):
                        overlapping_rna = rna
                        break

                # Exclude CDS if it overlaps with RNA (unless allowed)
                if overlapping_rna and not allow_cds_rna_overlap:
                    excluded_cds.append(
                        {
                            "contig": contig_id,
                            "cds_coords": f"{cds.start}..{cds.end}",
                            "strand": cds.strand,
                            "overlapping_rna_type": overlapping_rna.feature_type,
                            "original_id": cds.original_id,
                        }
                    )
                    print(
                        f"Excluding CDS which overlaps existing RNA ({overlapping_rna.feature_type}) "
                        f"at {contig_id}:{cds.start}..{cds.end} on {cds.strand} strand"
                    )
                else:
                    contig_features.append(cds)

        # Sort features by start position (prokka rule)
        contig_features.sort(key=lambda f: f.start)
        integrated_features[contig_id] = contig_features

    return dict(integrated_features), excluded_cds


def assign_locus_tags(
    integrated_features: Dict[str, List[Feature]],
    locus_tag_prefix,
    increment: int = 1,
) -> Tuple[Dict[str, List[Feature]], Dict[str, str]]:
    """
    Assign locus tags following prokka rules.
    Returns features with assigned IDs and mapping of new IDs to original prodigal IDs.
    """
    id_mapping = {}
    feature_counter = 0

    # Process all contigs in order
    for contig_id in sorted(integrated_features.keys()):
        features = integrated_features[contig_id]

        for feature in features:
            # Only assign locus tags to CDS and RNA features (prokka rule)
            if feature.feature_type in ["CDS", "tRNA", "tmRNA", "rRNA", "misc_RNA"]:
                feature_counter += 1
                locus_tag = f"{locus_tag_prefix}_{feature_counter * increment:05d}"

                # Add locus_tag and ID to attributes
                feature.attributes["ID"] = locus_tag

                # Track mapping for original prodigal IDs
                if feature.original_id:
                    id_mapping[locus_tag] = feature.original_id



    return integrated_features, id_mapping


def write_gff3(
    integrated_features: Dict[str, List[Feature]],
    prefix: str,
    contig_lengths: Optional[Dict[str, int]] = None,
):
    """
    Write integrated features to GFF3 format following prokka style.
    """
    output_file = prefix + '_merged.gff'
    with open(output_file, "w") as f:
        f.write("##gff-version 3\n")

        # Write sequence regions if lengths provided
        if contig_lengths:
            for contig_id in sorted(integrated_features.keys()):
                if contig_id in contig_lengths:
                    f.write(
                        f"##sequence-region {contig_id} 1 {contig_lengths[contig_id]}\n"
                    )

        # Write features
        for contig_id in sorted(integrated_features.keys()):
            features = integrated_features[contig_id]

            for feature in features:
                # Format attributes
                attr_strings = []

                attr_strings.append(f"ID={feature.attributes['ID']}")

                # Add all other attributes (excluding ID to avoid duplication)
                for key, value in feature.attributes.items():
                    if key != "ID" and key != "ori_id":  # Skip ID since we already added it first
                        attr_strings.append(f"{key}={value}")

                attributes_str = ";".join(attr_strings)

                # Write GFF3 line
                f.write(
                    f"{feature.contig_id}\t{feature.source}\t{feature.feature_type}\t"
                    f"{feature.start}\t{feature.end}\t.\t{feature.strand}\t0\t{attributes_str}\n"
                )



def write_faa(prodigal_faa, id_mapping, prefix):
    """
    Renaming proteins fasta file to be in line with merged GFF file
    """
    reversed_mapping = {v: k for k, v in id_mapping.items()}
    output_file = prefix + '_renamed.faa'
    with open(output_file, "w") as output:
        with gzip.open(prodigal_faa, "rt") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                # Removing '*' from stop codon
                clean_seq = str(record.seq).upper().replace('*','') 

                # Example of prodigal header: 
                #contig_1_5111 # 5328136 # 5328669 # -1 # ID=1_5111;partial=00;start_type=ATG;rbs_motif=GGAGG;rbs_spacer=5-10bp;gc_cont=0.515
                faa_desc_list = (str(record.description)).split(' # ')
                protein_id = faa_desc_list[0]
                start = faa_desc_list[1]
                end = faa_desc_list[2]
                strand_char = faa_desc_list[3]
                strand = "-" if strand_char == "-1" else "+"

                long_id = '_'.join([
                    protein_id,
                    start,
                    end,
                    strand
                ])

                if long_id in reversed_mapping:
                    output.write('>' + reversed_mapping[long_id] + "\n")
                    output.write(clean_seq + "\n" )


def write_mapping_file(id_mapping: Dict[str, str], excluded_cds: List[Dict]):
    """
    Write mapping file to track original prodigal CDS names.
    """
    with open("names.map", "w") as f:
        f.write("# Mapping of prokka locus_tags to original prodigal IDs\n")
        f.write("locus_tag\toriginal_prodigal_id\tstatus\n")

        # Write included CDS mappings
        for locus_tag, original_id in sorted(id_mapping.items()):
            f.write(f"{locus_tag}\t{original_id}\tincluded\n")

        # Write excluded CDS
        if excluded_cds:
            f.write("\n# Excluded CDS (overlapped with RNA features)\n")
            f.write(
                "# contig\tcoordinates\tstrand\toverlapping_rna_type\toriginal_prodigal_id\n"
            )
            for excluded in excluded_cds:
                f.write(
                    f"# {excluded['contig']}\t{excluded['cds_coords']}\t{excluded['strand']}\t"
                    f"{excluded['overlapping_rna_type']}\t{excluded['original_id']}\n"
                )


def main():
    parser = argparse.ArgumentParser(
        description="A merger that itegrates prodigal CDSs prediction and tRNAs predicted by aragorn. Merging and QC follows PROKKA rules"
    )
    parser.add_argument(
        "--prodigal_gff",
        type=str,
        help="Result of Prodigal in gff format",
        required=True,
    )
    parser.add_argument(
        "--aragorn",
        type=str,
        help="Result of Aragorn in tabular format",
        required=True,
    )
    parser.add_argument(
        "--prodigal_faa",
        type=str,
        help="Result of Prodigal proteins in fasta format (faa file)",
        required=True,
    )
    parser.add_argument(
        "--prefix",
        type=str,
        help="Output and proteins prefix to be used",
        required=True,
    )
    args = parser.parse_args()

    # Parse hardcoded PROKKA options
    allow_overlap = False
    increment = 1

    print("=== PROKKA-STYLE INTEGRATION ===")
    print(f"Prodigal gff file: {args.prodigal_gff}")
    print(f"Prodigal faa file: {args.prodigal_faa}")
    print(f"Aragorn tbl file: {args.aragorn}")
    print(f"Locus tag and outputs prefix: {args.prefix}")
    print(f"Allow CDS-RNA overlap: {allow_overlap}")
    print()

    # Parse input files
    print("Parsing Aragorn output...")
    trna_features = parse_aragorn_output(args.aragorn)
    total_trnas = sum(len(features) for features in trna_features.values())
    print(
        f"Found {total_trnas} tRNA/tmRNA features across {len(trna_features)} contigs"
    )

    print("\nParsing Prodigal output...")
    cds_features = parse_prodigal_output(args.prodigal_gff)
    total_cds = sum(len(features) for features in cds_features.values())
    print(f"Found {total_cds} CDS features across {len(cds_features)} contigs")

    # Integrate features
    print("\nIntegrating features (checking for overlaps)...")
    integrated_features, excluded_cds = integrate_features(
        cds_features, trna_features, allow_overlap
    )

    total_integrated = sum(
        len(features) for features in integrated_features.values()
    )
    print(f"Integrated {total_integrated} features")
    print(f"Excluded {len(excluded_cds)} CDS due to RNA overlaps")

    # Assign locus tags
    print(f"\nAssigning locus tags with prefix '{args.prefix}'...")
    integrated_features, id_mapping = assign_locus_tags(
        integrated_features, args.prefix, increment
    )

    # Write outputs
    print(f"\nWriting GFF3 output to: {args.prefix}_merged.gff")
    write_gff3(integrated_features, args.prefix)

    print(f"\nWriting renamed faa output to: {args.prefix}_renamed.faa")
    write_faa(args.prodigal_faa, id_mapping, args.prefix)

    print(f"Reporting ID mapping:")
    write_mapping_file(id_mapping, excluded_cds)

    # Summary
    print("\n=== INTEGRATION SUMMARY ===")
    print(f"Total features integrated: {total_integrated}")
    print(f"  - tRNA/tmRNA: {total_trnas}")
    print(f"  - CDS: {total_cds - len(excluded_cds)}")
    print(f"Excluded CDS (RNA overlaps): {len(excluded_cds)}")
    print(f"Locus tags assigned: {len(id_mapping)}")

    # Per-contig summary
    print(f"\nPer-contig breakdown:")
    for contig_id in sorted(integrated_features.keys()):
        features = integrated_features[contig_id]
        cds_count = sum(1 for f in features if f.feature_type == "CDS")
        rna_count = sum(
            1 for f in features if f.feature_type in ["tRNA", "tmRNA", "rRNA"]
        )
        print(
            f"  {contig_id}: {len(features)} features ({cds_count} CDS, {rna_count} RNA)"
        )

    print(f"\nIntegration completed successfully!")


if __name__ == "__main__":
    main()
