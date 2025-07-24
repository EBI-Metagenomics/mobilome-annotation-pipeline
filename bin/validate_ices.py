#!/usr/bin/env python3
"""
ICE Validation Script
Based on validation rules from icefinder2 single.py script
Filters ICE predictions to retain only high-quality candidates
"""
import sys
import csv
import argparse
from collections import defaultdict


def validate_ice_entry(row):
    """
    Validate a single ICE entry based on icefinder2 rules
    Args:
        row (dict): Dictionary containing ICE data with expected headers
    Returns:
        tuple: (is_valid, validation_messages)
    """
    validation_messages = []
    is_valid = True

    try:
        # Extract and convert numeric values
        refined_start = int(row["refined_start"])
        refined_end = int(row["refined_end"])
        refined_length = int(row["refined_length"])
        num_trnas = int(row["num_trnas"])
        gc_content = float(row["gc_content"])

        # Extract DR coordinates (handle empty values)
        dr1_start = (
            int(row["dr1_start"])
            if row["dr1_start"]
            and row["dr1_start"] != "NA"
            and row["dr1_start"].strip()
            else None
        )
        dr1_end = (
            int(row["dr1_end"])
            if row["dr1_end"] and row["dr1_end"] != "NA" and row["dr1_end"].strip()
            else None
        )
        dr2_start = (
            int(row["dr2_start"])
            if row["dr2_start"]
            and row["dr2_start"] != "NA"
            and row["dr2_start"].strip()
            else None
        )
        dr2_end = (
            int(row["dr2_end"])
            if row["dr2_end"] and row["dr2_end"] != "NA" and row["dr2_end"].strip()
            else None
        )

    except (ValueError, TypeError) as e:
        validation_messages.append(f"Invalid numeric data: {e}")
        return False, validation_messages

    # 1. Size validation: ICEs must be between 5kb and 500kb
    if refined_length < 5000:
        validation_messages.append(f"ICE too small: {refined_length} bp < 5000 bp")
        is_valid = False
    elif refined_length > 500000:
        validation_messages.append(f"ICE too large: {refined_length} bp > 500000 bp")
        is_valid = False

    # 2. tRNA count validation: Must have at least 1 tRNA
    if num_trnas == 0:
        validation_messages.append(f"Insufficient tRNAs: {num_trnas}")
        is_valid = False

    # 3. Boundary consistency validation
    if refined_start >= refined_end:
        validation_messages.append(
            f"Invalid boundaries: start ({refined_start}) >= end ({refined_end})"
        )
        is_valid = False

    # 4. Length consistency validation
    calculated_length = refined_end - refined_start + 1
    if abs(calculated_length - refined_length) > 10:  # Allow small discrepancies
        validation_messages.append(
            f"Length inconsistency: calculated {calculated_length} vs reported {refined_length}"
        )
        is_valid = False

    # 5. GC content validation (reasonable range)
    if gc_content < 0.0 or gc_content > 1.0:
        validation_messages.append(f"Invalid GC content: {gc_content}")
        is_valid = False
    elif gc_content < 0.20 or gc_content > 0.80:
        validation_messages.append(
            f"Unusual GC content: {gc_content} (outside 20-80% range)"
        )
        # Note: This is a warning, not a failure

    # 6. Direct Repeat validation (if present)
    if all(dr is not None for dr in [dr1_start, dr1_end, dr2_start, dr2_end]):
        # DR coordinate consistency
        if dr1_start >= dr1_end:
            validation_messages.append(
                f"Invalid DR1 coordinates: {dr1_start} >= {dr1_end}"
            )
            is_valid = False
        if dr2_start >= dr2_end:
            validation_messages.append(
                f"Invalid DR2 coordinates: {dr2_start} >= {dr2_end}"
            )
            is_valid = False

        # DR length validation (typical DR length range)
        dr1_length = dr1_end - dr1_start + 1
        dr2_length = dr2_end - dr2_start + 1
        if dr1_length < 10 or dr1_length > 1000:
            validation_messages.append(f"DR1 length out of range: {dr1_length} bp")
            is_valid = False
        if dr2_length < 10 or dr2_length > 1000:
            validation_messages.append(f"DR2 length out of range: {dr2_length} bp")
            is_valid = False

        # DR position relative to ICE boundaries
        if dr1_end > refined_start or dr2_start < refined_end:
            validation_messages.append("DRs should flank the ICE region")
            is_valid = False

    return is_valid, validation_messages


def parse_macsyfinder_results(
    macsyfinder_file, completeness_threshold=0.7, verbose=False
):
    """
    Parse MacSyFinder results to extract gene composition for each system
    Args:
        macsyfinder_file: Path to all_systems.tsv file
        completeness_threshold: Minimum system completeness to include (default: 0.7)
        verbose: Print debug information
    """
    system_genes = defaultdict(list)
    system_gene_names = defaultdict(list)
    all_systems = set()
    complete_systems = set()

    try:
        with open(macsyfinder_file, "r") as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith("#") or not line.strip():
                    continue

                fields = line.strip().split("\t")
                if len(fields) >= 9:  # Make sure we have enough fields
                    sys_id = fields[5]  # system ID
                    try:
                        sys_wholeness = float(fields[8])  # system completeness
                    except ValueError:
                        if verbose and fields[8] != "sys_wholeness":  # Skip header
                            print(
                                f"Warning: Invalid wholeness value '{fields[8]}' on line {line_num}",
                                file=sys.stderr,
                            )
                        continue

                    all_systems.add(sys_id)

                    if verbose:
                        print(
                            f"Line {line_num}: sys_id={sys_id}, sys_wholeness={sys_wholeness}",
                            file=sys.stderr,
                        )

                    # Use completeness threshold instead of strict 1.0
                    if sys_wholeness >= completeness_threshold:
                        complete_systems.add(sys_id)

                        hit_id = fields[1]  # gene ID
                        gene_name = fields[2]  # gene name

                        system_genes[sys_id].append(hit_id)
                        system_gene_names[sys_id].append(gene_name)

    except FileNotFoundError:
        print(
            f"Warning: MacSyFinder file '{macsyfinder_file}' not found. Skipping gene-based filtering.",
            file=sys.stderr,
        )
        return {}, {}, set(), set()
    except Exception as e:
        print(
            f"Warning: Error parsing MacSyFinder file: {e}. Skipping gene-based filtering.",
            file=sys.stderr,
        )
        return {}, {}, set(), set()

    if verbose:
        print(
            f"Found {len(all_systems)} total systems, {len(complete_systems)} systems >= {completeness_threshold} completeness",
            file=sys.stderr,
        )
        print(f"Complete systems: {complete_systems}", file=sys.stderr)

    return dict(system_genes), dict(system_gene_names), all_systems, complete_systems


def check_ime_gene_requirements(gene_names):
    """
    Check if IME has required gene composition based on icefinder2 logic
    Args:
        gene_names: List of gene names for the IME system
    Returns:
        bool: True if IME meets gene requirements
    """
    # Integration genes mapping from icefinder2
    integrase_genes = {
        "Phage_integrase",
        "UPF0236",
        "Recombinase",
        "rve",
        "TIGR02224",
        "TIGR02249",
        "TIGR02225",
        "PB001819",
    }

    has_mobilization = False
    has_integration = False

    for gene in gene_names:
        # Check for mobilization genes
        if "Relaxase_" in gene or "T4SS_MOB" in gene:
            has_mobilization = True

        # Check for integration genes
        if gene in integrase_genes:
            has_integration = True

    # IME must have both mobilization and integration genes
    return has_mobilization and has_integration


def map_ice_to_macsyfinder_system(ice_system_id, macsyfinder_systems, verbose=False):
    """
    Map ICE system ID to MacSyFinder system ID
    ICE format: IME_contig_1_1, T4SS_typeFA_contig_1_2
    MacSyFinder format: NC_000964_IME_1, NC_000964_T4SS_typeFA_2
    """
    # Direct match
    if ice_system_id in macsyfinder_systems:
        return ice_system_id

    # Extract type and number from ICE system ID
    # IME_contig_1_1 -> IME, 1
    # T4SS_typeFA_contig_1_2 -> T4SS_typeFA, 2
    ice_parts = ice_system_id.split("_")

    if len(ice_parts) >= 4:
        if ice_parts[0] == "IME":
            ice_type = "IME"
            ice_number = ice_parts[-1]  # Last part is usually the number
        elif ice_parts[0] == "T4SS":
            ice_type = f"{ice_parts[0]}_{ice_parts[1]}"  # T4SS_typeFA
            ice_number = ice_parts[-1]  # Last part is usually the number
        else:
            ice_type = ice_parts[0]
            ice_number = ice_parts[-1]

        # Look for MacSyFinder systems with matching type and number
        for mac_sys_id in macsyfinder_systems:
            if ice_type in mac_sys_id and ice_number in mac_sys_id:
                if verbose:
                    print(
                        f"Mapped {ice_system_id} -> {mac_sys_id} (type: {ice_type}, number: {ice_number})",
                        file=sys.stderr,
                    )
                return mac_sys_id

    # Fallback: try partial matching
    for mac_sys_id in macsyfinder_systems:
        # Check if ICE system ID is contained in MacSyFinder system ID
        if ice_system_id in mac_sys_id:
            if verbose:
                print(
                    f"Mapped {ice_system_id} -> {mac_sys_id} (ICE in MacSyFinder)",
                    file=sys.stderr,
                )
            return mac_sys_id

        # Check if MacSyFinder system ID is contained in ICE system ID
        if mac_sys_id in ice_system_id:
            if verbose:
                print(
                    f"Mapped {ice_system_id} -> {mac_sys_id} (MacSyFinder in ICE)",
                    file=sys.stderr,
                )
            return mac_sys_id

    if verbose:
        print(f"No mapping found for {ice_system_id}", file=sys.stderr)
        print(
            f"Available MacSyFinder systems: {list(macsyfinder_systems)}",
            file=sys.stderr,
        )
    return None


def filter_overlapping_ices_icefinder2(
    valid_ices, macsyfinder_file, completeness_threshold=0.7, verbose=False
):
    """
    Filter overlapping ICEs based on actual icefinder2 subset logic using MacSyFinder results
    """
    # Parse MacSyFinder results
    (
        system_genes,
        system_gene_names,
        all_systems,
        complete_systems,
    ) = parse_macsyfinder_results(macsyfinder_file, completeness_threshold, verbose)

    if not system_genes:
        if verbose:
            print(
                "No MacSyFinder data found, skipping gene-based filtering",
                file=sys.stderr,
            )
        return valid_ices, []

    # Separate systems by type following icefinder2 logic
    ice_systems = {}  # ICE and T4SS systems
    ime_systems = {}  # IME systems
    aice_systems = {}  # AICE systems

    # Map ICE predictions to MacSyFinder systems
    ice_to_system = {}
    unmapped_systems = []

    if verbose:
        print(
            f"Mapping {len(valid_ices)} ICE predictions to MacSyFinder systems...",
            file=sys.stderr,
        )

    for ice in valid_ices:
        system_id = ice["system_id"]
        ice_type = ice["type"]

        # Find corresponding MacSyFinder system
        mac_sys_id = map_ice_to_macsyfinder_system(system_id, complete_systems, verbose)

        if mac_sys_id:
            ice_to_system[system_id] = mac_sys_id

            if "IME" in ice_type:
                ime_systems[system_id] = system_genes[mac_sys_id]
            elif "AICE" in ice_type:
                aice_systems[system_id] = system_genes[mac_sys_id]
            else:
                ice_systems[system_id] = system_genes[mac_sys_id]

            if verbose:
                print(
                    f"Mapped {system_id} ({ice_type}) -> {mac_sys_id} with {len(system_genes[mac_sys_id])} genes",
                    file=sys.stderr,
                )
        else:
            unmapped_systems.append(system_id)
            if verbose:
                print(
                    f"Could not map {system_id} to any MacSyFinder system",
                    file=sys.stderr,
                )

    if verbose:
        print(f"Successfully mapped: {len(ice_to_system)} systems", file=sys.stderr)
        print(f"Unmapped systems: {unmapped_systems}", file=sys.stderr)
        print(f"ICE systems: {list(ice_systems.keys())}", file=sys.stderr)
        print(f"IME systems: {list(ime_systems.keys())}", file=sys.stderr)

    # Apply icefinder2 IME gene requirements filtering
    valid_ime_systems = []
    rejected_overlaps = []

    for ime_id, ime_genes in ime_systems.items():
        if ime_id in ice_to_system:
            mac_sys_id = ice_to_system[ime_id]
            gene_names = system_gene_names.get(mac_sys_id, [])

            if verbose:
                print(
                    f"Checking IME {ime_id} gene requirements: {gene_names}",
                    file=sys.stderr,
                )

            if check_ime_gene_requirements(gene_names):
                valid_ime_systems.append(ime_id)
                if verbose:
                    print(f"IME {ime_id} passed gene requirements", file=sys.stderr)
            else:
                # Find the ICE entry to reject
                for ice in valid_ices:
                    if ice["system_id"] == ime_id:
                        rejected_overlaps.append(
                            {
                                "row_data": ice,
                                "rejection_reasons": "IME lacks required gene composition (mobilization + integration genes)",
                            }
                        )
                        if verbose:
                            print(
                                f"IME {ime_id} rejected: lacks required gene composition",
                                file=sys.stderr,
                            )
                        break

    # Apply icefinder2 subset filtering logic
    # IMEs are rejected if their gene set is a subset of any ICE's gene set
    final_ime_systems = []

    if verbose:
        print(
            f"Applying subset filtering to {len(valid_ime_systems)} valid IMEs against {len(ice_systems)} ICEs",
            file=sys.stderr,
        )

    for ime_id in valid_ime_systems:
        ime_genes = set(ime_systems[ime_id])
        is_subset = False

        if verbose:
            print(f"Checking IME {ime_id} genes: {ime_genes}", file=sys.stderr)

        for ice_id, ice_genes in ice_systems.items():
            ice_gene_set = set(ice_genes)

            if verbose:
                print(f"  Against ICE {ice_id} genes: {ice_gene_set}", file=sys.stderr)

            # Check if IME genes are a subset of ICE genes
            if ime_genes.issubset(ice_gene_set):
                is_subset = True
                if verbose:
                    print(f"  IME {ime_id} is subset of ICE {ice_id}", file=sys.stderr)

                # Find the ICE entry to reject
                for ice in valid_ices:
                    if ice["system_id"] == ime_id:
                        rejected_overlaps.append(
                            {
                                "row_data": ice,
                                "rejection_reasons": f"IME gene set is subset of {ice_id} (icefinder2 subset rule)",
                            }
                        )
                        break
                break
            else:
                if verbose:
                    print(
                        f"  IME {ime_id} is NOT subset of ICE {ice_id}", file=sys.stderr
                    )

        if not is_subset:
            final_ime_systems.append(ime_id)
            if verbose:
                print(f"IME {ime_id} retained (not a subset)", file=sys.stderr)

    # Collect final filtered ICEs
    filtered_ices = []
    rejected_system_ids = {
        rejected["row_data"]["system_id"] for rejected in rejected_overlaps
    }

    for ice in valid_ices:
        if ice["system_id"] not in rejected_system_ids:
            filtered_ices.append(ice)

    if verbose:
        print(
            f"Final result: {len(filtered_ices)} ICEs retained, {len(rejected_overlaps)} rejected",
            file=sys.stderr,
        )

    return filtered_ices, rejected_overlaps


def main():
    parser = argparse.ArgumentParser(
        description="Validate ICE predictions based on icefinder2 rules"
    )
    parser.add_argument("input_file", help="Input CSV/TSV file with ICE predictions")
    parser.add_argument(
        "-m",
        "--macsyfinder",
        help="MacSyFinder all_systems.tsv file for gene-based filtering",
    )
    parser.add_argument(
        "-o", "--output", help="Output file for validated ICEs (default: stdout)"
    )
    parser.add_argument(
        "-r", "--rejected", help="Output file for rejected ICEs with reasons"
    )
    parser.add_argument(
        "-d", "--delimiter", default="\t", help="Field delimiter (default: tab)"
    )
    parser.add_argument(
        "-v", "--verbose", action="store_true", help="Print validation statistics"
    )
    parser.add_argument(
        "--no-gene-filter", action="store_true", help="Skip gene-based subset filtering"
    )
    parser.add_argument(
        "--completeness-threshold",
        type=float,
        default=0.7,
        help="Minimum system completeness threshold (default: 0.7)",
    )

    args = parser.parse_args()

    # Expected headers
    expected_headers = [
        "system_id",
        "contig",
        "type",
        "original_start",
        "original_end",
        "refined_start",
        "refined_end",
        "refined_length",
        "flank_start",
        "flank_end",
        "dr1_start",
        "dr1_end",
        "dr2_start",
        "dr2_end",
        "num_trnas",
        "gc_content",
    ]

    valid_ices = []
    rejected_ices = []
    total_count = 0

    try:
        with open(args.input_file, "r", newline="") as infile:
            # Read first line to check for headers
            first_line = infile.readline().strip()
            infile.seek(0)

            # Check if first line contains expected headers
            first_line_fields = first_line.split(args.delimiter)
            has_header = any(
                header in first_line_fields for header in expected_headers[:3]
            )

            if has_header:
                reader = csv.DictReader(infile, delimiter=args.delimiter)
                # Validate that all expected headers are present
                if not all(header in reader.fieldnames for header in expected_headers):
                    missing_headers = [
                        h for h in expected_headers if h not in reader.fieldnames
                    ]
                    print(
                        f"Error: Missing required headers: {missing_headers}",
                        file=sys.stderr,
                    )
                    print(f"Found headers: {reader.fieldnames}", file=sys.stderr)
                    sys.exit(1)
            else:
                # No header, use expected headers as fieldnames
                reader = csv.DictReader(
                    infile, fieldnames=expected_headers, delimiter=args.delimiter
                )

            for row in reader:
                total_count += 1
                is_valid, messages = validate_ice_entry(row)

                if is_valid:
                    valid_ices.append(row)
                else:
                    rejected_ices.append(
                        {"row_data": row, "rejection_reasons": "; ".join(messages)}
                    )

    except FileNotFoundError:
        print(f"Error: Input file '{args.input_file}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading input file: {e}", file=sys.stderr)
        sys.exit(1)

    # Apply gene-based subset filtering if MacSyFinder file is provided
    if args.macsyfinder and not args.no_gene_filter and len(valid_ices) > 1:
        if args.verbose:
            print(
                f"Applying icefinder2 gene-based subset filtering to {len(valid_ices)} valid ICEs...",
                file=sys.stderr,
            )

        filtered_ices, rejected_overlaps = filter_overlapping_ices_icefinder2(
            valid_ices, args.macsyfinder, args.completeness_threshold, args.verbose
        )

        # Add gene-based rejections to main rejected list
        rejected_ices.extend(rejected_overlaps)
        valid_ices = filtered_ices

        if args.verbose:
            print(
                f"After gene-based filtering: {len(valid_ices)} ICEs retained, {len(rejected_overlaps)} rejected due to gene subset rules",
                file=sys.stderr,
            )
    elif not args.macsyfinder and args.verbose:
        print(
            "Warning: No MacSyFinder file provided. Skipping gene-based subset filtering.",
            file=sys.stderr,
        )

    # Write ONLY validated and filtered ICEs to the main output
    output_file = open(args.output, "w", newline="") if args.output else sys.stdout
    try:
        writer = csv.DictWriter(
            output_file, fieldnames=expected_headers, delimiter=args.delimiter
        )
        writer.writeheader()
        for valid_ice in valid_ices:
            writer.writerow(valid_ice)
    finally:
        if args.output:
            output_file.close()

    # Write rejected ICEs if requested
    if args.rejected and rejected_ices:
        with open(args.rejected, "w", newline="") as rejected_file:
            fieldnames = expected_headers + ["rejection_reasons"]
            writer = csv.DictWriter(
                rejected_file, fieldnames=fieldnames, delimiter=args.delimiter
            )
            writer.writeheader()
            for rejected in rejected_ices:
                row_with_reasons = rejected["row_data"].copy()
                row_with_reasons["rejection_reasons"] = rejected["rejection_reasons"]
                writer.writerow(row_with_reasons)

    # Print statistics if verbose
    if args.verbose:
        print(f"Validation Summary:", file=sys.stderr)
        print(f"  Total ICEs processed: {total_count}", file=sys.stderr)
        print(
            f"  Valid ICEs: {len(valid_ices)} ({len(valid_ices)/total_count*100:.1f}%)",
            file=sys.stderr,
        )
        print(
            f"  Rejected ICEs: {len(rejected_ices)} ({len(rejected_ices)/total_count*100:.1f}%)",
            file=sys.stderr,
        )

        if rejected_ices:
            print(f"  Rejection reasons breakdown:", file=sys.stderr)
            reason_counts = {}
            for rejected in rejected_ices:
                for reason in rejected["rejection_reasons"].split("; "):
                    reason_type = reason.split(":")[0]
                    if (
                        "subset" in reason.lower()
                        or "gene composition" in reason.lower()
                    ):
                        reason_type = "Gene-based filtering"
                    reason_counts[reason_type] = reason_counts.get(reason_type, 0) + 1

            for reason, count in sorted(
                reason_counts.items(), key=lambda x: x[1], reverse=True
            ):
                print(f"    {reason}: {count}", file=sys.stderr)


if __name__ == "__main__":
    main()
