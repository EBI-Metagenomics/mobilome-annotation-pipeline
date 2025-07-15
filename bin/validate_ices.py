#!/usr/bin/env python3
"""
ICE Validation Script
Based on validation rules from icefinder2 single.py script
Filters ICE predictions to retain only high-quality candidates
"""
import sys
import csv
import argparse

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
        refined_start = int(row['refined_start'])
        refined_end = int(row['refined_end'])
        refined_length = int(row['refined_length'])
        num_trnas = int(row['num_trnas'])
        gc_content = float(row['gc_content'])
        
        # Extract DR coordinates (handle empty values)
        dr1_start = int(row['dr1_start']) if row['dr1_start'] and row['dr1_start'] != 'NA' and row['dr1_start'].strip() else None
        dr1_end = int(row['dr1_end']) if row['dr1_end'] and row['dr1_end'] != 'NA' and row['dr1_end'].strip() else None
        dr2_start = int(row['dr2_start']) if row['dr2_start'] and row['dr2_start'] != 'NA' and row['dr2_start'].strip() else None
        dr2_end = int(row['dr2_end']) if row['dr2_end'] and row['dr2_end'] != 'NA' and row['dr2_end'].strip() else None
        
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
        validation_messages.append(f"Invalid boundaries: start ({refined_start}) >= end ({refined_end})")
        is_valid = False
    
    # 4. Length consistency validation
    calculated_length = refined_end - refined_start + 1
    if abs(calculated_length - refined_length) > 10:  # Allow small discrepancies
        validation_messages.append(f"Length inconsistency: calculated {calculated_length} vs reported {refined_length}")
        is_valid = False
    
    # 5. GC content validation (reasonable range)
    if gc_content < 0.0 or gc_content > 1.0:
        validation_messages.append(f"Invalid GC content: {gc_content}")
        is_valid = False
    elif gc_content < 0.20 or gc_content > 0.80:
        validation_messages.append(f"Unusual GC content: {gc_content} (outside 20-80% range)")
        # Note: This is a warning, not a failure
    
    # 6. Direct Repeat validation (if present)
    if all(dr is not None for dr in [dr1_start, dr1_end, dr2_start, dr2_end]):
        # DR coordinate consistency
        if dr1_start >= dr1_end:
            validation_messages.append(f"Invalid DR1 coordinates: {dr1_start} >= {dr1_end}")
            is_valid = False
        if dr2_start >= dr2_end:
            validation_messages.append(f"Invalid DR2 coordinates: {dr2_start} >= {dr2_end}")
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

def main():
    parser = argparse.ArgumentParser(description='Validate ICE predictions based on icefinder2 rules')
    parser.add_argument('input_file', help='Input CSV/TSV file with ICE predictions')
    parser.add_argument('-o', '--output', help='Output file for validated ICEs (default: stdout)')
    parser.add_argument('-r', '--rejected', help='Output file for rejected ICEs with reasons')
    parser.add_argument('-d', '--delimiter', default='\t', help='Field delimiter (default: tab)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print validation statistics')
    
    args = parser.parse_args()
    
    # Expected headers
    expected_headers = [
        'system_id', 'contig', 'type', 'original_start', 'original_end',
        'refined_start', 'refined_end', 'refined_length', 'flank_start', 'flank_end',
        'dr1_start', 'dr1_end', 'dr2_start', 'dr2_end', 'num_trnas', 'gc_content'
    ]
    
    valid_ices = []
    rejected_ices = []
    total_count = 0
    
    try:
        with open(args.input_file, 'r', newline='') as infile:
            # Read first line to check for headers
            first_line = infile.readline().strip()
            infile.seek(0)
            
            # Check if first line contains expected headers
            first_line_fields = first_line.split(args.delimiter)
            has_header = any(header in first_line_fields for header in expected_headers[:3])  # Check first few headers
            
            if has_header:
                reader = csv.DictReader(infile, delimiter=args.delimiter)
                # Validate that all expected headers are present
                if not all(header in reader.fieldnames for header in expected_headers):
                    missing_headers = [h for h in expected_headers if h not in reader.fieldnames]
                    print(f"Error: Missing required headers: {missing_headers}", file=sys.stderr)
                    print(f"Found headers: {reader.fieldnames}", file=sys.stderr)
                    sys.exit(1)
            else:
                # No header, use expected headers as fieldnames
                reader = csv.DictReader(infile, fieldnames=expected_headers, delimiter=args.delimiter)
            
            for row in reader:
                total_count += 1
                is_valid, messages = validate_ice_entry(row)
                
                if is_valid:
                    valid_ices.append(row)
                else:
                    rejected_ices.append({
                        'row_data': row,
                        'rejection_reasons': '; '.join(messages)
                    })
                    
    except FileNotFoundError:
        print(f"Error: Input file '{args.input_file}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error reading input file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Write ONLY validated ICEs to the main output
    output_file = open(args.output, 'w', newline='') if args.output else sys.stdout
    try:
        writer = csv.DictWriter(output_file, fieldnames=expected_headers, delimiter=args.delimiter)
        writer.writeheader()
        # FIXED: Only write valid_ices, not all entries
        for valid_ice in valid_ices:
            writer.writerow(valid_ice)
    finally:
        if args.output:
            output_file.close()
    
    # Write rejected ICEs if requested
    if args.rejected and rejected_ices:
        with open(args.rejected, 'w', newline='') as rejected_file:
            fieldnames = expected_headers + ['rejection_reasons']
            writer = csv.DictWriter(rejected_file, fieldnames=fieldnames, delimiter=args.delimiter)
            writer.writeheader()
            for rejected in rejected_ices:
                row_with_reasons = rejected['row_data'].copy()
                row_with_reasons['rejection_reasons'] = rejected['rejection_reasons']
                writer.writerow(row_with_reasons)
    
    # Print statistics if verbose
    if args.verbose:
        print(f"Validation Summary:", file=sys.stderr)
        print(f"  Total ICEs processed: {total_count}", file=sys.stderr)
        print(f"  Valid ICEs: {len(valid_ices)} ({len(valid_ices)/total_count*100:.1f}%)", file=sys.stderr)
        print(f"  Rejected ICEs: {len(rejected_ices)} ({len(rejected_ices)/total_count*100:.1f}%)", file=sys.stderr)
        
        if rejected_ices:
            print(f"  Common rejection reasons:", file=sys.stderr)
            reason_counts = {}
            for rejected in rejected_ices:
                for reason in rejected['rejection_reasons'].split('; '):
                    reason_type = reason.split(':')[0]
                    reason_counts[reason_type] = reason_counts.get(reason_type, 0) + 1
            
            for reason, count in sorted(reason_counts.items(), key=lambda x: x[1], reverse=True):
                print(f"    {reason}: {count}", file=sys.stderr)

if __name__ == '__main__':
    main()
