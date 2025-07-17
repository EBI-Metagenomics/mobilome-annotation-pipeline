#!/usr/bin/env python3
"""
ICE Boundaries Refinement Script
Adapted from icefinder2 single.py methodology for use with direct tRNA predictions on ICE regions
"""

import sys
import csv
import argparse
from collections import defaultdict

def parse_macsyfinder_results(macsyfinder_file):
    """Parse MacSyFinder results file to extract ICE predictions with coordinates"""
    ice_predictions = {}
    
    try:
        with open(macsyfinder_file, 'r') as f:
            # Skip header if present
            first_line = f.readline().strip()
            if not first_line.startswith('system_id'):
                f.seek(0)  # Reset if no header
            
            for line in f:
                if line.strip() and not line.startswith('system_id'):
                    fields = line.strip().split('\t')
                    if len(fields) >= 11:
                        system_id = fields[0]
                        contig = fields[1]
                        ice_type = fields[2]
                        start = int(fields[3])
                        end = int(fields[4])
                        length = int(fields[5])
                        gc_content = float(fields[6])
                        num_genes = int(fields[7])
                        flank_start = int(fields[8])
                        flank_end = int(fields[9])
                        contig_length = int(fields[10])
                        
                        ice_predictions[system_id] = {
                            'system_id': system_id,
                            'contig': contig,
                            'type': ice_type,
                            'original_start': start,
                            'original_end': end,
                            'original_length': length,
                            'gc_content': gc_content,
                            'num_genes': num_genes,
                            'flank_start': flank_start,
                            'flank_end': flank_end,
                            'contig_length': contig_length
                        }
    except FileNotFoundError:
        print(f"Error: MacSyFinder file '{macsyfinder_file}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error parsing MacSyFinder file: {e}", file=sys.stderr)
        sys.exit(1)
    
    return ice_predictions

def parse_trna_gff(trna_file):
    """Parse tRNA predictions in GFF format"""
    trna_predictions = defaultdict(list)
    
    try:
        with open(trna_file, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    fields = line.strip().split('\t')
                    if len(fields) >= 9 and fields[2] == 'tRNA':
                        seqname = fields[0]
                        start = int(fields[3])
                        end = int(fields[4])
                        strand = fields[6]
                        
                        # Extract product from attributes
                        attributes = fields[8]
                        product = 'tRNA'
                        if 'product=' in attributes:
                            product = attributes.split('product=')[1].split(';')[0]
                        
                        trna_predictions[seqname].append({
                            'start': start,
                            'end': end,
                            'product': product
                        })
    except FileNotFoundError:
        print(f"Error: tRNA file '{trna_file}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error parsing tRNA file: {e}", file=sys.stderr)
        sys.exit(1)
    
    return trna_predictions

def parse_direct_repeats(dr_file):
    """Parse direct repeats from vmatch output"""
    direct_repeats = defaultdict(list)

    try:
        with open(dr_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                # Seqname is the contig ID
                seqname = fields[0]
                start1 = int(fields[1])
                end1 = int(fields[2])
                start2 = int(fields[3])
                end2 = int(fields[4])
                direct_repeats[seqname].append({
                    'start1': start1,
                    'end1': end1,
                    'start2': start2,
                    'end2': end2,
                    'length1': end1 - start1 + 1,
                    'length2': end2 - start2 + 1,
                    'distance': start2 - end1
                })
    except FileNotFoundError:
        print(f"Error: Direct repeats file '{dr_file}' not found", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"Error parsing direct repeats file: {e}", file=sys.stderr)
        sys.exit(1)

    return direct_repeats


def find_most_distal_trna(ice_prediction, trna_list):
    """
    Find the most distal tRNA relative to the original ICE boundaries.
    The most distal tRNA should be closer to flank_start and/or flank_end and outside the original region.
    
    Args:
        ice_prediction: Dictionary containing ICE prediction data
        trna_list: List of tRNA predictions with coordinates relative to extended region
    
    Returns:
        Dictionary with distal tRNA information and gap analysis
    """
    original_start = ice_prediction['original_start']
    original_end = ice_prediction['original_end']
    flank_start = ice_prediction['flank_start']
    flank_end = ice_prediction['flank_end']
    contig_length = ice_prediction['contig_length']

    # Convert tRNA coordinates from extended region to original contig coordinates
    # tRNA coordinates are relative to the extended region starting at flank_start
    corrected_trnas = []
    for trna in trna_list:
        corrected_start = trna['start'] + flank_start - 1  # Convert to contig coordinates
        corrected_end = trna['end'] + flank_start - 1
        corrected_trnas.append({
            'corrected_start': corrected_start,
            'corrected_end': corrected_end,
            'product': trna['product']
        })
    
    # Create extended boundaries (equivalent to nfICEnum, neICEnum in single.py)
    extended_start = flank_start
    extended_end = flank_end
    
    # Build position list for gap analysis (following single.py ICEtagnum logic)
    position_list = [extended_start, extended_end]
    
    # Add tRNA positions within extended region (outside original ICE boundaries)
    distal_trnas = []
    for trna in corrected_trnas:
        trna_pos = trna['corrected_start']  # Use start position as representative
        
        # Only consider tRNAs outside the original ICE boundaries
        if (trna_pos < original_start or trna_pos > original_end) and \
           (extended_start <= trna_pos <= extended_end):
            position_list.append(trna_pos)
            distal_trnas.append(trna)
    
    position_list.sort()
    
    # Find largest gap between consecutive positions (following single.py find_max_distance)
    gap_start = None
    gap_end = None
    gap_location = None  # 'start', 'end', or 'internal'
    
    if len(position_list) >= 2:
        max_gap = 0
        for i in range(len(position_list) - 1):
            gap = position_list[i + 1] - position_list[i]
            if gap > max_gap:
                max_gap = gap
                gap_start = position_list[i]
                gap_end = position_list[i + 1]
        
        # Determine gap location relative to original ICE (following single.py logic)
        if gap_end == extended_end:
            gap_location = 'end'  # Gap is at the END boundary
        elif gap_start == extended_start:
            gap_location = 'start'  # Gap is at the START boundary
        else:
            gap_location = 'internal'  # Gap is internal
    
    # Find the most distal tRNAs (those closest to flank boundaries)
    most_distal_trnas = []
    if distal_trnas:
        # Find tRNAs closest to flank_start (upstream distal)
        upstream_distal = min(distal_trnas, key=lambda t: abs(t['corrected_start'] - flank_start))
        if upstream_distal['corrected_start'] < original_start:
            most_distal_trnas.append(('upstream', upstream_distal))
        
        # Find tRNAs closest to flank_end (downstream distal)
        downstream_distal = min(distal_trnas, key=lambda t: abs(t['corrected_start'] - flank_end))
        if downstream_distal['corrected_start'] > original_end:
            most_distal_trnas.append(('downstream', downstream_distal))
    
    return {
        'gap_start': gap_start,
        'gap_end': gap_end,
        'gap_location': gap_location,
        'most_distal_trnas': most_distal_trnas,
        'all_distal_trnas': distal_trnas,
        'position_list': position_list
    }


def find_dr_flanking_ice_and_trna(ice_prediction, distal_trna_coordinates, direct_repeats):
    """
    Parse the list of direct repeat pairs to find those flanking the ICE + distal tRNAs.
    
    Args:
        ice_prediction: Dictionary containing ICE prediction data
        distal_trna_coordinates: Result from find_most_distal_trna function
        direct_repeats: List of direct repeat pairs
    
    Returns:
        Dictionary with selected DR pair and refined boundaries
    """

    original_start = ice_prediction['original_start']
    original_end = ice_prediction['original_end']
    flank_start = ice_prediction['flank_start']
    flank_end = ice_prediction['flank_end']
    
    gap_start = distal_trna_coordinates['gap_start']
    gap_end = distal_trna_coordinates['gap_end']
    gap_location = distal_trna_coordinates['gap_location']
    all_distal_trnas = distal_trna_coordinates['all_distal_trnas']
    
    # Initialize refined boundaries based on gap analysis
    refined_start = original_start
    refined_end = original_end
    
    # Apply gap-based boundary refinement (following single.py logic)
    if gap_location == 'end' and gap_start is not None:
        # Gap at END - move start boundary to gap start
        refined_start = gap_start
    elif gap_location == 'start' and gap_end is not None:
        # Gap at START - move end boundary to gap end
        refined_end = gap_end
    
    # Find suitable DR pairs (following single.py criteria)
    suitable_drs = []
    for dr in direct_repeats:
        # Filter by distance (5kb - 500kb total span)
        total_span = dr['end2'] - dr['start1']
        if total_span > 500000 or total_span < 5000:
            continue
        
        # Check if DR pair flanks the refined ICE + distal tRNAs
        left_dr_start = dr['start1']
        left_dr_end = dr['end1']
        right_dr_start = dr['start2']
        right_dr_end = dr['end2']

        # Ensure proper flanking (left DR before refined region, right DR after)
        if not (left_dr_end <= refined_start and right_dr_start >= refined_end):
            continue
        
        # Count tRNAs between DR pair (following single.py checktrna logic)
        trnas_between = 0
        for trna in all_distal_trnas:
            trna_start = trna['corrected_start']
            trna_end = trna['corrected_end']
            
            # Check if tRNA is between the inner edges of the DR pair
            if left_dr_end <= trna_start <= right_dr_start and left_dr_end <= trna_end <= right_dr_start:
                trnas_between += 1
        
        # Require at least 2 tRNAs between DRs (following single.py requirement)
        if trnas_between >= 1:
            # Additional check based on gap location (following single.py positioning logic)
            valid_positioning = True
            
            if gap_location == 'end':
                # Gap at END - check if left DR is within refined start region
                if refined_start - 1000 <= left_dr_start <= refined_start + 1000:
                    valid_positioning = True
            elif gap_location == 'start':
                # Gap at START - check if right DR is within refined end region  
                if refined_end - 1000 <= right_dr_end <= refined_end + 1000:
                    valid_positioning = True
            else:
                # No specific gap - general flanking is sufficient
                valid_positioning = True
        
            if valid_positioning:
                suitable_drs.append({
                    'dr': dr,
                    'trnas_between': trnas_between,
                    'total_span': total_span,
                    'score': trnas_between  # Score based on tRNA count
                })
    
    # Select best DR pair
    selected_dr = None
    dr1_start = dr1_end = dr2_start = dr2_end = None
    final_refined_start = refined_start
    final_refined_end = refined_end
    
    if suitable_drs:
        # Sort by score (number of tRNAs between DRs), then by smaller total span
        suitable_drs.sort(key=lambda x: (-x['score'], x['total_span']))
        best_dr_info = suitable_drs[0]
        selected_dr = best_dr_info['dr']
        
        # Set DR coordinates (attL and attR)
        dr1_start = selected_dr['start1']  # attL start
        dr1_end = selected_dr['end1']      # attL end
        dr2_start = selected_dr['start2']  # attR start
        dr2_end = selected_dr['end2']      # attR end
        
        # Update final boundaries to be between the inner edges of DRs
        final_refined_start = dr1_end + 1
        final_refined_end = dr2_start - 1
    
    # Ensure boundaries are within flanking regions
    final_refined_start = max(flank_start, final_refined_start)
    final_refined_end = min(flank_end, final_refined_end)
    final_refined_length = max(0, final_refined_end - final_refined_start + 1)
    
    return {
        'refined_start': final_refined_start,
        'refined_end': final_refined_end,
        'refined_length': final_refined_length,
        'dr1_start': dr1_start if dr1_start is not None else 'NA',
        'dr1_end': dr1_end if dr1_end is not None else 'NA',
        'dr2_start': dr2_start if dr2_start is not None else 'NA',
        'dr2_end': dr2_end if dr2_end is not None else 'NA',
        'selected_dr': selected_dr,
        'suitable_drs_count': len(suitable_drs),
        'gap_location': gap_location
    }




def refine_ice_boundaries(ice_prediction, trna_list, dr_list, verbose=False):
    """
    Refine ICE boundaries using tRNAs and direct repeats
    Adapted from single.py merge_tRNA function to work with genomic coordinates
    """
    system_id = ice_prediction['system_id']
    
    # Find most distal tRNAs and gap analysis
    distal_trna_coordinates = find_most_distal_trna(ice_prediction, trna_list)
    
    # Find DR pairs flanking ICE + distal tRNAs
    refined_boundaries = find_dr_flanking_ice_and_trna(ice_prediction, distal_trna_coordinates, dr_list)
    
    if verbose:
        print(f"System {system_id}: Gap location: {distal_trna_coordinates['gap_location']}", file=sys.stderr)
        print(f"System {system_id}: Found {len(distal_trna_coordinates['all_distal_trnas'])} distal tRNAs", file=sys.stderr)
        print(f"System {system_id}: Found {refined_boundaries['suitable_drs_count']} suitable DR pairs", file=sys.stderr)
        print(f"System {system_id}: Refined boundaries: {refined_boundaries['refined_start']}-{refined_boundaries['refined_end']}", file=sys.stderr)
    
    return {
        'refined_start': refined_boundaries['refined_start'],
        'refined_end': refined_boundaries['refined_end'],
        'refined_length': refined_boundaries['refined_length'],
        'dr1_start': refined_boundaries['dr1_start'],
        'dr1_end': refined_boundaries['dr1_end'],
        'dr2_start': refined_boundaries['dr2_start'],
        'dr2_end': refined_boundaries['dr2_end'],
        'num_trnas': len(distal_trna_coordinates['all_distal_trnas']),
        'gap_location': refined_boundaries['gap_location']
    }


def main():
    parser = argparse.ArgumentParser(description='Refine ICE boundaries using tRNAs and direct repeats')
    parser.add_argument('--macsyfinder', help='MacSyFinder results file with ICE coordinates')
    parser.add_argument('--trna', help='tRNA predictions in GFF format')
    parser.add_argument('--drs', help='Direct repeats in TSV format')
    parser.add_argument('-o', '--output', help='Output file for refined ICE boundaries (default: stdout)')
    parser.add_argument('-v', '--verbose', action='store_true', help='Print detailed processing information')
    
    args = parser.parse_args()
    
    if args.verbose:
        print("Parsing input files...", file=sys.stderr)
    
    # Parse input files
    ice_predictions = parse_macsyfinder_results(args.macsyfinder)
    trna_predictions = parse_trna_gff(args.trna)
    direct_repeats = parse_direct_repeats(args.drs)
    
    if args.verbose:
        print(f"Found {len(ice_predictions)} ICE predictions", file=sys.stderr)
        print(f"Found tRNA predictions for {len(trna_predictions)} sequences", file=sys.stderr)
        print(f"Found direct repeats for {len(direct_repeats)} sequences", file=sys.stderr)
    
    # Process each ICE prediction
    refined_results = []
    
    for system_id, ice_pred in ice_predictions.items():
        if args.verbose:
            print(f"Processing system {system_id}...", file=sys.stderr)
        
        # Find corresponding tRNAs and DRs
        flanked_seqname = f"{system_id}_with_flanks"
        contig_id = ice_pred['contig']
        
        trnas = trna_predictions.get(flanked_seqname, [])
        drs = direct_repeats.get(contig_id, [])
        
        # Also try without _with_flanks suffix
        if not trnas:
            trnas = trna_predictions.get(system_id, [])
        if not drs:
            drs = direct_repeats.get(system_id, [])
        
        if args.verbose:
            print(f"  Found {len(trnas)} tRNAs and {len(drs)} DR pairs", file=sys.stderr)
        
        # Refine boundaries
        refinement = refine_ice_boundaries(ice_pred, trnas, drs, args.verbose)
        
        # Combine original prediction with refinement results
        result = {
            'system_id': system_id,
            'contig': ice_pred['contig'],
            'type': ice_pred['type'],
            'original_start': ice_pred['original_start'],
            'original_end': ice_pred['original_end'],
            'refined_start': refinement['refined_start'],
            'refined_end': refinement['refined_end'],
            'refined_length': refinement['refined_length'],
            'flank_start': ice_pred['flank_start'],
            'flank_end': ice_pred['flank_end'],
            'dr1_start': refinement['dr1_start'],
            'dr1_end': refinement['dr1_end'],
            'dr2_start': refinement['dr2_start'],
            'dr2_end': refinement['dr2_end'],
            'num_trnas': refinement['num_trnas'],
            'gc_content': ice_pred['gc_content']
        }
        
        refined_results.append(result)
        
        if args.verbose:
            print(f"  Refined: {ice_pred['original_start']}-{ice_pred['original_end']} -> "
                  f"{refinement['refined_start']}-{refinement['refined_end']} "
                  f"({refinement['num_trnas']} tRNAs)", file=sys.stderr)
    
    # Write output
    output_file = open(args.output, 'w', newline='') if args.output else sys.stdout
    
    try:
        fieldnames = [
            'system_id', 'contig', 'type', 'original_start', 'original_end',
            'refined_start', 'refined_end', 'refined_length', 'flank_start', 'flank_end',
            'dr1_start', 'dr1_end', 'dr2_start', 'dr2_end', 'num_trnas', 'gc_content'
        ]
        
        writer = csv.DictWriter(output_file, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        
        for result in refined_results:
            writer.writerow(result)
    
    finally:
        if args.output:
            output_file.close()
    
    if args.verbose:
        print(f"Processed {len(refined_results)} ICE predictions", file=sys.stderr)
        print("Refinement complete!", file=sys.stderr)

if __name__ == '__main__':
    main()
