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
                            'strand': strand,
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
                if line.strip():
                    fields = line.strip().split('\t')
                    if len(fields) >= 5:
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


def find_most_distal_trna( ice_prediction, trna_list ):
    # Use the coordinates of the original prediction (original_start and original_end) and the coordinates of the extended region 5kb upstreand and downstream the original prediction (corresponding to flank_start and flank_end) to find the most distal RNA. The most distal tRNA should be the closer to the flank_start and/or the flank_end and outside the original region. Consider that RNAs were predicted on the fasta file of the extended region, meaning that the coordinates of the tRNA have to be corrected based on the flank start



def find_dr_flanking_ice_and_trna( ice_prediction, distal_trna_coordinates, direct_repeats):
    # Parse the list of pair of direct repeats coordinates flanking the ICE + distal RNA(s). Consider that direct repeats were predicted on the fasta file of the extended region, meaning that the coordinates of the tRNA have to be corrected based on the flank start 


def refine_ice_boundaries(ice_prediction, trna_list, dr_list, verbose=False):
    """
    Refine ICE boundaries using tRNAs and direct repeats
    Adapted from single.py merge_tRNA function to work with genomic coordinates
    """

    distal_trna_coordinates = find_most_distal_trna(ice_prediction, trna_list)


    refined_boundaries = find_dr_flanking_ice_and_trna( ice_prediction, distal_trna_coordinates, direct_repeats )



    '''
    system_id = ice_prediction['system_id']
    original_start = ice_prediction['original_start']
    original_end = ice_prediction['original_end']
    flank_start = ice_prediction['flank_start']  # This is the start of the extended region
    flank_end = ice_prediction['flank_end']
    contig_length = ice_prediction['contig_length']


    # Initialize refined boundaries
    refined_start = original_start
    refined_end = original_end
    dr1_start = dr1_end = dr2_start = dr2_end = None
    num_trnas = len(trna_list)
    

    # Create position list for gap analysis (adapted from ICEtagnum logic)
    # Include ICE boundaries and tRNA positions
    position_list = [original_start, original_end]
    
    # Add tRNA positions (using start positions as representative points)
    for trna in trna_list:
        position_list.append(trna['start'])
    
    # Sort positions
    position_list = sorted(set(position_list))  # Remove duplicates and sort
    
    if verbose:
        print(f"System {system_id}: Found {num_trnas} tRNAs", file=sys.stderr)
        print(f"Position list: {position_list}", file=sys.stderr)
    
    # Find largest gap between positions (following single.py gap logic)
    if len(position_list) >= 2:
        gap_start, gap_end = find_largest_gap_in_positions(position_list)
        
        if gap_start is not None and gap_end is not None:
            # The largest gap indicates where the ICE likely ends/begins
            # Following the original logic: the gap represents the ICE boundaries
            if gap_start >= original_start and gap_end <= original_end:
                # Gap is within the ICE - this suggests the ICE should be split
                # Use the gap to refine boundaries
                refined_start = gap_start
                refined_end = gap_end
            elif gap_start < original_start:
                # Gap before ICE start - extend start to gap end
                refined_start = min(gap_end, original_start)
            elif gap_end > original_end:
                # Gap after ICE end - extend end to gap start  
                refined_end = max(gap_start, original_end)
            
            if verbose:
                print(f"System {system_id}: Largest gap between {gap_start} and {gap_end}", file=sys.stderr)
    
    # Process direct repeats (following single.py DR selection logic)
    suitable_drs = []
    for dr in dr_list:
        dr_distance = dr['distance']
        
        # Filter DRs based on distance criteria (from single.py)
        if dr_distance > 500000 or dr_distance < 5000:
            continue
        
        # Check if DR appropriately flanks the refined ICE region
        dr_start1 = dr['start1']
        dr_end1 = dr['end1']
        dr_start2 = dr['start2']
        dr_end2 = dr['end2']
        
        # Ensure DRs flank the ICE (one before, one after)
        if (dr_end1 <= refined_start and dr_start2 >= refined_end) or \
           (dr_end2 <= refined_start and dr_start1 >= refined_end):
            
            # Count tRNAs between the DR pair (checktrna logic from single.py)
            trnas_between = 0
            dr_region_start = min(dr_start1, dr_start2)
            dr_region_end = max(dr_end1, dr_end2)
            
            for trna in trna_list:
                if dr_region_start <= trna['start'] <= dr_region_end:
                    trnas_between += 1
            
            # Require at least 2 tRNAs between DRs (following single.py logic)
            if trnas_between >= 2:
                suitable_drs.append({
                    'dr': dr,
                    'trnas_between': trnas_between,
                    'score': trnas_between  # Score based on tRNA count
                })
    
    # Select best DR pair if available
    if suitable_drs:
        # Sort by score (number of tRNAs between DRs)
        suitable_drs.sort(key=lambda x: x['score'], reverse=True)
        best_dr = suitable_drs[0]['dr']
        
        # Update boundaries based on selected DR pair
        if best_dr['end1'] <= refined_start:
            # DR1 is upstream, DR2 is downstream
            dr1_start = best_dr['start1']
            dr1_end = best_dr['end1']
            dr2_start = best_dr['start2']
            dr2_end = best_dr['end2']
            refined_start = dr1_end + 1
            refined_end = dr2_start - 1
        else:
            # DR2 is upstream, DR1 is downstream
            dr1_start = best_dr['start2']
            dr1_end = best_dr['end2']
            dr2_start = best_dr['start1']
            dr2_end = best_dr['end1']
            refined_start = dr1_end + 1
            refined_end = dr2_start - 1
        
        if verbose:
            print(f"System {system_id}: Selected DR pair with {suitable_drs[0]['score']} tRNAs between", file=sys.stderr)
    
    # Ensure refined boundaries are within flanking regions
    refined_start = max(flank_start, refined_start)
    refined_end = min(flank_end, refined_end)
    refined_length = max(0, refined_end - refined_start + 1)
    
    return {
        'refined_start': refined_start,
        'refined_end': refined_end,
        'refined_length': refined_length,
        'dr1_start': dr1_start if dr1_start is not None else 'NA',
        'dr1_end': dr1_end if dr1_end is not None else 'NA',
        'dr2_start': dr2_start if dr2_start is not None else 'NA',
        'dr2_end': dr2_end if dr2_end is not None else 'NA',
        'num_trnas': num_trnas
    }


    '''


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
        
        trnas = trna_predictions.get(flanked_seqname, [])
        drs = direct_repeats.get(flanked_seqname, [])
        
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
