#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import argparse
from Bio import SeqIO
from Bio.SeqUtils import GC

def parse_macsyfinder_results(macsyfinder_file):
    """Parse MacsyFinder results file"""
    ice_dict = {}
    with open(macsyfinder_file, 'r') as f:
        header = f.readline().strip().split('\t')
        for line in f:
            fields = line.strip().split('\t')
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
            
            ice_dict[system_id] = {
                'contig': contig,
                'type': ice_type,
                'start': start,
                'end': end,
                'length': length,
                'gc_content': gc_content,
                'num_genes': num_genes,
                'flank_start': flank_start,
                'flank_end': flank_end,
                'contig_length': contig_length
            }
    return ice_dict

def parse_trna_gff(trna_file):
    """Parse tRNA GFF file"""
    trna_dict = {}
    with open(trna_file, 'r') as f:
        for line in f:
            if line.startswith('#'):
                continue
            fields = line.strip().split('\t')
            if len(fields) < 9:
                continue
            
            sequence_name = fields[0]
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]
            
            # Extract tRNA type from attributes
            trna_type = "tRNA"
            if 'product=' in attributes:
                trna_type = attributes.split('product=')[1].split(';')[0]
            
            if sequence_name not in trna_dict:
                trna_dict[sequence_name] = []
            
            trna_dict[sequence_name].append({
                'start': start,
                'end': end,
                'strand': strand,
                'type': trna_type
            })
    
    return trna_dict

def parse_direct_repeats(dr_file):
    """Parse direct repeats TSV file"""
    dr_list = []
    with open(dr_file, 'r') as f:
        for line in f:
            fields = line.strip().split('\t')
            if len(fields) >= 5:
                sequence_name = fields[0]
                start1 = int(fields[1])
                end1 = int(fields[2])
                start2 = int(fields[3])
                end2 = int(fields[4])
                
                dr_list.append({
                    'sequence_name': sequence_name,
                    'start1': start1,
                    'end1': end1,
                    'start2': start2,
                    'end2': end2
                })
    
    return dr_list

def find_max_distance(numbers):
    """Find the maximum distance between consecutive numbers"""
    if len(numbers) < 2:
        return None, None
    
    max_distance = -1
    max_distance_index = -1
    
    for i in range(len(numbers) - 1):
        distance = abs(numbers[i] - numbers[i+1])
        if distance > max_distance:
            max_distance = distance
            max_distance_index = i
    
    if max_distance_index == -1:
        return None, None
    
    return numbers[max_distance_index], numbers[max_distance_index + 1]

def refine_ice_boundaries(ice_dict, trna_dict, dr_list):
    """Refine ICE boundaries using tRNA and direct repeats"""
    refined_ices = []
    
    for system_id, ice_info in ice_dict.items():
        # Get the extended sequence name
        extended_name = f"{system_id}_with_flanks"
        
        # Initial boundaries from MacsyFinder
        ice_start = ice_info['start']
        ice_end = ice_info['end']
        flank_start = ice_info['flank_start']
        flank_end = ice_info['flank_end']
        
        # Initialize refined boundaries
        refined_start = ice_start
        refined_end = ice_end
        dr1_start = ""
        dr1_end = ""
        dr2_start = ""
        dr2_end = ""
        
        # Check for tRNAs in the extended region
        trna_positions = []
        trna_info = []
        
        if extended_name in trna_dict:
            for trna in trna_dict[extended_name]:
                trna_positions.append(trna['start'])
                trna_positions.append(trna['end'])
                trna_info.append(trna)
        
        # Create position list for boundary refinement
        boundary_positions = [flank_start, flank_end]
        boundary_positions.extend(trna_positions)
        boundary_positions.sort()
        
        # Find maximum distance gap (potential ICE boundaries)
        if len(boundary_positions) > 2:
            gap_start, gap_end = find_max_distance(boundary_positions)
            if gap_start is not None and gap_end is not None:
                # Check for direct repeats that could define boundaries
                for dr in dr_list:
                    if dr['sequence_name'] == extended_name:
                        dr_span = dr['end2'] - dr['start1']
                        
                        # Filter DRs by size (5kb to 500kb as in reference)
                        if dr_span < 5000 or dr_span > 500000:
                            continue
                        
                        # Check if DR encompasses the ICE region
                        if (dr['start1'] <= ice_start <= dr['end2'] and 
                            dr['start1'] <= ice_end <= dr['end2']):
                            
                            # Count tRNAs within DR region
                            trna_count = 0
                            for trna in trna_info:
                                if (dr['start1'] <= trna['start'] <= dr['end2'] and
                                    dr['start1'] <= trna['end'] <= dr['end2']):
                                    trna_count += 1
                            
                            # If sufficient tRNAs, use DR boundaries
                            if trna_count >= 2:
                                refined_start = dr['start1']
                                refined_end = dr['end2']
                                dr1_start = str(dr['start1'])
                                dr1_end = str(dr['end1'])
                                dr2_start = str(dr['start2'])
                                dr2_end = str(dr['end2'])
                                break
        
        # Prepare output information
        refined_ice = {
            'system_id': system_id,
            'contig': ice_info['contig'],
            'type': ice_info['type'],
            'original_start': ice_start,
            'original_end': ice_end,
            'refined_start': refined_start,
            'refined_end': refined_end,
            'refined_length': refined_end - refined_start + 1,
            'flank_start': flank_start,
            'flank_end': flank_end,
            'dr1_start': dr1_start,
            'dr1_end': dr1_end,
            'dr2_start': dr2_start,
            'dr2_end': dr2_end,
            'num_trnas': len(trna_info),
            'gc_content': ice_info['gc_content']
        }
        
        refined_ices.append(refined_ice)
    
    return refined_ices

def write_output(refined_ices, output_file):
    """Write refined ICE results to output file"""
    header = [
        'system_id', 'contig', 'type', 'original_start', 'original_end',
        'refined_start', 'refined_end', 'refined_length', 'flank_start', 'flank_end',
        'dr1_start', 'dr1_end', 'dr2_start', 'dr2_end', 'num_trnas', 'gc_content'
    ]
    
    with open(output_file, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for ice in refined_ices:
            row = [str(ice[col]) for col in header]
            f.write('\t'.join(row) + '\n')

def main():
    parser = argparse.ArgumentParser(description='Refine ICE boundaries using tRNA and direct repeats')
    parser.add_argument('--macsyfinder', required=True, help='MacsyFinder results file')
    parser.add_argument('--trna', required=True, help='tRNA GFF file')
    parser.add_argument('--direct-repeats', required=True, help='Direct repeats TSV file')
    parser.add_argument('--output', required=True, help='Output file')
    
    args = parser.parse_args()
    
    # Parse input files
    print("Parsing MacsyFinder results...")
    ice_dict = parse_macsyfinder_results(args.macsyfinder)
    
    print("Parsing tRNA predictions...")
    trna_dict = parse_trna_gff(args.trna)
    
    print("Parsing direct repeats...")
    dr_list = parse_direct_repeats(args.direct_repeats)
    
    # Refine ICE boundaries
    print("Refining ICE boundaries...")
    refined_ices = refine_ice_boundaries(ice_dict, trna_dict, dr_list)
    
    # Write output
    print(f"Writing results to {args.output}...")
    write_output(refined_ices, args.output)
    
    print(f"Processed {len(refined_ices)} ICE elements")

if __name__ == "__main__":
    main()
