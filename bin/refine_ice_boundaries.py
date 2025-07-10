#!/usr/bin/env python3

import pandas as pd
import sys
import argparse
from collections import defaultdict

def refine_ice_boundaries(coordinates_file, trna_file, repeats_file, output_file, log_file, 
                         min_ice_size=10000, max_ice_size=500000):
    '''Refine ICE boundaries using tRNA and direct repeat information'''
    
    # Read input files
    try:
        ice_coords = pd.read_csv(coordinates_file, sep='\t')
        print(f"Loaded {len(ice_coords)} ICE coordinates")
    except Exception as e:
        print(f"Error reading ICE coordinates file: {e}", file=sys.stderr)
        return 0
    
    try:
        # Try to read tRNA summary file (derived from GFF3 filename)
        if trna_file.endswith('.gff3'):
            trna_summary_file = trna_file.replace('.gff3', '_summary.tsv')
        else:
            trna_summary_file = trna_file
        trna_data = pd.read_csv(trna_summary_file, sep='\t')
        print(f"Loaded {len(trna_data)} tRNA annotations")
    except Exception as e:
        print(f"Warning: No tRNA data available - {e}")
        trna_data = pd.DataFrame()
    
    try:
        repeat_data = pd.read_csv(repeats_file, sep='\t')
        print(f"Loaded {len(repeat_data)} direct repeat entries")
    except Exception as e:
        print(f"Warning: No direct repeat data available - {e}")
        repeat_data = pd.DataFrame()
    
    refined_boundaries = []
    log_entries = []
    
    for _, ice_region in ice_coords.iterrows():
        system_id = ice_region['system_id']
        original_start = ice_region['ice_start']
        original_end = ice_region['ice_end']
        extended_start = ice_region['extended_start']
        extended_end = ice_region['extended_end']
        contig = ice_region['contig']
        
        log_entries.append(f"Processing ICE system: {system_id}")
        log_entries.append(f"Original boundaries: {original_start}-{original_end}")
        
        # Find tRNAs in the extended region
        sequence_name = f"{system_id}_extended"
        region_trnas = trna_data[trna_data['sequence_name'] == sequence_name] if not trna_data.empty else pd.DataFrame()
        
        # Find direct repeats for this region
        region_repeats = repeat_data[repeat_data['sequence_id'] == sequence_name] if not repeat_data.empty else pd.DataFrame()
        
        refined_start = original_start
        refined_end = original_end
        refinement_method = "original"
        
        # Strategy 1: Use direct repeats if available
        if not region_repeats.empty:
            best_repeat = region_repeats.iloc[0]  # Already sorted by length
            
            # Adjust coordinates back to original contig coordinates
            dr_start1 = extended_start + best_repeat['start_pos1'] - 1
            dr_end1 = extended_start + best_repeat['end_pos1'] - 1
            dr_start2 = extended_start + best_repeat['start_pos2'] - 1
            dr_end2 = extended_start + best_repeat['end_pos2'] - 1
            
            # Use the outer boundaries of the direct repeats
            refined_start = min(dr_start1, dr_start2)
            refined_end = max(dr_end1, dr_end2)
            refinement_method = "direct_repeats"
            
            log_entries.append(f"Found direct repeats: {best_repeat['repeat_length']}bp")
            log_entries.append(f"DR boundaries: {refined_start}-{refined_end}")
        
        # Strategy 2: Use tRNA boundaries if no good direct repeats
        elif not region_trnas.empty and len(region_trnas) >= 2:
            # Find tRNAs that could flank the ICE
            trna_positions = []
            for _, trna in region_trnas.iterrows():
                # Convert back to original contig coordinates
                trna_pos = extended_start + trna['start'] - 1
                trna_positions.append(trna_pos)
            
            trna_positions.sort()
            
            # Use outermost tRNAs as boundaries
            if len(trna_positions) >= 2:
                refined_start = trna_positions[0]
                refined_end = trna_positions[-1]
                refinement_method = "trna_flanked"
                
                log_entries.append(f"Found {len(region_trnas)} tRNAs")
                log_entries.append(f"tRNA boundaries: {refined_start}-{refined_end}")
        
        # Strategy 3: Extend by fixed amount if no other evidence
        else:
            extension = 2000  # 2kb extension
            refined_start = max(1, original_start - extension)
            refined_end = original_end + extension
            refinement_method = "fixed_extension"
            
            log_entries.append("No tRNA or DR evidence, using fixed extension")
        
        # Ensure boundaries are within reasonable limits
        ice_length = refined_end - refined_start + 1
        if ice_length < min_ice_size:
            refined_start = max(1, original_start - 1000)
            refined_end = original_end + 1000
            refinement_method += "_size_adjusted"
            log_entries.append("Adjusted boundaries due to minimum size constraint")
        elif ice_length > max_ice_size:
            center = (original_start + original_end) // 2
            half_max = max_ice_size // 2
            refined_start = max(1, center - half_max)
            refined_end = center + half_max
            refinement_method += "_size_limited"
            log_entries.append("Limited boundaries due to maximum size constraint")
        
        refined_boundaries.append({
            'system_id': system_id,
            'contig': contig,
            'original_start': original_start,
            'original_end': original_end,
            'refined_start': refined_start,
            'refined_end': refined_end,
            'ice_length': refined_end - refined_start + 1,
            'refinement_method': refinement_method,
            'num_trnas': len(region_trnas),
            'num_direct_repeats': len(region_repeats),
            'confidence': 'high' if refinement_method.startswith('direct_repeats') else 
                        'medium' if refinement_method.startswith('trna') else 'low'
        })
        
        log_entries.append(f"Final boundaries: {refined_start}-{refined_end} ({refinement_method})")
        log_entries.append("---")
    
    # Write output files
    try:
        if refined_boundaries:
            df = pd.DataFrame(refined_boundaries)
            df.to_csv(output_file, sep='\t', index=False)
        else:
            with open(output_file, 'w') as f:
                f.write('system_id\tcontig\toriginal_start\toriginal_end\trefined_start\trefined_end\tice_length\trefinement_method\tnum_trnas\tnum_direct_repeats\tconfidence\n')
    except Exception as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Write log file
    try:
        with open(log_file, 'w') as f:
            f.write('\n'.join(log_entries))
    except Exception as e:
        print(f"Error writing log file: {e}", file=sys.stderr)
        sys.exit(1)
    
    return len(refined_boundaries)

def main():
    parser = argparse.ArgumentParser(description='Refine ICE boundaries using tRNA and direct repeat information')
    parser.add_argument('coordinates_file', help='Input ICE coordinates TSV file')
    parser.add_argument('trna_file', help='Input tRNA annotations file (GFF3 or summary TSV)')
    parser.add_argument('repeats_file', help='Input direct repeats TSV file')
    parser.add_argument('output_file', help='Output refined boundaries TSV file')
    parser.add_argument('log_file', help='Output log file')
    parser.add_argument('--min-ice-size', type=int, default=10000,
                       help='Minimum ICE size (default: 10000)')
    parser.add_argument('--max-ice-size', type=int, default=500000,
                       help='Maximum ICE size (default: 500000)')
    
    args = parser.parse_args()
    
    try:
        num_refined = refine_ice_boundaries(
            args.coordinates_file, 
            args.trna_file, 
            args.repeats_file,
            args.output_file, 
            args.log_file,
            args.min_ice_size,
            args.max_ice_size
        )
        print(f"Successfully refined boundaries for {num_refined} ICE systems")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
