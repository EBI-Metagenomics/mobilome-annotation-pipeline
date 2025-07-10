#!/usr/bin/env python3

import pandas as pd
import sys
import argparse

def validate_ice_elements(boundaries_file, trna_file, output_file, summary_file, 
                         min_ice_size=10000, max_ice_size=500000, min_trna_count=1):
    '''Validate ICE elements and generate final results'''
    
    # Read input files
    try:
        boundaries = pd.read_csv(boundaries_file, sep='\t')
        print(f"Loaded {len(boundaries)} ICE boundaries")
    except Exception as e:
        print(f"Error reading boundaries file: {e}", file=sys.stderr)
        return 0
    
    try:
        # Try to read tRNA summary file (derived from GFF3 filename)
        if trna_file.endswith('.gff3'):
            trna_summary_file = trna_file.replace('.gff3', '_summary.tsv')
        else:
            trna_summary_file = trna_file
        trna_summary = pd.read_csv(trna_summary_file, sep='\t')
        print(f"Loaded {len(trna_summary)} tRNA annotations")
    except Exception as e:
        print(f"Warning: No tRNA summary available - {e}")
        trna_summary = pd.DataFrame()
    
    validated_elements = []
    validation_stats = {
        'total_elements': len(boundaries),
        'high_confidence': 0,
        'medium_confidence': 0,
        'low_confidence': 0,
        'failed': 0,
        'with_trnas': 0,
        'with_direct_repeats': 0,
        'size_distribution': {}
    }
    
    for _, ice in boundaries.iterrows():
        system_id = ice['system_id']
        ice_length = ice['ice_length']
        confidence = ice['confidence']
        num_trnas = ice['num_trnas']
        num_drs = ice['num_direct_repeats']
        
        # Validation criteria
        is_valid = True
        validation_notes = []
        
        # Size validation
        if ice_length < min_ice_size:
            validation_notes.append(f"Below minimum size ({ice_length} < {min_ice_size})")
            is_valid = False
        elif ice_length > max_ice_size:
            validation_notes.append(f"Above maximum size ({ice_length} > {max_ice_size})")
            is_valid = False
        
        # tRNA validation
        if num_trnas < min_trna_count:
            validation_notes.append(f"Insufficient tRNAs ({num_trnas} < {min_trna_count})")
            if confidence == 'high':
                confidence = 'medium'  # Downgrade confidence
        
        # Update statistics
        if is_valid:
            validation_stats[f'{confidence}_confidence'] += 1
        else:
            validation_stats['failed'] += 1
            
        if num_trnas > 0:
            validation_stats['with_trnas'] += 1
        if num_drs > 0:
            validation_stats['with_direct_repeats'] += 1
        
        # Size distribution
        size_category = 'small' if ice_length < 20000 else 'medium' if ice_length < 100000 else 'large'
        validation_stats['size_distribution'][size_category] = validation_stats['size_distribution'].get(size_category, 0) + 1
        
        if is_valid:
            validated_elements.append({
                'system_id': system_id,
                'contig': ice['contig'],
                'start': ice['refined_start'],
                'end': ice['refined_end'],
                'length': ice_length,
                'confidence': confidence,
                'refinement_method': ice['refinement_method'],
                'num_trnas': num_trnas,
                'num_direct_repeats': num_drs,
                'validation_status': 'PASS',
                'notes': '; '.join(validation_notes) if validation_notes else 'All criteria met'
            })
        else:
            validated_elements.append({
                'system_id': system_id,
                'contig': ice['contig'],
                'start': ice['refined_start'],
                'end': ice['refined_end'],
                'length': ice_length,
                'confidence': 'failed',
                'refinement_method': ice['refinement_method'],
                'num_trnas': num_trnas,
                'num_direct_repeats': num_drs,
                'validation_status': 'FAIL',
                'notes': '; '.join(validation_notes)
            })
    
    # Write validated elements
    try:
        if validated_elements:
            df = pd.DataFrame(validated_elements)
            df.to_csv(output_file, sep='\t', index=False)
        else:
            with open(output_file, 'w') as f:
                f.write('system_id\tcontig\tstart\tend\tlength\tconfidence\trefinement_method\tnum_trnas\tnum_direct_repeats\tvalidation_status\tnotes\n')
    except Exception as e:
        print(f"Error writing output file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Write summary report
    try:
        with open(summary_file, 'w') as f:
            f.write("ICE Boundary Detection Summary Report\n")
            f.write("=" * 40 + "\n\n")
            
            f.write(f"Total ICE elements processed: {validation_stats['total_elements']}\n")
            f.write(f"High confidence elements: {validation_stats['high_confidence']}\n")
            f.write(f"Medium confidence elements: {validation_stats['medium_confidence']}\n")
            f.write(f"Low confidence elements: {validation_stats['low_confidence']}\n")
            f.write(f"Failed validation: {validation_stats['failed']}\n")
            f.write(f"Elements with tRNAs: {validation_stats['with_trnas']}\n")
            f.write(f"Elements with direct repeats: {validation_stats['with_direct_repeats']}\n\n")
            
            f.write("Size Distribution:\n")
            for size_cat, count in validation_stats['size_distribution'].items():
                f.write(f"  {size_cat}: {count}\n")
            
            passed_count = len([e for e in validated_elements if e['validation_status'] == 'PASS'])
            failed_count = len([e for e in validated_elements if e['validation_status'] == 'FAIL'])
            
            f.write(f"\nValidation passed: {passed_count}\n")
            f.write(f"Validation failed: {failed_count}\n")
            
            # Additional statistics
            f.write(f"\nSuccess rate: {passed_count/len(validated_elements)*100:.1f}%\n")
            
            # Method breakdown
            f.write("\nRefinement method breakdown:\n")
            method_counts = {}
            for element in validated_elements:
                method = element['refinement_method']
                method_counts[method] = method_counts.get(method, 0) + 1
            
            for method, count in sorted(method_counts.items()):
                f.write(f"  {method}: {count}\n")
                
    except Exception as e:
        print(f"Error writing summary file: {e}", file=sys.stderr)
        sys.exit(1)
    
    return len([e for e in validated_elements if e['validation_status'] == 'PASS'])

def main():
    parser = argparse.ArgumentParser(description='Validate ICE elements and generate final results')
    parser.add_argument('boundaries_file', help='Input refined boundaries TSV file')
    parser.add_argument('trna_file', help='Input tRNA annotations file (GFF3 or summary TSV)')
    parser.add_argument('output_file', help='Output validated ICE elements TSV file')
    parser.add_argument('summary_file', help='Output summary report file')
    parser.add_argument('--min-ice-size', type=int, default=10000,
                       help='Minimum ICE size (default: 10000)')
    parser.add_argument('--max-ice-size', type=int, default=500000,
                       help='Maximum ICE size (default: 500000)')
    parser.add_argument('--min-trna-count', type=int, default=1,
                       help='Minimum tRNA count (default: 1)')
    
    args = parser.parse_args()
    
    try:
        num_valid = validate_ice_elements(
            args.boundaries_file,
            args.trna_file,
            args.output_file,
            args.summary_file,
            args.min_ice_size,
            args.max_ice_size,
            args.min_trna_count
        )
        print(f"Successfully validated {num_valid} ICE elements")
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
