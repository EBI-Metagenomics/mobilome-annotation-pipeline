#!/usr/bin/env python3

import json
import pandas as pd
from Bio import SeqIO
from Bio.SeqUtils import GC
import sys
import argparse

def calculate_gc_content(genome_file, start, end):
    '''Calculate GC content for a genomic region'''
    try:
        record = SeqIO.read(genome_file, "fasta")
        sequence = record.seq[int(start)-1:int(end)]
        return round(GC(sequence), 2)
    except Exception as e:
        print(f"Warning: Could not calculate GC content: {e}", file=sys.stderr)
        return 0.0

def calculate_gc_profile(genome_file, start, end, window_size=1000, step_size=500):
    '''Calculate GC content profile across ICE region'''
    try:
        record = SeqIO.read(genome_file, "fasta")
        sequence = record.seq[int(start)-1:int(end)]
        
        gc_profile = []
        positions = []
        
        for i in range(0, len(sequence) - window_size + 1, step_size):
            window = sequence[i:i+window_size]
            gc_content = GC(window)
            gc_profile.append(round(gc_content, 2))
            positions.append(int(start) + i + window_size//2)
        
        return positions, gc_profile
    except Exception as e:
        print(f"Warning: Could not calculate GC profile: {e}", file=sys.stderr)
        return [], []

def parse_macsyfinder_features(gff3_file):
    '''Parse MacSyFinder GFF3 to extract ICE gene features'''
    features = {}
    
    try:
        with open(gff3_file, 'r') as f:
            for line in f:
                if line.startswith('#'):
                    continue
                
                fields = line.strip().split('\t')
                if len(fields) < 9:
                    continue
                
                attributes = fields[8]
                attr_dict = {}
                for attr in attributes.split(';'):
                    if '=' in attr:
                        key, value = attr.split('=', 1)
                        attr_dict[key] = value
                
                if 'system' in attr_dict and 'Name' in attr_dict:
                    system_id = attr_dict['system']
                    gene_name = attr_dict['Name']
                    start = int(fields[3])
                    end = int(fields[4])
                    
                    if system_id not in features:
                        features[system_id] = []
                    
                    features[system_id].append({
                        'gene_name': gene_name,
                        'start': start,
                        'end': end,
                        'feature_type': fields[2],
                        'strand': fields[6]
                    })
    except Exception as e:
        print(f"Warning: Could not parse MacSyFinder GFF3: {e}", file=sys.stderr)
    
    return features

def parse_trna_data(trna_file):
    '''Parse tRNA annotation data'''
    trna_data = {}
    
    try:
        if trna_file.endswith('.gff3'):
            with open(trna_file, 'r') as f:
                for line in f:
                    if line.startswith('#'):
                        continue
                    
                    fields = line.strip().split('\t')
                    if len(fields) >= 9 and fields[2] == 'tRNA':
                        sequence_name = fields[0]
                        start = int(fields[3])
                        end = int(fields[4])
                        strand = fields[6]
                        
                        attributes = fields[8]
                        trna_type = 'tRNA'
                        anticodon = 'unknown'
                        
                        for attr in attributes.split(';'):
                            if 'anticodon=' in attr:
                                anticodon = attr.split('=')[1]
                            elif 'Name=' in attr:
                                trna_type = attr.split('=')[1]
                        
                        if sequence_name not in trna_data:
                            trna_data[sequence_name] = []
                        
                        trna_data[sequence_name].append({
                            'start': start,
                            'end': end,
                            'strand': strand,
                            'type': trna_type,
                            'anticodon': anticodon
                        })
        else:
            # Handle TSV format
            df = pd.read_csv(trna_file, sep='\t')
            for _, row in df.iterrows():
                seq_name = row['sequence_name']
                if seq_name not in trna_data:
                    trna_data[seq_name] = []
                
                trna_data[seq_name].append({
                    'start': row['start'],
                    'end': row['end'],
                    'strand': row['strand'],
                    'type': row['trna_type'],
                    'anticodon': row.get('anticodon', 'unknown')
                })
    except Exception as e:
        print(f"Warning: Could not parse tRNA data: {e}", file=sys.stderr)
    
    return trna_data

def parse_direct_repeats(dr_file):
    '''Parse direct repeat data'''
    dr_data = []
    
    try:
        df = pd.read_csv(dr_file, sep='\t')
        for _, row in df.iterrows():
            dr_data.append({
                'sequence_id': row['sequence_id'],
                'repeat_length': row['repeat_length'],
                'start_pos1': row['start_pos1'],
                'end_pos1': row['end_pos1'],
                'start_pos2': row['start_pos2'],
                'end_pos2': row['end_pos2'],
                'distance': row['distance']
            })
    except Exception as e:
        print(f"Warning: Could not parse direct repeat data: {e}", file=sys.stderr)
    
    return dr_data

def parse_orit_results(orit_file):
    '''Parse oriT detection results'''
    orit_data = {}
    
    try:
        df = pd.read_csv(orit_file, sep='\t')
        for _, row in df.iterrows():
            ice_id = row['ice_id']
            orit_data[ice_id] = {
                'orit_id': row['orit_id'],
                'identity_percent': row['identity_percent'],
                'coverage': row['coverage'],
                'h_value': row['h_value'],
                'query_start': row['query_start'],
                'query_end': row['query_end'],
                'evalue': row['evalue']
            }
    except Exception as e:
        print(f"Warning: Could not parse oriT results: {e}", file=sys.stderr)
    
    return orit_data

def integrate_ice_information(validated_ice_elements, macsyfinder_gff3, trna_annotations, 
                             direct_repeats, orit_results, genome_fasta, prefix,
                             gc_window_size=1000, gc_step_size=500):
    '''Main function to integrate all ICE information (equivalent to get_map())'''
    
    # Read input files
    try:
        ice_elements = pd.read_csv(validated_ice_elements, sep='\t')
        print(f"Loaded {len(ice_elements)} validated ICE elements")
    except Exception as e:
        print(f"Error reading validated ICE elements: {e}", file=sys.stderr)
        return
    
    macsyfinder_features = parse_macsyfinder_features(macsyfinder_gff3)
    trna_data = parse_trna_data(trna_annotations)
    dr_data = parse_direct_repeats(direct_repeats)
    orit_data = parse_orit_results(orit_results)
    
    integrated_info = {}
    characteristics_data = []
    
    for _, ice in ice_elements.iterrows():
        system_id = ice['system_id']
        contig = ice['contig']
        start = ice['start']
        end = ice['end']
        length = ice['length']
        
        # Calculate GC content
        gc_content = calculate_gc_content(genome_fasta, start, end)
        gc_positions, gc_profile = calculate_gc_profile(genome_fasta, start, end, 
                                                       gc_window_size, gc_step_size)
        
        # Get MacSyFinder features
        ice_features = macsyfinder_features.get(system_id, [])
        
        # Get associated tRNAs
        extended_seq_name = f"{system_id}_extended"
        associated_trnas = trna_data.get(extended_seq_name, [])
        
        # Get direct repeats
        ice_drs = [dr for dr in dr_data if dr['sequence_id'] == extended_seq_name]
        
        # Get oriT information
        ice_seq_name = f"{system_id}_ice"
        orit_info = orit_data.get(ice_seq_name, {})
        
        # Determine flanking tRNA
        closest_trna = None
        if associated_trnas:
            # Find closest tRNA to ICE boundaries
            min_distance = float('inf')
            for trna in associated_trnas:
                distance = min(abs(trna['start'] - start), abs(trna['end'] - end))
                if distance < min_distance:
                    min_distance = distance
                    closest_trna = trna
        
        # Determine direct repeat boundaries
        dr_boundaries = None
        if ice_drs:
            # Use the best direct repeat (first one, as they're sorted by quality)
            best_dr = ice_drs[0]
            dr_boundaries = {
                'attL_start': best_dr['start_pos1'],
                'attL_end': best_dr['end_pos1'],
                'attR_start': best_dr['start_pos2'],
                'attR_end': best_dr['end_pos2'],
                'repeat_length': best_dr['repeat_length']
            }
        
        # Compile integrated information
        integrated_info[system_id] = {
            'basic_info': {
                'system_id': system_id,
                'contig': contig,
                'start': start,
                'end': end,
                'length': length,
                'gc_content': gc_content
            },
            'features': ice_features,
            'trna_info': {
                'associated_trnas': associated_trnas,
                'closest_trna': closest_trna,
                'trna_count': len(associated_trnas)
            },
            'direct_repeats': {
                'boundaries': dr_boundaries,
                'all_repeats': ice_drs
            },
            'orit_info': orit_info,
            'gc_profile': {
                'positions': gc_positions,
                'values': gc_profile
            }
        }
        
        # Prepare characteristics data for TSV output
        characteristics_data.append({
            'system_id': system_id,
            'contig': contig,
            'start': start,
            'end': end,
            'length': length,
            'gc_content': gc_content,
            'num_features': len(ice_features),
            'num_trnas': len(associated_trnas),
            'num_direct_repeats': len(ice_drs),
            'has_orit': bool(orit_info),
            'orit_identity': orit_info.get('identity_percent', 0),
            'orit_coverage': orit_info.get('coverage', 0),
            'closest_trna_type': closest_trna['type'] if closest_trna else 'None',
            'dr_repeat_length': dr_boundaries['repeat_length'] if dr_boundaries else 0
        })
    
    # Write integrated information to JSON
    try:
        with open(f'{prefix}_integrated_ice_info.json', 'w') as f:
            json.dump(integrated_info, f, indent=2)
    except Exception as e:
        print(f"Error writing JSON file: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Write characteristics to TSV
    try:
        if characteristics_data:
            df_char = pd.DataFrame(characteristics_data)
            df_char.to_csv(f'{prefix}_ice_characteristics.tsv', sep='\t', index=False)
        else:
            # Create empty file with headers
            with open(f'{prefix}_ice_characteristics.tsv', 'w') as f:
                f.write('system_id\tcontig\tstart\tend\tlength\tgc_content\tnum_features\tnum_trnas\tnum_direct_repeats\thas_orit\torit_identity\torit_coverage\tclosest_trna_type\tdr_repeat_length\n')
    except Exception as e:
        print(f"Error writing TSV file: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Successfully integrated information for {len(integrated_info)} ICE systems")

def main():
    parser = argparse.ArgumentParser(description='Integrate ICE information from multiple sources')
    parser.add_argument('validated_ice_elements', help='Input validated ICE elements TSV file')
    parser.add_argument('macsyfinder_gff3', help='Input MacSyFinder GFF3 file')
    parser.add_argument('trna_annotations', help='Input tRNA annotations file')
    parser.add_argument('direct_repeats', help='Input direct repeats TSV file')
    parser.add_argument('orit_results', help='Input oriT results TSV file')
    parser.add_argument('genome_fasta', help='Input genome FASTA file')
    parser.add_argument('prefix', help='Output file prefix')
    parser.add_argument('--gc-window-size', type=int, default=1000,
                       help='GC content window size (default: 1000)')
    parser.add_argument('--gc-step-size', type=int, default=500,
                       help='GC content step size (default: 500)')
    
    args = parser.parse_args()
    
    try:
        integrate_ice_information(
            args.validated_ice_elements,
            args.macsyfinder_gff3,
            args.trna_annotations,
            args.direct_repeats,
            args.orit_results,
            args.genome_fasta,
            args.prefix,
            args.gc_window_size,
            args.gc_step_size
        )
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
