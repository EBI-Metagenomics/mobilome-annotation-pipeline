#!/usr/bin/env python3

import json
import pandas as pd
import sys
import argparse

def classify_ice_features(features):
    '''Classify ICE features based on MacSyFinder results'''
    
    feature_categories = {
        'integrase': ['Phage_integrase', 'UPF0236', 'Recombinase', 'rve', 
                     'TIGR02224', 'TIGR02249', 'TIGR02225', 'PB001819'],
        'relaxase': ['Relaxase_', 'T4SS_MOB'],
        't4ss': ['FATA_', 'FA_', 'T4SS_'],
        't4cp': ['t4cp', 'tcpA'],
        'replication': ['RepSAv2', 'DUF3631', 'Prim-Pol'],
        'transfer': ['FtsK_SpoIIIE']
    }
    
    classified_features = {
        'integrase': [],
        'relaxase': [],
        't4ss': [],
        't4cp': [],
        'replication': [],
        'transfer': [],
        'other': []
    }
    
    for feature in features:
        gene_name = feature['gene_name']
        classified = False
        
        for category, patterns in feature_categories.items():
            for pattern in patterns:
                if pattern in gene_name:
                    classified_features[category].append({
                        'gene_name': gene_name,
                        'start': feature['start'],
                        'end': feature['end'],
                        'strand': feature['strand']
                    })
                    classified = True
                    break
            if classified:
                break
        
        if not classified:
            classified_features['other'].append({
                'gene_name': gene_name,
                'start': feature['start'],
                'end': feature['end'],
                'strand': feature['strand']
            })
    
    return classified_features

def determine_ice_type(classified_features, orit_info, trna_info):
    '''Determine ICE type based on feature composition'''
    
    has_integrase = len(classified_features['integrase']) > 0
    has_relaxase = len(classified_features['relaxase']) > 0
    has_t4ss = len(classified_features['t4ss']) > 0
    has_t4cp = len(classified_features['t4cp']) > 0
    has_orit = bool(orit_info)
    has_trna = trna_info['trna_count'] > 0
    
    # Classification logic based on icefinder2 script
    if has_integrase and has_relaxase and (has_t4ss or has_t4cp):
        if has_orit:
            ice_type = 'T4SS-type ICE'
            confidence = 'high'
        else:
            ice_type = 'T4SS-type ICE'
            confidence = 'medium'
    elif has_integrase and has_relaxase:
        ice_type = 'IME'  # Integrative and Mobilizable Element
        confidence = 'medium' if has_trna else 'low'
    elif has_integrase and not has_relaxase:
        ice_type = 'Integrative Element'
        confidence = 'low'
    elif has_relaxase and not has_integrase:
        ice_type = 'Mobilizable Element'
        confidence = 'low'
    else:
        ice_type = 'Unknown'
        confidence = 'very_low'
    
    # Check for AICE (Antibiotic resistance ICE) - would need additional AR gene detection
    # For now, we'll classify based on the presence of specific patterns
    if 'AICE' in str(classified_features):
        ice_type = 'AICE'
    
    return ice_type, confidence

def determine_mobilization_type(relaxase_features):
    '''Determine mobilization type based on relaxase genes'''
    
    mob_types = []
    for feature in relaxase_features:
        gene_name = feature['gene_name']
        
        # Extract MOB type from gene name
        if 'MOB' in gene_name:
            parts = gene_name.split('_')
            for part in parts:
                if part.startswith('MOB'):
                    mob_types.append(part)
        elif 'Relaxase_' in gene_name:
            mob_type = gene_name.replace('Relaxase_', '')
            mob_types.append(mob_type)
    
    return list(set(mob_types))  # Remove duplicates

def determine_t4ss_type(t4ss_features, t4cp_features):
    '''Determine T4SS (Type IV Secretion System) type'''
    
    t4ss_types = []
    
    # Analyze T4SS genes
    for feature in t4ss_features:
        gene_name = feature['gene_name']
        if '_' in gene_name:
            parts = gene_name.split('_')
            if len(parts) > 1:
                t4ss_types.append(parts[-1])
    
    # Analyze T4CP genes
    for feature in t4cp_features:
        gene_name = feature['gene_name']
        if '_' in gene_name:
            parts = gene_name.split('_')
            if len(parts) > 1:
                t4ss_types.append(parts[-1])
    
    return list(set(t4ss_types))  # Remove duplicates

def classify_ice_elements(integrated_ice_info, prefix):
    '''Main function for ICE classification and typing'''
    
    # Read integrated ICE information
    try:
        with open(integrated_ice_info, 'r') as f:
            ice_data = json.load(f)
        print(f"Loaded integrated information for {len(ice_data)} ICE systems")
    except Exception as e:
        print(f"Error reading integrated ICE information: {e}", file=sys.stderr)
        return
    
    classified_elements = []
    typing_log = []
    
    typing_log.append("ICE Classification and Typing Report")
    typing_log.append("=" * 50)
    typing_log.append("")
    
    for system_id, ice_info in ice_data.items():
        typing_log.append(f"Processing ICE system: {system_id}")
        
        # Extract information
        basic_info = ice_info['basic_info']
        features = ice_info['features']
        orit_info = ice_info['orit_info']
        trna_info = ice_info['trna_info']
        
        # Classify features
        classified_features = classify_ice_features(features)
        
        # Determine ICE type
        ice_type, confidence = determine_ice_type(classified_features, orit_info, trna_info)
        
        # Determine mobilization type
        mob_types = determine_mobilization_type(classified_features['relaxase'])
        
        # Determine T4SS type
        t4ss_types = determine_t4ss_type(classified_features['t4ss'], classified_features['t4cp'])
        
        # Count features by category
        feature_counts = {category: len(features) for category, features in classified_features.items()}
        
        typing_log.append(f"  ICE Type: {ice_type} (confidence: {confidence})")
        typing_log.append(f"  Feature counts: {feature_counts}")
        typing_log.append(f"  Mobilization types: {mob_types}")
        typing_log.append(f"  T4SS types: {t4ss_types}")
        typing_log.append(f"  Has oriT: {bool(orit_info)}")
        typing_log.append(f"  Associated tRNAs: {trna_info['trna_count']}")
        typing_log.append("")
        
        # Compile classification results
        classified_elements.append({
            'system_id': system_id,
            'contig': basic_info['contig'],
            'start': basic_info['start'],
            'end': basic_info['end'],
            'length': basic_info['length'],
            'ice_type': ice_type,
            'confidence': confidence,
            'num_integrase': feature_counts['integrase'],
            'num_relaxase': feature_counts['relaxase'],
            'num_t4ss': feature_counts['t4ss'],
            'num_t4cp': feature_counts['t4cp'],
            'num_replication': feature_counts['replication'],
            'num_transfer': feature_counts['transfer'],
            'num_other': feature_counts['other'],
            'mobilization_types': ';'.join(mob_types) if mob_types else 'None',
            't4ss_types': ';'.join(t4ss_types) if t4ss_types else 'None',
            'has_orit': bool(orit_info),
            'orit_identity': orit_info.get('identity_percent', 0) if orit_info else 0,
            'num_trnas': trna_info['trna_count'],
            'closest_trna': trna_info['closest_trna']['type'] if trna_info['closest_trna'] else 'None'
        })
    
    # Write classification results
    try:
        if classified_elements:
            df = pd.DataFrame(classified_elements)
            df.to_csv(f'{prefix}_classified_ice_elements.tsv', sep='\t', index=False)
        else:
            # Create empty file with headers
            headers = ['system_id', 'contig', 'start', 'end', 'length', 'ice_type', 'confidence',
                      'num_integrase', 'num_relaxase', 'num_t4ss', 'num_t4cp', 'num_replication',
                      'num_transfer', 'num_other', 'mobilization_types', 't4ss_types', 'has_orit',
                      'orit_identity', 'num_trnas', 'closest_trna']
            with open(f'{prefix}_classified_ice_elements.tsv', 'w') as f:
                f.write('\t'.join(headers) + '\n')
    except Exception as e:
        print(f"Error writing classification results: {e}", file=sys.stderr)
        sys.exit(1)
    
    # Write typing report
    try:
        with open(f'{prefix}_ice_typing_report.txt', 'w') as f:
            f.write('\n'.join(typing_log))
    except Exception as e:
        print(f"Error writing typing report: {e}", file=sys.stderr)
        sys.exit(1)
    
    print(f"Successfully classified {len(classified_elements)} ICE elements")
    
    # Summary statistics
    type_counts = {}
    for element in classified_elements:
        ice_type = element['ice_type']
        type_counts[ice_type] = type_counts.get(ice_type, 0) + 1
    
    print("ICE type distribution:")
    for ice_type, count in type_counts.items():
        print(f"  {ice_type}: {count}")

def main():
    parser = argparse.ArgumentParser(description='Classify and type ICE elements based on integrated information')
    parser.add_argument('integrated_ice_info', help='Input integrated ICE information JSON file')
    parser.add_argument('prefix', help='Output file prefix')
    
    args = parser.parse_args()
    
    try:
        classify_ice_elements(args.integrated_ice_info, args.prefix)
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
