#!/usr/bin/env python3

import json
import pandas as pd
import sys
import argparse

def assess_ice_quality(classified_elements, integrated_info, 
                      min_size=10000, 
                      max_size=500000,
                      min_trna_count=1,
                      min_confidence_score=0.5,
                      require_integrase=True,
                      require_mobilization=False):
    '''Assess quality of ICE elements based on multiple criteria'''
    
    quality_criteria = {
        'size_check': 'ICE size within acceptable range',
        'integrase_check': 'Presence of integrase gene',
        'mobilization_check': 'Presence of mobilization machinery',
        'trna_check': 'Sufficient tRNA associations',
        'confidence_check': 'Minimum confidence score met',
        'structural_check': 'Structural integrity (boundaries, repeats)',
        'completeness_check': 'Feature completeness assessment'
    }
    
    passed_elements = []
    failed_elements = []
    quality_stats = {
        'total_elements': len(classified_elements),
        'passed_elements': 0,
        'failed_elements': 0,
        'failure_reasons': {},
        'quality_distribution': {}
    }
    
    assessment_log = []
    assessment_log.append("ICE Quality Assessment Report")
    assessment_log.append("=" * 40)
    assessment_log.append("")
    assessment_log.append("Quality Criteria:")
    for criterion, description in quality_criteria.items():
        assessment_log.append(f"  {criterion}: {description}")
    assessment_log.append("")
    
    for _, ice in classified_elements.iterrows():
        system_id = ice['system_id']
        assessment_log.append(f"Assessing ICE: {system_id}")
        
        # Initialize quality assessment
        quality_score = 0.0
        max_possible_score = 0.0
        passed_criteria = []
        failed_criteria = []
        failure_reasons = []
        
        # Criterion 1: Size check (weight: 1.0)
        max_possible_score += 1.0
        if min_size <= ice['length'] <= max_size:
            quality_score += 1.0
            passed_criteria.append('size_check')
            assessment_log.append(f"  ✓ Size check: {ice['length']} bp (within {min_size}-{max_size} bp)")
        else:
            failed_criteria.append('size_check')
            failure_reasons.append(f"Size out of range: {ice['length']} bp")
            assessment_log.append(f"  ✗ Size check: {ice['length']} bp (outside {min_size}-{max_size} bp)")
        
        # Criterion 2: Integrase check (weight: 2.0)
        max_possible_score += 2.0
        if ice['num_integrase'] > 0:
            quality_score += 2.0
            passed_criteria.append('integrase_check')
            assessment_log.append(f"  ✓ Integrase check: {ice['num_integrase']} integrase(s) found")
        else:
            failed_criteria.append('integrase_check')
            if require_integrase:
                failure_reasons.append("No integrase gene found (required)")
            assessment_log.append(f"  ✗ Integrase check: No integrase genes found")
        
        # Criterion 3: Mobilization check (weight: 1.5)
        max_possible_score += 1.5
        has_mobilization = ice['num_relaxase'] > 0 or ice['num_t4ss'] > 0 or ice['num_t4cp'] > 0
        if has_mobilization:
            quality_score += 1.5
            passed_criteria.append('mobilization_check')
            assessment_log.append(f"  ✓ Mobilization check: Relaxase({ice['num_relaxase']}), T4SS({ice['num_t4ss']}), T4CP({ice['num_t4cp']})")
        else:
            failed_criteria.append('mobilization_check')
            if require_mobilization:
                failure_reasons.append("No mobilization machinery found (required)")
            assessment_log.append(f"  ✗ Mobilization check: No mobilization genes found")
        
        # Criterion 4: tRNA check (weight: 1.0)
        max_possible_score += 1.0
        if ice['num_trnas'] >= min_trna_count:
            quality_score += 1.0
            passed_criteria.append('trna_check')
            assessment_log.append(f"  ✓ tRNA check: {ice['num_trnas']} tRNAs (≥{min_trna_count} required)")
        else:
            failed_criteria.append('trna_check')
            failure_reasons.append(f"Insufficient tRNAs: {ice['num_trnas']} < {min_trna_count}")
            assessment_log.append(f"  ✗ tRNA check: {ice['num_trnas']} tRNAs (<{min_trna_count} required)")
        
        # Criterion 5: Confidence check (weight: 1.0)
        max_possible_score += 1.0
        confidence_mapping = {'high': 1.0, 'medium': 0.7, 'low': 0.4, 'very_low': 0.1, 'failed': 0.0}
        confidence_score = confidence_mapping.get(ice['confidence'], 0.0)
        
        if confidence_score >= min_confidence_score:
            quality_score += confidence_score
            passed_criteria.append('confidence_check')
            assessment_log.append(f"  ✓ Confidence check: {ice['confidence']} (score: {confidence_score})")
        else:
            failed_criteria.append('confidence_check')
            failure_reasons.append(f"Low confidence: {ice['confidence']} (score: {confidence_score})")
            assessment_log.append(f"  ✗ Confidence check: {ice['confidence']} (score: {confidence_score} < {min_confidence_score})")
        
        # Criterion 6: Structural check (weight: 0.5)
        max_possible_score += 0.5
        has_orit = ice['has_orit']
        if has_orit or ice['num_trnas'] > 0:
            quality_score += 0.5
            passed_criteria.append('structural_check')
            assessment_log.append(f"  ✓ Structural check: oriT({has_orit}) or tRNAs({ice['num_trnas']})")
        else:
            failed_criteria.append('structural_check')
            assessment_log.append(f"  ✗ Structural check: No oriT or tRNA evidence")
        
        # Criterion 7: Completeness check (weight: 0.5)
        max_possible_score += 0.5
        total_features = (ice['num_integrase'] + ice['num_relaxase'] + 
                        ice['num_t4ss'] + ice['num_t4cp'] + 
                        ice['num_replication'] + ice['num_transfer'])
        if total_features >= 3:  # Minimum feature count for completeness
            quality_score += 0.5
            passed_criteria.append('completeness_check')
            assessment_log.append(f"  ✓ Completeness check: {total_features} total features")
        else:
            failed_criteria.append('completeness_check')
            assessment_log.append(f"  ✗ Completeness check: {total_features} features (minimum 3)")
        
        # Calculate final quality score
        final_quality_score = quality_score / max_possible_score if max_possible_score > 0 else 0.0
        
        # Determine pass/fail status
        is_high_quality = (
            final_quality_score >= min_confidence_score and
            (not require_integrase or ice['num_integrase'] > 0) and
            (not require_mobilization or has_mobilization) and
            ice['num_trnas'] >= min_trna_count and
            min_size <= ice['length'] <= max_size
        )
        
        # Compile assessment results
        assessment_result = {
            'system_id': system_id,
            'contig': ice['contig'],
            'start': ice['start'],
            'end': ice['end'],
            'length': ice['length'],
            'ice_type': ice['ice_type'],
            'confidence': ice['confidence'],
            'quality_score': round(final_quality_score, 3),
            'passed_criteria': ';'.join(passed_criteria),
            'failed_criteria': ';'.join(failed_criteria),
            'failure_reasons': ';'.join(failure_reasons),
            'is_high_quality': is_high_quality,
            'assessment_status': 'PASS' if is_high_quality else 'FAIL',
            'num_integrase': ice['num_integrase'],
            'num_relaxase': ice['num_relaxase'],
            'num_t4ss': ice['num_t4ss'],
            'num_t4cp': ice['num_t4cp'],
            'num_trnas': ice['num_trnas'],
            'has_orit': ice['has_orit'],
            'orit_identity': ice['orit_identity'],
            'mobilization_types': ice['mobilization_types'],
            't4ss_types': ice['t4ss_types']
        }
        
        if is_high_quality:
            passed_elements.append(assessment_result)
            quality_stats['passed_elements'] += 1
            assessment_log.append(f"  → PASSED (Quality score: {final_quality_score:.3f})")
        else:
            failed_elements.append(assessment_result)
            quality_stats['failed_elements'] += 1
            assessment_log.append(f"  → FAILED (Quality score: {final_quality_score:.3f})")
            
            # Track failure reasons
            for reason in failure_reasons:
                quality_stats['failure_reasons'][reason] = quality_stats['failure_reasons'].get(reason, 0) + 1
        
        # Track quality distribution
        quality_bin = 'high' if final_quality_score >= 0.8 else 'medium' if final_quality_score >= 0.5 else 'low'
        quality_stats['quality_distribution'][quality_bin] = quality_stats['quality_distribution'].get(quality_bin, 0) + 1
        
        assessment_log.append("")
    
    return passed_elements, failed_elements, quality_stats, assessment_log

def main():
    parser = argparse.ArgumentParser(description='Assess quality of ICE elements based on multiple criteria')
    parser.add_argument('classified_ice_elements', help='Input classified ICE elements TSV file')
    parser.add_argument('integrated_ice_info', help='Input integrated ICE information JSON file')
    parser.add_argument('prefix', help='Output file prefix')
    parser.add_argument('--min-size', type=int, default=10000,
                       help='Minimum ICE size (default: 10000)')
    parser.add_argument('--max-size', type=int, default=500000,
                       help='Maximum ICE size (default: 500000)')
    parser.add_argument('--min-trna-count', type=int, default=1,
                       help='Minimum tRNA count (default: 1)')
    parser.add_argument('--min-confidence-score', type=float, default=0.5,
                       help='Minimum confidence score (default: 0.5)')
    parser.add_argument('--require-integrase', action='store_true', default=True,
                       help='Require integrase gene (default: True)')
    parser.add_argument('--require-mobilization', action='store_true', default=False,
                       help='Require mobilization machinery (default: False)')
    
    args = parser.parse_args()
    
    try:
        # Read input files
        classified_elements = pd.read_csv(args.classified_ice_elements, sep='\t')
        print(f"Loaded {len(classified_elements)} classified ICE elements")
        
        with open(args.integrated_ice_info, 'r') as f:
            integrated_info = json.load(f)
        print(f"Loaded integrated information for {len(integrated_info)} ICE systems")
        
        # Perform quality assessment
        passed_elements, failed_elements, quality_stats, assessment_log = assess_ice_quality(
            classified_elements, integrated_info,
            args.min_size, args.max_size, args.min_trna_count,
            args.min_confidence_score, args.require_integrase, args.require_mobilization
        )
        
        # Write passed elements
        if passed_elements:
            df_passed = pd.DataFrame(passed_elements)
            df_passed.to_csv(f'{args.prefix}_quality_assessed_ice.tsv', sep='\t', index=False)
        else:
            # Create empty file with headers
            headers = ['system_id', 'contig', 'start', 'end', 'length', 'ice_type', 'confidence',
                      'quality_score', 'passed_criteria', 'failed_criteria', 'failure_reasons',
                      'is_high_quality', 'assessment_status', 'num_integrase', 'num_relaxase',
                      'num_t4ss', 'num_t4cp', 'num_trnas', 'has_orit', 'orit_identity',
                      'mobilization_types', 't4ss_types']
            with open(f'{args.prefix}_quality_assessed_ice.tsv', 'w') as f:
                f.write('\t'.join(headers) + '\n')
        
        # Write failed elements
        if failed_elements:
            df_failed = pd.DataFrame(failed_elements)
            df_failed.to_csv(f'{args.prefix}_failed_ice_elements.tsv', sep='\t', index=False)
        else:
            # Create empty file with headers
            headers = ['system_id', 'contig', 'start', 'end', 'length', 'ice_type', 'confidence',
                      'quality_score', 'passed_criteria', 'failed_criteria', 'failure_reasons',
                      'is_high_quality', 'assessment_status', 'num_integrase', 'num_relaxase',
                      'num_t4ss', 'num_t4cp', 'num_trnas', 'has_orit', 'orit_identity',
                      'mobilization_types', 't4ss_types']
            with open(f'{args.prefix}_failed_ice_elements.tsv', 'w') as f:
                f.write('\t'.join(headers) + '\n')
        
        # Generate quality report
        report_lines = assessment_log + [
            "",
            "Quality Assessment Summary",
            "=" * 30,
            f"Total ICE elements assessed: {quality_stats['total_elements']}",
            f"Passed quality assessment: {quality_stats['passed_elements']}",
            f"Failed quality assessment: {quality_stats['failed_elements']}",
            f"Pass rate: {quality_stats['passed_elements']/quality_stats['total_elements']*100:.1f}%" if quality_stats['total_elements'] > 0 else "Pass rate: 0%",
            "",
            "Quality Distribution:",
        ]
        
        for quality_level, count in quality_stats['quality_distribution'].items():
            report_lines.append(f"  {quality_level}: {count}")
        
        if quality_stats['failure_reasons']:
            report_lines.extend([
                "",
                "Common Failure Reasons:",
            ])
            for reason, count in quality_stats['failure_reasons'].items():
                report_lines.append(f"  {reason}: {count}")
        
        # Write quality report
        with open(f'{args.prefix}_quality_report.txt', 'w') as f:
            f.write('\n'.join(report_lines))
        
        print(f"Quality assessment completed:")
        print(f"  Passed: {quality_stats['passed_elements']}")
        print(f"  Failed: {quality_stats['failed_elements']}")
        print(f"  Total: {quality_stats['total_elements']}")
        
    except Exception as e:
        print(f"Error: {e}", file=sys.stderr)
        sys.exit(1)

if __name__ == "__main__":
    main()
