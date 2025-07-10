#!/usr/bin/env python3

import os
import sys
import argparse
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
import pandas as pd

def run_blast_search(query_file, database, output_file, evalue='0.01', word_size='11'):
    '''Run BLAST search against oriT database'''
    
    blast_cmd = [
        'blastn',
        '-db', database,
        '-query', query_file,
        '-evalue', evalue,
        '-word_size', word_size,
        '-outfmt', '6 std qlen slen',
        '-num_alignments', '1',
        '-out', output_file
    ]
    
    try:
        result = subprocess.run(blast_cmd, capture_output=True, text=True, check=True)
        print(f"BLAST search completed successfully")
        return True
    except subprocess.CalledProcessError as e:
        print(f"BLAST search failed: {e}", file=sys.stderr)
        print(f"STDERR: {e.stderr}", file=sys.stderr)
        return False

def parse_blast_results(blast_file, identity_threshold=49.0, coverage_threshold=0.49):
    '''Parse BLAST results and identify significant oriT matches'''
    
    orit_results = []
    
    if not os.path.exists(blast_file) or os.path.getsize(blast_file) == 0:
        print("No BLAST results found or file is empty")
        return orit_results
    
    try:
        with open(blast_file, 'r') as f:
            for line in f:
                if line.strip():
                    fields = line.strip().split('\t')
                    if len(fields) >= 14:
                        query_id = fields[0]
                        subject_id = fields[1]
                        identity = float(fields[2])
                        alignment_length = int(fields[3])
                        mismatches = int(fields[4])
                        gap_opens = int(fields[5])
                        query_start = int(fields[6])
                        query_end = int(fields[7])
                        subject_start = int(fields[8])
                        subject_end = int(fields[9])
                        evalue = float(fields[10])
                        bit_score = float(fields[11])
                        query_length = int(fields[12])
                        subject_length = int(fields[13])
                        
                        # Calculate coverage and combined score (H-value from original script)
                        coverage = alignment_length / subject_length
                        h_value = coverage * identity
                        
                        # Apply thresholds (based on original script: hvalue > 0.49)
                        if h_value > coverage_threshold:
                            orit_results.append({
                                'ice_id': query_id,
                                'orit_id': subject_id,
                                'identity_percent': identity,
                                'alignment_length': alignment_length,
                                'mismatches': mismatches,
                                'gap_opens': gap_opens,
                                'query_start': query_start,
                                'query_end': query_end,
                                'subject_start': subject_start,
                                'subject_end': subject_end,
                                'evalue': evalue,
                                'bit_score': bit_score,
                                'query_length': query_length,
                                'subject_length': subject_length,
                                'coverage': coverage,
                                'h_value': h_value,
                                'orit_detected': 'Yes'
                            })
    except Exception as e:
        print(f"Error parsing BLAST results: {e}", file=sys.stderr)
    
    return orit_results

def extract_orit_sequences(ice_sequences_file, orit_results):
    '''Extract oriT sequences from ICE sequences based on BLAST results'''
    
    orit_sequences = []
    
    if not orit_results:
        return orit_sequences
    
    try:
        # Parse ICE sequences
        ice_records = {}
        for record in SeqIO.parse(ice_sequences_file, "fasta"):
            ice_records[record.id] = record
        
        for result in orit_results:
            ice_id = result['ice_id']
            query_start = result['query_start']
            query_end = result['query_end']
            
            if ice_id in ice_records:
                ice_record = ice_records[ice_id]
                
                # Extract oriT sequence (adjust for 0-based indexing)
                orit_start = min(query_start, query_end) - 1
                orit_end = max(query_start, query_end)
                orit_seq = ice_record.seq[orit_start:orit_end]
                
                # Create oriT sequence record
                orit_record = SeqRecord(
                    orit_seq,
                    id=f"{ice_id}_oriT",
                    description=f"oriT sequence from {ice_id} ({query_start}-{query_end}) | Match: {result['orit_id']} | Identity: {result['identity_percent']:.1f}% | H-value: {result['h_value']:.3f}"
                )
                
                orit_sequences.append(orit_record)
    except Exception as e:
        print(f"Error extracting oriT sequences: {e}", file=sys.stderr)
    
    return orit_sequences

def main():
    parser = argparse.ArgumentParser(description='Detect oriT sequences in ICE elements using BLAST')
    parser.add_argument('ice_sequences', help='Input ICE sequences FASTA file')
    parser.add_argument('orit_database', help='oriT reference database (BLAST format)')
    parser.add_argument('prefix', help='Output file prefix')
    parser.add_argument('--evalue', default='0.01', help='BLAST E-value threshold (default: 0.01)')
    parser.add_argument('--word-size', default='11', help='BLAST word size (default: 11)')
    parser.add_argument('--identity-threshold', type=float, default=49.0, 
                       help='Identity threshold (default: 49.0)')
    parser.add_argument('--coverage-threshold', type=float, default=0.49,
                       help='Coverage threshold for H-value (default: 0.49)')
    
    args = parser.parse_args()
    
    blast_output = f"{args.prefix}_blast_output.txt"
    results_file = f"{args.prefix}_orit_results.tsv"
    sequences_file = f"{args.prefix}_orit_sequences.fasta"
    
    print(f"Starting oriT detection for: {args.ice_sequences}")
    print(f"Using oriT database: {args.orit_database}")
    
    # Run BLAST search
    blast_success = run_blast_search(args.ice_sequences, args.orit_database, blast_output, 
                                   args.evalue, args.word_size)
    
    if not blast_success:
        print("BLAST search failed, creating empty results")
        # Create empty results file
        with open(results_file, 'w') as f:
            f.write('ice_id\torit_id\tidentity_percent\talignment_length\tmismatches\tgap_opens\tquery_start\tquery_end\tsubject_start\tsubject_end\tevalue\tbit_score\tquery_length\tsubject_length\tcoverage\th_value\torit_detected\n')
        
        # Create empty BLAST output if it doesn't exist
        if not os.path.exists(blast_output):
            with open(blast_output, 'w') as f:
                f.write('# No BLAST results - search failed\n')
        
        return
    
    # Parse BLAST results
    orit_results = parse_blast_results(blast_output, args.identity_threshold, args.coverage_threshold)
    
    print(f"Found {len(orit_results)} significant oriT matches")
    
    # Write results to TSV file
    try:
        if orit_results:
            df = pd.DataFrame(orit_results)
            df.to_csv(results_file, sep='\t', index=False)
            
            # Extract oriT sequences
            orit_sequences = extract_orit_sequences(args.ice_sequences, orit_results)
            
            if orit_sequences:
                SeqIO.write(orit_sequences, sequences_file, "fasta")
                print(f"Extracted {len(orit_sequences)} oriT sequences")
            else:
                print("No oriT sequences could be extracted")
        else:
            # Create empty results file with headers
            with open(results_file, 'w') as f:
                f.write('ice_id\torit_id\tidentity_percent\talignment_length\tmismatches\tgap_opens\tquery_start\tquery_end\tsubject_start\tsubject_end\tevalue\tbit_score\tquery_length\tsubject_length\tcoverage\th_value\torit_detected\n')
            print("No significant oriT matches found")
    except Exception as e:
        print(f"Error writing results: {e}", file=sys.stderr)
        sys.exit(1)
    
    print("oriT detection completed")

if __name__ == "__main__":
    main()
