#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import sys
import json
import argparse
from collections import defaultdict
from Bio import SeqIO
from Bio.SeqUtils import gc_fraction
import pandas as pd

class MacsyfinderProcessor:
    """
    Process macsyfinder output for boundary delineation analysis
    """
    
    def __init__(self, all_systems_file, genome_file, output_dir, annotation_file=None):
        self.all_systems_file = all_systems_file
        self.genome_file = genome_file
        self.output_dir = output_dir
        self.annotation_file = annotation_file
        
        # Create output directory if it doesn't exist
        os.makedirs(output_dir, exist_ok=True)
        
        # Initialize data structures
        self.systems_data = []
        self.filtered_systems = {}
        self.boundary_data = {}
        self.annotation_dict = {}
        self.genome_records = {}  # Store all contigs
        
        # Load genome records
        self.load_genome_records()
        
    def load_genome_records(self):
        """Load all genome records (contigs) into memory"""
        try:
            print("Loading genome records...")
            for record in SeqIO.parse(self.genome_file, "fasta"):
                self.genome_records[record.id] = record
            print(f"Loaded {len(self.genome_records)} contigs")
        except Exception as e:
            print(f"Error loading genome records: {e}")
            sys.exit(1)
    
    def load_annotation_data(self):
        """Load annotation data from GFF/GBK file if provided"""
        if not self.annotation_file:
            return
            
        if self.annotation_file.endswith('.gff'):
            self._parse_gff()
        elif self.annotation_file.endswith('.gbk') or self.annotation_file.endswith('.gb'):
            self._parse_genbank()
    
    def _parse_gff(self):
        """Parse GFF annotation file"""
        with open(self.annotation_file, 'r') as f:
            for line in f:
                if line.startswith('#') or not line.strip():
                    continue
                    
                parts = line.strip().split('\t')
                if len(parts) >= 9:
                    contig_id = parts[0]  # First column is contig/sequence ID
                    feature_type = parts[2]
                    start = int(parts[3])
                    end = int(parts[4])
                    strand = parts[6]
                    attributes = parts[8]
                    
                    # Extract gene ID and product
                    gene_id = self._extract_attribute(attributes, 'ID')
                    product = self._extract_attribute(attributes, 'product', 'hypothetical protein')
                    
                    if gene_id:
                        self.annotation_dict[gene_id] = {
                            'contig': contig_id,
                            'start': start,
                            'end': end,
                            'strand': strand,
                            'type': feature_type,
                            'product': product
                        }
    
    def _parse_genbank(self):
        """Parse GenBank annotation file"""
        for record in SeqIO.parse(self.annotation_file, "genbank"):
            contig_id = record.id
            for feature in record.features:
                if feature.type in ['CDS', 'tRNA', 'rRNA']:
                    start = int(feature.location.start) + 1  # Convert to 1-based
                    end = int(feature.location.end)
                    strand = '+' if feature.location.strand == 1 else '-'
                    
                    # Extract gene information
                    gene_id = feature.qualifiers.get('locus_tag', [''])[0]
                    product = feature.qualifiers.get('product', ['hypothetical protein'])[0]
                    
                    if gene_id:
                        self.annotation_dict[gene_id] = {
                            'contig': contig_id,
                            'start': start,
                            'end': end,
                            'strand': strand,
                            'type': feature.type,
                            'product': product
                        }
    
    def _extract_attribute(self, attributes, key, default=''):
        """Extract attribute value from GFF attributes string"""
        for attr in attributes.split(';'):
            if '=' in attr:
                k, v = attr.split('=', 1)
                if k.strip() == key:
                    return v.strip()
        return default
   

    def parse_macsyfinder_output(self):
        """Parse the all_systems.tsv file from macsyfinder and group by contig"""
        try:
            # Read the file, skipping comment lines that start with #
            with open(self.all_systems_file, 'r') as f:
                lines = [line for line in f.readlines() if line.strip()]

            # Find the header line (first non-comment line)
            header_idx = 0
            for i, line in enumerate(lines):
                if not line.startswith('#') and line.strip():
                    header_idx = i
                    break
    
            # Read the data starting from the header
            df = pd.read_csv(self.all_systems_file, sep='\t', skiprows=header_idx)
    
            # Group by sys_id first, then by contig within each system
            system_groups = df.groupby('sys_id')
    
            for sys_id, system_group in system_groups:
                # Extract system type from model_fqn
                model_fqn = system_group.iloc[0]['model_fqn']
                system_type = model_fqn.split('/')[-1]
            
                # Group genes by contig within this system
                genes_with_contigs = []
                genes_without_contigs = []
            
                for _, row in system_group.iterrows():
                    gene_info = {
                        'gene_name': row['gene_name'],
                        'hit_id': row['hit_id'],
                        'hit_pos': row['hit_pos'],
                        'model_fqn': row['model_fqn'],
                        'hit_status': row['hit_status'],
                        'hit_score': row['hit_score'],
                        'profile_cov': row['hit_profile_cov'],
                        'seq_cov': row['hit_seq_cov']
                    }
                
                    # Add contig information if available from annotation
                    if row['hit_id'] in self.annotation_dict:
                        contig_id = self.annotation_dict[row['hit_id']]['contig']
                        gene_info['contig'] = contig_id
                        genes_with_contigs.append(gene_info)
                    else:
                        genes_without_contigs.append(gene_info)
            
                # Group genes by contig
                contig_groups = defaultdict(list)
                for gene in genes_with_contigs:
                    contig_groups[gene['contig']].append(gene)
            
                # Create separate systems for each contig
                for contig_id, contig_genes in contig_groups.items():
                    # Calculate boundaries for this contig's genes
                    start_pos = min(gene['hit_pos'] for gene in contig_genes)
                    end_pos = max(gene['hit_pos'] for gene in contig_genes)
                
                    contig_system_info = {
                        'original_system_name': sys_id,
                        'system_name': f"{sys_id}_{contig_id}",
                        'system_type': system_type,
                        'genes': contig_genes,
                        'start_pos': start_pos,
                        'end_pos': end_pos,
                        'replicon': system_group.iloc[0]['replicon'],
                        'contig': contig_id,
                        'wholeness': len(contig_genes) / len(system_group),  # Recalculate wholeness for this contig
                        'original_wholeness': system_group.iloc[0]['sys_wholeness']
                    }
                
                    self.systems_data.append(contig_system_info)
            
                # Handle genes without contig information (if any)
                if genes_without_contigs:
                    print(f"Warning: {len(genes_without_contigs)} genes from system {sys_id} have no contig information")
                    # You could optionally create a system for these genes too
        
        except Exception as e:
            print(f"Error parsing macsyfinder output: {e}")
            print(f"Please check the format of {self.all_systems_file}")
            sys.exit(1)


    def filter_systems(self):
        """Filter systems based on criteria similar to original ICE_filter function"""

        # Define system type mappings (same as ICEfinder2)
        integrase_genes = ['Phage_integrase', 'UPF0236', 'Recombinase', 'rve', 
                      'TIGR02224', 'TIGR02249', 'TIGR02225', 'PB001819']
        relaxase_genes = ['Relaxase_', 'T4SS_MOB']
        t4ss_genes = ['FATA_', 'FA_', 'T4SS_']

        for system in self.systems_data:
            # Extract gene names from the system
            system_genes = [gene['gene_name'] for gene in system['genes']]
    
            # Check for essential components
            has_integrase = any(gene in integrase_genes for gene in system_genes)
            has_relaxase = any(any(rg in gene for rg in relaxase_genes) for gene in system_genes)
            has_t4ss = any(any(t4g in gene for t4g in t4ss_genes) for gene in system_genes)
    
            # Use the system_type from macsyfinder output
            system_type = system['system_type']
    
            # Apply filtering criteria based on system type (following ICEfinder2 logic)
            if system_type == 'IME':
                if has_integrase and has_relaxase:
                    pass  # Keep IME systems with integrase and relaxase
                else:
                    continue  # Skip incomplete IMEs
            elif system_type == 'AICE':
                pass  # Keep all AICE systems
            elif system_type == 'ICE':
                if has_integrase:
                    pass  # Keep ICE systems with integrase
                else:
                    continue  # Skip ICE systems without integrase
            else:
                # For other system types, require integrase
                if not has_integrase:
                    continue
    
            # Filter based on wholeness - but now it's per-contig wholeness
            # You might want to adjust this threshold since it's now calculated per contig
            min_wholeness = 0.5  # Lower threshold since we're looking at contig fragments
            if system['wholeness'] < min_wholeness:
                continue
        
            # Filter based on minimum number of genes per contig
            min_genes_per_contig = 3  # Require at least 3 genes per contig
            if len(system['genes']) < min_genes_per_contig:
                continue
            
            system['type'] = system_type
    
            # Create unique system ID
            system_id = f"{system_type}_{system['contig']}_{len(self.filtered_systems)+1}"
            self.filtered_systems[system_id] = system
        
            print(f"Kept system: {system_id} with {len(system['genes'])} genes on contig {system['contig']}")


    def calculate_boundaries(self, flank_size=5):
        """Calculate system boundaries with flanking regions"""
        
        for system_id, system in self.filtered_systems.items():
            # Get gene positions from annotation if available
            gene_positions = []
            system_contig = None
            
            for gene in system['genes']:
                hit_id = gene['hit_id']
                if hit_id in self.annotation_dict:
                    ann = self.annotation_dict[hit_id]
                    gene_positions.append({
                        'gene_id': hit_id,
                        'start': ann['start'],
                        'end': ann['end'],
                        'strand': ann['strand'],
                        'product': ann['product'],
                        'contig': ann['contig']
                    })
                    
                    # Set system contig (assuming all genes are on the same contig)
                    if system_contig is None:
                        system_contig = ann['contig']
            
            if gene_positions and system_contig:
                # Sort by position
                gene_positions.sort(key=lambda x: x['start'])
                
                # Calculate system boundaries
                system_start = min(pos['start'] for pos in gene_positions)
                system_end = max(pos['end'] for pos in gene_positions)
                
                # Get contig length to avoid going beyond contig boundaries
                contig_length = len(self.genome_records[system_contig].seq) if system_contig in self.genome_records else float('inf')
                
                # Add flanking regions (respecting contig boundaries)
                flank_start = max(1, system_start - flank_size * 1000)  # 5kb flank
                flank_end = min(contig_length, system_end + flank_size * 1000)
                
                self.boundary_data[system_id] = {
                    'contig': system_contig,
                    'system_start': system_start,
                    'system_end': system_end,
                    'flank_start': flank_start,
                    'flank_end': flank_end,
                    'genes': gene_positions,
                    'type': system['type'],
                    'length': system_end - system_start + 1,
                    'contig_length': contig_length
                }
    
    def calculate_gc_content(self, contig_id, start, end):
        """Calculate GC content for a genomic region on a specific contig"""
        try:
            if contig_id not in self.genome_records:
                print(f"Warning: Contig {contig_id} not found in genome file")
                return 0.0
                
            record = self.genome_records[contig_id]
            sequence = record.seq[start-1:end]
            return round(gc_fraction(sequence), 2)
        except Exception as e:
            print(f"Error calculating GC content for {contig_id}:{start}-{end}: {e}")
            return 0.0
    
    def generate_boundary_output(self):
        """Generate boundary delineation output files"""
        
        # Summary file
        summary_data = []
        
        for system_id, boundary in self.boundary_data.items():
            gc_content = self.calculate_gc_content(
                boundary['contig'], 
                boundary['system_start'], 
                boundary['system_end']
            )
            
            summary_entry = {
                'system_id': system_id,
                'contig': boundary['contig'],
                'type': boundary['type'],
                'start': boundary['system_start'],
                'end': boundary['system_end'],
                'length': boundary['length'],
                'gc_content': gc_content,
                'num_genes': len(boundary['genes']),
                'flank_start': boundary['flank_start'],
                'flank_end': boundary['flank_end'],
                'contig_length': boundary['contig_length']
            }
            summary_data.append(summary_entry)
        
        # Write summary file
        summary_file = os.path.join(self.output_dir, 'boundary_summary.json')
        with open(summary_file, 'w') as f:
            json.dump(summary_data, f, indent=2)
        
        # Write detailed boundary files for each system
        for system_id, boundary in self.boundary_data.items():
            detail_file = os.path.join(self.output_dir, f'{system_id}_boundaries.json')
            with open(detail_file, 'w') as f:
                json.dump(boundary, f, indent=2, default=str)
        
        # Write TSV format for easy parsing
        tsv_file = os.path.join(self.output_dir, 'boundaries.tsv')
        with open(tsv_file, 'w') as f:
            f.write('system_id\tcontig\ttype\tstart\tend\tlength\tgc_content\tnum_genes\tflank_start\tflank_end\tcontig_length\n')
            for entry in summary_data:
                f.write(f"{entry['system_id']}\t{entry['contig']}\t{entry['type']}\t{entry['start']}\t"
                       f"{entry['end']}\t{entry['length']}\t{entry['gc_content']}\t"
                       f"{entry['num_genes']}\t{entry['flank_start']}\t{entry['flank_end']}\t{entry['contig_length']}\n")

    def extract_sequences(self):
        """Extract sequences for each identified system into single FASTA files"""
        try:
            # Single file for all system sequences (core systems only)
            systems_file = os.path.join(self.output_dir, 'all_systems.fasta')
            
            # Single file for all systems with flanking regions
            flanks_file = os.path.join(self.output_dir, 'all_systems_with_flanks.fasta')
        
            with open(systems_file, 'w') as sys_f, open(flanks_file, 'w') as flank_f:
                for system_id, boundary in self.boundary_data.items():
                    contig_id = boundary['contig']
                    
                    if contig_id not in self.genome_records:
                        print(f"Warning: Contig {contig_id} not found, skipping system {system_id}")
                        continue
                    
                    record = self.genome_records[contig_id]
                    
                    # Extract core system sequence
                    system_seq = record.seq[boundary['system_start']-1:boundary['system_end']]
                
                    # Write core system to systems file
                    sys_header = (f'{system_id}_system contig={contig_id} type={boundary["type"]} '
                             f'coords={boundary["system_start"]}-{boundary["system_end"]} '
                             f'length={boundary["length"]}bp genes={len(boundary["genes"])}')
                
                    sys_f.write(f'>{sys_header}\n')
                    sys_f.write(str(system_seq) + '\n')
                
                    # Extract system with flanking regions
                    flank_seq = record.seq[boundary['flank_start']-1:boundary['flank_end']]
                    flank_length = boundary['flank_end'] - boundary['flank_start'] + 1
                
                    # Calculate actual flank sizes
                    left_flank = boundary['system_start'] - boundary['flank_start']
                    right_flank = boundary['flank_end'] - boundary['system_end']
                
                    # Write system with flanks to flanks file
                    flank_header = (f'{system_id}_with_flanks contig={contig_id} type={boundary["type"]} '
                               f'coords={boundary["flank_start"]}-{boundary["flank_end"]} '
                               f'total_length={flank_length}bp '
                               f'system_length={boundary["length"]}bp '
                               f'left_flank={left_flank}bp right_flank={right_flank}bp '
                               f'genes={len(boundary["genes"])}')
                
                    flank_f.write(f'>{flank_header}\n')
                    flank_f.write(str(flank_seq) + '\n')
                
            print(f"All system sequences written to: {systems_file}")
            print(f"All systems with flanks written to: {flanks_file}")
                    
        except Exception as e:
            print(f"Error extracting sequences: {e}")

    def run_analysis(self):
        """Run the complete boundary delineation analysis"""
        print("Loading annotation data...")
        self.load_annotation_data()
        
        print("Parsing macsyfinder output...")
        self.parse_macsyfinder_output()
        
        print("Filtering systems...")
        self.filter_systems()
        
        print("Calculating boundaries...")
        self.calculate_boundaries()
        
        print("Generating output files...")
        self.generate_boundary_output()
        
        print("Extracting sequences...")
        self.extract_sequences()
        
        print(f"Analysis complete. Found {len(self.filtered_systems)} systems.")
        print(f"Results written to: {self.output_dir}")

def main():
    parser = argparse.ArgumentParser(description='Process macsyfinder output for boundary delineation')
    parser.add_argument('--all-systems', required=True, 
                       help='Path to all_systems.tsv file from macsyfinder')
    parser.add_argument('--genome', required=True,
                       help='Path to genome FASTA file')
    parser.add_argument('--annotation', 
                       help='Path to annotation file (GFF or GenBank format)')
    parser.add_argument('--output-dir', required=True,
                       help='Output directory for results')
    parser.add_argument('--flank-size', type=int, default=5,
                       help='Flanking region size in kb (default: 5)')
    
    args = parser.parse_args()
    
    # Validate input files
    if not os.path.exists(args.all_systems):
        print(f"Error: all_systems.tsv file not found: {args.all_systems}")
        sys.exit(1)
    
    if not os.path.exists(args.genome):
        print(f"Error: Genome file not found: {args.genome}")
        sys.exit(1)
    
    if args.annotation and not os.path.exists(args.annotation):
        print(f"Error: Annotation file not found: {args.annotation}")
        sys.exit(1)
    
    # Initialize processor
    processor = MacsyfinderProcessor(
        all_systems_file=args.all_systems,
        genome_file=args.genome,
        output_dir=args.output_dir,
        annotation_file=args.annotation
    )
    
    # Run analysis
    processor.run_analysis()

if __name__ == "__main__":
    main()
