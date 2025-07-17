#!/usr/bin/env python3

import argparse
import sys
from collections import defaultdict, Counter
import logging
import re

def setup_logging(verbose=False):
    """Setup logging configuration"""
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(
        level=level,
        format='%(asctime)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )

def validate_contig_consistency(ice_predictions, trna_coordinates, direct_repeats, gene_positions=None):
    """
    Validate that contig names are consistent across all input files
    and report data availability per contig
    """
    # Get all contig names from each data source
    ice_contigs = set(ice['seqname'] for ice in ice_predictions)
    trna_contigs = set(trna_coordinates.keys())
    dr_contigs = set(direct_repeats.keys())
    gene_contigs = set(gene_positions.keys()) if gene_positions else set()
    
    # Report data availability
    all_contigs = ice_contigs | trna_contigs | dr_contigs | gene_contigs
    
    logging.info(f"Data availability summary:")
    logging.info(f"  Total unique contigs: {len(all_contigs)}")
    logging.info(f"  ICE predictions: {len(ice_contigs)} contigs")
    logging.info(f"  tRNA coordinates: {len(trna_contigs)} contigs")
    logging.info(f"  Direct repeats: {len(dr_contigs)} contigs")
    if gene_positions:
        logging.info(f"  Gene annotations: {len(gene_contigs)} contigs")
    
    # Check for missing data per contig
    missing_data_report = []
    
    for contig in sorted(all_contigs):
        missing = []
        if contig not in ice_contigs:
            missing.append("ICE")
        if contig not in trna_contigs:
            missing.append("tRNA")
        if contig not in dr_contigs:
            missing.append("DR")
        if gene_positions and contig not in gene_contigs:
            missing.append("genes")
        
        if missing:
            missing_data_report.append(f"  {contig}: missing {', '.join(missing)}")
    
    if missing_data_report:
        logging.warning(f"Contigs with missing data:")
        for report in missing_data_report[:10]:  # Show first 10
            logging.warning(report)
        if len(missing_data_report) > 10:
            logging.warning(f"  ... and {len(missing_data_report) - 10} more contigs")
    
    # Return contigs that have at least ICE predictions
    processable_contigs = ice_contigs
    logging.info(f"Will process {len(processable_contigs)} contigs with ICE predictions")
    
    return processable_contigs, {
        'ice_contigs': ice_contigs,
        'trna_contigs': trna_contigs,
        'dr_contigs': dr_contigs,
        'gene_contigs': gene_contigs,
        'missing_data': missing_data_report
    }

def parse_ice_predictions(ice_file):
    """Parse ICE predictions from MacSyFinder TSV format"""
    ice_predictions = []
    contig_counts = Counter()
    
    try:
        with open(ice_file, 'r') as f:
            # Read and skip header
            header = next(f).strip()
            logging.debug(f"ICE file header: {header}")
            
            for line_num, line in enumerate(f, 2):  # Start from line 2
                if line.startswith('#') or not line.strip():
                    continue
                    
                fields = line.strip().split('\t')
                if len(fields) >= 5:  # MacSyFinder format has more fields
                    system_id = fields[0]
                    contig = fields[1]
                    ice_type = fields[2]
                    start = int(fields[3])
                    end = int(fields[4])
                    length = int(fields[5]) if len(fields) > 5 else end - start + 1
                    
                    contig_counts[contig] += 1
                    ice_predictions.append({
                        'seqname': contig,  # Use contig name as sequence name
                        'original_start': start,
                        'original_end': end,
                        'original_length': length,
                        'system_id': system_id,
                        'ice_type': ice_type,
                        'line_number': line_num
                    })
                else:
                    logging.warning(f"Skipping malformed line {line_num} in ICE predictions: insufficient fields")
    
    except FileNotFoundError:
        logging.error(f"ICE predictions file not found: {ice_file}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error parsing ICE predictions: {e}")
        sys.exit(1)
    
    # Report ICE distribution per contig
    logging.info(f"ICE predictions per contig:")
    for contig, count in contig_counts.most_common(10):
        logging.info(f"  {contig}: {count} ICEs")
    if len(contig_counts) > 10:
        logging.info(f"  ... and {len(contig_counts) - 10} more contigs")
    
    return ice_predictions

def parse_trna_gff(trna_gff_file):
    """
    Parse tRNA coordinates from aragorn GFF format
    Expected format: contig_1	aragorn	tRNA	18551	18626	.	-	.	ID=tRNA2;product=tRNA-Lys
    """
    trna_coordinates = defaultdict(list)
    contig_counts = Counter()
    
    try:
        with open(trna_gff_file, 'r') as f:
            for line in f:
                fields = line.strip().split('\t')
                if len(fields) == 9 and 'tRNA' in fields[2]:
                    seqname = fields[0]
                    start = int(fields[3])
                    end = int(fields[4])
                    strand = fields[6]
                    attributes = fields[8]

                    for att in attributes.split(';'):
                        att_key, att_value = att.split('=')
                        if att_key == 'ID':
                            trna_id = att_value
                        if att_key == 'product':
                            product = att_value
     
                    trna_coordinates[seqname].append({
                        'start': start,
                        'end': end,
                        'length': end - start + 1,
                        'strand': strand,
                        'id': trna_id,
                        'product': product,
                    })
                    contig_counts[seqname] += 1
                
                    logging.debug(f"Added tRNA: {seqname}:{trna_id} ({start}-{end}) {product}")
    
    except FileNotFoundError:
        logging.error(f"tRNA GFF file not found: {trna_gff_file}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error parsing tRNA GFF file: {e}")
        sys.exit(1)
    
    # Report tRNA distribution per contig
    logging.info(f"tRNA coordinates per contig:")
    for contig, count in contig_counts.most_common(10):
        logging.info(f"  {contig}: {count} tRNAs")
    if len(contig_counts) > 10:
        logging.info(f"  ... and {len(contig_counts) - 10} more contigs")
    
    return trna_coordinates

def parse_direct_repeats(dr_file):
    """
    Parse direct repeats from TSV format with length validation
    Expected format: contig_1	15291	15305	955	969
    """
    direct_repeats = defaultdict(list)
    contig_counts = Counter()
    filtered_counts = Counter()
    
    try:
        with open(dr_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                if line.startswith('#') or not line.strip():
                    continue
                
                fields = line.strip().split('\t')
                seqname = fields[0]
                start1 = int(fields[1])
                end1 = int(fields[2])
                start2 = int(fields[3])
                end2 = int(fields[4])
                    
                length1 = end1 - start1 + 1
                length2 = end2 - start2 + 1
                total_span = end2 - start1 + 1
                    
                # Track all DRs before filtering
                contig_counts[seqname] += 1
                    
                # Apply filtering rules from single.py
                if length1 < 15 or length2 < 15:
                    filtered_counts['too_short'] += 1
                    logging.debug(f"Skipping DR shorter than 15bp: {length1}, {length2}")
                    continue
                    
                if total_span > 500000:
                    filtered_counts['too_long'] += 1
                    logging.debug(f"Skipping DR with span > 500kb: {total_span}")
                    continue
                    
                if total_span < 5000:
                    filtered_counts['too_small'] += 1
                    logging.debug(f"Skipping DR with span < 5kb: {total_span}")
                    continue
                    
                direct_repeats[seqname].append({
                    'start1': start1,
                    'end1': end1,
                    'start2': start2,
                    'end2': end2,
                    'length1': length1,
                    'length2': length2,
                    'total_span': total_span,
                    'inner_distance': start2 - end1 - 1,
                    'line_number': line_num
                })
            else:
                logging.warning(f"Skipping malformed line {line_num} in direct repeats: insufficient fields")
    
    except FileNotFoundError:
        logging.error(f"Direct repeats file not found: {dr_file}")
        sys.exit(1)
    except Exception as e:
        logging.error(f"Error parsing direct repeats: {e}")
        sys.exit(1)
    
    # Report filtering statistics
    total_drs = sum(contig_counts.values())
    kept_drs = sum(len(drs) for drs in direct_repeats.values())
    
    logging.info(f"Direct repeat filtering summary:")
    logging.info(f"  Total DRs found: {total_drs}")
    logging.info(f"  DRs kept after filtering: {kept_drs}")
    logging.info(f"  Filtered out: {total_drs - kept_drs}")
    for reason, count in filtered_counts.items():
        logging.info(f"    {reason}: {count}")
    
    # Report DR distribution per contig
    logging.info(f"Valid direct repeats per contig:")
    dr_per_contig = {contig: len(drs) for contig, drs in direct_repeats.items()}
    for contig, count in Counter(dr_per_contig).most_common(10):
        logging.info(f"  {contig}: {count} DRs")
    if len(dr_per_contig) > 10:
        logging.info(f"  ... and {len(dr_per_contig) - 10} more contigs")
    
    return direct_repeats


def parse_gff_file(gff_file):
    """
    Parse GFF/GTF file to extract gene positions with contig tracking
    Returns gene positions organized by sequence and gene ID
    """
    gene_positions = defaultdict(dict)
    contig_counts = Counter()
    feature_counts = Counter()
    
    if not gff_file:
        logging.info("No GFF file provided, skipping gene boundary validation")
        return gene_positions
    
    try:
        with open(gff_file, 'r') as f:
            for line in f:
                line_l = line.rstrip().split("\t")
                # Annotation lines have exactly 9 columns
                if len(line_l) == 9:
                    (
                        seqname,
                        seq_source,
                        feature_type,
                        start,
                        end,
                        score,
                        strand,
                        phase,
                        attributes,
                    ) = line.rstrip().split("\t")
                
                    if feature_type == 'CDS':
                        for att in attributes.split(';'):
                            att_key, att_value = att.split('=')
                            if att_key == 'ID':
                                gene_id = att_value
                    
                        # Store gene position
                        gene_positions[seqname][gene_id] = {
                            'start': int(start),
                            'end': int(end),
                        }
                        contig_counts[seqname] += 1
                
                        logging.debug(f"Added gene: {seqname}:{gene_id} ({start}-{end})")
        
        # Report parsing statistics
        total_genes = sum(contig_counts.values())
        logging.info(f"GFF parsing summary:")
        logging.info(f"  Total genes parsed: {total_genes}")
        logging.info(f"  Contigs with genes: {len(gene_positions)}")
        
        # Report genes per contig
        logging.info(f"Genes per contig:")
        for contig, count in contig_counts.most_common(10):
            logging.info(f"  {contig}: {count} genes")
        if len(contig_counts) > 10:
            logging.info(f"  ... and {len(contig_counts) - 10} more contigs")
        
    except FileNotFoundError:
        logging.warning(f"GFF file not found: {gff_file}. Skipping gene boundary validation.")
        return defaultdict(dict)
    except Exception as e:
        logging.warning(f"Error parsing GFF file: {e}. Skipping gene boundary validation.")
        return defaultdict(dict)
    
    return gene_positions

def find_genes_near_coordinates(gene_positions, seqname, start, end, window=5):
    """
    Find genes near given coordinates (within ±window genes)
    Following single.py logic of ±5 genes from ICE boundaries
    """
    if seqname not in gene_positions:
        return []
    
    genes = gene_positions[seqname]
    
    # Get all genes sorted by start position
    sorted_genes = sorted(genes.items(), key=lambda x: x[1]['start'])
    
    # Find genes that overlap or are near the ICE region
    nearby_genes = []
    
    for gene_id, gene_info in sorted_genes:
        gene_start = gene_info['start']
        gene_end = gene_info['end']
        
        # Check if gene overlaps with ICE region or is nearby
        if (gene_end >= start - 10000 and gene_start <= end + 10000):  # 10kb window
            nearby_genes.append((gene_id, gene_info))
    
    return nearby_genes

def validate_dr_gene_overlap(dr, ice_prediction, gene_positions):
    """
    Check if DR positions overlap with ICE gene boundaries
    Following single.py logic for gene boundary validation
    """
    seqname = ice_prediction['seqname']
    ice_start = ice_prediction['original_start']
    ice_end = ice_prediction['original_end']
    
    if seqname not in gene_positions:
        logging.debug(f"No gene positions available for {seqname}")
        return True  # Skip validation if no gene data
    
    # Find genes near the ICE region
    nearby_genes = find_genes_near_coordinates(gene_positions, seqname, ice_start, ice_end)
    
    if not nearby_genes:
        logging.debug(f"No genes found near ICE region {seqname}:{ice_start}-{ice_end}")
        return True  # Skip validation if no nearby genes
    
    # Check if DR boundaries overlap with any nearby gene boundaries
    for gene_id, gene_info in nearby_genes:
        gene_start = gene_info['start']
        gene_end = gene_info['end']
        
        # Scenario A: DR start within gene boundaries
        if gene_start <= dr['start1'] <= gene_end:
            logging.debug(f"DR start {dr['start1']} within gene {gene_id} boundaries ({gene_start}-{gene_end})")
            return True
        
        # Scenario B: DR end within gene boundaries  
        if gene_start <= dr['end2'] <= gene_end:
            logging.debug(f"DR end {dr['end2']} within gene {gene_id} boundaries ({gene_start}-{gene_end})")
            return True
        
        # Additional check: gene boundaries within DR span
        if dr['start1'] <= gene_start <= dr['end2'] or dr['start1'] <= gene_end <= dr['end2']:
            logging.debug(f"Gene {gene_id} boundaries within DR span")
            return True
    
    logging.debug(f"DR {dr['start1']}-{dr['end2']} doesn't overlap with gene boundaries")
    return False

def count_trnas_in_span(trna_list, start, end):
    """Count tRNAs within a given genomic span"""
    count = 0
    trnas_in_span = []
    
    for trna in trna_list:
        # Check if tRNA is completely within the span
        if start <= trna['start'] <= end and start <= trna['end'] <= end:
            count += 1
            trnas_in_span.append(trna)
    
    return count, trnas_in_span

def find_max_gap_between_trnas(trnas_in_span):
    """
    Find the maximum gap between tRNAs to identify ICE insertion site
    Following single.py gap analysis logic
    """
    if len(trnas_in_span) < 2:
        return None, 0
    
    # Sort tRNAs by start position
    sorted_trnas = sorted(trnas_in_span, key=lambda x: x['start'])
    
    max_gap = 0
    gap_location = None
    
    for i in range(len(sorted_trnas) - 1):
        current_end = sorted_trnas[i]['end']
        next_start = sorted_trnas[i + 1]['start']
        gap = next_start - current_end - 1
        
        if gap > max_gap:
            max_gap = gap
            gap_location = {
                'gap_start': current_end + 1,
                'gap_end': next_start - 1,
                'gap_size': gap,
                'upstream_trna': sorted_trnas[i],
                'downstream_trna': sorted_trnas[i + 1]
            }
    
    return gap_location, max_gap

def find_dr_flanking_ice_and_trna(ice_prediction, trna_coordinates, direct_repeats, gene_positions=None):
    """
    Find direct repeats that flank ICE and contain sufficient tRNAs
    Implements single.py logic with ≥2 tRNA requirement
    """
    seqname = ice_prediction['seqname']
    ice_start = ice_prediction['original_start']
    ice_end = ice_prediction['original_end']
    
    if seqname not in trna_coordinates or seqname not in direct_repeats:
        logging.debug(f"No tRNA or DR data for sequence: {seqname}")
        return None
    
    trna_list = trna_coordinates[seqname]
    dr_list = direct_repeats[seqname]
    
    best_dr = None
    best_score = 0
    best_trnas = []
    best_gap = None
    
    logging.debug(f"Evaluating {len(dr_list)} direct repeats for ICE {seqname}:{ice_start}-{ice_end}")
    
    for dr in dr_list:
        # Check if DR spans the ICE region
        if dr['start1'] > ice_start or dr['end2'] < ice_end:
            logging.debug(f"DR {dr['start1']}-{dr['end2']} doesn't span ICE {ice_start}-{ice_end}")
            continue
        
        # Validate gene boundary overlap if gene positions available
        if gene_positions and not validate_dr_gene_overlap(dr, ice_prediction, gene_positions):
            logging.debug(f"DR doesn't meet gene boundary overlap criteria")
            continue
        
        # Count tRNAs within DR span - CRITICAL: Must be ≥2 (not ≥1)
        trnas_between, trnas_in_span = count_trnas_in_span(trna_list, dr['start1'], dr['end2'])
        
        if trnas_between < 2:  # Following single.py requirement
            logging.debug(f"DR has only {trnas_between} tRNAs, need ≥2")
            continue
        
        # Find the largest gap between tRNAs
        gap_location, max_gap = find_max_gap_between_trnas(trnas_in_span)
        
        # Score this DR (prefer more tRNAs and larger gaps)
        score = trnas_between * 1000 + max_gap
        
        if score > best_score:
            best_score = score
            best_dr = dr
            best_trnas = trnas_in_span
            best_gap = gap_location
            
        logging.debug(f"DR {dr['start1']}-{dr['end2']}: {trnas_between} tRNAs, max gap: {max_gap}, score: {score}")
    
    if best_dr is None:
        logging.debug(f"No suitable DR found for ICE {seqname}:{ice_start}-{ice_end}")
        return None
    
    return {
        'selected_dr': best_dr,
        'trnas_in_span': best_trnas,
        'gap_location': best_gap,
        'num_trnas': len(best_trnas),
        'score': best_score
    }

def refine_ice_boundaries(ice_prediction, trna_coordinates, direct_repeats, gene_positions=None, verbose=False):
    """
    Refine ICE boundaries using direct repeats and tRNA analysis
    Implements single.py methodology with proper fallback logic
    """
    seqname = ice_prediction['seqname']
    original_start = ice_prediction['original_start']
    original_end = ice_prediction['original_end']
    original_length = ice_prediction['original_length']
    
    logging.debug(f"Refining boundaries for ICE {seqname}:{original_start}-{original_end}")
    
    # Find the best DR flanking the ICE
    dr_result = find_dr_flanking_ice_and_trna(ice_prediction, trna_coordinates, direct_repeats, gene_positions)
    
    if dr_result is None:
        # Fallback logic: use original boundaries (like single.py)
        logging.debug(f"No suitable DR found for {seqname}, using original boundaries")
        return {
            'seqname': seqname,
            'refined_start': original_start,
            'refined_end': original_end,
            'refined_length': original_length,
            'dr1_start': original_start,  # attL site
            'dr1_end': original_start,
            'dr2_start': original_end,    # attR site  
            'dr2_end': original_end,
            'num_trnas': len(trna_coordinates.get(seqname, [])),
            'gap_location': 'none',
            'refinement_method': 'fallback_original'
        }
    
    # Use DR boundaries for refined ICE coordinates
    selected_dr = dr_result['selected_dr']
    gap_info = dr_result['gap_location']
    
    refined_start = selected_dr['start1']
    refined_end = selected_dr['end2']
    refined_length = refined_end - refined_start + 1
    
    # Determine gap location description
    gap_desc = 'none'
    if gap_info:
        gap_desc = f"{gap_info['gap_start']}-{gap_info['gap_end']}"
    
    logging.debug(f"Refined {seqname}: {refined_start}-{refined_end} (length: {refined_length})")
    logging.debug(f"Found {dr_result['num_trnas']} tRNAs, largest gap: {gap_desc}")
    
    return {
        'seqname': seqname,
        'refined_start': refined_start,
        'refined_end': refined_end,
        'refined_length': refined_length,
        'dr1_start': selected_dr['start1'],  # attL site (left DR start)
        'dr1_end': selected_dr['end1'],      # Left DR end
        'dr2_start': selected_dr['start2'],  # Right DR start
        'dr2_end': selected_dr['end2'],      # attR site (right DR end)
        'num_trnas': dr_result['num_trnas'],
        'gap_location': gap_desc,
        'refinement_method': 'direct_repeat'
    }

def write_refined_boundaries(refined_results, output_file):
    """Write refined ICE boundaries to output file with contig organization"""
    try:
        with open(output_file, 'w') as f:
            # Write header
            header = [
                'seqname', 'refined_start', 'refined_end', 'refined_length',
                'dr1_start', 'dr1_end', 'dr2_start', 'dr2_end',
                'num_trnas', 'gap_location', 'refinement_method'
            ]
            f.write('\t'.join(header) + '\n')
            
            # Sort results by contig name and then by start position
            sorted_results = sorted(refined_results, key=lambda x: (x['seqname'], x['refined_start']))
            
            # Write results
            for result in sorted_results:
                row = [
                    result['seqname'],
                    str(result['refined_start']),
                    str(result['refined_end']),
                    str(result['refined_length']),
                    str(result['dr1_start']),
                    str(result['dr1_end']),
                    str(result['dr2_start']),
                    str(result['dr2_end']),
                    str(result['num_trnas']),
                    result['gap_location'],
                    result['refinement_method']
                ]
                f.write('\t'.join(row) + '\n')
                
        logging.info(f"Results written to: {output_file}")
        
    except Exception as e:
        logging.error(f"Error writing output file: {e}")
        sys.exit(1)

def generate_summary_report(refined_results, data_summary):
    """Generate a comprehensive summary report of the refinement process"""
    
    # Overall statistics
    total_ices = len(refined_results)
    successful_refinements = sum(1 for r in refined_results if r['refinement_method'] == 'direct_repeat')
    
    # Per-contig statistics
    contig_stats = defaultdict(lambda: {'total': 0, 'refined': 0, 'fallback': 0})
    
    for result in refined_results:
        contig = result['seqname']
        contig_stats[contig]['total'] += 1
        if result['refinement_method'] == 'direct_repeat':
            contig_stats[contig]['refined'] += 1
        else:
            contig_stats[contig]['fallback'] += 1
    
    logging.info("="*60)
    logging.info("REFINEMENT SUMMARY REPORT")
    logging.info("="*60)
    
    logging.info(f"Overall Statistics:")
    logging.info(f"  Total ICEs processed: {total_ices}")
    logging.info(f"  Successfully refined: {successful_refinements}")
    logging.info(f"  Used fallback boundaries: {total_ices - successful_refinements}")
    logging.info(f"  Success rate: {successful_refinements/total_ices*100:.1f}%")
    
    logging.info(f"\nPer-contig refinement success:")
    for contig in sorted(contig_stats.keys())[:10]:  # Show top 10
        stats = contig_stats[contig]
        success_rate = stats['refined'] / stats['total'] * 100 if stats['total'] > 0 else 0
        logging.info(f"  {contig}: {stats['refined']}/{stats['total']} ({success_rate:.1f}%)")
    
    if len(contig_stats) > 10:
        logging.info(f"  ... and {len(contig_stats) - 10} more contigs")
    
    # Data availability impact
    logging.info(f"\nData availability impact:")
    missing_trna = len([c for c in data_summary['ice_contigs'] if c not in data_summary['trna_contigs']])
    missing_dr = len([c for c in data_summary['ice_contigs'] if c not in data_summary['dr_contigs']])
    missing_genes = len([c for c in data_summary['ice_contigs'] if c not in data_summary['gene_contigs']]) if data_summary['gene_contigs'] else 0
    
    logging.info(f"  ICE contigs missing tRNA data: {missing_trna}")
    logging.info(f"  ICE contigs missing DR data: {missing_dr}")
    if data_summary['gene_contigs']:
        logging.info(f"  ICE contigs missing gene data: {missing_genes}")
    
    logging.info("="*60)

def main():
    parser = argparse.ArgumentParser(
        description='Refine ICE boundaries using direct repeats and tRNA analysis (single.py compatible)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s -i ice_predictions.tsv -t trna_coords.gff -d direct_repeats.tsv -o refined_boundaries.tsv
  %(prog)s -i ice_predictions.tsv -t trna_coords.gff -d direct_repeats.tsv -g annotation.gff -o refined_boundaries.tsv -v

Input file formats:
  ICE predictions (MacSyFinder): system_id<TAB>contig<TAB>type<TAB>start<TAB>end<TAB>length...
  tRNA coordinates (GFF): contig<TAB>aragorn<TAB>tRNA<TAB>start<TAB>end<TAB>.<TAB>strand<TAB>.<TAB>attributes
  Direct repeats (TSV): contig<TAB>start1<TAB>end1<TAB>start2<TAB>end2
  GFF annotation: Standard GFF3/GTF format
        """
    )
    
    parser.add_argument('-i', '--ice-predictions', required=True,
                        help='Input file with ICE predictions from MacSyFinder (TSV format)')
    parser.add_argument('-t', '--trna-coordinates', required=True,
                        help='Input file with tRNA coordinates from aragorn (GFF format)')
    parser.add_argument('-d', '--direct-repeats', required=True,
                        help='Input file with direct repeats (TSV format)')
    parser.add_argument('-g', '--gff-file', 
                        help='Optional: GFF/GTF annotation file for gene boundary validation')
    parser.add_argument('-o', '--output', required=True,
                        help='Output file for refined ICE boundaries (TSV format)')
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='Enable verbose logging')
    parser.add_argument('--validate-contigs', action='store_true',
                        help='Perform strict contig name validation across input files')
    
    args = parser.parse_args()
    
    # Setup logging
    setup_logging(args.verbose)
    
    logging.info("Starting ICE boundary refinement (single.py compatible)")
    logging.info(f"ICE predictions: {args.ice_predictions}")
    logging.info(f"tRNA coordinates: {args.trna_coordinates}")
    logging.info(f"Direct repeats: {args.direct_repeats}")
    if args.gff_file:
        logging.info(f"GFF annotation: {args.gff_file}")
    logging.info(f"Output file: {args.output}")
    
    # Parse input files with enhanced contig tracking
    logging.info("\nParsing input files...")
    ice_predictions = parse_ice_predictions(args.ice_predictions)
    trna_coordinates = parse_trna_gff(args.trna_coordinates)  # Updated to use GFF parser
    direct_repeats = parse_direct_repeats(args.direct_repeats)
    
    gene_positions = None
    if args.gff_file:
        gene_positions = parse_gff_file(args.gff_file)
    
    # Validate contig consistency and data availability
    logging.info("\nValidating data consistency across contigs...")
    processable_contigs, data_summary = validate_contig_consistency(
        ice_predictions, trna_coordinates, direct_repeats, gene_positions
    )
    
    # Optional strict validation
    if args.validate_contigs:
        ice_contigs = data_summary['ice_contigs']
        trna_contigs = data_summary['trna_contigs']
        dr_contigs = data_summary['dr_contigs']
        
        missing_essential_data = []
        for contig in ice_contigs:
            if contig not in trna_contigs:
                missing_essential_data.append(f"{contig}: missing tRNA data")
            if contig not in dr_contigs:
                missing_essential_data.append(f"{contig}: missing DR data")
        
        if missing_essential_data:
            logging.error("Strict validation failed. Missing essential data:")
            for issue in missing_essential_data[:10]:
                logging.error(f"  {issue}")
            if len(missing_essential_data) > 10:
                logging.error(f"  ... and {len(missing_essential_data) - 10} more issues")
            logging.error("Use --validate-contigs=false to proceed anyway")
            sys.exit(1)
    
    # Process each ICE prediction with enhanced error handling
    logging.info(f"\nProcessing {len(ice_predictions)} ICE predictions...")
    refined_results = []
    successful_refinements = 0
    failed_processing = []
    
    # Group ICEs by contig for better progress reporting
    ices_by_contig = defaultdict(list)
    for ice in ice_predictions:
        ices_by_contig[ice['seqname']].append(ice)
    
    processed_contigs = 0
    total_contigs = len(ices_by_contig)
    
    for contig, contig_ices in ices_by_contig.items():
        processed_contigs += 1
        logging.info(f"Processing contig {processed_contigs}/{total_contigs}: {contig} ({len(contig_ices)} ICEs)")
        
        contig_successes = 0
        
        for ice_prediction in contig_ices:
            try:
                result = refine_ice_boundaries(
                    ice_prediction, 
                    trna_coordinates, 
                    direct_repeats, 
                    gene_positions, 
                    args.verbose
                )
                refined_results.append(result)
                
                if result['refinement_method'] == 'direct_repeat':
                    successful_refinements += 1
                    contig_successes += 1
                    
            except Exception as e:
                error_msg = f"Error processing ICE {ice_prediction['seqname']}:{ice_prediction['original_start']}-{ice_prediction['original_end']}: {e}"
                logging.error(error_msg)
                failed_processing.append(error_msg)
                continue
        
        if contig_ices:
            success_rate = contig_successes / len(contig_ices) * 100
            logging.info(f"  Contig {contig}: {contig_successes}/{len(contig_ices)} refined ({success_rate:.1f}%)")
    
    # Write results with enhanced organization
    logging.info(f"\nWriting results...")
    write_refined_boundaries(refined_results, args.output)
    
    # Generate comprehensive summary report
    generate_summary_report(refined_results, data_summary)
    
    # Report any processing failures
    if failed_processing:
        logging.warning(f"\nProcessing failures ({len(failed_processing)}):")
        for failure in failed_processing[:5]:  # Show first 5
            logging.warning(f"  {failure}")
        if len(failed_processing) > 5:
            logging.warning(f"  ... and {len(failed_processing) - 5} more failures")
    
    # Final success check
    if successful_refinements == 0:
        logging.warning("No ICE boundaries were successfully refined!")
        logging.warning("Check that:")
        logging.warning("  1. Direct repeats span the ICE regions")
        logging.warning("  2. At least 2 tRNAs are present within DR spans")
        logging.warning("  3. DR lengths are ≥15bp and spans are 5kb-500kb")
        logging.warning("  4. Contig names match across all input files")
    
    logging.info(f"\nProcessing complete. Check {args.output} for results.")

if __name__ == '__main__':
    main()

