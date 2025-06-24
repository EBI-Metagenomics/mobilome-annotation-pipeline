#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright 2025 EMBL - European Bioinformatics Institute
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


import argparse
import sys
from collections import defaultdict, namedtuple, Counter
from typing import List, Tuple, Dict, Optional
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
import multiprocessing as mp
from concurrent.futures import ProcessPoolExecutor, as_completed
from functools import partial

# Data structures
ContigInfo = namedtuple('ContigInfo', ['name', 'sequence', 'length'])
Window = namedtuple('Window', ['start', 'end', 'gc_content', 'kmer_profile', 'gc_score', 'kmer_score', 'combined_score'])
Outlier = namedtuple('Outlier', ['contig', 'start', 'end', 'gc_score', 'kmer_score', 'combined_score', 'confidence'])
Repeat = namedtuple('Repeat', ['start', 'end', 'sequence', 'type'])
Prediction = namedtuple('Prediction', ['outlier', 'left_repeat', 'right_repeat'])

class FastKmerAnalyzer:
    # k-mer analyzer using 3-mer and 4-mer
    def __init__(self, k_sizes: List[int] = [3, 4]):
        self.k_sizes = k_sizes
        # pre-compute all possible k-mers
        self.all_kmers = {}
        for k in k_sizes:
            bases = ['A', 'T', 'G', 'C']
            self.all_kmers[k] = [''.join(kmer) for kmer in self._product(bases, k)]
    
    def _product(self, bases, k):
        # itertools.product for small k
        if k == 1:
            return [[b] for b in bases]
        elif k == 2:
            return [[b1, b2] for b1 in bases for b2 in bases]
        elif k == 3:
            return [[b1, b2, b3] for b1 in bases for b2 in bases for b3 in bases]
        elif k == 4:
            return [[b1, b2, b3, b4] for b1 in bases for b2 in bases for b3 in bases for b4 in bases]
    
    def calculate_kmer_profile(self, sequence: str) -> Dict[str, float]:
        # calculate the k-mer profile
        sequence = sequence.upper()
        combined_profile = {}
        
        for k in self.k_sizes:
            kmer_counts = Counter()
            
            # counting k-mers
            seq_len = len(sequence)
            for i in range(seq_len - k + 1):
                kmer = sequence[i:i + k]
                if len(kmer) == k and all(base in 'ATGC' for base in kmer):
                    kmer_counts[kmer] += 1
            
            # Normalise frequencies
            total_kmers = sum(kmer_counts.values())
            if total_kmers > 0:
                for kmer in self.all_kmers[k]:
                    key = f"{k}mer_{''.join(kmer)}"
                    combined_profile[key] = kmer_counts.get(''.join(kmer), 0) / total_kmers
            else:
                for kmer in self.all_kmers[k]:
                    key = f"{k}mer_{''.join(kmer)}"
                    combined_profile[key] = 0.0
        
        return combined_profile
    
    def profile_to_vector(self, profile: Dict[str, float]) -> np.ndarray:
        # convert profile to numpy vector for distance calculations
        return np.array(list(profile.values()))

class FastCompositionAnalyzer: 
    def __init__(self, window_size, step_size):
        self.window_size = window_size
        self.step_size = step_size
        self.kmer_analyzer = FastKmerAnalyzer()
    
    def calculate_gc_content(self, sequence: str) -> float:
        sequence = sequence.upper()
        gc_count = sequence.count('G') + sequence.count('C')
        return gc_count / len(sequence) if len(sequence) > 0 else 0.0
    
    def sliding_window_analysis(self, contig: ContigInfo) -> List[Window]:
        # sliding window analysis
        windows = []
        sequence = contig.sequence
        seq_len = len(sequence)
        
        for start in range(0, seq_len - self.window_size + 1, self.step_size):
            end = start + self.window_size
            window_seq = sequence[start:end]
            
            # calculate features
            gc_content = self.calculate_gc_content(window_seq)
            kmer_profile = self.kmer_analyzer.calculate_kmer_profile(window_seq)
            
            windows.append(Window(
                start, end, gc_content, kmer_profile, 0.0, 0.0, 0.0
            ))
        
        return windows
    
    def calculate_fast_outlier_scores(self, windows: List[Window]) -> List[Window]:
        if len(windows) < 3:
            return windows
        
        # calculate basic statistics
        gc_values = np.array([w.gc_content for w in windows])
        mean_gc = np.mean(gc_values)
        std_gc = np.std(gc_values)
        
        # k-mer composition analysis
        kmer_vectors = np.array([self.kmer_analyzer.profile_to_vector(w.kmer_profile) for w in windows])
        mean_kmer_vector = np.mean(kmer_vectors, axis=0)
        
        # calculate distances 
        kmer_distances = np.linalg.norm(kmer_vectors - mean_kmer_vector, axis=1)
        mean_kmer_dist = np.mean(kmer_distances)
        std_kmer_dist = np.std(kmer_distances)
        
        # calculate scores
        scored_windows = []
        for i, window in enumerate(windows):
            gc_z_score = abs(window.gc_content - mean_gc) / std_gc if std_gc > 0 else 0.0
            kmer_z_score = abs(kmer_distances[i] - mean_kmer_dist) / std_kmer_dist if std_kmer_dist > 0 else 0.0
            
            # combined score
            combined_score = max(gc_z_score, kmer_z_score)
            
            scored_windows.append(Window(
                window.start, window.end, window.gc_content, window.kmer_profile,
                gc_z_score, kmer_z_score, combined_score
            ))
        
        return scored_windows

class FastOutlierDetector:
    def __init__(self, score_threshold, min_region_size: int = 5000):
        self.score_threshold = score_threshold
        self.min_region_size = min_region_size
    
    def filter_significant_outliers(self, windows: List[Window]) -> List[Window]:
        return [w for w in windows if w.combined_score >= self.score_threshold]
    
    def merge_overlapping_outliers(self, outliers: List[Outlier]) -> List[Outlier]:
        if not outliers:
            return []
        
        # sort by start position
        sorted_outliers = sorted(outliers, key=lambda x: x.start)
        merged = [sorted_outliers[0]]
        
        for current in sorted_outliers[1:]:
            last_merged = merged[-1]
            
            # check for proximity (within 2kb)
            if current.start <= last_merged.end + 2000:
                # Merge regions
                new_end = max(last_merged.end, current.end)
                new_gc_score = max(last_merged.gc_score, current.gc_score)
                new_kmer_score = max(last_merged.kmer_score, current.kmer_score)
                new_combined_score = max(last_merged.combined_score, current.combined_score)
                new_confidence = max(last_merged.confidence, current.confidence)
                
                merged[-1] = Outlier(
                    last_merged.contig, last_merged.start, new_end,
                    new_gc_score, new_kmer_score, new_combined_score, new_confidence
                )
            else:
                merged.append(current)
        
        # filter by minimum size
        return [outlier for outlier in merged if (outlier.end - outlier.start) >= self.min_region_size]

class FastRepeatFinder:
    # Finding direct and reverse repeats flanking the outliers
    # This is extra evidence of insertion in the chromosome
    # Only exact alignments are considered

    def __init__(self, flank_window: int = 500, min_repeat_length: int = 10, max_repeat_length: int = 100):
        self.flank_window = flank_window
        self.min_repeat_length = min_repeat_length
        self.max_repeat_length = max_repeat_length
    
    def find_exact_matches(self, seq1, seq2, is_inverted: bool) -> List[Tuple[int, int, int, str]]:
        matches = []
        if is_inverted:
            seq2 = str(Seq(seq2).reverse_complement())

        # Parsing all possible lengths 
        for length in list(range(self.max_repeat_length, self.min_repeat_length, -1)):
            if length > len(seq1) or length > len(seq2):
                continue
                
            for i in range(len(seq1) - length + 1):
                subseq1 = seq1[i:i + length]
                
                # search in seq2
                pos = seq2.find(subseq1)
                if pos != -1:
                    matches.append((i, pos, length, subseq1))

        return matches
    
    def extract_flanking_regions(self, sequence: str, outlier: Outlier) -> Tuple[str, str]:
        seq_len = len(sequence)
        
        # Left flanking region
        left_start = max(0, outlier.start - self.flank_window)
        left_end = outlier.start
        left_flank = sequence[left_start:left_end]
        
        # Right flanking region  
        right_start = outlier.end
        right_end = min(seq_len, outlier.end + self.flank_window)
        right_flank = sequence[right_start:right_end]
        
        return left_flank, right_flank
    
    def find_flanking_repeats(self, sequence: str, outlier: Outlier) -> Optional[Tuple[Repeat, Repeat]]:
        left_flank, right_flank = self.extract_flanking_regions(sequence, outlier)
        
        # Try direct repeats first
        direct_matches = self.find_exact_matches(left_flank, right_flank, is_inverted=False)
        
        # Try inverted repeats if no direct matches
        inverted_matches = []
        if not direct_matches:
            inverted_matches = self.find_exact_matches(left_flank, right_flank, is_inverted=True)
        
        # Selecting the longest match
        best_match = None
        repeat_type = None
        
        if direct_matches:
            best_match = max(direct_matches, key=lambda x: x[2])
            repeat_type = "direct"
        elif inverted_matches:
            best_match = max(inverted_matches, key=lambda x: x[2])
            repeat_type = "inverted"
        
        if not best_match:
            return None
        
        # calculate absolute coordinates
        left_start = max(0, outlier.start - self.flank_window)
        right_start = outlier.end
        
        left_repeat_start = left_start + best_match[0]
        left_repeat_end = left_repeat_start + best_match[2]
        
        right_repeat_start = right_start + best_match[1]
        right_repeat_end = right_repeat_start + best_match[2]
        
        left_repeat = Repeat(left_repeat_start, left_repeat_end, best_match[3], repeat_type)
        right_repeat = Repeat(right_repeat_start, right_repeat_end, best_match[3], repeat_type)
       
        return left_repeat, right_repeat

def process_contig_fast(contig_info: Tuple[ContigInfo, Dict]) -> Tuple[List[Prediction], str]:
    contig, params = contig_info
    
    try:
        composition_analyzer = FastCompositionAnalyzer(
            params['window_size'], 
            params['step_size']
        )
        outlier_detector = FastOutlierDetector(
            params['score_threshold']
        )
        repeat_finder = FastRepeatFinder()
        
        # processing contig
        windows = composition_analyzer.sliding_window_analysis(contig)
        
        if len(windows) < 3:
            return [], contig.name
        
        scored_windows = composition_analyzer.calculate_fast_outlier_scores(windows)
        significant_windows = outlier_detector.filter_significant_outliers(scored_windows)
        
        if not significant_windows:
            return [], contig.name
        
        # Saving outliers
        outliers = []
        for w in significant_windows:
            confidence = min(w.combined_score / 3.0, 1.0)  # Simple confidence
            outlier = Outlier(
                contig.name, w.start, w.end, w.gc_score, w.kmer_score, 
                w.combined_score, confidence
            )
            outliers.append(outlier)
        
        # merge overlapping outliers and find repeats
        merged_outliers = outlier_detector.merge_overlapping_outliers(outliers)
        
        predictions = []
        for outlier in merged_outliers:
            repeat_pair = repeat_finder.find_flanking_repeats(contig.sequence, outlier)
            
            if repeat_pair:
                left_repeat, right_repeat = repeat_pair
                corrected_outlier = Outlier(
                    outlier.contig, left_repeat.start, right_repeat.end,
                    outlier.gc_score, outlier.kmer_score, outlier.combined_score, outlier.confidence
                )
                predictions.append(Prediction(corrected_outlier, left_repeat, right_repeat))
        
        return predictions, contig.name
        
    except Exception as e:
        print(f"Error processing contig {contig.name}: {e}", file=sys.stderr)
        return [], contig.name

class FastGenomicOutlierPipeline:
    
    def __init__(self, min_contig_length, window_size, 
                 step_size, score_threshold, n_threads):
        self.min_contig_length = min_contig_length
        self.window_size = window_size
        self.step_size = step_size
        self.score_threshold = score_threshold
        
        # Determine number of threads
        if n_threads is None:
            self.n_threads = min(mp.cpu_count(), 4)  # Conservative for memory
        else:
            self.n_threads = min(n_threads, mp.cpu_count())
    
    def load_and_filter_contigs(self, fasta_file: str) -> List[ContigInfo]:
        # Reading fasta file and filter contigs by min lenght
        contigs = []
        
        try:
            for record in SeqIO.parse(fasta_file, "fasta"):
                if len(record.seq) >= self.min_contig_length:
                    contigs.append(ContigInfo(record.id, str(record.seq), len(record.seq)))
            
            print(f"Loaded {len(contigs)} contigs >= {self.min_contig_length} bp", file=sys.stderr)
            return contigs
            
        except Exception as e:
            print(f"Error loading FASTA file: {e}", file=sys.stderr)
            sys.exit(1)
    
    def process_contigs_parallel(self, contigs: List[ContigInfo]) -> List[Prediction]:
        # processing multiple contigs in parallel
        all_predictions = []
        total_contigs = len(contigs)

        print(f"Processing {total_contigs} contigs using {self.n_threads} threads...", file=sys.stderr)
        
        # Prepare parameters
        params = {
            'window_size': self.window_size,
            'step_size': self.step_size,
            'score_threshold': self.score_threshold
        }
        
        # Prepare input for multiprocessing
        contig_params = [(contig, params) for contig in contigs]
        
        if self.n_threads == 1 or len(contigs) == 1:
            # Single-threaded processing
            for contig_param in contig_params:
                predictions, _ = process_contig_fast(contig_param)
                all_predictions.extend(predictions)
        else:
            # Multi-threaded processing
            with ProcessPoolExecutor(max_workers=self.n_threads) as executor:
                # Submit all jobs
                future_to_contig = {
                    executor.submit(process_contig_fast, contig_param): contig_param[0] 
                    for contig_param in contig_params
                }
                
                # Collect results
                processed = 0
                for future in as_completed(future_to_contig):
                    contig = future_to_contig[future]
                    #try:
                    predictions, contig_name = future.result()
                    all_predictions.extend(predictions)
                    #    processed += 1
                    #   if processed % 10 == 0:
                    #      print(f"  Processed {processed}/{total_contigs} contigs", file=sys.stderr)
                    #except Exception as e:
                    #    print(f"Error processing contig {contig.name}: {e}", file=sys.stderr)
        
        return all_predictions
    
    def run_pipeline(self, fasta_file: str, output_file: str):
        print("Starting Fast Compositional Outlier Detection Pipeline...", file=sys.stderr)
        print(f"Score threshold: {self.score_threshold}", file=sys.stderr)
        print(f"Threads: {self.n_threads}", file=sys.stderr)
        
        # Load and filter contigs
        contigs = self.load_and_filter_contigs(fasta_file)
        
        # Process contigs in parallel
        all_predictions = self.process_contigs_parallel(contigs)
        
        print(f"Total predictions found: {len(all_predictions)}", file=sys.stderr)
        
        # Sort predictions by confidence score (descending)
        all_predictions.sort(key=lambda p: p.outlier.confidence, reverse=True)
        
        # Write output
        self.write_output(all_predictions, output_file)
        print(f"Results written to {output_file}", file=sys.stderr)
    
    def write_output(self, predictions: List[Prediction], output_file: str):
        # Write predictions to BED format file
        with open(output_file, 'w') as f:
            # Write header
            f.write(f"# Fast compositional outlier analysis\n")
            f.write(f"# Score threshold: {self.score_threshold}\n")
            f.write(f"# Format: contig\tstart\tend\tinfo\n")
            
            for i, pred in enumerate(predictions, 1):
                outlier = pred.outlier
                left_repeat = pred.left_repeat
                right_repeat = pred.right_repeat

                right_repeat_sequence = str(right_repeat.sequence)
                if right_repeat.type == 'inverted':
                    right_repeat_sequence = Seq(right_repeat_sequence).reverse_complement()

                # Create prediction name
                pred_name = f"{outlier.contig}_{i:03d}"
                outlier_coords = f"{left_repeat.start}-{right_repeat.end}"
                
                # Score information
                score_info = (f"gc:{outlier.gc_score:.2f}|"
                            f"kmer:{outlier.kmer_score:.2f}|"
                            f"combined:{outlier.combined_score:.2f}")
                
                confidence_info = f"conf:{outlier.confidence:.3f}"
                
                # Write left repeat
                left_info = (f"{pred_name}:{outlier_coords}|"
                           f"{left_repeat.type}_1|{left_repeat.sequence}|"
                           f"{score_info}|{confidence_info}")
                f.write(f"{outlier.contig}\t{left_repeat.start}\t{left_repeat.end}\t{left_info}\n")
                
                # Write right repeat
                right_info = (f"{pred_name}:{outlier_coords}|"
                            f"{right_repeat.type}_2|{right_repeat_sequence}|"
                            f"{score_info}|{confidence_info}")
                f.write(f"{outlier.contig}\t{right_repeat.start}\t{right_repeat.end}\t{right_info}\n")

def main():
    parser = argparse.ArgumentParser(
        description="Fast compositional outlier detection in long contigs",
        formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument(
        "--input_fasta",
        help="Input FASTA file with contigs"
    )
    parser.add_argument(
        "--output_bed",
        help="Output BED file with predictions"
    )
    parser.add_argument(
        "--score-threshold",
        type=float,
        default=2.0,
        help="Z-score threshold for outlier detection (default: 2.0)"
    )
    parser.add_argument(
        "--min-contig-length",
        type=int,
        default=100000,
        help="Minimum contig length to consider (default: 100000)"
    )
    parser.add_argument(
        "--window-size",
        type=int,
        default=5000,
        help="Window size for compositional analysis (default: 5000)"
    )
    parser.add_argument(
        "--threads",
        type=int,
        help="Number of threads to use (default: auto-detect, max 4)"
    )
    args = parser.parse_args()
    
    # Running the pipeline
    pipeline = FastGenomicOutlierPipeline(
        min_contig_length=int(args.min_contig_length),
        window_size=int(args.window_size),
        step_size=int(args.window_size/2),
        score_threshold=float(args.score_threshold),
        n_threads=int(args.threads)
    )
    
    pipeline.run_pipeline(args.input_fasta, args.output_bed)

if __name__ == "__main__":
    main()
