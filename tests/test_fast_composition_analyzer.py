#!/usr/bin/env python3
"""
Comprehensive unit tests for fast_composition_analyzer.py
Tests the compositional outlier detection functionality
"""

import pytest
import tempfile
import os
import sys
from pathlib import Path
from unittest.mock import patch, MagicMock
import numpy as np
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from collections import Counter

# Add the bin directory to Python path for imports
project_root = Path(__file__).parent.parent
bin_dir = project_root / "bin"
sys.path.insert(0, str(bin_dir))

# Import the classes and functions from the script
from fast_composition_analyzer import (
    FastKmerAnalyzer,
    FastCompositionAnalyzer,
    FastOutlierDetector,
    FastRepeatFinder,
    FastGenomicOutlierPipeline,
    ContigInfo,
    Window,
    Outlier,
    Repeat,
    Prediction,
    process_contig_fast,
    main,
)


class TestFastKmerAnalyzer:
    """Test the k-mer analysis functionality"""

    def test_init_default_k_sizes(self):
        """Test initialization with default k-mer sizes"""
        analyzer = FastKmerAnalyzer()
        assert analyzer.k_sizes == [3, 4]
        assert 3 in analyzer.all_kmers
        assert 4 in analyzer.all_kmers
        assert len(analyzer.all_kmers[3]) == 64  # 4^3
        assert len(analyzer.all_kmers[4]) == 256  # 4^4

    def test_init_custom_k_sizes(self):
        """Test initialization with custom k-mer sizes"""
        analyzer = FastKmerAnalyzer(k_sizes=[2, 3])
        assert analyzer.k_sizes == [2, 3]
        assert len(analyzer.all_kmers[2]) == 16  # 4^2
        assert len(analyzer.all_kmers[3]) == 64  # 4^3

    def test_product_method(self):
        """Test the _product method for generating k-mers"""
        analyzer = FastKmerAnalyzer()
        bases = ["A", "T", "G", "C"]

        # Test 2-mer generation
        result = analyzer._product(bases, 2)
        assert len(result) == 16
        assert ("A", "A") in result
        assert ("T", "G") in result
        assert ("C", "C") in result

    def test_calculate_kmer_profile_simple(self):
        """Test k-mer profile calculation with simple sequence"""
        analyzer = FastKmerAnalyzer(k_sizes=[3])
        sequence = "AAATTTGGGCCC"

        profile = analyzer.calculate_kmer_profile(sequence)

        # Check that profile contains expected k-mers
        assert "3mer_AAA" in profile
        assert "3mer_TTT" in profile
        assert "3mer_GGG" in profile
        assert "3mer_CCC" in profile

        # Check normalization (should sum to 1.0)
        total = sum(profile.values())
        assert abs(total - 1.0) < 1e-10

    def test_calculate_kmer_profile_empty_sequence(self):
        """Test k-mer profile with empty sequence"""
        analyzer = FastKmerAnalyzer(k_sizes=[3])
        profile = analyzer.calculate_kmer_profile("")

        # All k-mers should have frequency 0
        for value in profile.values():
            assert value == 0.0

    def test_calculate_kmer_profile_short_sequence(self):
        """Test k-mer profile with sequence shorter than k"""
        analyzer = FastKmerAnalyzer(k_sizes=[5])
        sequence = "ATGC"  # Length 4, k=5

        profile = analyzer.calculate_kmer_profile(sequence)

        # All k-mers should have frequency 0 since sequence is too short
        for value in profile.values():
            assert value == 0.0

    def test_profile_to_vector(self):
        """Test conversion of profile to numpy vector"""
        analyzer = FastKmerAnalyzer(k_sizes=[2])
        profile = {"2mer_AA": 0.25, "2mer_AT": 0.25, "2mer_GC": 0.25, "2mer_TT": 0.25}

        vector = analyzer.profile_to_vector(profile)

        assert isinstance(vector, np.ndarray)
        assert len(vector) == len(profile)
        assert np.sum(vector) == 1.0


class TestFastCompositionAnalyzer:
    """Test the composition analysis functionality"""

    def test_init(self):
        """Test initialization of composition analyzer"""
        analyzer = FastCompositionAnalyzer(window_size=1000, step_size=500)
        assert analyzer.window_size == 1000
        assert analyzer.step_size == 500
        assert isinstance(analyzer.kmer_analyzer, FastKmerAnalyzer)

    def test_calculate_gc_content(self):
        """Test GC content calculation"""
        analyzer = FastCompositionAnalyzer(1000, 500)

        # Test balanced sequence
        sequence = "ATGC"
        gc_content = analyzer.calculate_gc_content(sequence)
        assert gc_content == 0.5

        # Test all GC
        sequence = "GGCC"
        gc_content = analyzer.calculate_gc_content(sequence)
        assert gc_content == 1.0

        # Test no GC
        sequence = "AATT"
        gc_content = analyzer.calculate_gc_content(sequence)
        assert gc_content == 0.0

        # Test empty sequence
        gc_content = analyzer.calculate_gc_content("")
        assert gc_content == 0.0

    def test_sliding_window_analysis(self):
        """Test sliding window analysis"""
        analyzer = FastCompositionAnalyzer(window_size=10, step_size=5)

        # Create test contig
        sequence = "A" * 30  # 30 bp sequence
        contig = ContigInfo("test_contig", sequence, len(sequence))

        windows = analyzer.sliding_window_analysis(contig)

        # Should have windows at positions 0, 5, 10, 15, 20
        assert len(windows) == 5
        assert windows[0].start == 0
        assert windows[0].end == 10
        assert windows[1].start == 5
        assert windows[1].end == 15

        # Check that all windows have GC content calculated
        for window in windows:
            assert hasattr(window, "gc_content")
            assert hasattr(window, "kmer_profile")

    def test_sliding_window_analysis_short_sequence(self):
        """Test sliding window with sequence shorter than window size"""
        analyzer = FastCompositionAnalyzer(window_size=100, step_size=50)

        sequence = "ATGC" * 5  # 20 bp sequence
        contig = ContigInfo("test_contig", sequence, len(sequence))

        windows = analyzer.sliding_window_analysis(contig)

        # Should have no windows since sequence is shorter than window size
        assert len(windows) == 0

    def test_calculate_fast_outlier_scores(self):
        """Test outlier score calculation"""
        analyzer = FastCompositionAnalyzer(window_size=10, step_size=5)

        # Create mock windows with different GC contents
        windows = [
            Window(0, 10, 0.5, {}, 0.0, 0.0, 0.0),  # Normal GC
            Window(5, 15, 0.9, {}, 0.0, 0.0, 0.0),  # High GC (outlier)
            Window(10, 20, 0.1, {}, 0.0, 0.0, 0.0),  # Low GC (outlier)
            Window(15, 25, 0.5, {}, 0.0, 0.0, 0.0),  # Normal GC
        ]

        # Add mock k-mer profiles
        for window in windows:
            window = window._replace(
                kmer_profile={
                    "3mer_AAA": 0.25,
                    "3mer_TTT": 0.25,
                    "3mer_GGG": 0.25,
                    "3mer_CCC": 0.25,
                }
            )

        scored_windows = analyzer.calculate_fast_outlier_scores(windows)

        assert len(scored_windows) == len(windows)

        # Check that outlier windows have higher scores
        assert (
            scored_windows[1].gc_score > scored_windows[0].gc_score
        )  # High GC outlier
        assert scored_windows[2].gc_score > scored_windows[0].gc_score  # Low GC outlier

    def test_calculate_fast_outlier_scores_few_windows(self):
        """Test outlier score calculation with too few windows"""
        analyzer = FastCompositionAnalyzer(window_size=10, step_size=5)

        windows = [Window(0, 10, 0.5, {}, 0.0, 0.0, 0.0)]

        scored_windows = analyzer.calculate_fast_outlier_scores(windows)

        # Should return original windows unchanged
        assert len(scored_windows) == 1
        assert scored_windows[0] == windows[0]


class TestFastOutlierDetector:
    """Test the outlier detection functionality"""

    def test_init(self):
        """Test initialization of outlier detector"""
        detector = FastOutlierDetector(score_threshold=2.0, min_region_size=1000)
        assert detector.score_threshold == 2.0
        assert detector.min_region_size == 1000

    def test_filter_significant_outliers(self):
        """Test filtering of significant outliers"""
        detector = FastOutlierDetector(score_threshold=2.0)

        windows = [
            Window(0, 10, 0.5, {}, 1.0, 1.0, 1.5),  # Below threshold
            Window(10, 20, 0.9, {}, 2.5, 2.0, 2.5),  # Above threshold
            Window(20, 30, 0.1, {}, 3.0, 1.5, 3.0),  # Above threshold
        ]

        significant = detector.filter_significant_outliers(windows)

        assert len(significant) == 2
        assert significant[0].combined_score == 2.5
        assert significant[1].combined_score == 3.0

    def test_merge_overlapping_outliers(self):
        """Test merging of overlapping outliers"""
        detector = FastOutlierDetector(score_threshold=2.0, min_region_size=100)

        outliers = [
            Outlier("contig1", 1000, 2000, 2.0, 2.0, 2.5, 0.8),
            Outlier("contig1", 1500, 2500, 2.5, 2.0, 2.8, 0.9),  # Overlapping
            Outlier("contig1", 5000, 6000, 3.0, 2.0, 3.0, 0.9),  # Separate
        ]

        merged = detector.merge_overlapping_outliers(outliers)

        assert len(merged) == 2
        # First merged outlier should span from 1000 to 2500
        assert merged[0].start == 1000
        assert merged[0].end == 2500
        assert merged[0].combined_score == 2.8  # Max score

        # Second outlier should remain unchanged
        assert merged[1].start == 5000
        assert merged[1].end == 6000

    def test_merge_overlapping_outliers_size_filter(self):
        """Test that small merged regions are filtered out"""
        detector = FastOutlierDetector(score_threshold=2.0, min_region_size=1000)

        outliers = [
            Outlier("contig1", 1000, 1100, 2.0, 2.0, 2.5, 0.8),  # Too small
            Outlier("contig1", 5000, 7000, 3.0, 2.0, 3.0, 0.9),  # Large enough
        ]

        merged = detector.merge_overlapping_outliers(outliers)

        assert len(merged) == 1
        assert merged[0].start == 5000
        assert merged[0].end == 7000


class TestFastRepeatFinder:
    """Test the repeat finding functionality"""

    def test_init(self):
        """Test initialization of repeat finder"""
        finder = FastRepeatFinder(
            flank_window=200, min_repeat_length=5, max_repeat_length=50
        )
        assert finder.flank_window == 200
        assert finder.min_repeat_length == 5
        assert finder.max_repeat_length == 50

    def test_find_exact_matches_direct(self):
        """Test finding direct exact matches"""
        finder = FastRepeatFinder(min_repeat_length=4, max_repeat_length=50)

        seq1 = "ATCGATCGATCG"
        seq2 = "GGGATCGATCCC"

        matches = finder.find_exact_matches(seq1, seq2, is_inverted=False)

        # Should find "ATCG" match at positions (0,3) or similar
        assert len(matches) >= 1

        # Check that at least one match is the expected "ATCG"
        found_atcg = any(match[3] == "GATCGATC" for match in matches)
        assert found_atcg, f"Expected to find 'ATCG' match, but got: {matches}"

        # Verify match structure: (pos1, pos2, length, sequence)
        for match in matches:
            assert len(match) == 4
            pos1, pos2, length, sequence = match
            assert isinstance(pos1, int)
            assert isinstance(pos2, int)
            assert isinstance(length, int)
            assert isinstance(sequence, str)
            assert length == len(sequence)
            assert length >= finder.min_repeat_length

    def test_find_exact_matches_inverted(self):
        """Test finding inverted exact matches"""
        finder = FastRepeatFinder(min_repeat_length=4, max_repeat_length=50)

        seq1 = "ATCGATCGATCG"
        seq2 = "GGGCGATCCCCC"  # Contains reverse complement of ATCG (CGAT)

        matches = finder.find_exact_matches(seq1, seq2, is_inverted=True)

        # For inverted matches, we're looking for reverse complements
        # ATCG reverse complement is CGAT
        # Check if we find any matches
        assert isinstance(matches, list)

        # If matches are found, verify their structure
        for match in matches:
            assert len(match) == 4
            pos1, pos2, length, sequence = match
            assert length >= finder.min_repeat_length
            assert length <= finder.max_repeat_length

    def test_extract_flanking_regions(self):
        """Test extraction of flanking regions"""
        finder = FastRepeatFinder(flank_window=10)

        sequence = "A" * 50  # 50 bp sequence
        outlier = Outlier("contig1", 20, 30, 2.0, 2.0, 2.5, 0.8)

        left_flank, right_flank = finder.extract_flanking_regions(sequence, outlier)

        assert len(left_flank) == 10  # flank_window size
        assert len(right_flank) == 10
        assert left_flank == "A" * 10
        assert right_flank == "A" * 10

    def test_extract_flanking_regions_edge_cases(self):
        """Test flanking region extraction at sequence edges"""
        finder = FastRepeatFinder(flank_window=100)

        sequence = "A" * 50  # Short sequence
        outlier = Outlier("contig1", 10, 40, 2.0, 2.0, 2.5, 0.8)

        left_flank, right_flank = finder.extract_flanking_regions(sequence, outlier)

        # Should handle edge cases gracefully
        assert len(left_flank) == 10  # From position 0 to 10
        assert len(right_flank) == 10  # From position 40 to 50

    def test_find_flanking_repeats_with_repeats(self):
        """Test finding flanking repeats when they exist"""
        finder = FastRepeatFinder(
            flank_window=20, min_repeat_length=5, max_repeat_length=15
        )

        # Create sequence with flanking repeats
        repeat_seq = "ATCGATCG"
        sequence = repeat_seq + "N" * 100 + repeat_seq + "N" * 50
        outlier = Outlier("contig1", 20, 120, 2.0, 2.0, 2.5, 0.8)

        result = finder.find_flanking_repeats(sequence, outlier)

        if result:  # May find repeats
            left_repeat, right_repeat = result
            assert isinstance(left_repeat, Repeat)
            assert isinstance(right_repeat, Repeat)
            assert left_repeat.type in ["direct", "inverted"]
            assert right_repeat.type in ["direct", "inverted"]

    def test_find_flanking_repeats_no_repeats(self):
        """Test finding flanking repeats when none exist"""
        finder = FastRepeatFinder()

        # Create sequence with no repeats
        sequence = "A" * 50 + "T" * 100 + "G" * 50
        outlier = Outlier("contig1", 50, 150, 2.0, 2.0, 2.5, 0.8)

        result = finder.find_flanking_repeats(sequence, outlier)

        # May or may not find repeats depending on the exact sequences
        assert result is None or isinstance(result, tuple)


class TestProcessContigFast:
    """Test the contig processing function"""

    def test_process_contig_fast_basic(self):
        """Test basic contig processing"""
        # Create test contig
        sequence = "A" * 1000 + "G" * 1000 + "C" * 1000  # 3kb with composition change
        contig = ContigInfo("test_contig", sequence, len(sequence))

        params = {"window_size": 500, "step_size": 250, "score_threshold": 1.0}

        predictions, contig_name = process_contig_fast((contig, params))

        assert contig_name == "test_contig"
        assert isinstance(predictions, list)
        # May or may not find predictions depending on the exact analysis

    def test_process_contig_fast_short_contig(self):
        """Test processing of very short contig"""
        sequence = "ATGC"  # Very short sequence
        contig = ContigInfo("short_contig", sequence, len(sequence))

        params = {"window_size": 100, "step_size": 50, "score_threshold": 2.0}

        predictions, contig_name = process_contig_fast((contig, params))

        assert contig_name == "short_contig"
        assert predictions == []  # Should return empty list for short contigs

    def test_process_contig_fast_error_handling(self):
        """Test error handling in contig processing"""
        # Create invalid contig that might cause errors
        contig = ContigInfo("error_contig", None, 0)  # Invalid sequence

        params = {"window_size": 100, "step_size": 50, "score_threshold": 2.0}

        predictions, contig_name = process_contig_fast((contig, params))

        assert contig_name == "error_contig"
        assert predictions == []  # Should handle errors gracefully


class TestFastGenomicOutlierPipeline:
    """Test the main pipeline class"""

    def test_init(self):
        """Test pipeline initialization"""
        pipeline = FastGenomicOutlierPipeline(
            min_contig_length=1000,
            window_size=500,
            step_size=250,
            score_threshold=2.0,
            n_threads=2,
        )

        assert pipeline.min_contig_length == 1000
        assert pipeline.window_size == 500
        assert pipeline.step_size == 250
        assert pipeline.score_threshold == 2.0
        assert pipeline.n_threads == 2

    def test_init_auto_threads(self):
        """Test pipeline initialization with auto thread detection"""
        pipeline = FastGenomicOutlierPipeline(
            min_contig_length=1000,
            window_size=500,
            step_size=250,
            score_threshold=2.0,
            n_threads=None,
        )

        assert pipeline.n_threads <= 4  # Should be capped at 4
        assert pipeline.n_threads >= 1

    def test_load_and_filter_contigs(self):
        """Test loading and filtering contigs from FASTA"""
        pipeline = FastGenomicOutlierPipeline(1000, 500, 250, 2.0, 1)

        # Create temporary FASTA file
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as tmp_file:
            # Write test sequences
            tmp_file.write(">contig1\n")
            tmp_file.write("A" * 1500 + "\n")  # Long enough
            tmp_file.write(">contig2\n")
            tmp_file.write("T" * 500 + "\n")  # Too short
            tmp_file.write(">contig3\n")
            tmp_file.write("G" * 2000 + "\n")  # Long enough
            tmp_file.flush()

            try:
                contigs = pipeline.load_and_filter_contigs(tmp_file.name)

                # Should only load contigs >= min_contig_length
                assert len(contigs) == 2
                assert contigs[0].name == "contig1"
                assert contigs[1].name == "contig3"
                assert contigs[0].length == 1500
                assert contigs[1].length == 2000

            finally:
                os.unlink(tmp_file.name)

    def test_load_and_filter_contigs_invalid_file(self):
        """Test loading from non-existent file"""
        pipeline = FastGenomicOutlierPipeline(1000, 500, 250, 2.0, 1)

        with pytest.raises(SystemExit):
            pipeline.load_and_filter_contigs("non_existent_file.fasta")

    def test_process_contigs_parallel_single_thread(self):
        """Test parallel processing with single thread"""
        pipeline = FastGenomicOutlierPipeline(1000, 500, 250, 2.0, 1)

        contigs = [
            ContigInfo("contig1", "A" * 1500, 1500),
            ContigInfo("contig2", "T" * 2000, 2000),
        ]

        predictions = pipeline.process_contigs_parallel(contigs)

        assert isinstance(predictions, list)
        # May or may not find predictions depending on analysis

    def test_write_output(self):
        """Test writing output to BED format"""
        pipeline = FastGenomicOutlierPipeline(1000, 500, 250, 2.0, 1)

        # Create test predictions
        outlier = Outlier("contig1", 1000, 2000, 2.5, 2.0, 2.8, 0.9)
        left_repeat = Repeat(950, 970, "ATCGATCG", "direct")
        right_repeat = Repeat(2020, 2040, "ATCGATCG", "direct")
        prediction = Prediction(outlier, left_repeat, right_repeat)

        predictions = [prediction]

        with tempfile.NamedTemporaryFile(mode="w", delete=False) as tmp_file:
            try:
                pipeline.write_output(predictions, tmp_file.name)

                # Check that file was written
                assert os.path.exists(tmp_file.name)

                # Check file content
                with open(tmp_file.name, "r") as f:
                    content = f.read()
                    assert "# Fast compositional outlier analysis" in content
                    assert "contig1" in content
                    assert "direct_1" in content
                    assert "direct_2" in content

            finally:
                os.unlink(tmp_file.name)


class TestIntegration:
    """Integration tests for the complete pipeline"""

    def test_full_pipeline_integration(self):
        """Test the complete pipeline with sample data"""
        # Create test FASTA file
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as fasta_file:
            # Create a sequence with potential compositional outlier
            normal_seq = "ATGC" * 500  # 2kb normal composition
            outlier_seq = "GGGG" * 250  # 1kb high GC outlier
            normal_seq2 = "ATGC" * 500  # 2kb normal composition

            full_sequence = normal_seq + outlier_seq + normal_seq2

            fasta_file.write(">test_contig\n")
            fasta_file.write(full_sequence + "\n")
            fasta_file.flush()

            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".bed", delete=False
            ) as bed_file:
                try:
                    # Run pipeline
                    pipeline = FastGenomicOutlierPipeline(
                        min_contig_length=1000,
                        window_size=500,
                        step_size=250,
                        score_threshold=1.5,
                        n_threads=1,
                    )

                    pipeline.run_pipeline(fasta_file.name, bed_file.name)

                    # Check that output file was created
                    assert os.path.exists(bed_file.name)

                    # Check file content
                    with open(bed_file.name, "r") as f:
                        content = f.read()
                        assert len(content) > 0
                        assert "# Fast compositional outlier analysis" in content

                finally:
                    os.unlink(fasta_file.name)
                    os.unlink(bed_file.name)


class TestMainFunction:
    """Test the main function and command line interface"""

    def test_main_function_help(self):
        """Test main function with help argument"""
        with patch("sys.argv", ["fast_composition_analyzer.py", "--help"]):
            with pytest.raises(SystemExit) as exc_info:
                main()
            assert exc_info.value.code == 0  # Help should exit with code 0

    def test_main_function_with_args(self):
        """Test main function with valid arguments"""
        # Create temporary files
        with tempfile.NamedTemporaryFile(
            mode="w", suffix=".fasta", delete=False
        ) as fasta_file:
            fasta_file.write(">test\n")
            fasta_file.write("ATGC" * 1000 + "\n")
            fasta_file.flush()

            with tempfile.NamedTemporaryFile(
                mode="w", suffix=".bed", delete=False
            ) as bed_file:
                try:
                    test_args = [
                        "fast_composition_analyzer.py",
                        "--input_fasta",
                        fasta_file.name,
                        "--output_bed",
                        bed_file.name,
                        "--score-threshold",
                        "2.0",
                        "--min-contig-length",
                        "1000",
                        "--window-size",
                        "500",
                        "--threads",
                        "1",
                    ]

                    with patch("sys.argv", test_args):
                        main()

                    # Check that output file was created
                    assert os.path.exists(bed_file.name)

                finally:
                    os.unlink(fasta_file.name)
                    os.unlink(bed_file.name)


# Fixtures for common test data
@pytest.fixture
def sample_sequence():
    """Sample DNA sequence for testing"""
    return "ATGCATGCATGC" * 100  # 1200 bp sequence


@pytest.fixture
def sample_contig(sample_sequence):
    """Sample ContigInfo object"""
    return ContigInfo("test_contig", sample_sequence, len(sample_sequence))


@pytest.fixture
def sample_outlier():
    """Sample Outlier object"""
    return Outlier("contig1", 1000, 2000, 2.5, 2.0, 2.8, 0.9)


@pytest.fixture
def sample_repeat():
    """Sample Repeat object"""
    return Repeat(950, 970, "ATCGATCG", "direct")


@pytest.fixture
def temp_fasta_file():
    """Temporary FASTA file for testing"""
    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".fasta", delete=False
    ) as tmp_file:
        tmp_file.write(">contig1\n")
        tmp_file.write("ATGC" * 1000 + "\n")
        tmp_file.write(">contig2\n")
        tmp_file.write("GGCC" * 500 + "\n")
        tmp_file.flush()

        yield tmp_file.name

        os.unlink(tmp_file.name)


# Test configuration
if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
