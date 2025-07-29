#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import pytest
import tempfile
import os
from unittest.mock import patch, mock_open, MagicMock
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import sys

# Add the script directory to path for importing
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

# Import the functions from the script
from ice_boundary_refinement import (
    parse_blast_uniprot,
    parse_merged_gff,
    get_DR,
    ICE_filter,
    get_feat,
    getnum,
    zill,
    find_max_distance,
    pos_tag,
    merge_tRNA,
    get_ICE,
    get_map,
    get_sequence,
    get_gc,
    main
)


class TestParseBlastUniprot:
    """Test parse_blast_uniprot function"""
    
    def test_parse_blast_uniprot_valid_input(self):
        """Test parsing valid UniProt BLAST results"""
        test_data = "query1\tsubject1\tproduct1\n" \
                   "query2\tsubject2\tproduct2\n" \
                   "query3\tsubject3\tproduct3\n"
        
        with patch("builtins.open", mock_open(read_data=test_data)):
            result = parse_blast_uniprot("test_file.txt")
        
        expected = {
            "query1": "product1",
            "query2": "product2", 
            "query3": "product3"
        }
        assert result == expected
    
    def test_parse_blast_uniprot_empty_file(self):
        """Test parsing empty file"""
        with patch("builtins.open", mock_open(read_data="")):
            result = parse_blast_uniprot("empty_file.txt")
        
        assert result == {}
    
    def test_parse_blast_uniprot_single_line(self):
        """Test parsing single line"""
        test_data = "single_query\tsingle_subject\tsingle_product\n"
        
        with patch("builtins.open", mock_open(read_data=test_data)):
            result = parse_blast_uniprot("single_line.txt")
        
        expected = {"single_query": "single_product"}
        assert result == expected


class TestParseMergedGff:
    """Test parse_merged_gff function"""
    
    def test_parse_merged_gff_with_trna_and_cds(self):
        """Test parsing GFF with tRNA and CDS entries"""
        test_gff = """contig_1\tprodigal\ttRNA\t100\t175\t.\t+\t.\tID=CP003200_00017;product=tRNA-Glu
contig_1\tprodigal\tCDS\t200\t500\t.\t+\t.\tID=CP003200_00018;ori_id=contig_1_17
contig_1\tprodigal\tCDS\t600\t900\t.\t-\t.\tID=CP003200_00019;ori_id=contig_1_18"""
        
        uniprot_dict = {"contig_1_17": "hypothetical protein", "contig_1_18": "integrase"}
        
        with patch("builtins.open", mock_open(read_data=test_gff)):
            names_map, trnadict, posdict, header, totalnum_dict, locusdict, prots_contigs = \
                parse_merged_gff("test.gff", uniprot_dict)
        
        assert "contig_1_17" in names_map
        assert "CP003200_00017" in trnadict
        assert "CP003200_00018" in posdict
        assert totalnum_dict["contig_1"] == 3
        assert "CP003200_00018" in locusdict
    
    def test_parse_merged_gff_empty_file(self):
        """Test parsing empty GFF file"""
        with patch("builtins.open", mock_open(read_data="")):
            result = parse_merged_gff("empty.gff", {})
        
        names_map, trnadict, posdict, header, totalnum_dict, locusdict, prots_contigs = result
        assert names_map == {}
        assert trnadict == {}
        assert posdict == {}
        assert totalnum_dict == {}
        assert locusdict == {}
        assert prots_contigs == {}


class TestGetDR:
    """Test get_DR function"""
    
    def test_get_dr_valid_input(self):
        """Test parsing valid direct repeat file"""
        test_data = "contig_1 100 120 500 520\n" \
                   "contig_1 200 220 600 620\n" \
                   "contig_2 300 320 700 720\n"
        
        with patch("builtins.open", mock_open(read_data=test_data)):
            result = get_DR("test_dr.tsv")
        
        expected = {
            "contig_1": ["100|120|500|520", "200|220|600|620"],
            "contig_2": ["300|320|700|720"]
        }
        assert result == expected
    
    def test_get_dr_empty_file(self):
        """Test parsing empty DR file"""
        with patch("builtins.open", mock_open(read_data="")):
            result = get_DR("empty_dr.tsv")
        
        assert result == {}


class TestICEFilter:
    """Test ICE_filter function"""
    
    def test_ice_filter_valid_input(self):
        """Test filtering valid MacSyFinder output"""
        test_data = """replicon_name	hit_id	gene_name	hit_pos	hit_sequence_length	used_in	system_id	system_occ	system_wholeness
Chromosome	prot_001	Phage_integrase	1	300	ICE_1	ICE_1	1	1.0
Chromosome	prot_002	T4SS_MOB_F	2	250	ICE_1	ICE_1	1	1.0
Chromosome	prot_003	Relaxase_MOB_F	3	200	IME_1	IME_1	1	1.0"""
        
        with patch("builtins.open", mock_open(read_data=test_data)):
            concat_array, filtered_lines = ICE_filter("test_macsyfinder.tsv")
        
        assert len(concat_array) > 0
        assert len(filtered_lines) > 0
    
    def test_ice_filter_no_valid_systems(self):
        """Test filtering with no valid systems"""
        test_data = """replicon_name	hit_id	gene_name	hit_pos	hit_sequence_length	used_in	system_id	system_occ	system_wholeness
Chromosome	prot_001	Phage_integrase	1	300	ICE_1	ICE_1	0	0.5"""
        
        with patch("builtins.open", mock_open(read_data=test_data)):
            concat_array, filtered_lines = ICE_filter("test_macsyfinder.tsv")
        
        assert len(concat_array) == 0
        assert len(filtered_lines) == 0


class TestGetFeat:
    """Test get_feat function"""
    
    def test_get_feat_integrase(self):
        """Test feature mapping for integrase"""
        assert get_feat("Phage_integrase") == "Integrase@Phage_integrase"
        assert get_feat("Recombinase") == "Integrase@Recombinase"
        assert get_feat("rve") == "Integrase@rve"
    
    def test_get_feat_relaxase(self):
        """Test feature mapping for relaxase"""
        assert get_feat("T4SS_MOB_F") == "Relaxase@MOB"
        assert get_feat("Relaxase_MOB_P") == "Relaxase@MOB_P"
    
    def test_get_feat_t4cp(self):
        """Test feature mapping for T4CP"""
        assert get_feat("t4cp_F") == "T4CP@F"
        assert get_feat("tcpA_F") == "T4CP@F"
    
    def test_get_feat_t4ss(self):
        """Test feature mapping for T4SS"""
        assert get_feat("FATA_F") == "T4SS@F"
        assert get_feat("FA_F") == "T4SS@F"
        assert get_feat("T4SS_virB4") == "T4SS@virB4"
    
    def test_get_feat_replication(self):
        """Test feature mapping for replication proteins"""
        assert get_feat("RepSAv2") == "Rep@RepSAv2"
        assert get_feat("DUF3631") == "Rep@DUF3631"
        assert get_feat("Prim-Pol") == "Rep@Prim-Pol"


class TestUtilityFunctions:
    """Test utility functions"""
    
    def test_getnum(self):
        """Test getnum function"""
        assert getnum("contig_00001") == 1
        assert getnum("contig_00123") == 123
        assert getnum("contig_01000") == 1000
    
    def test_zill(self):
        """Test zill function"""
        assert zill("contig", 1) == "contig_00001"
        assert zill("contig", 123) == "contig_00123"
        assert zill("contig", 1000) == "contig_01000"
    
    def test_find_max_distance(self):
        """Test find_max_distance function"""
        assert find_max_distance([1, 3, 5, 15, 17]) == (5, 15)
        assert find_max_distance([1, 2, 3, 4]) == (1, 2)
        assert find_max_distance([10]) is None
        assert find_max_distance([]) is None
    
    def test_pos_tag(self):
        """Test pos_tag function"""
        posdict = {
            "contig_00001": [100, 200, "+", "protein1"],
            "contig_00002": [300, 400, "-", "protein2"],
            "contig_00003": [500, 600, "+", "protein3"]
        }
        
        # Test start direction
        tICE, tfinal = pos_tag(150, posdict, 2, 2, 10, "s")
        assert tICE == 1
        assert tfinal == 1  # max(1, 1-5)
        
        # Test end direction
        tICE, tfinal = pos_tag(350, posdict, 2, 2, 10, "e")
        assert tICE == 2
        assert tfinal == 7  # min(10, 2+5)


class TestSequenceFunctions:
    """Test sequence-related functions"""
    
    @patch('ice_boundary_refinement.SeqIO.parse')
    def test_get_sequence(self, mock_parse):
        """Test get_sequence function"""
        # Mock SeqRecord
        mock_record = MagicMock()
        mock_record.id = "contig1"
        mock_record.seq = Seq("ATCGATCGATCGATCG")
        mock_parse.return_value = [mock_record]
        
        result = get_sequence("test.fasta", "contig1", 2, 8)
        assert result == "CGATCG"
    
    def test_get_gc(self):
        """Test get_gc function"""
        assert get_gc("ATCG") == "0.50"
        assert get_gc("AAAA") == "0.00"
        assert get_gc("CCCC") == "1.00"
        assert get_gc("") == "0.00"
        assert get_gc("atcg") == "0.50"  # Test case insensitivity


class TestIntegrationFunctions:
    """Test integration functions that combine multiple components"""
    
    def test_merge_trna_basic(self):
        """Test merge_tRNA function with basic input"""
        ice_dict = {
                'CP003200_00010': 'Integrase@Phage_integrase',
                'CP003200_00011': 'T4SS@T_virB1',
                'CP003200_00012': 'T4SS@T_virB2',
            }

        dr_dict = {"contig_1": ["1000|1020|5000|5020"]}
        trnadict = {"CP003200_00009": [3433680, 3433710, "+", "tRNA-Ala"]}
        posdict = {
                'CP003200_00009': [3433680, 3433710, "+", "tRNA-Ala"],
                'CP003200_00010': [3433717, 3434979, '+', 'Prophage integrase IntA'],
                'CP003200_00011': [3468340, 3469050, '+', ''],
                'CP003200_00012': [3469050, 3469343, '+', 'Type IV secretion system protein PtlA']
                }
        totalnum_dict = {"contig_1": 15}
        locusdict = {"CP003200_00010": "CP003200_00010"}
        prots_contigs = {"CP003200_00010": "contig_1"}
        
        listgff = [trnadict, posdict, "CP003200", totalnum_dict, locusdict]
        
        result = merge_tRNA("ICE1", ice_dict, dr_dict, listgff, prots_contigs)
        
        assert len(result) == 13  # Should return tuple with 13 elements
        assert result[0] == "contig_1"  # contig name


class TestMainFunction:
    """Test main function and argument parsing"""
    
    @patch('ice_boundary_refinement.get_map')
    @patch('ice_boundary_refinement.get_ICE')
    @patch('ice_boundary_refinement.get_DR')
    @patch('ice_boundary_refinement.parse_merged_gff')
    @patch('ice_boundary_refinement.parse_blast_uniprot')
    @patch('argparse.ArgumentParser.parse_args')
    def test_main_function_flow(self, mock_args, mock_parse_blast, mock_parse_gff, 
                               mock_get_dr, mock_get_ice, mock_get_map):
        """Test main function execution flow"""
        # Mock arguments
        mock_args.return_value = MagicMock(
            assembly="test.fasta",
            gff_file="test.gff",
            macsyfinder_out="test_macsyfinder.tsv",
            uniprot_annot="test_uniprot.txt",
            drs_tsv="test_drs.tsv",
            prefix="test_output"
        )
        
        # Mock return values
        mock_parse_blast.return_value = {}
        mock_parse_gff.return_value = ({}, {}, {}, "header", {}, {}, {})
        mock_get_dr.return_value = {}
        mock_get_ice.return_value = ({}, {}, {}, "header", {})
        
        # Run main function
        main()
        
        # Verify all functions were called
        mock_parse_blast.assert_called_once()
        mock_parse_gff.assert_called_once()
        mock_get_dr.assert_called_once()
        mock_get_ice.assert_called_once()
        mock_get_map.assert_called_once()


class TestErrorHandling:
    """Test error handling and edge cases"""
    
    def test_parse_blast_uniprot_file_not_found(self):
        """Test handling of missing file"""
        with patch("builtins.open", side_effect=FileNotFoundError):
            with pytest.raises(FileNotFoundError):
                parse_blast_uniprot("nonexistent_file.txt")
    
    def test_get_dr_malformed_line(self):
        """Test handling of malformed DR file"""
        test_data = "contig1 100 120\n"  # Missing columns
        
        with patch("builtins.open", mock_open(read_data=test_data)):
            with pytest.raises(ValueError):
                get_DR("malformed_dr.tsv")
    
    def test_getnum_invalid_format(self):
        """Test getnum with invalid gene ID format"""
        with pytest.raises((ValueError, IndexError)):
            getnum("invalid_gene_id")


class TestFileIO:
    """Test file I/O operations"""
    
    def test_output_file_creation(self):
        """Test that output files are created correctly"""
        with tempfile.TemporaryDirectory() as temp_dir:
            prefix = os.path.join(temp_dir, "test")
            
            # Mock the required data structures
            drs_ice_dict = {
                "ICE1": [
                    "contig_1", 
                    '3433540', 
                    '3433556', 
                    '3495689', 
                    '3495705', 
                    10, 
                    20, 
                    10, 
                    20, 
                    [
                        [
                            3433480, 
                            3433555, 
                            '+', 
                            'tRNA-Asn'
                        ]
                    ], {
                        'CP003200_00001': 'CP003200_00001', 
                        'CP003200_00002': 'CP003200_00002', 
                        'CP003200_00003': 'CP003200_00003', 
                        'CP003200_00004': 'CP003200_00004', 
                        'CP003200_00005': 'CP003200_00005', 
                        'CP003200_00006': 'CP003200_00006', 
                        'CP003200_00007': 'CP003200_00007', 
                        'CP003200_00008': 'CP003200_00008', 
                        'CP003200_00009': 'CP003200_00009', 
                        'CP003200_00010': 'CP003200_00010', 
                        'CP003200_00011': 'CP003200_00011', 
                        'CP003200_00012': 'CP003200_00012', 
                        'CP003200_00013': 'CP003200_00013', 
                        'CP003200_00014': 'CP003200_00014', 
                        'CP003200_00015': 'CP003200_00015', 
                        'CP003200_00016': 'CP003200_00016', 
                        'CP003200_00017': 'CP003200_00017', 
                        'CP003200_00018': 'CP003200_00018', 
                        'CP003200_00019': 'CP003200_00019', 
                        'CP003200_00020': 'CP003200_00020', 
                        'CP003200_00021': 'CP003200_00021', 
                        'CP003200_00022': 'CP003200_00022', 
                        'CP003200_00023': 'CP003200_00023', 
                        'CP003200_00024': 'CP003200_00024', 
                        'CP003200_00025': 'CP003200_00025'
                    }]
            }

            genes_icedict = {
                'ICE1': {
                    'CP003200_00010': 'Integrase@Phage_integrase', 
                    'CP003200_00011': 'T4SS@T_virB1', 
                    'CP003200_00012': 'T4SS@T_virB2', 
                    'CP003200_00013': 'T4SS@virb4', 
                    'CP003200_00015': 'T4SS@T_virB6', 
                    'CP003200_00016': 'T4SS@T_virB8', 
                    'CP003200_00017': 'T4SS@T_virB9', 
                    'CP003200_00018': 'T4SS@T_virB10', 
                    'CP003200_00019': 'T4SS@T_virB11', 
                    'CP003200_00020': 'T4CP@t4cp2', 
                    'CP003200_00021': 'Relaxase@MOBC'
                }
            }
            
            posdict = {
                'CP003200_00010': [3433717, 3434979, '+', 'Prophage integrase IntA'],
                'CP003200_00011': [3468340, 3469050, '+', ''],
                'CP003200_00012': [3469050, 3469343, '+', 'Type IV secretion system protein PtlA'],
                'CP003200_00013': [3469356, 3472094, '+', 'Type IV secretion system protein virB4'],
                'CP003200_00014': [3473073, 3474146, '+', ''],
                'CP003200_00015': [3474368, 3475051, '+', 'Type IV secretion system protein PtlE'],
                'CP003200_00016': [3475051, 3475956, '+', 'Type IV secretion system protein virB9'],
                'CP003200_00017': [3476000, 3477250, '+', ''],
                'CP003200_00018': [3477240, 3478265, '+', 'Type IV secretion system protein VirB11'],
                'CP003200_00019': [3480058, 3481917, '+', ''],
                'CP003200_00020': [3481927, 3482673, '+', '']
            }

            header = "CP003200"

            infodict = {'ICE1': {'mob': ['MOBC'], 'mpf': ['typeT']}}

            assembly = "test.fasta"
            
            with patch('ice_boundary_refinement.get_sequence', return_value="ATCG"):
                with patch('ice_boundary_refinement.get_gc', return_value="0.50"):
                    get_map(drs_ice_dict, genes_icedict, posdict, header, 
                           infodict, assembly, prefix)
            
            # Check if output files were created
            assert os.path.exists(f"{prefix}_ices.tsv")
            assert os.path.exists(f"{prefix}_ice_genes.tsv")


if __name__ == "__main__":
    pytest.main([__file__, "-v"])

