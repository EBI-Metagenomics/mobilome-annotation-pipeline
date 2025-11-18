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

from textwrap import dedent as _

from gff_mapping import (
    sort_gff_file,
    mobilome_parser,
    gff_updater,
)


def test_sort_gff_with_fasta(tmp_path):
    """
    Test GFF sorting with FASTA section preserved at the end.

    This test ensures that:
    - GFF entries are sorted by contig and position
    - Headers remain at the top
    - FASTA section (##FASTA and sequences) is preserved at the end
    """
    test_file = tmp_path / "test_fasta.gff"
    gff_content = _(
        """
        ##gff-version 3
        ctg123\t.\texon\t5000\t5500\t.\t+\t.\tID=exon00004
        ctg123\t.\texon\t1300\t1500\t.\t+\t.\tID=exon00001
        ctg123\t.\texon\t3000\t3902\t.\t+\t.\tID=exon00003
        ##FASTA
        >ctg123
        CTTCTGGGCGTACCCGATTCTCGGAGAACTTGCCGCACCATTCCGCCTTG
        TGTTCATTGCTGCCTGCATGTTCATTGTCTACCTCGGCTACGTGTGGCTA
    """
    )
    test_file.write_text(gff_content)

    sort_gff_file(str(test_file))

    lines = test_file.read_text().splitlines(keepends=True)

    # Check header is first
    assert lines[0].strip() == "##gff-version 3"
    # Check entries are sorted by position
    assert "1300" in lines[1]
    assert "3000" in lines[2]
    assert "5000" in lines[3]
    # Check FASTA section is at the end
    assert lines[4].strip() == "##FASTA"
    assert lines[5].strip() == ">ctg123"
    assert "CTTCTGGGCGTACCCGATTCTCGGAGAACTTGCCGCACCATTCCGCCTTG" in lines[6]


def test_parse_mobilome_annotations(tmp_path):
    """
    Test parsing mobilome annotations from different prediction tools.

    Verifies that the parser correctly identifies and organizes:
    - Mobile genetic elements from different tools (geNomad, ISEScan, ICEfinder)
    - Contig grouping
    - Element type classification
    """
    test_file = tmp_path / "mobilome.gff"
    gff_content = """##gff-version 3
contig1\tgeNomad\tvirus\t1000\t2000\t.\t+\t.\tID=virus001
contig1\tISEScan\tIS\t3000\t4000\t.\t+\t.\tID=is001
contig2\tICEfinder\tICE\t5000\t6000\t.\t-\t.\tID=ice001
"""
    test_file.write_text(gff_content)

    proteins_annot, mobilome_annot, mges_dict, mob_types = mobilome_parser(
        str(test_file)
    )

    # Check mobilome_annot
    assert "contig1" in mobilome_annot
    assert "contig2" in mobilome_annot
    assert len(mobilome_annot["contig1"]) == 2
    assert len(mobilome_annot["contig2"]) == 1

    # Check mges_dict
    assert "contig1" in mges_dict
    assert (1000, 2000) in mges_dict["contig1"]
    assert (3000, 4000) in mges_dict["contig1"]

    # Check mob_types
    assert mob_types[("contig1", 1000, 2000)] == "virus"
    assert mob_types[("contig1", 3000, 4000)] == "IS"
    assert mob_types[("contig2", 5000, 6000)] == "ICE"


def test_parse_protein_annotations(tmp_path):
    """
    Test parsing protein annotations with VIRify extra attributes.

    Ensures that viphog and viphog_taxonomy annotations are correctly
    extracted from protein features.
    """
    test_file = tmp_path / "proteins.gff"
    gff_content = """##gff-version 3
contig1\tProdigal\tCDS\t1000\t2000\t.\t+\t0\tID=prot001;viphog=VOG0001;viphog_taxonomy=Viruses
contig1\tProdigal\tCDS\t3000\t4000\t.\t-\t0\tID=prot002;product=hypothetical
"""
    test_file.write_text(gff_content)

    proteins_annot, mobilome_annot, mges_dict, mob_types = mobilome_parser(
        str(test_file)
    )

    # Check proteins_annot
    key1 = ("contig1", "1000", "2000", "+")
    assert key1 in proteins_annot
    assert "viphog=VOG0001" in proteins_annot[key1]
    assert "viphog_taxonomy=Viruses" in proteins_annot[key1]

    # Protein without extra annotations should not be in proteins_annot
    key2 = ("contig1", "3000", "4000", "-")
    assert key2 not in proteins_annot


def test_update_with_mobilome(tmp_path):
    """
    Test the main GFF update workflow.

    Verifies that:
    - Mobilome elements are added to all output files
    - Passenger proteins are identified and labeled
    - Three output files are created (clean, extra, full)
    """
    user_gff = tmp_path / "user.gff"
    user_content = """##gff-version 3
contig1\tProdigal\tCDS\t1500\t1800\t.\t+\t0\tID=prot001
contig1\tProdigal\tCDS\t5000\t5500\t.\t+\t0\tID=prot002
"""
    user_gff.write_text(user_content)

    # Create mobilome annotations
    mobilome_annot = {
        "contig1": ["contig1\tgeNomad\tvirus\t1000\t2000\t.\t+\t.\tID=virus001"]
    }
    mges_dict = {"contig1": [(1000, 2000)]}
    mob_types = {("contig1", 1000, 2000): "virus"}
    proteins_annot = {}

    output_prefix = tmp_path / "output"
    gff_updater(
        str(user_gff),
        str(output_prefix),
        proteins_annot,
        mobilome_annot,
        mges_dict,
        mob_types,
    )

    # Check that output files were created
    assert (tmp_path / "output_user_mobilome_clean.gff").exists()
    assert (tmp_path / "output_user_mobilome_extra.gff").exists()
    assert (tmp_path / "output_user_mobilome_full.gff").exists()

    # Check clean file contains mobilome annotation and passenger protein
    clean_content = (tmp_path / "output_user_mobilome_clean.gff").read_text()

    assert "geNomad" in clean_content
    assert "virus" in clean_content
    assert "mge_location=virus" in clean_content


def test_passenger_protein_detection(tmp_path):
    """
    Test passenger protein detection based on 75% coverage threshold.

    Tests the core logic that identifies proteins within mobile elements:
    - Protein from 1500-1800 (length 300)
    - MGE from 1000-2000
    - Overlap: 300bp, coverage: 300/300 = 100% > 75% threshold
    - Therefore prot001 should be marked as passenger
    - prot002 (no overlap) should NOT be marked as passenger
    """
    user_gff = tmp_path / "user.gff"
    user_content = """##gff-version 3
contig1\tProdigal\tCDS\t1500\t1800\t.\t+\t0\tID=prot001
contig1\tProdigal\tCDS\t5000\t5500\t.\t+\t0\tID=prot002
"""
    user_gff.write_text(user_content)

    mobilome_annot = {
        "contig1": ["contig1\tgeNomad\tvirus\t1000\t2000\t.\t+\t.\tID=virus001"]
    }
    mges_dict = {"contig1": [(1000, 2000)]}
    mob_types = {("contig1", 1000, 2000): "virus"}
    proteins_annot = {}

    output_prefix = tmp_path / "output"
    gff_updater(
        str(user_gff),
        str(output_prefix),
        proteins_annot,
        mobilome_annot,
        mges_dict,
        mob_types,
    )

    lines = (tmp_path / "output_user_mobilome_clean.gff").read_text().splitlines()

    # Check that prot001 is marked as passenger (has mge_location)
    prot001_line = [l for l in lines if "prot001" in l][0]
    assert "mge_location=virus" in prot001_line

    # Check that prot002 is NOT marked as passenger
    prot002_line = [l for l in lines if "prot002" in l][0]
    assert "mge_location" not in prot002_line
