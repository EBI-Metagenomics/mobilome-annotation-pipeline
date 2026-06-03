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

import gzip
import os
import sys
import subprocess
from pathlib import Path

import pytest

FIXTURES = Path(__file__).parent.parent / "fixtures" / "integrator"
BIN_DIR = Path(__file__).parent.parent.parent / "bin"
INTEGRATOR = BIN_DIR / "mge_integrator.py"


def _run_integrator(extra_args: list[str], tmp_path: Path) -> subprocess.CompletedProcess:
    env = os.environ.copy()
    pythonpath = str(BIN_DIR)
    if "PYTHONPATH" in env:
        pythonpath = pythonpath + ":" + env["PYTHONPATH"]
    env["PYTHONPATH"] = pythonpath

    cmd = [
        sys.executable, str(INTEGRATOR),
        "--gff_file", str(FIXTURES / "test_merged.gff"),
        "--map", str(FIXTURES / "test_contigID.map"),
        "--iss_tsv", str(FIXTURES / "test_1kb_contigs.fasta.tsv"),
        "--inf_tsv", str(FIXTURES / "contig_dummy.summary"),
        "--inf_gbks", str(FIXTURES / "contig_dummy.gbk"),
        "--icf_tsv",
        "--geno_out", str(FIXTURES / "test_5kb_contigs_virus_summary.tsv"),
        "--geno_plas", str(FIXTURES / "test_5kb_contigs_plasmid_summary.tsv"),
        "--comp_bed",
        "--checkv_genomad", str(FIXTURES / "quality_summary.tsv"),
        "--prefix", str(tmp_path / "test"),
    ] + extra_args
    return subprocess.run(cmd, capture_output=True, text=True, env=env)


def test_read_checkv_quality():
    """read_checkv_quality correctly parses the fixture quality_summary.tsv."""
    from mge_integrator import read_checkv_quality

    quality = read_checkv_quality(str(FIXTURES / "quality_summary.tsv"))

    assert set(quality.keys()) == {
        "contig_1", "contig_2", "contig_3", "contig_4",
        "contig_5", "contig_6", "contig_7",
    }

    contig3 = quality["contig_3"]
    assert "checkv_quality=Medium-quality" in contig3
    assert "checkv_miuvig_quality=Genome-fragment" in contig3
    assert "checkv_provirus=Yes" in contig3
    assert "checkv_viral_genes=21" in contig3
    assert "checkv_kmer_freq=1.0" in contig3

    contig1 = quality["contig_1"]
    assert "checkv_quality=Low-quality" in contig1
    assert "checkv_provirus=No" in contig1
    assert "checkv_viral_genes=2" in contig1


def test_mge_integrator_creates_output(tmp_path):
    """mge_integrator.py produces the expected output files without errors."""
    result = _run_integrator([], tmp_path)
    assert result.returncode == 0, result.stderr

    assert (tmp_path / "test_mobilome.gff.gz").exists()
    assert (tmp_path / "test_overlap_report.txt").exists()
    assert (tmp_path / "test_discarded_mge.txt").exists()


def test_mge_integrator_gff_matches_fixture(tmp_path):
    """mge_integrator.py GFF output matches the checked-in fixture."""
    result = _run_integrator([], tmp_path)
    assert result.returncode == 0, result.stderr

    with gzip.open(tmp_path / "test_mobilome.gff.gz", "rt") as fh:
        actual = fh.read()
    with gzip.open(FIXTURES / "expected/test_mobilome.gff.gz", "rt") as fh:
        expected = fh.read()

    assert actual == expected


def test_mge_integrator_all_viruses_present(tmp_path):
    """Quality-passing contigs from the virus summary appear as viral_sequence entries.

    MGYG000518629_154 (contig_2, 0 viral genes, Low-quality) is intentionally
    absent — it does not pass the quality_decision filter.

    Every contig from the virus summary appears as a viral_sequence in the output.
    """
    result = _run_integrator([], tmp_path)
    assert result.returncode == 0, result.stderr

    with gzip.open(tmp_path / "test_mobilome.gff.gz", "rt") as fh:
        content = fh.read()

    passing_contigs = [
        "MGYG000518629_154",
    ]
    expected_contigs = [
        "MGYG000535607_62",
        "MGYG000518621_235",
        "MGYG000518644_131",
        "MGYG000535607_58",
        "MGYG000535607_34",
        "MGYG000518629_136",
    ]
    for contig_id in passing_contigs:
        assert contig_id in content, f"{contig_id} missing from mobilome GFF"

    assert "MGYG000518629_154" not in content, "Low-quality contig (0 viral genes) should be filtered"

    for contig_id in expected_contigs:
        assert contig_id in content, f"{contig_id} missing from mobilome GFF"


def test_mge_integrator_checkv_attributes_in_gff(tmp_path):
    """CheckV quality attributes are embedded in viral_sequence GFF entries."""
    result = _run_integrator([], tmp_path)
    assert result.returncode == 0, result.stderr

    with gzip.open(tmp_path / "test_mobilome.gff.gz", "rt") as fh:
        lines = fh.readlines()

    viral_lines = [l for l in lines if "\tviral_sequence\t" in l]
    assert viral_lines, "No viral_sequence lines found"

    for line in viral_lines:
        assert "checkv_quality=" in line
        assert "checkv_provirus=" in line
        assert "checkv_miuvig_quality=" in line
        assert "checkv_viral_genes=" in line
        assert "checkv_kmer_freq=" in line


def test_mge_integrator_provirus_contig(tmp_path):
    """contig_3 (provirus=Yes, Medium-quality) is correctly annotated in the output."""
    result = _run_integrator([], tmp_path)
    assert result.returncode == 0, result.stderr

    with gzip.open(tmp_path / "test_mobilome.gff.gz", "rt") as fh:
        content = fh.read()

    contig3_line = next(
        l for l in content.splitlines() if "MGYG000535607_34" in l and "\tviral_sequence\t" in l
    )
    assert "checkv_provirus=Yes" in contig3_line
    assert "checkv_quality=Medium-quality" in contig3_line


def test_mge_integrator_empty_plasmid_summary(tmp_path):
    """Pipeline runs cleanly when the plasmid summary has no entries."""
    result = _run_integrator([], tmp_path)
    assert result.returncode == 0, result.stderr

    with gzip.open(tmp_path / "test_mobilome.gff.gz", "rt") as fh:
        content = fh.read()

    assert "\tplasmid\t" not in content


def test_mge_integrator_no_overlap(tmp_path):
    """Overlap report is created and contains only a header when no overlaps exist."""
    result = _run_integrator([], tmp_path)
    assert result.returncode == 0, result.stderr

    overlap_report = (tmp_path / "test_overlap_report.txt").read_text()
    lines = [l for l in overlap_report.splitlines() if l.strip()]
    assert lines[0].startswith("contig\t")
    assert len(lines) == 1, f"Unexpected overlap entries: {lines[1:]}"
