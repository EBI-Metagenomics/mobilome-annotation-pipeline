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

import pytest

from map_tools.genomad_parser import genomad_viral

# Header + columns mirror a real geNomad virus_summary.tsv. Only seq_name (0),
# length (1), coordinates (3), virus_score (6) and taxonomy (10) are consumed.
HEADER = (
    "seq_name\tlength\ttopology\tcoordinates\tn_genes\tgenetic_code\t"
    "virus_score\tfdr\tn_hallmarks\tmarker_enrichment\ttaxonomy"
)
TAXONOMY = "Viruses;Duplodnaviria;Heunggongvirae;Uroviricota;Caudoviricetes;;"

# A geNomad provirus call: seq_name carries the |provirus_<start>_<end> suffix
# and coordinates hold the sub-contig range. CheckV keys on the full seq_name.
PROVIRUS = (
    "contig_5|provirus_47020_74455\t27436\tProvirus\t47020-74455\t20\t11\t"
    f"0.9996\tNA\t4\t22.0675\t{TAXONOMY}"
)
# A whole-contig viral sequence (no provirus suffix; coordinates NA).
VIRAL = f"contig_14\t5874\tNo terminal repeats\tNA\t17\t11\t0.9991\tNA\t0\t5.8542\t{TAXONOMY}"
# Below the 0.8 virus_score cutoff: skipped before any CheckV lookup.
LOW_SCORE = f"contig_35\t14310\tNo terminal repeats\tNA\t21\t11\t0.7500\tNA\t0\t1.7170\t{TAXONOMY}"


def _checkv(quality, viral_genes, kmer_freq=1.0):
    """Build a CheckV attribute string as read_checkv_quality emits it."""
    return (
        f"checkv_kmer_freq={kmer_freq};checkv_miuvig_quality=Genome-fragment;"
        f"checkv_provirus=Yes;checkv_quality={quality};checkv_viral_genes={viral_genes}"
    )


def _write_summary(tmp_path, rows):
    path = tmp_path / "virus_summary.tsv"
    path.write_text(HEADER + "\n" + "\n".join(rows) + "\n")
    return str(path)


def test_provirus_quality_lookup_uses_full_seqname(tmp_path):
    """A provirus is matched against its CheckV entry by the full seq_name,
    stored on the parent contig at the provirus coordinates, with CheckV
    attributes embedded."""
    geno_out = _write_summary(tmp_path, [PROVIRUS, VIRAL])
    quality = {
        "contig_5|provirus_47020_74455": _checkv("Medium-quality", 21),
        "contig_14": _checkv("High-quality", 5),
    }

    mge_data = genomad_viral(geno_out, {}, quality)

    contig, description, coord = mge_data["vir1_1"]
    assert contig == "contig_5"                       # parent contig, suffix stripped
    assert coord == (47020, 74455)                    # provirus sub-contig range
    assert "mobile_element_type=prophage" in description
    assert "checkv_quality=Medium-quality" in description
    assert "checkv_provirus=Yes" in description

    contig2, description2, coord2 = mge_data["vir1_2"]
    assert contig2 == "contig_14"
    assert coord2 == (1, 5874)
    assert "mobile_element_type=viral_sequence" in description2


def test_provirus_requires_full_seqname_key(tmp_path):
    """Keying CheckV by the stripped contig name must NOT match a provirus —
    guards against regressing to quality.get(contig)."""
    geno_out = _write_summary(tmp_path, [PROVIRUS])
    quality = {"contig_5": _checkv("High-quality", 5)}  # stripped key (wrong)

    with pytest.raises(SystemExit) as excinfo:
        genomad_viral(geno_out, {}, quality)
    assert str(excinfo.value) == "No checkV values for record contig_5|provirus_47020_74455"


def test_missing_checkv_entry_exits(tmp_path):
    """A kept record with no CheckV entry aborts with a clear message."""
    geno_out = _write_summary(tmp_path, [VIRAL])

    with pytest.raises(SystemExit) as excinfo:
        genomad_viral(geno_out, {}, {})
    assert str(excinfo.value) == "No checkV values for record contig_14"


def test_low_quality_record_is_retained_with_attrs(tmp_path):
    """A record passing the geNomad score cutoff is kept regardless of its
    CheckV quality tier; the CheckV values are embedded in the attributes so the
    user can judge them. We no longer filter on CheckV (database bias)."""
    geno_out = _write_summary(tmp_path, [VIRAL])
    quality = {"contig_14": _checkv("Low-quality", 0)}

    mge_data = genomad_viral(geno_out, {}, quality)

    assert "vir1_1" in mge_data
    _, description, _ = mge_data["vir1_1"]
    assert "checkv_quality=Low-quality" in description
    assert "checkv_viral_genes=0" in description


def test_low_score_record_skips_checkv_lookup(tmp_path):
    """Records below the 0.8 virus_score cutoff are skipped before the CheckV
    lookup, so a missing entry does not raise."""
    geno_out = _write_summary(tmp_path, [LOW_SCORE])

    mge_data = genomad_viral(geno_out, {}, {})
    assert mge_data == {}
