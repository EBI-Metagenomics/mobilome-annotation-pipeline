#!/usr/bin/env python

import sys
from pathlib import Path

import pytest

# Add the script directory to the path
script_dir = Path(__file__).parent.parent / "resources" / "usr" / "bin"
sys.path.insert(0, str(script_dir))

# Import functions from the script
from pathofact2_integrator import (  # noqa: E402
    _add_annotation_to_protein,
    _has_gff_record_9cols,
    _pred_support_has_data_rows,
    parse_cdd,
    parse_gff,
    parse_ips,
    parse_pathofact_support,
    validate_inputs,
)


@pytest.fixture
def tmp_files(tmp_path):
    """Create temporary test files with realistic data."""

    # PathOFact2 support file data
    pred_support_data = {
        "header": [
            "sequence_id",
            "detection_method",
            "support_value_type",
            "support_value",
            "vfdb_hit",
        ],
        "rows": [
            ["WP_000031783.1", "blastp", "evalue", "6.5e-248", "VFG046459"],
            ["WP_000242755.1", "blastp", "evalue", "2.58e-94", "VFG042734"],
            ["WP_001273183.1", "blastp", "evalue", "0.0", "VFG049144"],
            ["WP_000242755.1", "pathofact2_tox", "probability", "1.0", ""],
            ["WP_001273183.1", "pathofact2_tox", "probability", "1.0", ""],
            ["WP_000163771.1", "pathofact2_vf", "probability", "1.0", ""],
            ["WP_001273183.1", "pathofact2_vf", "probability", "0.9977883", ""],
        ],
    }

    # CDD annotation file data
    cdd_annot_data = {
        "header": [
            "query",
            "hit_type",
            "pssm_id",
            "from",
            "to",
            "evalue",
            "bitscore",
            "accession",
            "short_name",
            "incomplete",
            "superfamily_pssm_id",
        ],
        "rows": [
            [
                "WP_000031783.1",
                "Specific",
                "206671",
                "11",
                "203",
                "4.85688e-133",
                "377.694",
                "cd01884",
                "EF_Tu",
                "-",
                "476819",
            ],
            [
                "WP_000031783.1",
                "Specific",
                "294006",
                "300",
                "389",
                "1.02934e-65",
                "202.357",
                "cd03707",
                "EFTU_III",
                "-",
                "470675",
            ],
            [
                "WP_000031783.1",
                "Specific",
                "293898",
                "211",
                "297",
                "1.80452e-52",
                "168.466",
                "cd03697",
                "EFTU_II",
                "-",
                "445922",
            ],
            [
                "WP_000242755.1",
                "Specific",
                "237999",
                "8",
                "117",
                "3.94612e-35",
                "118.968",
                "cd00038",
                "CAP_ED",
                "-",
                "469590",
            ],
            [
                "WP_000242755.1",
                "Specific",
                "238044",
                "140",
                "207",
                "1.57469e-15",
                "66.9205",
                "cd00092",
                "HTH_CRP",
                "-",
                "481199",
            ],
        ],
    }

    # InterProScan annotation file data (no header)
    ips_annot_data = {
        "rows": [
            [
                "WP_000242755.1",
                "9a58754a5c561c0784b75360a962f090",
                "210",
                "CDD",
                "cd00092",
                "HTH_CRP",
                "140",
                "207",
                "1.14133E-15",
                "T",
                "15-01-2026",
                "-",
                "-",
                "-",
                "-",
            ],
            [
                "WP_000242755.1",
                "9a58754a5c561c0784b75360a962f090",
                "210",
                "CDD",
                "cd00038",
                "CAP_ED",
                "8",
                "117",
                "2.86014E-35",
                "T",
                "15-01-2026",
                "IPR000595",
                "Cyclic nucleotide-binding domain",
                "-",
                "Reactome:R-BTA-1296072",
            ],
            [
                "WP_000031783.1",
                "3ba5ec971eda5dc5f5c1c093a5df4e11",
                "394",
                "CDD",
                "cd03697",
                "EFTU_II",
                "211",
                "297",
                "1.30792E-52",
                "T",
                "15-01-2026",
                "IPR033720",
                "Elongation factor Tu, domain 2",
                "-",
                "Reactome:R-HSA-5389840",
            ],
            [
                "WP_000031783.1",
                "3ba5ec971eda5dc5f5c1c093a5df4e11",
                "394",
                "CDD",
                "cd01884",
                "EF_Tu",
                "11",
                "203",
                "0.0",
                "T",
                "15-01-2026",
                "IPR041709",
                "Elongation factor Tu (EF-Tu), GTP-binding domain",
                "-",
                "Reactome:R-HSA-5389840",
            ],
            [
                "WP_000031783.1",
                "3ba5ec971eda5dc5f5c1c093a5df4e11",
                "394",
                "CDD",
                "cd03707",
                "EFTU_III",
                "300",
                "389",
                "7.46066E-66",
                "T",
                "15-01-2026",
                "-",
                "-",
                "-",
                "-",
            ],
        ]
    }

    # GFF file data
    gff_data = {
        "header": [
            "##gff-version 3",
            "##sequence-region Contig001 1 100000",
            "##sequence-region Contig002 1 80000",
        ],
        "rows": [
            [
                "Contig001",
                "Prodigal",
                "CDS",
                "51000",
                "51837",
                ".",
                "+",
                "0",
                "ID=WP_000259031.1;Name=WP_000259031.1;product=hypothetical protein",
            ],
            [
                "Contig001",
                "Prodigal",
                "CDS",
                "87000",
                "87630",
                ".",
                "+",
                "0",
                "ID=WP_000242755.1;Name=WP_000242755.1;product=hypothetical protein",
            ],
            [
                "Contig002",
                "Prodigal",
                "CDS",
                "3000",
                "3744",
                ".",
                "+",
                "0",
                "ID=WP_000163771.1;Name=WP_000163771.1;product=hypothetical protein",
            ],
            [
                "Contig002",
                "Prodigal",
                "CDS",
                "38000",
                "39182",
                ".",
                "-",
                "0",
                "ID=WP_000031783.1;Name=WP_000031783.1;product=hypothetical protein",
            ],
            [
                "Contig003",
                "Prodigal",
                "CDS",
                "12000",
                "15102",
                ".",
                "+",
                "0",
                "ID=WP_001273183.1;Name=WP_001273183.1;product=hypothetical protein",
            ],
        ],
    }

    # Create pred_support.tsv
    pred_support = tmp_path / "pred_support.tsv"
    pred_support_content = "\t".join(pred_support_data["header"]) + "\n"
    pred_support_content += "\n".join(
        ["\t".join(row) for row in pred_support_data["rows"]]
    )
    pred_support.write_text(pred_support_content + "\n")

    # Create CDD annotation file
    cdd_annot = tmp_path / "cdd_annot.tsv"
    cdd_annot_content = "\t".join(cdd_annot_data["header"]) + "\n"
    cdd_annot_content += "\n".join(["\t".join(row) for row in cdd_annot_data["rows"]])
    cdd_annot.write_text(cdd_annot_content + "\n")

    # Create IPS annotation file (no header)
    ips_annot = tmp_path / "ips_annot.tsv"
    ips_annot_content = "\n".join(["\t".join(row) for row in ips_annot_data["rows"]])
    ips_annot.write_text(ips_annot_content + "\n")

    # Create GFF file
    gff_file = tmp_path / "input.gff"
    gff_content = "\n".join(gff_data["header"]) + "\n"
    gff_content += "\n".join(["\t".join(row) for row in gff_data["rows"]])
    gff_file.write_text(gff_content + "\n")

    # Create empty file
    empty_file = tmp_path / "empty.tsv"
    empty_file.write_text("")

    # Create header-only pred_support file
    header_only = tmp_path / "header_only.tsv"
    header_only.write_text(
        "sequence_id\tdetection_method\tsupport_value_type\tsupport_value\tvfdb_hit\n"
    )

    # Create GFF file with only comments (no 9-column records)
    gff_comments_only = tmp_path / "gff_comments_only.gff"
    gff_comments_only.write_text(
        "##gff-version 3\n##sequence-region Contig001 1 100000\n"
    )

    # Create GFF file with invalid format (not 9 columns)
    gff_invalid = tmp_path / "gff_invalid.gff"
    gff_invalid.write_text("##gff-version 3\nContig001\tProdigal\tCDS\n")

    # Create pred_support with blank lines
    pred_support_blank_lines = tmp_path / "pred_support_blank.tsv"
    pred_support_blank_content = "\t".join(pred_support_data["header"]) + "\n"
    pred_support_blank_content += "\n"  # blank line
    pred_support_blank_content += "\t\t\t\t\n"  # empty columns
    pred_support_blank_lines.write_text(pred_support_blank_content)

    return {
        "pred_support": str(pred_support),
        "cdd_annot": str(cdd_annot),
        "ips_annot": str(ips_annot),
        "gff_file": str(gff_file),
        "empty_file": str(empty_file),
        "header_only": str(header_only),
        "gff_comments_only": str(gff_comments_only),
        "gff_invalid": str(gff_invalid),
        "pred_support_blank_lines": str(pred_support_blank_lines),
        "tmp_path": tmp_path,
    }


class TestHelperFunctions:
    """Test private helper functions."""

    def test_pred_support_has_data_rows_valid(self, tmp_files):
        """Test that pred_support with data rows returns True."""
        result = _pred_support_has_data_rows(tmp_files["pred_support"])
        assert result is True

    def test_pred_support_has_data_rows_empty(self, tmp_files):
        """Test that empty file returns False."""
        result = _pred_support_has_data_rows(tmp_files["empty_file"])
        assert result is False

    def test_pred_support_has_data_rows_header_only(self, tmp_files):
        """Test that header-only file returns False."""
        result = _pred_support_has_data_rows(tmp_files["header_only"])
        assert result is False

    def test_pred_support_has_data_rows_blank_lines(self, tmp_files):
        """Test that file with only blank lines returns False."""
        result = _pred_support_has_data_rows(tmp_files["pred_support_blank_lines"])
        assert result is False

    def test_has_gff_record_9cols_valid(self, tmp_files):
        """Test that valid GFF with 9 columns returns True."""
        result = _has_gff_record_9cols(tmp_files["gff_file"])
        assert result is True

    def test_has_gff_record_9cols_comments_only(self, tmp_files):
        """Test that GFF with only comments returns False."""
        result = _has_gff_record_9cols(tmp_files["gff_comments_only"])
        assert result is False

    def test_has_gff_record_9cols_invalid_format(self, tmp_files):
        """Test that GFF with invalid format returns False."""
        result = _has_gff_record_9cols(tmp_files["gff_invalid"])
        assert result is False

    def test_add_annotation_to_protein_new_protein(self):
        """Test adding annotation to a new protein."""
        pathofact_attrs = {"protein1": "vfdb=VFG123;blastp_eval=1e-10"}
        proteins_annotation = {}

        _add_annotation_to_protein(
            "protein1", "cd12345", "Domain_Name", pathofact_attrs, proteins_annotation
        )

        assert "protein1" in proteins_annotation
        assert (
            proteins_annotation["protein1"]
            == "vfdb=VFG123;blastp_eval=1e-10;cdd=cd12345:Domain_Name"
        )

    def test_add_annotation_to_protein_existing_protein(self):
        """Test adding annotation to existing protein with CDD annotations."""
        pathofact_attrs = {"protein1": "vfdb=VFG123;blastp_eval=1e-10"}
        proteins_annotation = {
            "protein1": "vfdb=VFG123;blastp_eval=1e-10;cdd=cd11111:Domain1"
        }

        _add_annotation_to_protein(
            "protein1", "cd22222", "Domain2", pathofact_attrs, proteins_annotation
        )

        assert "cd11111:Domain1" in proteins_annotation["protein1"]
        assert "cd22222:Domain2" in proteins_annotation["protein1"]
        assert proteins_annotation["protein1"].count("cdd=") == 1
        assert "," in proteins_annotation["protein1"]

    def test_add_annotation_to_protein_not_in_pathofact(self):
        """Test that proteins not in pathofact_attrs are not annotated."""
        pathofact_attrs = {"protein1": "vfdb=VFG123;blastp_eval=1e-10"}
        proteins_annotation = {}

        _add_annotation_to_protein(
            "protein2", "cd12345", "Domain_Name", pathofact_attrs, proteins_annotation
        )

        assert "protein2" not in proteins_annotation


class TestValidateInputs:
    """Test input validation function with new validation logic."""

    def test_valid_inputs_all_files(self, tmp_files):
        """Test validation passes for all valid files."""
        to_validate = {
            "gff": tmp_files["gff_file"],
            "annotation": tmp_files["cdd_annot"],
            "pred_support": tmp_files["pred_support"],
        }
        result = validate_inputs(to_validate)
        assert len(result) == 0, "Valid files should pass validation"

    def test_missing_file(self, tmp_files):
        """Test validation fails for non-existent files."""
        to_validate = {
            "gff": "/nonexistent/file.gff",
            "annotation": tmp_files["cdd_annot"],
            "pred_support": tmp_files["pred_support"],
        }
        result = validate_inputs(to_validate)
        assert len(result) == 1
        assert "File not found" in result["gff"]

    def test_pred_support_empty_file(self, tmp_files):
        """Test validation fails for empty pred_support file."""
        to_validate = {
            "gff": tmp_files["gff_file"],
            "annotation": tmp_files["cdd_annot"],
            "pred_support": tmp_files["empty_file"],
        }
        result = validate_inputs(to_validate)
        assert len(result) == 1
        assert "pred_support" in result
        assert "No predicted proteins" in result["pred_support"]

    def test_pred_support_header_only(self, tmp_files):
        """Test validation fails for header-only pred_support file."""
        to_validate = {
            "gff": tmp_files["gff_file"],
            "annotation": tmp_files["cdd_annot"],
            "pred_support": tmp_files["header_only"],
        }
        result = validate_inputs(to_validate)
        assert len(result) == 1
        assert "pred_support" in result
        assert "No predicted proteins" in result["pred_support"]

    def test_gff_no_valid_records(self, tmp_files):
        """Test validation fails for GFF with no valid 9-column records."""
        to_validate = {
            "gff": tmp_files["gff_comments_only"],
            "annotation": tmp_files["cdd_annot"],
            "pred_support": tmp_files["pred_support"],
        }
        result = validate_inputs(to_validate)
        assert len(result) == 1
        assert "gff" in result
        assert "No valid GFF records" in result["gff"]

    def test_gff_not_validated_when_pred_support_empty(self, tmp_files):
        """Test that GFF validation is skipped when pred_support is empty."""
        # Even if GFF is invalid, validation should stop at pred_support
        to_validate = {
            "gff": tmp_files["gff_comments_only"],
            "annotation": tmp_files["cdd_annot"],
            "pred_support": tmp_files["empty_file"],
        }
        result = validate_inputs(to_validate)
        # Should only report pred_support error, not GFF error
        assert "pred_support" in result
        assert "gff" not in result

    def test_none_input(self):
        """Test validation handles None inputs."""
        to_validate = {
            "gff": None,
            "annotation": "some_file.tsv",
            "pred_support": "support.tsv",
        }
        result = validate_inputs(to_validate)
        assert len(result) >= 1
        assert "gff" in result
        assert result["gff"] == "No input provided"

    def test_empty_string_input(self):
        """Test validation handles empty string inputs."""
        to_validate = {
            "gff": "",
            "annotation": "some_file.tsv",
            "pred_support": "support.tsv",
        }
        result = validate_inputs(to_validate)
        assert len(result) >= 1
        assert "gff" in result
        assert result["gff"] == "No input provided"

    def test_annotation_file_can_be_empty(self, tmp_files):
        """Test that annotation file can be empty (no content validation)."""
        to_validate = {
            "gff": tmp_files["gff_file"],
            "annotation": tmp_files["empty_file"],
            "pred_support": tmp_files["pred_support"],
        }
        result = validate_inputs(to_validate)
        # Annotation file should pass even if empty
        assert len(result) == 0

    def test_multiple_invalid_files(self, tmp_files):
        """Test validation reports multiple errors."""
        to_validate = {
            "gff": tmp_files["gff_comments_only"],
            "annotation": tmp_files["cdd_annot"],
            "pred_support": tmp_files["header_only"],
        }
        result = validate_inputs(to_validate)
        # Only pred_support should be reported (stops early)
        assert "pred_support" in result


class TestParsePathofactSupport:
    """Test parsing of PathOFact2 support file."""

    def test_parse_support_file(self, tmp_files):
        """Test parsing support file with mixed detection methods."""
        result = parse_pathofact_support(tmp_files["pred_support"])

        # Check all proteins are present
        assert "WP_000031783.1" in result
        assert "WP_000242755.1" in result
        assert "WP_001273183.1" in result
        assert "WP_000163771.1" in result

        # Check blastp annotation format
        assert "vfdb=VFG046459" in result["WP_000031783.1"]
        assert "blastp_eval=6.5e-248" in result["WP_000031783.1"]

        # Check pathofact2 annotation format
        assert "pathofact2_tox_prob=1.0" in result["WP_000242755.1"]
        assert "pathofact2_vf_prob=1.0" in result["WP_000163771.1"]

        # Check multiple annotations are concatenated
        assert "vfdb=VFG042734" in result["WP_000242755.1"]
        assert "pathofact2_tox_prob=1.0" in result["WP_000242755.1"]

        # WP_001273183.1 has all three: blastp, tox, vf
        wp_001273183_annot = result["WP_001273183.1"]
        assert "vfdb=VFG049144" in wp_001273183_annot
        assert "pathofact2_tox_prob=1.0" in wp_001273183_annot
        assert "pathofact2_vf_prob=0.9977883" in wp_001273183_annot

    def test_parse_support_file_semicolon_separator(self, tmp_files):
        """Test that annotations are separated by semicolons."""
        result = parse_pathofact_support(tmp_files["pred_support"])

        # Check semicolon separation
        assert ";" in result["WP_000031783.1"]
        assert ";" in result["WP_000242755.1"]
        assert ";" in result["WP_001273183.1"]


class TestParseCDD:
    """Test parsing of CDD annotation file."""

    def test_parse_cdd_basic(self, tmp_files):
        """Test basic CDD parsing."""
        pathofact_attrs = parse_pathofact_support(tmp_files["pred_support"])
        result = parse_cdd(pathofact_attrs, tmp_files["cdd_annot"])

        # Check proteins with both pathofact and CDD annotations
        assert "WP_000031783.1" in result
        assert "WP_000242755.1" in result

        # Check CDD format: accession:short_name
        assert "cd01884:EF_Tu" in result["WP_000031783.1"]
        assert "cd00038:CAP_ED" in result["WP_000242755.1"]
        assert "cd00092:HTH_CRP" in result["WP_000242755.1"]

    def test_parse_cdd_multiple_hits(self, tmp_files):
        """Test CDD with multiple hits per protein."""
        pathofact_attrs = parse_pathofact_support(tmp_files["pred_support"])
        result = parse_cdd(pathofact_attrs, tmp_files["cdd_annot"])

        # WP_000031783.1 has 3 CDD hits
        wp_031783_annot = result["WP_000031783.1"]
        assert "cd01884:EF_Tu" in wp_031783_annot
        assert "cd03707:EFTU_III" in wp_031783_annot
        assert "cd03697:EFTU_II" in wp_031783_annot

        # Check they're comma-separated
        assert wp_031783_annot.count(",") == 2

    def test_parse_cdd_preserves_pathofact(self, tmp_files):
        """Test that CDD parsing preserves PathOFact annotations."""
        pathofact_attrs = parse_pathofact_support(tmp_files["pred_support"])
        result = parse_cdd(pathofact_attrs, tmp_files["cdd_annot"])

        # Check pathofact annotations are still present
        assert "vfdb=VFG046459" in result["WP_000031783.1"]
        assert "blastp_eval=6.5e-248" in result["WP_000031783.1"]
        assert ";cdd=" in result["WP_000031783.1"]

    def test_parse_cdd_only_annotates_pathofact_proteins(self, tmp_files):
        """Test that only proteins in pathofact_attrs get CDD annotations."""
        # Create minimal pathofact_attrs with only one protein
        pathofact_attrs = {"WP_000242755.1": "vfdb=VFG042734;blastp_eval=2.58e-94"}
        result = parse_cdd(pathofact_attrs, tmp_files["cdd_annot"])

        # Only WP_000242755.1 should be in result
        assert "WP_000242755.1" in result
        assert "WP_000031783.1" not in result

    def test_parse_cdd_uses_helper_function(self, tmp_files):
        """Test that parse_cdd uses _add_annotation_to_protein helper."""
        pathofact_attrs = {"WP_000031783.1": "vfdb=VFG046459;blastp_eval=6.5e-248"}
        result = parse_cdd(pathofact_attrs, tmp_files["cdd_annot"])

        # Should have CDD annotations added via helper
        assert "WP_000031783.1" in result
        assert "cdd=" in result["WP_000031783.1"]
        assert "cd01884:EF_Tu" in result["WP_000031783.1"]


class TestParseIPS:
    """Test parsing of InterProScan annotation file."""

    def test_parse_ips_basic(self, tmp_files):
        """Test basic IPS parsing."""
        pathofact_attrs = parse_pathofact_support(tmp_files["pred_support"])
        result = parse_ips(pathofact_attrs, tmp_files["ips_annot"])

        # Check proteins with both pathofact and IPS CDD annotations
        assert "WP_000031783.1" in result
        assert "WP_000242755.1" in result

    def test_parse_ips_filters_cdd_only(self, tmp_files):
        """Test that IPS parser only extracts CDD annotations."""
        pathofact_attrs = parse_pathofact_support(tmp_files["pred_support"])
        result = parse_ips(pathofact_attrs, tmp_files["ips_annot"])

        # Check CDD annotations are present
        assert "cd00092:HTH_CRP" in result["WP_000242755.1"]
        assert "cd00038:CAP_ED" in result["WP_000242755.1"]

        # Check format matches: accession:description
        assert "cd03697:EFTU_II" in result["WP_000031783.1"]
        assert "cd01884:EF_Tu" in result["WP_000031783.1"]
        assert "cd03707:EFTU_III" in result["WP_000031783.1"]

    def test_parse_ips_multiple_hits(self, tmp_files):
        """Test IPS with multiple CDD hits per protein."""
        pathofact_attrs = parse_pathofact_support(tmp_files["pred_support"])
        result = parse_ips(pathofact_attrs, tmp_files["ips_annot"])

        # WP_000031783.1 has 3 CDD annotations
        wp_031783_annot = result["WP_000031783.1"]
        assert wp_031783_annot.count(",") == 2  # comma-separated

        # WP_000242755.1 has 2 CDD annotations
        wp_000242755_annot = result["WP_000242755.1"]
        assert wp_000242755_annot.count(",") == 1

    def test_parse_ips_uses_helper_function(self, tmp_files):
        """Test that parse_ips uses _add_annotation_to_protein helper."""
        pathofact_attrs = {"WP_000242755.1": "vfdb=VFG042734;blastp_eval=2.58e-94"}
        result = parse_ips(pathofact_attrs, tmp_files["ips_annot"])

        # Should have CDD annotations added via helper
        assert "WP_000242755.1" in result
        assert "cdd=" in result["WP_000242755.1"]
        assert "cd00038:CAP_ED" in result["WP_000242755.1"]


class TestParseGFF:
    """Test GFF parsing and annotation integration."""

    def test_parse_gff_integration(self, tmp_files):
        """Test GFF parsing integrates annotations correctly."""
        pathofact_attrs = parse_pathofact_support(tmp_files["pred_support"])
        proteins_annotation = parse_cdd(pathofact_attrs, tmp_files["cdd_annot"])

        output_file = tmp_files["tmp_path"] / "output.gff"
        parse_gff(tmp_files["gff_file"], str(output_file), proteins_annotation)

        # Read output
        content = output_file.read_text()
        lines = content.strip().split("\n")

        # Check header
        assert lines[0] == "##gff-version 3"

        # Check annotated lines contain pathofact + CDD annotations
        annotated_lines = [line for line in lines if "WP_000242755.1" in line]
        assert len(annotated_lines) == 1
        assert "vfdb=VFG042734" in annotated_lines[0]
        assert "cd00038:CAP_ED" in annotated_lines[0]
        assert "cd00092:HTH_CRP" in annotated_lines[0]

    def test_parse_gff_only_annotated_proteins(self, tmp_files):
        """Test that only proteins with annotations are in output."""
        pathofact_attrs = parse_pathofact_support(tmp_files["pred_support"])
        proteins_annotation = parse_cdd(pathofact_attrs, tmp_files["cdd_annot"])

        output_file = tmp_files["tmp_path"] / "output.gff"
        parse_gff(tmp_files["gff_file"], str(output_file), proteins_annotation)

        content = output_file.read_text()

        # Proteins with annotations should be present
        assert "WP_000242755.1" in content
        assert "WP_000031783.1" in content

        # Protein without PathOFact prediction should NOT be in output
        assert "WP_000259031.1" not in content

    def test_parse_gff_preserves_coordinates(self, tmp_files):
        """Test that GFF coordinates are preserved correctly."""
        pathofact_attrs = parse_pathofact_support(tmp_files["pred_support"])
        proteins_annotation = parse_cdd(pathofact_attrs, tmp_files["cdd_annot"])

        output_file = tmp_files["tmp_path"] / "output.gff"
        parse_gff(tmp_files["gff_file"], str(output_file), proteins_annotation)

        content = output_file.read_text()
        lines = [line for line in content.split("\n") if "WP_000242755.1" in line]

        assert len(lines) == 1
        fields = lines[0].split("\t")

        # Check coordinates match input
        assert fields[0] == "Contig001"
        assert fields[3] == "87000"
        assert fields[4] == "87630"
        assert fields[6] == "+"

    def test_parse_gff_attribute_format(self, tmp_files):
        """Test that attributes are formatted correctly."""
        pathofact_attrs = {
            "WP_000242755.1": "vfdb=VFG042734;blastp_eval=2.58e-94;cdd=cd00038:CAP_ED"
        }

        output_file = tmp_files["tmp_path"] / "output.gff"
        parse_gff(tmp_files["gff_file"], str(output_file), pathofact_attrs)

        content = output_file.read_text()
        lines = [line for line in content.split("\n") if "WP_000242755.1" in line]

        assert len(lines) == 1
        attributes = lines[0].split("\t")[8]

        # Check attribute starts with ID
        assert attributes.startswith("ID=WP_000242755.1;")

        # Check semicolon-separated format
        assert "vfdb=VFG042734" in attributes
        assert "blastp_eval=2.58e-94" in attributes
        assert "cdd=cd00038:CAP_ED" in attributes

    def test_parse_gff_handles_comments(self, tmp_files):
        """Test that GFF parser handles comment lines correctly."""
        pathofact_attrs = parse_pathofact_support(tmp_files["pred_support"])
        proteins_annotation = parse_cdd(pathofact_attrs, tmp_files["cdd_annot"])

        output_file = tmp_files["tmp_path"] / "output.gff"
        parse_gff(tmp_files["gff_file"], str(output_file), proteins_annotation)

        content = output_file.read_text()
        lines = content.strip().split("\n")

        # Only first line should be the header we write
        assert lines[0] == "##gff-version 3"

        # All other lines should be data lines (9 columns)
        for line in lines[1:]:
            assert len(line.split("\t")) == 9


class TestIntegration:
    """Integration tests for complete workflow."""

    def test_full_workflow_cdd(self, tmp_files):
        """Test complete workflow with CDD annotations."""
        # Step 1: Validate inputs
        to_validate = {
            "gff": tmp_files["gff_file"],
            "annotation": tmp_files["cdd_annot"],
            "pred_support": tmp_files["pred_support"],
        }
        validation_result = validate_inputs(to_validate)
        assert len(validation_result) == 0

        # Step 2: Parse PathOFact support
        pathofact_attrs = parse_pathofact_support(tmp_files["pred_support"])
        assert len(pathofact_attrs) == 4

        # Step 3: Parse CDD annotations
        proteins_annotation = parse_cdd(pathofact_attrs, tmp_files["cdd_annot"])
        assert len(proteins_annotation) == 2  # Only 2 proteins have CDD hits

        # Step 4: Parse GFF and write output
        output_file = tmp_files["tmp_path"] / "final_output.gff"
        parse_gff(tmp_files["gff_file"], str(output_file), proteins_annotation)

        # Verify output
        content = output_file.read_text()
        assert "##gff-version 3" in content
        assert "WP_000242755.1" in content
        assert "WP_000031783.1" in content
        assert "vfdb=" in content
        assert "cdd=" in content

    def test_full_workflow_ips(self, tmp_files):
        """Test complete workflow with IPS annotations."""
        # Step 1: Validate inputs
        to_validate = {
            "gff": tmp_files["gff_file"],
            "annotation": tmp_files["ips_annot"],
            "pred_support": tmp_files["pred_support"],
        }
        validation_result = validate_inputs(to_validate)
        assert len(validation_result) == 0

        # Step 2: Parse PathOFact support
        pathofact_attrs = parse_pathofact_support(tmp_files["pred_support"])

        # Step 3: Parse IPS annotations
        proteins_annotation = parse_ips(pathofact_attrs, tmp_files["ips_annot"])
        assert len(proteins_annotation) == 2

        # Step 4: Parse GFF and write output
        output_file = tmp_files["tmp_path"] / "final_output_ips.gff"
        parse_gff(tmp_files["gff_file"], str(output_file), proteins_annotation)

        # Verify output
        content = output_file.read_text()
        assert "##gff-version 3" in content
        assert "WP_000242755.1" in content
        assert "WP_000031783.1" in content

    def test_workflow_with_invalid_pred_support(self, tmp_files):
        """Test workflow fails gracefully with invalid pred_support."""
        to_validate = {
            "gff": tmp_files["gff_file"],
            "annotation": tmp_files["cdd_annot"],
            "pred_support": tmp_files["header_only"],
        }
        validation_result = validate_inputs(to_validate)

        # Should fail at validation
        assert len(validation_result) > 0
        assert "pred_support" in validation_result


class TestFullWorkflowIPS:
    """Test complete workflow with IPS annotations."""

    def test_workflow_ips_complete(self, tmp_files):
        """Test full workflow using IPS annotation file."""
        # Step 1: Parse PathOFact support
        pathofact_attrs = parse_pathofact_support(tmp_files["pred_support"])
        assert len(pathofact_attrs) > 0

        # Step 2: Parse IPS annotations
        proteins_annotation = parse_ips(pathofact_attrs, tmp_files["ips_annot"])
        assert len(proteins_annotation) == 2

        # Step 3: Parse GFF and write output
        output_file = tmp_files["tmp_path"] / "final_output_ips.gff"
        parse_gff(tmp_files["gff_file"], str(output_file), proteins_annotation)

        # Verify output
        content = output_file.read_text()
        assert "##gff-version 3" in content
        assert "WP_000242755.1" in content
        assert "WP_000031783.1" in content


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
