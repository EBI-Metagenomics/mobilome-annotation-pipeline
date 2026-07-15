import gzip
from unittest.mock import patch

import pytest

from pathofact2_report import (
    assign_mge_types,
    build_rows,
    calculate_overlap_length,
    join_or_dash,
    merge_unique_preserving_order,
    parse_amr_gff,
    parse_bgc_gff,
    parse_interproscan_signalp,
    parse_mobilome_gff,
    parse_pathofact2_gff,
    path_is_missing_or_empty,
    remap_contigs,
    split_csv_value,
)

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

MODULE = "pathofact2_report._iter_gff_rows"


def gff_row(contig, feature_type, start, end, attrs):
    raw = f"{contig}\t.\t{feature_type}\t{start}\t{end}\t.\t+\t0\t{attrs}"
    return (raw, contig, feature_type, start, end, _parse(attrs))


def _parse(attr_str):
    result = {}
    for part in attr_str.split(";"):
        if "=" in part:
            k, v = part.split("=", 1)
            result[k] = v
    return result


# ---------------------------------------------------------------------------
# path_is_missing_or_empty
# ---------------------------------------------------------------------------


class TestPathIsMissingOrEmpty:
    def test_none(self):
        assert path_is_missing_or_empty(None) is True

    def test_nonexistent_path(self, tmp_path):
        assert path_is_missing_or_empty(tmp_path / "ghost.gff") is True

    def test_empty_file(self, tmp_path):
        f = tmp_path / "empty.gff"
        f.write_text("")
        assert path_is_missing_or_empty(f) is True

    def test_populated_file(self, tmp_path):
        f = tmp_path / "data.gff"
        f.write_text("##gff-version 3\n")
        assert path_is_missing_or_empty(f) is False


# ---------------------------------------------------------------------------
# split_csv_value
# ---------------------------------------------------------------------------


class TestSplitCsvValue:
    def test_empty_string(self):
        assert split_csv_value("") == []

    def test_dash(self):
        assert split_csv_value("-") == []

    def test_single_value(self):
        assert split_csv_value("RiPP") == ["RiPP"]

    def test_multiple_values(self):
        assert split_csv_value("sanntis,gecco,antismash") == ["sanntis", "gecco", "antismash"]

    def test_strips_spaces(self):
        assert split_csv_value("sanntis, gecco , antismash") == ["sanntis", "gecco", "antismash"]


# ---------------------------------------------------------------------------
# merge_unique_preserving_order
# ---------------------------------------------------------------------------


class TestMergeUniquePreservingOrder:
    def test_no_overlap(self):
        assert merge_unique_preserving_order(["a"], ["b", "c"]) == ["a", "b", "c"]

    def test_duplicates_dropped(self):
        assert merge_unique_preserving_order(["a", "b"], ["b", "c"]) == ["a", "b", "c"]

    def test_empty_new_values(self):
        assert merge_unique_preserving_order(["a"], []) == ["a"]

    def test_empty_existing(self):
        assert merge_unique_preserving_order([], ["x", "y"]) == ["x", "y"]

    def test_order_preserved(self):
        result = merge_unique_preserving_order(["z", "a"], ["b", "a"])
        assert result == ["z", "a", "b"]


# ---------------------------------------------------------------------------
# join_or_dash
# ---------------------------------------------------------------------------


class TestJoinOrDash:
    def test_non_empty(self):
        assert join_or_dash(["sanntis", "gecco"]) == "sanntis,gecco"

    def test_empty_list(self):
        assert join_or_dash([]) == "-"

    def test_single_item(self):
        assert join_or_dash(["sanntis"]) == "sanntis"


# ---------------------------------------------------------------------------
# calculate_overlap_length
# ---------------------------------------------------------------------------


class TestCalculateOverlapLength:
    def test_no_overlap_before(self):
        assert calculate_overlap_length(100, 200, 300, 400) == 0

    def test_no_overlap_after(self):
        assert calculate_overlap_length(300, 400, 100, 200) == 0

    def test_partial_overlap(self):
        assert calculate_overlap_length(100, 200, 150, 250) == 51

    def test_full_containment(self):
        assert calculate_overlap_length(100, 300, 150, 200) == 51

    def test_exact_match(self):
        assert calculate_overlap_length(100, 200, 100, 200) == 101

    def test_single_base_touch(self):
        assert calculate_overlap_length(100, 200, 200, 300) == 1


# ---------------------------------------------------------------------------
# parse_pathofact2_gff
# ---------------------------------------------------------------------------


class TestParsePathofact2GFF:
    def test_typical_cds(self):
        rows = [
            gff_row(
                "ctg1", "CDS", 100, 300,
                "ID=prot1;vfdb=VFG001;blastp_eval=1e-10;pathofact2_tox_prob=0.9;pathofact2_vf_prob=0.8;cdd=cd001",
            )
        ]
        with patch(MODULE, return_value=iter(rows)):
            coords, data = parse_pathofact2_gff("fake.gff")
        assert coords == {"prot1": ("ctg1", 100, 300)}
        assert data["prot1"] == ("VFG001", "1e-10", "0.9", "0.8", "cd001")

    def test_non_cds_rows_skipped(self):
        rows = [
            gff_row("ctg1", "gene", 100, 300, "ID=gene1"),
            gff_row("ctg1", "CDS", 100, 300, "ID=prot1;vfdb=VFG001;blastp_eval=1e-5;pathofact2_tox_prob=0.5;pathofact2_vf_prob=0.5;cdd=cd001"),
        ]
        with patch(MODULE, return_value=iter(rows)):
            coords, _ = parse_pathofact2_gff("fake.gff")
        assert "gene1" not in coords
        assert "prot1" in coords

    def test_missing_id_skipped(self):
        rows = [gff_row("ctg1", "CDS", 100, 300, "vfdb=VFG001")]
        with patch(MODULE, return_value=iter(rows)):
            coords, data = parse_pathofact2_gff("fake.gff")
        assert len(coords) == 0

    def test_missing_attributes_default_to_dash(self):
        rows = [gff_row("ctg1", "CDS", 100, 300, "ID=prot1")]
        with patch(MODULE, return_value=iter(rows)):
            _, data = parse_pathofact2_gff("fake.gff")
        assert data["prot1"] == ("-", "-", "-", "-", "-")

    def test_empty_gff(self):
        with patch(MODULE, return_value=iter([])):
            coords, data = parse_pathofact2_gff("fake.gff")
        assert coords == {}
        assert data == {}


# ---------------------------------------------------------------------------
# parse_amr_gff
# ---------------------------------------------------------------------------


class TestParseAmrGFF:
    def test_typical_cds(self):
        rows = [gff_row("ctg1", "CDS", 100, 300, "ID=prot1;drug_class=beta-lactam;amr_tool=AMRFinder;amr_tool_ident=99.5")]
        prots_coords = {}
        with patch(MODULE, return_value=iter(rows)):
            data = parse_amr_gff("fake.gff", prots_coords)
        assert data["prot1"] == ("beta-lactam", "AMRFinder", "99.5")

    def test_extends_prots_coords_for_new_protein(self):
        rows = [gff_row("ctg1", "CDS", 50, 150, "ID=prot_new;drug_class=X;amr_tool=T;amr_tool_ident=90")]
        prots_coords = {}
        with patch(MODULE, return_value=iter(rows)):
            parse_amr_gff("fake.gff", prots_coords)
        assert "prot_new" in prots_coords
        assert prots_coords["prot_new"] == ("ctg1", 50, 150)

    def test_does_not_overwrite_existing_coords(self):
        rows = [gff_row("ctg2", "CDS", 999, 1999, "ID=prot1;drug_class=X;amr_tool=T;amr_tool_ident=90")]
        prots_coords = {"prot1": ("ctg1", 100, 300)}
        with patch(MODULE, return_value=iter(rows)):
            parse_amr_gff("fake.gff", prots_coords)
        assert prots_coords["prot1"] == ("ctg1", 100, 300)

    def test_missing_id_skipped(self):
        rows = [gff_row("ctg1", "CDS", 100, 300, "drug_class=X;amr_tool=T;amr_tool_ident=90")]
        prots_coords = {}
        with patch(MODULE, return_value=iter(rows)):
            data = parse_amr_gff("fake.gff", prots_coords)
        assert len(data) == 0

    def test_non_cds_skipped(self):
        rows = [gff_row("ctg1", "mRNA", 100, 300, "ID=prot1;drug_class=X;amr_tool=T;amr_tool_ident=90")]
        prots_coords = {}
        with patch(MODULE, return_value=iter(rows)):
            data = parse_amr_gff("fake.gff", prots_coords)
        assert len(data) == 0

    def test_missing_amr_attributes_default_to_dash(self):
        rows = [gff_row("ctg1", "CDS", 100, 300, "ID=prot1")]
        prots_coords = {}
        with patch(MODULE, return_value=iter(rows)):
            data = parse_amr_gff("fake.gff", prots_coords)
        assert data["prot1"] == ("-", "-", "-")


# ---------------------------------------------------------------------------
# parse_mobilome_gff
# ---------------------------------------------------------------------------


class TestParseMobilomeGFF:
    def test_typical_mge(self):
        rows = [gff_row("ctg1", "insertion_sequence", 100, 2000, "mobile_element_type=IS3")]
        with patch(MODULE, return_value=iter(rows)):
            data = parse_mobilome_gff("fake.gff")
        assert data == {"ctg1": [(100, 2000, "IS3")]}

    def test_ignored_feature_types_skipped(self):
        ignored = ["inverted_repeat_element", "attC_site", "direct_repeat", "terminal_inverted_repeat_element"]
        for ft in ignored:
            rows = [gff_row("ctg1", ft, 100, 200, "ID=x")]
            with patch(MODULE, return_value=iter(rows)):
                data = parse_mobilome_gff("fake.gff")
            assert data == {}, f"{ft} should be ignored"

    def test_falls_back_to_feature_type_when_no_attribute(self):
        rows = [gff_row("ctg1", "integron", 100, 500, "ID=x")]
        with patch(MODULE, return_value=iter(rows)):
            data = parse_mobilome_gff("fake.gff")
        assert data["ctg1"][0][2] == "integron"

    def test_multiple_mges_same_contig(self):
        rows = [
            gff_row("ctg1", "insertion_sequence", 100, 500, "mobile_element_type=IS1"),
            gff_row("ctg1", "insertion_sequence", 600, 1000, "mobile_element_type=IS3"),
        ]
        with patch(MODULE, return_value=iter(rows)):
            data = parse_mobilome_gff("fake.gff")
        assert len(data["ctg1"]) == 2

    def test_empty_gff(self):
        with patch(MODULE, return_value=iter([])):
            data = parse_mobilome_gff("fake.gff")
        assert data == {}


# ---------------------------------------------------------------------------
# parse_bgc_gff
# ---------------------------------------------------------------------------


class TestParseBGCGFF:
    def test_typical_bgc(self):
        rows = [gff_row("ctg1", "CDS", 100, 300, "ID=prot1;bgc_tools=sanntis;nearest_MiBIG_class=RiPP")]
        prots_coords = {"prot1": ("ctg1", 100, 300)}
        with patch(MODULE, return_value=iter(rows)):
            data = parse_bgc_gff("fake.gff", prots_coords)
        assert data["prot1"] == ("sanntis", "RiPP")

    def test_protein_not_in_seed_skipped(self):
        rows = [gff_row("ctg1", "CDS", 100, 300, "ID=prot_other;bgc_tools=gecco;nearest_MiBIG_class=Terpene")]
        prots_coords = {"prot1": ("ctg1", 100, 300)}
        with patch(MODULE, return_value=iter(rows)):
            data = parse_bgc_gff("fake.gff", prots_coords)
        assert len(data) == 0

    def test_missing_id_skipped(self):
        rows = [gff_row("ctg1", "CDS", 100, 300, "bgc_tools=gecco")]
        prots_coords = {}
        with patch(MODULE, return_value=iter(rows)):
            data = parse_bgc_gff("fake.gff", prots_coords)
        assert len(data) == 0

    def test_all_type_sources_merged(self):
        rows = [gff_row(
            "ctg1", "CDS", 100, 300,
            "ID=prot1;bgc_tools=sanntis,gecco;nearest_MiBIG_class=RiPP;antismash_product=lanthipeptide;gecco_bgc_type=Polyketide",
        )]
        prots_coords = {"prot1": ("ctg1", 100, 300)}
        with patch(MODULE, return_value=iter(rows)):
            data = parse_bgc_gff("fake.gff", prots_coords)
        tools, types = data["prot1"]
        assert "sanntis" in tools and "gecco" in tools
        assert "RiPP" in types and "lanthipeptide" in types and "Polyketide" in types

    def test_duplicate_types_deduplicated(self):
        rows = [
            gff_row("ctg1", "CDS", 100, 300, "ID=prot1;bgc_tools=sanntis;nearest_MiBIG_class=RiPP;antismash_product=RiPP"),
        ]
        prots_coords = {"prot1": ("ctg1", 100, 300)}
        with patch(MODULE, return_value=iter(rows)):
            data = parse_bgc_gff("fake.gff", prots_coords)
        _, types = data["prot1"]
        assert types.count("RiPP") == 1

    def test_missing_bgc_attributes_default_to_dash(self):
        rows = [gff_row("ctg1", "CDS", 100, 300, "ID=prot1")]
        prots_coords = {"prot1": ("ctg1", 100, 300)}
        with patch(MODULE, return_value=iter(rows)):
            data = parse_bgc_gff("fake.gff", prots_coords)
        assert data["prot1"] == ("-", "-")


# ---------------------------------------------------------------------------
# parse_interproscan_signalp  (real temp files — no toolkit needed)
# ---------------------------------------------------------------------------

IPS_COLS = ["protein_id", "md5", "length", "analysis", "accession", "description",
            "start", "end", "score", "status", "date", "ipr_acc", "ipr_desc", "-", "-"]


def _ips_line(*vals):
    return "\t".join(vals)


class TestParseInterproscanSignalp:
    def _write_ips(self, path, rows):
        path.write_text("\n".join(rows) + "\n")

    def test_typical_signalp_row(self, tmp_path):
        f = tmp_path / "ips.tsv"
        self._write_ips(f, [_ips_line("prot1", "md5", "200", "SignalP_EUK", "Sec/SPI", "", "", "", "", "", "", "", "", "", "")])
        data = parse_interproscan_signalp(f)
        assert data["prot1"] == "Sec/SPI"

    def test_non_signalp_rows_skipped(self, tmp_path):
        f = tmp_path / "ips.tsv"
        self._write_ips(f, [
            _ips_line("prot1", "md5", "200", "Pfam", "PF00001", "domain", "", "", "", "", "", "", "", "", ""),
            _ips_line("prot2", "md5", "200", "SignalP_EUK", "Sec/SPI", "", "", "", "", "", "", "", "", "", ""),
        ])
        data = parse_interproscan_signalp(f)
        assert "prot1" not in data
        assert "prot2" in data

    def test_empty_annotation_skipped(self, tmp_path):
        f = tmp_path / "ips.tsv"
        self._write_ips(f, [
            _ips_line("prot1", "md5", "200", "SignalP_EUK", "", "", "", "", "", "", "", "", "", "", ""),
            _ips_line("prot2", "md5", "200", "SignalP_EUK", "-", "", "", "", "", "", "", "", "", "", ""),
        ])
        data = parse_interproscan_signalp(f)
        assert "prot1" not in data
        assert "prot2" not in data

    def test_multiple_signalp_hits_merged(self, tmp_path):
        f = tmp_path / "ips.tsv"
        self._write_ips(f, [
            _ips_line("prot1", "md5", "200", "SignalP_EUK", "Sec/SPI", "", "", "", "", "", "", "", "", "", ""),
            _ips_line("prot1", "md5", "200", "SignalP_GRAM_NEGATIVE", "Sec/SPII", "", "", "", "", "", "", "", "", "", ""),
        ])
        data = parse_interproscan_signalp(f)
        assert "Sec/SPI" in data["prot1"] and "Sec/SPII" in data["prot1"]

    def test_gzip_input(self, tmp_path):
        f = tmp_path / "ips.tsv.gz"
        content = _ips_line("prot1", "md5", "200", "SignalP_EUK", "Sec/SPI", "", "", "", "", "", "", "", "", "", "") + "\n"
        with gzip.open(f, "wt") as fh:
            fh.write(content)
        data = parse_interproscan_signalp(f)
        assert data["prot1"] == "Sec/SPI"

    def test_empty_file(self, tmp_path):
        f = tmp_path / "empty.tsv"
        f.write_text("")
        data = parse_interproscan_signalp(f)
        assert data == {}


# ---------------------------------------------------------------------------
# assign_mge_types
# ---------------------------------------------------------------------------


class TestAssignMgeTypes:
    def test_full_overlap_assigned(self):
        prots = {"prot1": ("ctg1", 100, 200)}
        mges = {"ctg1": [(50, 250, "IS3")]}
        result = assign_mge_types(prots, mges)
        assert result["prot1"] == "IS3"

    def test_no_overlap_dash(self):
        prots = {"prot1": ("ctg1", 100, 200)}
        mges = {"ctg1": [(300, 500, "IS3")]}
        result = assign_mge_types(prots, mges)
        assert result["prot1"] == "-"

    def test_overlap_below_threshold_not_assigned(self):
        # protein length 101, overlap = 10 → fraction ~0.10 < 0.90
        prots = {"prot1": ("ctg1", 100, 200)}
        mges = {"ctg1": [(191, 300, "IS3")]}
        result = assign_mge_types(prots, mges)
        assert result["prot1"] == "-"

    def test_overlap_at_threshold_assigned(self):
        # protein 100-200 (length 101), mge 109-300 → overlap 100-200 = 92 → 92/101 = 0.91
        prots = {"prot1": ("ctg1", 100, 200)}
        mges = {"ctg1": [(109, 300, "IS3")]}
        result = assign_mge_types(prots, mges)
        assert result["prot1"] == "IS3"

    def test_multiple_mge_types_merged(self):
        prots = {"prot1": ("ctg1", 100, 200)}
        mges = {"ctg1": [(50, 250, "IS3"), (50, 250, "integron")]}
        result = assign_mge_types(prots, mges)
        assert "IS3" in result["prot1"] and "integron" in result["prot1"]

    def test_protein_on_different_contig(self):
        prots = {"prot1": ("ctg2", 100, 200)}
        mges = {"ctg1": [(50, 250, "IS3")]}
        result = assign_mge_types(prots, mges)
        assert result["prot1"] == "-"


# ---------------------------------------------------------------------------
# remap_contigs
# ---------------------------------------------------------------------------


class TestRemapContigs:
    def test_renamed_contigs_translated_to_original(self):
        prots = {"prot1": ("1", 100, 200), "prot2": ("2", 300, 400)}
        names_equiv = {"1": "NZ_real_1", "2": "NZ_real_2"}
        result = remap_contigs(prots, names_equiv)
        assert result == {
            "prot1": ("NZ_real_1", 100, 200),
            "prot2": ("NZ_real_2", 300, 400),
        }

    def test_unmapped_contigs_pass_through_unchanged(self):
        # Original names (e.g. user-provided proteins) are not map keys -> left untouched
        prots = {"prot1": ("NZ_real_1", 100, 200)}
        names_equiv = {"1": "NZ_real_1"}
        result = remap_contigs(prots, names_equiv)
        assert result == {"prot1": ("NZ_real_1", 100, 200)}

    def test_empty_map_is_noop(self):
        prots = {"prot1": ("1", 100, 200)}
        assert remap_contigs(prots, {}) == prots


# ---------------------------------------------------------------------------
# build_rows
# ---------------------------------------------------------------------------


class TestBuildRows:
    def _base(self):
        return {
            "prot1": ("ctg1", 100, 200),
            "prot2": ("ctg1", 300, 400),
        }

    def test_all_annotations_present(self):
        rows = build_rows(
            prots_coords={"prot1": ("ctg1", 100, 200)},
            pathofact_data={"prot1": ("VFG001", "1e-5", "0.9", "0.8", "cd001")},
            amr_data={"prot1": ("beta-lactam", "AMRFinder", "99.0")},
            protein_mge_map={"prot1": "IS3"},
            bgc_data={"prot1": ("sanntis", "RiPP")},
            signalp_data={"prot1": "Sec/SPI"},
        )
        assert len(rows) == 1
        r = rows[0]
        assert r["protein_id"] == "prot1"
        assert r["vfdb_hit"] == "VFG001"
        assert r["amr_drug_class"] == "beta-lactam"
        assert r["mge_type"] == "IS3"
        assert r["bgc_tools"] == "sanntis"
        assert r["signalP"] == "Sec/SPI"

    def test_missing_optional_annotations_default_to_dash(self):
        rows = build_rows(
            prots_coords={"prot1": ("ctg1", 100, 200)},
            pathofact_data={},
            amr_data={},
            protein_mge_map={},
            bgc_data={},
            signalp_data={},
        )
        r = rows[0]
        assert r["vfdb_hit"] == "-"
        assert r["amr_drug_class"] == "-"
        assert r["mge_type"] == "-"
        assert r["bgc_type"] == "-"
        assert r["signalP"] == "-"

    def test_rows_sorted_by_protein_id(self):
        rows = build_rows(
            prots_coords={"prot_b": ("ctg1", 300, 400), "prot_a": ("ctg1", 100, 200)},
            pathofact_data={}, amr_data={}, protein_mge_map={}, bgc_data={}, signalp_data={},
        )
        assert rows[0]["protein_id"] == "prot_a"
        assert rows[1]["protein_id"] == "prot_b"
