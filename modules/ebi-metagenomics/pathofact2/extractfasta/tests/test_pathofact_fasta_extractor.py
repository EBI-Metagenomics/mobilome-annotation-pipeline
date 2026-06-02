#!/usr/bin/env python

import sys
from pathlib import Path

import pytest

pytest.importorskip("Bio")

# Add the script directory to the path
script_dir = Path(__file__).parent.parent / "resources" / "usr" / "bin"
sys.path.insert(0, str(script_dir))

# Import functions from the script
from pathofact_fasta_extractor import (  # noqa: E402
    BestBlastHit,
    _file_is_readable,
    _fmt_prob,
    collect_detected_sequence_ids,
    parse_blastp_best_hits,
    parse_pathofact2_predictions_tsv,
    read_fasta_to_dict,
    write_detected_fasta,
    write_support_table,
)


@pytest.fixture
def tmp_files(tmp_path):
    """Create temporary test files with realistic data."""

    # FASTA file data
    fasta_data = {
        "WP_001273183.1": "MANFFIQRPIFAWVLAIILMMAGALAILQLPVAQYPTIAPPAVSVSANYPGADAQTVQDTVTQVIEQNMNGIDNLMYMSSTSDSAGSVTITLTFQSGTDPDIAQVQVQNKLQLATPLLPQEVQQQGISVEKSSSSYLMVAGFVSDNPGTTQDDISDYVASNVKDTLSRLNGVGDVQLFGAQYAMRIWLDADLLNKYKLTPVDVINQLKVQNDQIAAGQLGGTPALPGQQLNASIIAQTRLKNPEEFGKVTLRVNSDGSVVRLKDVARVELGGENYNVIARINGKPAAGLGIKLATGANALDTAKAIKAKLAELQPFFPQGMKVLYPYDTTPFVQLSIQEVVKTLFEAIILVFLVMYLFLQNMRATLIPTIAVPVVLLGTFAILAAFGYSINTLTMFGMVLAIGLLVDDAIVVVENVERVMMEDKLPPREATEKSMSQIQGALVGIAMVLSAVFIPMAFFGGSTGAIYRQFSITIVSAMALSVLVALILTPALCATLLKPVSAEHHENKGGFFGWFNTTFDHSVNHYTNSVGKILGSTGRYLLIYALIVAGMVVLFLRLPSSFLPEEDQGVFLTMIQLPAGATQERTQKVLDQVTDYYLKNEKANVESVFTVNGFSFSGQAQNAGMAFVSLKPWEERSGDENSAEAVIHRAKMELGKIRDGFVIPFNMPAIVELGTATGFDFELIDQAGLGHDALTQARNQLLGMAAQHPASLVSVRPNGLEDTAQFKLEVDQEKAQALGVSLSDINQTISTALGGTYVNDFIDRGRVKKVYVQADAKFRMLPEDVDKLYVRSANGEMVPFSAFTTSHWVYGSPRLERYNGLPSMEIQGEAAPGTSSGDAMALMENLASKLPAGIGYDWTGMSYQERLSGNQAPALVAISFVVVFLCLAALYESWSIPVSVMLVVPLGIVGVLLAATLFNQKNDVYFMVGLLTTIGLSAKNAILIVEFAKDLMEKEGKGVVEATLMAVRMRLRPILMTSLAFILGVLPLAISNGAGSGAQNAVGIGVMGGMVSATLLAIFFVPVFFVVIRRCFKG",
        "WP_000242755.1": "MVLGKPQTDPTLEWFLSHCHIHKYPSKSTLIHQGEKAETLYYIVKGSVAVLIKDEEGKEMILSYLNQGDFIGELGLFEEGQERSAWVRAKTACEVAEISYKKFRQLIQVNPDILMRLSAQMARRLQVTSEKVGNLAFLDVTGRIAQTLLNLAKQPDAMTHPDGMQIKITRQEIGQIVGCSRETVGRILKMLEDQNLISAHGKTIVVYGTR",
        "WP_000163771.1": "MTKLKLLALGVLIATSAGVAHAEGKFSLGAGVGVVEHPYKDYDTDVYPVPVINYEGDNFWFRGLGGGYYLWNDATDKLSITAYWSPLYFKAKDSGDHQMRHLDDRKSTMMAGLSYAHFTQYGYLRTTLAGDTLDNSNGIVWDMAWLYRYTNGGLTVTPGIGVQWNSENQNEYYYGVSRKESARSGLRGYNPNDSWSPYLELSASYNFLGDWSVYGTARYTRLSDEVTDSPMVDKSWTGLISTGITYKF",
        "WP_000031783.1": "MSKEKFERTKPHVNVGTIGHVDHGKTTLTAAITTVLAKTYGGAARAFDQIDNAPEEKARGITINTSHVEYDTPTRHYAHVDCPGHADYVKNMITGAAQMDGAILVVAATDGPMPQTREHILLGRQVGVPYIIVFLNKCDMVDDEELLELVEMEVRELLSQYDFPGDDTPIVRGSALKALEGDAEWEAKILELAGFLDSYIPEPERAIDKPFLLPIEDVFSISGRGTVVTGRVERGIIKVGEEVEIVGIKETQKSTCTGVEMFRKLLDEGRAGENVGVLLRGIKREEIERGQVLAKPGTIKPHTKFESEVYILSKDEGGRHTPFFKGYRPQFYFRTTDVTGTIELPEGVEMVMPGDNIKMVVTLIHPIAMDDGLRFAIREGGRTVGAGVVAKVLG",
        "WP_000259031.1": "MVTVFGILNLTEDSFFDESRRLDPAGAVTAAIEMLRVGSDVVDVGPAASHPDARPVSPADEIRRIAPLLDALSDQMHRVSIDSFQPETQRYALKRGVGYLNDIQGFPDPALYPDIAEADCRLVVMHSAQRDGIATRTGHLRPEDALDEIVRFFEARVSALRRSGVAADRLILDPGMGFFLSPAPETSLHVLSNLQKLKSALGLPLLVSVSRKSFLGATVGLPVKDLGPASLAAELHAIGNGADYVRTHAPGDLRSAITFSETLAKFRSRDARDRGLDHA",
    }

    # Toxins TSV data
    toxins_data = {
        "header": ["Sequence", "Prediction", "Probability"],
        "rows": [
            ["WP_000242755.1", "1", "1.00"],
            ["WP_001273183.1", "1", "1.00"],
        ],
    }

    # Virulence TSV data
    virulence_data = {
        "header": ["Sequence", "Prediction", "Probability"],
        "rows": [
            ["WP_000163771.1", "1", "1.0"],
            ["WP_001273183.1", "1", "0.9977883"],
        ],
    }

    # DIAMOND BLAST output data (tab-separated, 8 columns)
    diamond_data = {
        "rows": [
            [
                "WP_000242755.1",
                "VFG042734(gb|NP_249343)",
                "66.8",
                "202",
                "210",
                "214",
                "2.58e-94",
                "274",
            ],
            [
                "WP_000031783.1",
                "VFG046459(gb|WP_004287053)",
                "81.7",
                "393",
                "394",
                "394",
                "6.50e-248",
                "677",
            ],
            [
                "WP_000031783.1",
                "VFG046474(gb|WP_014714676)",
                "81.4",
                "393",
                "394",
                "394",
                "9.23e-248",
                "677",
            ],
            [
                "WP_001273183.1",
                "VFG049144(gb|WP_002892069)",
                "78.0",
                "1033",
                "1034",
                "1048",
                "0.0",
                "1541",
            ],
        ]
    }

    # Create FASTA file
    fasta_file = tmp_path / "proteins.faa"
    fasta_content = ""
    for seq_id, sequence in fasta_data.items():
        fasta_content += f">{seq_id}\n{sequence}\n"
    fasta_file.write_text(fasta_content)

    # Create toxins TSV
    toxins_file = tmp_path / "classifier_toxins.tsv"
    toxins_content = "\t".join(toxins_data["header"]) + "\n"
    toxins_content += "\n".join(["\t".join(row) for row in toxins_data["rows"]])
    toxins_file.write_text(toxins_content + "\n")

    # Create virulence TSV
    virulence_file = tmp_path / "classifier_virulence.tsv"
    virulence_content = "\t".join(virulence_data["header"]) + "\n"
    virulence_content += "\n".join(["\t".join(row) for row in virulence_data["rows"]])
    virulence_file.write_text(virulence_content + "\n")

    # Create DIAMOND BLAST output
    diamond_file = tmp_path / "diamond_vfdb.txt"
    diamond_content = "\n".join(["\t".join(row) for row in diamond_data["rows"]])
    diamond_file.write_text(diamond_content + "\n")

    # Create empty files for testing
    empty_file = tmp_path / "empty.txt"
    empty_file.write_text("")

    # Create invalid TSV files for error testing
    invalid_header_tsv = tmp_path / "invalid_header.tsv"
    invalid_header_tsv.write_text("WrongHeader\tBadColumn\tNoGood\ndata\tmore\tdata\n")

    malformed_row_tsv = tmp_path / "malformed_row.tsv"
    malformed_row_tsv.write_text(
        "Sequence\tPrediction\tProbability\nWP_12345\t1\n"
    )  # Missing column

    invalid_probability_tsv = tmp_path / "invalid_prob.tsv"
    invalid_probability_tsv.write_text(
        "Sequence\tPrediction\tProbability\nWP_12345\t1\tNOT_A_NUMBER\n"
    )

    return {
        "fasta_file": str(fasta_file),
        "toxins_file": str(toxins_file),
        "virulence_file": str(virulence_file),
        "diamond_file": str(diamond_file),
        "empty_file": str(empty_file),
        "invalid_header_tsv": str(invalid_header_tsv),
        "malformed_row_tsv": str(malformed_row_tsv),
        "invalid_probability_tsv": str(invalid_probability_tsv),
        "tmp_path": tmp_path,
        "fasta_data": fasta_data,
        "toxins_data": toxins_data,
        "virulence_data": virulence_data,
        "diamond_data": diamond_data,
    }


class TestHelperFunctions:
    """Test utility helper functions."""

    def test_file_is_readable_valid(self, tmp_files):
        """Test readable file returns True."""
        assert _file_is_readable(tmp_files["fasta_file"]) is True

    def test_file_is_readable_empty(self, tmp_files):
        """Test empty file returns False."""
        assert _file_is_readable(tmp_files["empty_file"]) is False

    def test_file_is_readable_nonexistent(self):
        """Test nonexistent file returns False."""
        assert _file_is_readable("/nonexistent/file.txt") is False

    def test_file_is_readable_none(self):
        """Test None input returns False."""
        assert _file_is_readable(None) is False

    def test_fmt_prob_default(self):
        """Test probability formatting with default precision."""
        assert _fmt_prob(0.9977883) == "0.997788"
        assert _fmt_prob(1.0) == "1.000000"
        assert _fmt_prob(0.6) == "0.600000"

    def test_fmt_prob_custom_precision(self):
        """Test probability formatting with custom precision."""
        assert _fmt_prob(0.9977883, ndigits=2) == "1.00"
        assert _fmt_prob(0.9977883, ndigits=4) == "0.9978"
        assert _fmt_prob(0.1234567, ndigits=2) == "0.12"
        assert _fmt_prob(0.1234567, ndigits=3) == "0.123"


class TestReadFasta:
    """Test FASTA reading functionality."""

    def test_read_fasta_basic(self, tmp_files):
        """Test reading valid FASTA file."""
        sequences = read_fasta_to_dict(tmp_files["fasta_file"])

        assert len(sequences) == 5
        assert "WP_001273183.1" in sequences
        assert "WP_000242755.1" in sequences
        assert sequences["WP_000242755.1"].id == "WP_000242755.1"

    def test_read_fasta_sequence_content(self, tmp_files):
        """Test that sequences are read correctly."""
        sequences = read_fasta_to_dict(tmp_files["fasta_file"])

        # Check a specific sequence
        assert str(sequences["WP_000242755.1"].seq).startswith("MVLGKPQTDPTLEWFL")
        assert len(str(sequences["WP_001273183.1"].seq)) == 1034

    def test_read_fasta_empty_file(self, tmp_files):
        """Test reading empty file returns empty dict."""
        sequences = read_fasta_to_dict(tmp_files["empty_file"])
        assert sequences == {}

    def test_read_fasta_nonexistent_file(self):
        """Test reading nonexistent file returns empty dict."""
        sequences = read_fasta_to_dict("/nonexistent/file.fasta")
        assert sequences == {}


class TestParsePathofact2Predictions:
    """Test PathOFact2 prediction TSV parsing."""

    def test_parse_toxins_basic(self, tmp_files):
        """Test parsing toxins TSV with threshold."""
        preds = parse_pathofact2_predictions_tsv(tmp_files["toxins_file"], 0.6)

        assert len(preds) == 2
        assert "WP_000242755.1" in preds
        assert "WP_001273183.1" in preds
        assert preds["WP_000242755.1"] == 1.0
        assert preds["WP_001273183.1"] == 1.0

    def test_parse_virulence_basic(self, tmp_files):
        """Test parsing virulence TSV with threshold."""
        preds = parse_pathofact2_predictions_tsv(tmp_files["virulence_file"], 0.9)

        assert len(preds) == 2
        assert "WP_000163771.1" in preds
        assert "WP_001273183.1" in preds
        assert preds["WP_000163771.1"] == 1.0
        assert preds["WP_001273183.1"] == pytest.approx(0.9977883)

    def test_parse_predictions_threshold_filtering(self, tmp_files):
        """Test that threshold filtering works correctly."""
        # Higher threshold should filter out lower probability
        preds_low = parse_pathofact2_predictions_tsv(tmp_files["virulence_file"], 0.9)
        preds_high = parse_pathofact2_predictions_tsv(
            tmp_files["virulence_file"], 0.999
        )

        assert len(preds_low) == 2  # Both pass 0.9 threshold
        assert len(preds_high) == 1  # Only 1.0 passes 0.999 threshold
        assert "WP_000163771.1" in preds_high

    def test_parse_predictions_empty_file(self, tmp_files):
        """Test parsing empty file returns empty dict."""
        preds = parse_pathofact2_predictions_tsv(tmp_files["empty_file"], 0.5)
        assert preds == {}

    def test_parse_predictions_invalid_header(self, tmp_files):
        """Test that invalid header raises ValueError."""
        with pytest.raises(ValueError) as exc_info:
            parse_pathofact2_predictions_tsv(tmp_files["invalid_header_tsv"], 0.5)
        assert "Unexpected header" in str(exc_info.value)

    def test_parse_predictions_malformed_row(self, tmp_files):
        """Test that malformed row raises ValueError."""
        with pytest.raises(ValueError) as exc_info:
            parse_pathofact2_predictions_tsv(tmp_files["malformed_row_tsv"], 0.5)
        assert "Malformed row" in str(exc_info.value)

    def test_parse_predictions_invalid_probability(self, tmp_files):
        """Test that invalid probability value raises ValueError."""
        with pytest.raises(ValueError) as exc_info:
            parse_pathofact2_predictions_tsv(tmp_files["invalid_probability_tsv"], 0.5)
        assert "Invalid Probability value" in str(exc_info.value)


class TestParseBlastpBestHits:
    """Test BLASTP output parsing."""

    def test_parse_blastp_basic(self, tmp_files):
        """Test parsing BLASTP output."""
        hits = parse_blastp_best_hits(tmp_files["diamond_file"])

        assert len(hits) == 3  # 3 unique query IDs
        assert "WP_000242755.1" in hits
        assert "WP_000031783.1" in hits
        assert "WP_001273183.1" in hits

    def test_parse_blastp_best_hit_selection(self, tmp_files):
        """Test that best hit is selected (highest bitscore)."""
        hits = parse_blastp_best_hits(tmp_files["diamond_file"])

        # WP_000031783.1 has multiple hits, best should be VFG046459 (bitscore 677)
        assert hits["WP_000031783.1"].sseqid == "VFG046459"
        assert hits["WP_000031783.1"].bitscore == 677.0
        assert hits["WP_000031783.1"].evalue == pytest.approx(6.50e-248)

        # WP_001273183.1 has best bitscore 1541
        assert hits["WP_001273183.1"].bitscore == 1541.0

    def test_parse_blastp_vfdb_id_extraction(self, tmp_files):
        """Test VFDB ID extraction (before parenthesis)."""
        hits = parse_blastp_best_hits(tmp_files["diamond_file"])

        # Should extract VFG042734 from "VFG042734(gb|NP_249343)"
        assert hits["WP_000242755.1"].sseqid == "VFG042734"

    def test_parse_blastp_empty_file(self, tmp_files):
        """Test parsing empty file returns empty dict."""
        hits = parse_blastp_best_hits(tmp_files["empty_file"])
        assert hits == {}


class TestCollectDetectedSequenceIds:
    """Test sequence ID collection from multiple sources."""

    def test_collect_all_sources(self):
        """Test collecting IDs from all three sources."""
        blast_hits = {"seq1": BestBlastHit("vfg1", 1e-10, 100)}
        tox_preds = {"seq2": 0.9, "seq3": 0.8}
        vf_preds = {"seq3": 0.95, "seq4": 0.85}

        detected = collect_detected_sequence_ids(blast_hits, tox_preds, vf_preds)

        assert len(detected) == 4
        assert "seq1" in detected
        assert "seq2" in detected
        assert "seq3" in detected  # In both tox and vf
        assert "seq4" in detected

    def test_collect_empty_sources(self):
        """Test with empty input dictionaries."""
        detected = collect_detected_sequence_ids({}, {}, {})
        assert len(detected) == 0

    def test_collect_single_source(self):
        """Test with only one source having data."""
        blast_hits = {"seq1": BestBlastHit("vfg1", 1e-10, 100)}
        detected = collect_detected_sequence_ids(blast_hits, {}, {})

        assert len(detected) == 1
        assert "seq1" in detected


class TestWriteDetectedFasta:
    """Test writing detected sequences to FASTA."""

    def test_write_detected_fasta_basic(self, tmp_files):
        """Test writing subset of sequences to FASTA."""
        sequences = read_fasta_to_dict(tmp_files["fasta_file"])
        detected_ids = {"WP_000242755.1", "WP_000031783.1"}

        output_file = tmp_files["tmp_path"] / "detected.fasta"
        write_detected_fasta(sequences, detected_ids, str(output_file))

        # Read output and verify
        output_seqs = read_fasta_to_dict(str(output_file))
        assert len(output_seqs) == 2
        assert "WP_000242755.1" in output_seqs
        assert "WP_000031783.1" in output_seqs

    def test_write_detected_fasta_missing_sequence(self, tmp_files):
        """Test that missing sequences raise ValueError."""
        sequences = read_fasta_to_dict(tmp_files["fasta_file"])
        detected_ids = {"WP_000242755.1", "WP_NONEXISTENT"}

        output_file = tmp_files["tmp_path"] / "detected.fasta"

        with pytest.raises(ValueError) as exc_info:
            write_detected_fasta(sequences, detected_ids, str(output_file))

        assert "FATAL ERROR" in str(exc_info.value)
        assert "WP_NONEXISTENT" in str(exc_info.value)

    def test_write_detected_fasta_multiple_missing(self, tmp_files):
        """Test error message with multiple missing sequences."""
        sequences = read_fasta_to_dict(tmp_files["fasta_file"])
        detected_ids = {"WP_MISSING1", "WP_MISSING2", "WP_MISSING3"}

        output_file = tmp_files["tmp_path"] / "detected.fasta"

        with pytest.raises(ValueError) as exc_info:
            write_detected_fasta(sequences, detected_ids, str(output_file))

        error_msg = str(exc_info.value)
        assert "3 predicted protein(s) not found" in error_msg
        assert "WP_MISSING1" in error_msg
        assert "WP_MISSING2" in error_msg


class TestWriteSupportTable:
    """Test writing support table."""

    def test_write_support_table_basic(self, tmp_files):
        """Test writing support table with all sources."""
        blast_hits = {
            "WP_000031783.1": BestBlastHit("VFG046459", 6.50e-248, 677.0),
            "WP_000242755.1": BestBlastHit("VFG042734", 2.58e-94, 274.0),
        }
        tox_preds = {"WP_000242755.1": 1.0, "WP_001273183.1": 1.0}
        vf_preds = {"WP_000163771.1": 1.0, "WP_001273183.1": 0.9977883}

        output_file = tmp_files["tmp_path"] / "support.tsv"
        write_support_table(blast_hits, tox_preds, vf_preds, str(output_file))

        # Read and verify output
        content = output_file.read_text()
        lines = content.strip().split("\n")

        # Check header
        assert (
            lines[0]
            == "sequence_id\tdetection_method\tsupport_value_type\tsupport_value\tvfdb_hit"
        )

        # Check that all entries are present
        assert any("WP_000031783.1\tblastp" in line for line in lines)
        assert any("WP_000242755.1\tpathofact2_tox" in line for line in lines)
        assert any("WP_001273183.1\tpathofact2_vf" in line for line in lines)
        assert any("VFG046459" in line for line in lines)

    def test_write_support_table_format(self, tmp_files):
        """Test support table format details."""
        blast_hits = {"seq1": BestBlastHit("VFG001", 1e-100, 500.0)}
        tox_preds = {"seq2": 0.95}
        vf_preds = {}

        output_file = tmp_files["tmp_path"] / "support.tsv"
        write_support_table(blast_hits, tox_preds, vf_preds, str(output_file))

        content = output_file.read_text()
        lines = content.strip().split("\n")

        # Check BLASTP line format
        blastp_line = [line for line in lines if "blastp" in line][0]
        fields = blastp_line.split("\t")
        assert fields[0] == "seq1"
        assert fields[1] == "blastp"
        assert fields[2] == "evalue"
        assert fields[4] == "VFG001"

        # Check PathOFact2 line format
        tox_line = [line for line in lines if "pathofact2_tox" in line][0]
        fields = tox_line.split("\t")
        assert fields[0] == "seq2"
        assert fields[1] == "pathofact2_tox"
        assert fields[2] == "probability"
        # For PathOFact2 predictions, vfdb_hit field is empty (no DIAMOND hit)
        # When split by tab, trailing empty fields may be omitted
        assert len(fields) in [
            4,
            5,
        ]  # Either 4 fields (trailing empty omitted) or 5 fields
        if len(fields) == 5:
            assert (
                fields[4] == ""
            )  # If present, should be emptyassert fields[3] == "0.950000"


class TestIntegration:
    """Integration tests for complete workflow."""

    def test_full_workflow(self, tmp_files):
        """Test complete workflow with all components."""
        # Step 1: Read FASTA
        sequences = read_fasta_to_dict(tmp_files["fasta_file"])
        assert len(sequences) == 5

        # Step 2: Parse BLASTP
        blast_hits = parse_blastp_best_hits(tmp_files["diamond_file"])
        assert len(blast_hits) == 3

        # Step 3: Parse PathOFact2 predictions
        tox_preds = parse_pathofact2_predictions_tsv(tmp_files["toxins_file"], 0.6)
        vf_preds = parse_pathofact2_predictions_tsv(tmp_files["virulence_file"], 0.9)
        assert len(tox_preds) == 2
        assert len(vf_preds) == 2

        # Step 4: Collect detected IDs
        detected_ids = collect_detected_sequence_ids(blast_hits, tox_preds, vf_preds)
        # Unique IDs: WP_000031783.1, WP_000242755.1, WP_001273183.1, WP_000163771.1
        assert len(detected_ids) == 4

        # Step 5: Write outputs
        output_fasta = tmp_files["tmp_path"] / "pathofact2.fasta"
        output_tsv = tmp_files["tmp_path"] / "support.tsv"

        write_detected_fasta(sequences, detected_ids, str(output_fasta))
        write_support_table(blast_hits, tox_preds, vf_preds, str(output_tsv))

        # Verify outputs exist and have content
        assert output_fasta.exists()
        assert output_tsv.exists()
        assert output_fasta.stat().st_size > 0
        assert output_tsv.stat().st_size > 0

        # Verify FASTA content
        output_seqs = read_fasta_to_dict(str(output_fasta))
        assert len(output_seqs) == 4

        # Verify TSV content
        tsv_content = output_tsv.read_text()
        assert "blastp" in tsv_content
        assert "pathofact2_tox" in tsv_content
        assert "pathofact2_vf" in tsv_content


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
