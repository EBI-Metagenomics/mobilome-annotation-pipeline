# Test dataset

This test is composed by two input fasta files for a positive and a negative tests.

The [positive test](https://raw.githubusercontent.com/EBI-Metagenomics/nf-modules/refs/heads/feature/amr_prediction/subworkflows/ebi-metagenomics/amr_annotation/tests/data/pos_minitest.faa) consist of a protein FASTA file containing eleven sequences with known annotation as AMGs.

The [negative test](https://raw.githubusercontent.com/EBI-Metagenomics/nf-modules/refs/heads/feature/amr_prediction/subworkflows/ebi-metagenomics/amr_annotation/tests/data/neg_minitest.faa) consist of a protein FASTA file containing one non-AMG sequence.

The CARD database for RGI is downloaded and formatted during test. AMRFinderplus databse is also downloaded o the fly. Deeparg database have been manually curated for testing purposes and aggressively trimmed down.
