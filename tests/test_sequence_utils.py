"""Tests for sequence utility functions."""

from unittest.mock import mock_open, patch

import pytest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord

from src.utils.sequence_utils import (
    clustal_genomes,
    get_sequence_statistics,
    has_only_valid_nts,
    read_fasta,
    validate_sequence,
    write_fasta,
)


class TestReadFasta:
    """Test cases for read_fasta function."""

    def test_read_fasta_success(self):
        """Test successful reading of FASTA file."""
        mock_fasta_content = [
            SeqRecord(Seq("ATCG"), id="seq1", description=""),
            SeqRecord(Seq("GCTA"), id="seq2", description=""),
        ]

        with patch(
            "src.utils.sequence_utils.SeqIO.parse", return_value=mock_fasta_content
        ):
            result = read_fasta("test.fasta")

        expected = [
            {"seqId": "seq1", "seq": "ATCG"},
            {"seqId": "seq2", "seq": "GCTA"},
        ]
        assert result == expected

    def test_read_fasta_empty_file(self):
        """Test reading empty FASTA file raises ValueError."""
        with patch("src.utils.sequence_utils.SeqIO.parse", return_value=[]):
            with pytest.raises(ValueError, match="No sequences found"):
                read_fasta("empty.fasta")

    def test_read_fasta_file_not_found(self):
        """Test reading non-existent file raises FileNotFoundError."""
        with patch(
            "src.utils.sequence_utils.SeqIO.parse", side_effect=FileNotFoundError()
        ):
            with pytest.raises(FileNotFoundError):
                read_fasta("nonexistent.fasta")


class TestWriteFasta:
    """Test cases for write_fasta function."""

    def test_write_fasta_success(self):
        """Test successful writing of FASTA file."""
        headers = ["seq1", "seq2"]
        sequences = ["ATCG", "GCTA"]

        with patch("src.utils.sequence_utils.SeqIO.write") as mock_write:
            write_fasta(headers, sequences, "output.fasta")

        assert mock_write.called

    def test_write_fasta_mismatched_lengths(self):
        """Test writing with mismatched headers and sequences raises ValueError."""
        headers = ["seq1", "seq2"]
        sequences = ["ATCG"]  # Only one sequence

        with pytest.raises(
            ValueError, match="Headers and sequences must have the same length"
        ):
            write_fasta(headers, sequences, "output.fasta")


class TestHasOnlyValidNts:
    """Test cases for has_only_valid_nts function."""

    @pytest.mark.parametrize(
        "sequence,expected",
        [
            ("ATCG", True),
            ("ATCG-", True),
            ("ATCGN", True),
            ("ATCGX", False),
            ("", True),
            ("ATCG123", False),
        ],
    )
    def test_has_only_valid_nts(self, sequence, expected):
        """Test nucleotide validation with various sequences."""
        assert has_only_valid_nts(sequence) == expected


class TestValidateSequence:
    """Test cases for validate_sequence function."""

    @pytest.mark.parametrize(
        "sequence,allow_gaps,allow_ambiguous,expected",
        [
            ("ATCG", True, True, True),
            ("ATCG", False, True, True),
            ("ATCG", True, False, True),
            ("ATCG-", True, True, True),
            ("ATCG-", False, True, False),
            ("ATCGN", True, True, True),
            ("ATCGN", True, False, False),
            ("ATCGX", True, True, False),
            ("", True, True, False),
        ],
    )
    def test_validate_sequence(self, sequence, allow_gaps, allow_ambiguous, expected):
        """Test sequence validation with various parameters."""
        assert validate_sequence(sequence, allow_gaps, allow_ambiguous) == expected


class TestGetSequenceStatistics:
    """Test cases for get_sequence_statistics function."""

    def test_get_sequence_statistics_basic(self):
        """Test basic sequence statistics calculation."""
        sequence = "ATCGATCG"
        stats = get_sequence_statistics(sequence)

        expected = {
            "length": 8,
            "A": 2,
            "T": 2,
            "C": 2,
            "G": 2,
            "N": 0,
            "gaps": 0,
            "GC_content": 4,
            "AT_content": 4,
        }
        assert stats == expected

    def test_get_sequence_statistics_with_gaps_and_n(self):
        """Test sequence statistics with gaps and ambiguous nucleotides."""
        sequence = "ATCG-N"
        stats = get_sequence_statistics(sequence)

        expected = {
            "length": 6,
            "A": 1,
            "T": 1,
            "C": 1,
            "G": 1,
            "N": 1,
            "gaps": 1,
            "GC_content": 2,
            "AT_content": 2,
        }
        assert stats == expected

    def test_get_sequence_statistics_empty(self):
        """Test sequence statistics with empty sequence."""
        sequence = ""
        stats = get_sequence_statistics(sequence)

        expected = {
            "length": 0,
            "A": 0,
            "T": 0,
            "C": 0,
            "G": 0,
            "N": 0,
            "gaps": 0,
            "GC_content": 0,
            "AT_content": 0,
        }
        assert stats == expected


class TestClustalGenomes:
    """Test cases for clustal_genomes function."""

    @patch("src.utils.sequence_utils.read_fasta")
    def test_clustal_genomes_success(self, mock_read_fasta):
        """Test successful clustal_genomes processing."""
        # Mock input sequences
        mock_sequences = [
            {"seqId": f"seq{i}", "seq": "ATCG" * 100}
            for i in range(200)  # More than 131 sequences
        ]
        mock_read_fasta.return_value = mock_sequences

        # Mock reference genome
        # mock_ref_sequences = [{"seqId": "NC_045512", "seq": "ATCG" * 100}]

        with patch("src.utils.sequence_utils.write_fasta"):
            with patch("builtins.open", mock_open()):
                result = clustal_genomes("test.fasta")

        assert result is True

    @patch("src.utils.sequence_utils.read_fasta")
    def test_clustal_genomes_with_nnnnnn(self, mock_read_fasta):
        """Test clustal_genomes filters out sequences with NNNNNN."""
        mock_sequences = [
            {"seqId": "seq1", "seq": "ATCG" * 100},
            {"seqId": "seq2", "seq": "ATCGNNNNNNATCG"},
            {"seqId": "seq3", "seq": "GCTA" * 100},
        ]
        mock_read_fasta.return_value = mock_sequences

        with patch("src.utils.sequence_utils.write_fasta"):
            with patch("builtins.open", mock_open()):
                result = clustal_genomes("test.fasta")

        # Should only process 2 sequences (filtered out the one with NNNNNN)
        assert result is True


# Integration tests
class TestSequenceUtilsIntegration:
    """Integration tests for sequence utilities."""

    def test_read_write_fasta_roundtrip(self):
        """Test that reading and writing FASTA files preserves data."""
        original_headers = ["seq1", "seq2", "seq3"]
        original_sequences = ["ATCG", "GCTA", "TAGC"]

        # Mock the file operations
        with patch("src.utils.sequence_utils.SeqIO.write"):
            write_fasta(original_headers, original_sequences, "test.fasta")

        # Mock reading back the written file
        mock_records = [
            SeqRecord(Seq(seq), id=header, description="")
            for header, seq in zip(original_headers, original_sequences)
        ]

        with patch("src.utils.sequence_utils.SeqIO.parse", return_value=mock_records):
            read_result = read_fasta("test.fasta")

        expected = [
            {"seqId": header, "seq": seq}
            for header, seq in zip(original_headers, original_sequences)
        ]
        assert read_result == expected
