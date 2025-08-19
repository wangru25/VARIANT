# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-06-23 14:27:00
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-19 17:40:18
Email: wang.rui@nyu.edu
FilePath: /VARIANT/src/core/reference_genome.py
Description:
'''


class ReferenceGenome:
    """
    A class to handle reference genome data and operations.

    This class loads and manages reference genome sequences from FASTA files,
    providing easy access to genome data for mutation analysis.

    Attributes:
        ref_name (str): The name/ID of the reference genome
        ref_seq (str): The reference genome sequence
        ref_file_path (str): Path to the reference genome FASTA file
    """

    def __init__(self, ref_file_path: str) -> None:
        """
        Initialize the ReferenceGenome with data from a FASTA file.

        Args:
            ref_file_path: Path to the reference genome FASTA file

        Raises:
            FileNotFoundError: If the reference file doesn't exist
            ValueError: If the FASTA file is empty or invalid
        """
        self.ref_file_path = ref_file_path
        self._load_reference_data()

    def _load_reference_data(self) -> None:
        """Load reference genome data from the FASTA file."""
        try:
            fasta_data = read_fasta(self.ref_file_path)
            if not fasta_data:
                raise ValueError(f"No valid sequences found in {self.ref_file_path}")

            self.ref_name = fasta_data[0]["seqId"]
            self.ref_seq = fasta_data[0]["seq"]

        except FileNotFoundError:
            raise FileNotFoundError(
                f"Reference genome file not found: {self.ref_file_path}"
            )
        except Exception as e:
            raise ValueError(f"Error loading reference genome: {e}")

    def get_sequence_length(self) -> int:
        """
        Get the length of the reference genome sequence.

        Returns:
            Length of the reference genome sequence
        """
        return len(self.ref_seq)

    def get_sequence_region(self, start: int, end: int) -> str:
        """
        Get a specific region of the reference genome sequence.

        Args:
            start: Start position (0-based)
            end: End position (exclusive)

        Returns:
            The sequence region from start to end

        Raises:
            ValueError: If start or end positions are invalid
        """
        if start < 0 or end > len(self.ref_seq) or start >= end:
            raise ValueError(
                f"Invalid region: start={start}, end={end}, length={len(self.ref_seq)}"
            )

        return self.ref_seq[start:end]

    def to_dict(self) -> Dict[str, Any]:
        """
        Convert the reference genome to a dictionary representation.

        Returns:
            Dictionary containing reference genome data
        """
        return {
            "ref_name": self.ref_name,
            "ref_seq": self.ref_seq,
            "ref_file_path": self.ref_file_path,
            "length": len(self.ref_seq),
        }

    def __str__(self) -> str:
        """String representation of the reference genome."""
        return f"ReferenceGenome(name='{self.ref_name}', length={len(self.ref_seq)})"

    def __repr__(self) -> str:
        """Detailed string representation of the reference genome."""
        return f"ReferenceGenome(ref_file_path='{self.ref_file_path}')"
