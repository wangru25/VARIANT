"""Sequence processing utilities for FASTA files and sequence operations."""

import math
from typing import Dict, List

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def read_fasta(fasta_file: str) -> List[Dict[str, str]]:
    """
    Read sequences from a FASTA file.

    Args:
        fasta_file: Path to the FASTA file

    Returns:
        List of dictionaries containing sequence ID and sequence data

    Raises:
        FileNotFoundError: If the FASTA file doesn't exist
        ValueError: If the FASTA file is empty or invalid
    """
    try:
        records = []
        for record in SeqIO.parse(fasta_file, "fasta"):
            records.append({"seqId": record.id, "seq": str(record.seq)})

        if not records:
            raise ValueError(f"No sequences found in {fasta_file}")

        return records

    except FileNotFoundError:
        raise FileNotFoundError(f"FASTA file not found: {fasta_file}")
    except Exception as e:
        raise ValueError(f"Error reading FASTA file {fasta_file}: {e}")


def write_fasta(headers: List[str], sequences: List[str], file_name: str) -> None:
    """
    Write sequences to a FASTA file.

    Args:
        headers: List of sequence headers/IDs
        sequences: List of sequences
        file_name: Output FASTA file name

    Raises:
        ValueError: If headers and sequences have different lengths
        IOError: If there's an error writing the file
    """
    if len(headers) != len(sequences):
        raise ValueError("Headers and sequences must have the same length")

    try:
        records = []
        for header, seq in zip(headers, sequences):
            records.append(SeqRecord(Seq(seq), id=header, description=""))

        SeqIO.write(records, file_name, "fasta")

    except Exception as e:
        raise OSError(f"Error writing FASTA file {file_name}: {e}")


def clustal_genomes(seq_file_name: str) -> bool:
    """
    Split genome sequences into smaller files for ClustalW processing.

    This function is specifically designed for SARS-CoV-2 genomes and splits
    large sequence files into smaller chunks of 130 sequences each, adding
    a reference sequence to each chunk.

    Args:
        seq_file_name: Path to the input sequence file

    Returns:
        True if processing was successful

    Raises:
        FileNotFoundError: If input files don't exist
        ValueError: If the sequence file is invalid
    """
    try:
        # Read and filter sequences
        records = read_fasta(seq_file_name)
        filtered_records = [
            record for record in records if "NNNNNN" not in record["seq"]
        ]

        # Calculate number of output files needed
        file_name = seq_file_name.split(".")[0]
        num_files = math.ceil(len(filtered_records) / 131)

        # Read reference genome
        ref_file = "../data/refs/NC_045512.fasta"
        ref_records = read_fasta(ref_file)
        ref_seq = ref_records[0]["seq"]
        ref_header = ref_records[0]["seqId"]

        # Create output files
        for i in range(num_files):
            file_name_i = f"{file_name}_{i}.fasta"
            start = i * 130
            end = min(start + 130, len(filtered_records))

            # Extract sequences for this chunk
            chunk_records = filtered_records[start:end]
            sequences = [record["seq"] for record in chunk_records]
            headers = [record["seqId"] for record in chunk_records]

            # Add reference sequence at the beginning
            sequences = [ref_seq] + sequences
            headers = [ref_header] + headers

            # Write the chunk to file
            write_fasta(headers, sequences, file_name_i)

        return True

    except Exception as e:
        raise ValueError(f"Error in clustal_genomes: {e}")


def has_only_valid_nts(sequence: str) -> bool:
    """
    Check if a sequence contains only valid nucleotides.

    Args:
        sequence: DNA sequence to validate

    Returns:
        True if sequence contains only valid nucleotides, False otherwise
    """
    valid_nucleotides = {"A", "T", "C", "G", "-", "N"}
    return all(nt in valid_nucleotides for nt in sequence)


def validate_sequence(
    sequence: str, allow_gaps: bool = True, allow_ambiguous: bool = True
) -> bool:
    """
    Validate a DNA sequence for common issues.

    Args:
        sequence: DNA sequence to validate
        allow_gaps: Whether to allow gap characters ('-')
        allow_ambiguous: Whether to allow ambiguous nucleotides ('N')

    Returns:
        True if sequence is valid, False otherwise
    """
    if not sequence:
        return False

    valid_nucleotides = {"A", "T", "C", "G"}

    if allow_gaps:
        valid_nucleotides.add("-")

    if allow_ambiguous:
        valid_nucleotides.add("N")

    return all(nt in valid_nucleotides for nt in sequence.upper())


def get_sequence_statistics(sequence: str) -> Dict[str, int]:
    """
    Calculate basic statistics for a DNA sequence.

    Args:
        sequence: DNA sequence to analyze

    Returns:
        Dictionary containing sequence statistics
    """
    stats = {
        "length": len(sequence),
        "A": sequence.upper().count("A"),
        "T": sequence.upper().count("T"),
        "C": sequence.upper().count("C"),
        "G": sequence.upper().count("G"),
        "N": sequence.upper().count("N"),
        "gaps": sequence.count("-"),
    }

    stats["GC_content"] = stats["G"] + stats["C"]
    stats["AT_content"] = stats["A"] + stats["T"]

    return stats
