# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-08-04 13:30:27
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-05 16:07:44
Email: wang.rui@nyu.edu
FilePath: /viralytics-mut/src/utils/sequence_utils.py
Description: 
'''
"""Sequence processing utilities for FASTA files and sequence operations."""

import re
from typing import Dict, List, Tuple, Union

from Bio import SeqIO


def read_fasta(file_path: str) -> List[Dict[str, str]]:
    """
    Read a FASTA file and return a list of dictionaries containing sequence data.

    Args:
        file_path (str): Path to the FASTA file

    Returns:
        List[Dict[str, str]]: List of dictionaries with 'seqId' and 'seq' keys
    """
    sequences = []
    try:
        with open(file_path, "r") as handle:
            for record in SeqIO.parse(handle, "fasta"):
                sequences.append({"seqId": record.id, "seq": str(record.seq)})
    except FileNotFoundError:
        print(f"Error: File {file_path} not found.")
        return []
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
        return []
    return sequences


def write_fasta(headers: List[str], sequences: List[str], file_name: str) -> None:
    """
    Write sequences to a FASTA file.

    Args:
        headers (List[str]): List of sequence headers
        sequences (List[str]): List of sequences
        file_name (str): Output file name
    """
    with open(file_name, "w") as f:
        for header, seq in zip(headers, sequences):
            f.write(f">{header}\n{seq}\n")


def has_only_valid_nts(sequence: str) -> bool:
    """
    Check if a sequence contains only valid nucleotides.

    Args:
        sequence (str): DNA sequence to check

    Returns:
        bool: True if sequence contains only valid nucleotides, False otherwise
    """
    valid_nts = set("ATCGN-")
    return all(nt in valid_nts for nt in sequence.upper())


def parse_gene_coordinates(coordinates: str) -> List[Tuple[int, int]]:
    """
    Parse gene coordinates, handling both simple ranges and join operations.
    
    Args:
        coordinates (str): Gene coordinates in format like "13442..16236" or "join(13442..13468,13468..16236)"
    
    Returns:
        List[Tuple[int, int]]: List of (start, end) coordinate pairs (1-based)
    
    Examples:
        >>> parse_gene_coordinates("13442..16236")
        [(13442, 16236)]
        >>> parse_gene_coordinates("join(13442..13468,13468..16236)")
        [(13442, 13468), (13468, 16236)]
    """
    if coordinates.startswith("join(") and coordinates.endswith(")"):
        # Handle join operation
        join_content = coordinates[5:-1]  # Remove "join(" and ")"
        segments = join_content.split(",")
        coordinate_pairs = []
        
        for segment in segments:
            if ".." in segment:
                start_str, end_str = segment.split("..")
                start_pos = int(start_str.strip())
                end_pos = int(end_str.strip())
                coordinate_pairs.append((start_pos, end_pos))
        
        return coordinate_pairs
    else:
        # Handle simple range
        if ".." in coordinates:
            start_str, end_str = coordinates.split("..")
            start_pos = int(start_str.strip())
            end_pos = int(end_str.strip())
            return [(start_pos, end_pos)]
        else:
            raise ValueError(f"Invalid coordinate format: {coordinates}")


def extract_gene_sequence_with_join(genome_sequence: str, coordinates: str) -> str:
    """
    Extract gene sequence from genome, handling both simple ranges and join operations.
    
    Args:
        genome_sequence (str): Complete genome sequence
        coordinates (str): Gene coordinates in format like "13442..16236" or "join(13442..13468,13468..16236)"
    
    Returns:
        str: Extracted gene sequence
    
    Examples:
        >>> extract_gene_sequence_with_join(genome_seq, "13442..16236")
        "ATCG..."
        >>> extract_gene_sequence_with_join(genome_seq, "join(13442..13468,13468..16236)")
        "ATCG..."  # Joined sequence from both segments
    """
    coordinate_pairs = parse_gene_coordinates(coordinates)
    joined_sequence = ""
    
    for start_pos, end_pos in coordinate_pairs:
        # Convert to 0-based coordinates for sequence extraction
        start_idx = start_pos - 1
        end_idx = end_pos
        segment_sequence = genome_sequence[start_idx:end_idx]
        joined_sequence += segment_sequence
    
    return joined_sequence


def create_genome_to_joined_mapping(coordinates: str) -> Dict[int, int]:
    """
    Create a mapping from genome positions to joined sequence positions.
    
    Args:
        coordinates (str): Gene coordinates in format like "13442..16236" or "join(13442..13468,13468..16236)"
    
    Returns:
        Dict[int, int]: Mapping from genome position to joined sequence position
    
    Examples:
        >>> mapping = create_genome_to_joined_mapping("join(13442..13468,13468..16236)")
        >>> mapping[13442]  # First position in first segment
        1
        >>> mapping[13468]  # First position in second segment
        27  # 27 nucleotides from first segment
    """
    coordinate_pairs = parse_gene_coordinates(coordinates)
    genome_to_joined = {}
    joined_pos = 1
    
    for start_pos, end_pos in coordinate_pairs:
        for genome_pos in range(start_pos, end_pos + 1):
            genome_to_joined[genome_pos] = joined_pos
            joined_pos += 1
    
    return genome_to_joined


def calculate_amino_acid_position_from_joined(joined_position: int) -> int:
    """
    Calculate amino acid position from joined sequence position.
    
    Args:
        joined_position (int): Position in joined sequence (1-based)
    
    Returns:
        int: Amino acid position (1-based)
    """
    return (joined_position - 1) // 3 + 1


def extract_genome_id(full_header: str, virus_name: str = "") -> str:
    """
    Extract genome ID from MSA sequence header based on virus type.
    
    Args:
        full_header: Full sequence header from MSA file
        virus_name: Name of the virus (e.g., "SARS-CoV-2", "ZaireEbola")
    
    Returns:
        Extracted genome ID
    
    Examples:
        >>> extract_genome_id("hCoV-19/Shanghai/SJTU-235817/2022|EPI_ISL_16327572|2022-12-07", "SARS-CoV-2")
        "EPI_ISL_16327572"
        >>> extract_genome_id("AB050936v1", "ZaireEbola")
        "AB050936v1"
        >>> extract_genome_id("NC_002549.1|Zaire", "ZaireEbola")
        "NC_002549.1|Zaire"
    """
    # For SARS-CoV-2, extract EPI_ISL ID from the middle part
    if virus_name.upper() in ["SARS-COV-2", "SARS_COV_2", "SARSCOV2"] and "|" in full_header:
        parts = full_header.split("|")
        # Look for EPI_ISL pattern
        for part in parts:
            if "EPI_ISL" in part:
                return part.strip()
    
    # For other viruses or if no specific pattern found, return the full header
    # This preserves the original behavior for non-SARS-CoV-2 viruses
    return full_header


def calculate_amino_acid_position_from_genome(genome_position: int, genome_to_joined_mapping: Dict[int, int]) -> int:
    """
    Calculate amino acid position from genome position using the mapping.
    
    Args:
        genome_position (int): Position in genome (1-based)
        genome_to_joined_mapping (Dict[int, int]): Mapping from genome to joined positions
    
    Returns:
        int: Amino acid position (1-based), or None if position not in mapping
    """
    if genome_position in genome_to_joined_mapping:
        joined_position = genome_to_joined_mapping[genome_position]
        return calculate_amino_acid_position_from_joined(joined_position)
    return None


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
    import math
    
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
