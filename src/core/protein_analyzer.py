# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-08-04 13:30:27
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-05 17:58:57
Email: wang.rui@nyu.edu
FilePath: /viralytics-mut/src/core/protein_analyzer.py
Description:
'''
# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-08-04 13:30:27
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-05 17:46:24
Email: wang.rui@nyu.edu
FilePath: /viralytics-mut/src/core/protein_analyzer.py
Description:
'''
"""
Author: Rui Wang
Date: 2023-11-01 09:21:48
LastModifiedBy: Rui Wang
LastEditTime: 2023-11-17 10:37:16
Email: rw3594@nyu.edu
FilePath: /MutParser/src/core/protein_analyzer.py
Description:
"""

from itertools import takewhile
from typing import Dict, List, Tuple

from Bio import SeqIO, pairwise2
from Bio.Data import CodonTable
from fuzzysearch import find_near_matches

from .reference_genome import ReferenceGenome


class Proteome:
    def __init__(self, proteome_dir: str):
        self.proteome_dir = proteome_dir


class ORFProcessor:
    """
    Class for translating DNA sequences and analyzing the effect of mutations on Open Reading Frames (ORFs).
    """

    def __init__(self, ref_genome_instance: ReferenceGenome) -> None:
        self.ref_seq = ref_genome_instance.ref_seq
        self.codon_table = CodonTable.unambiguous_dna_by_id[1]

    def translate_dna_3_frames(
        self, seq: str, codon_type: str = "reference"
    ) -> List[str]:
        """
        Translate a DNA sequence in 3 frames and return the resulting ORFs.

        Args:
            seq (str): The DNA sequence to translate.
            codon_type (str): The codon table type. Default is 'reference'.

        Returns:
            List[str]: A list of ORFs resulting from the translation.
        """
        protein_seqs = ["", "", ""]

        for i in range(len(seq) - 2):
            codon = seq[i : i + 3]
            protein_seqs[i % 3] += self.codon_table.forward_table.get(codon, "_")

        if codon_type == "reference":
            protein_seqs = self.get_good_orfs(protein_seqs)

        return [orf for orf in protein_seqs if orf]

    def translate_dna_3_frames_row(self, seq: str) -> List[str]:
        """
        Translate a DNA sequence in 3 frames for a specific row and return the resulting ORFs.

        Args:
            seq (str): The DNA sequence to translate.

        Returns:
            List[str]: A list of ORFs resulting from the translation.
        """
        protein_seqs = ["", "", ""]

        for i in range(len(seq)):
            codon = seq[i : i + 3]
            if len(codon) == 3:
                protein_seqs[i % 3] += self.codon_table.forward_table.get(codon, "_")

        return [orf for orf in protein_seqs if orf]

    def get_orfs_snps(self, snp_pos: int, mutated_to: str) -> List[List[str]]:
        """
        Get the ORFs before and after a Single Nucleotide Polymorphism (SNP).

        Args:
            snp_pos (int): The position of the SNP in genome order (starts from 1).
            mutated_to (str): The nucleotide that the original nucleotide is mutated to.

        Returns:
            List[List[str]]: A list containing two lists of ORFs, one before the mutation and one after.
        """
        pos = int(snp_pos) - 1

        # Increase window size to capture more downstream effects
        # Original: pos - 2 : pos + 28 (30 nucleotides)
        # New: pos - 5 : pos + 55 (60 nucleotides) to capture more amino acids
        seq_t = self.ref_seq[pos - 5 : pos + 55]
        orfs = self.translate_dna_3_frames(seq_t)

        seq_m = seq_t[:5] + mutated_to + seq_t[6:]
        m_orfs = self.translate_dna_3_frames(seq_m, codon_type="mutant")

        return [orfs, m_orfs]

    def get_orfs_dels(self, poss: str) -> Tuple[List[str], List[str]]:
        """
        Get ORFs before and after a deletion.

        Args:
            poss (str): The position range of the deletion in the genome in the format 'start:end'.

        Returns:
            Tuple[List[str], List[str]]: A tuple of two lists, ORFs before and after the deletion.
        """
        pos_a, pos_b = map(int, poss.split(":"))
        pos_a -= 1

        num_pre_nts = 15
        num_post_nts = 15
        seq_t = self.ref_seq[pos_a - num_pre_nts : pos_b + num_post_nts]
        seq_t_del = (
            self.ref_seq[pos_a - num_pre_nts : pos_a]
            + self.ref_seq[pos_b : pos_b + num_post_nts]
        )

        orfs = self.translate_dna_3_frames(seq_t, codon_type="reference")
        m_orfs = self.translate_dna_3_frames(seq_t_del, codon_type="mutant")

        return orfs, m_orfs

    def get_orfs_ins(self, pos: int, ins_seq: str) -> Tuple[List[str], List[str]]:
        """
        Get ORFs before and after an insertion.

        Args:
            pos (int): The position of the insertion in the genome.
            ins_seq (str): The sequence of the insertion.

        Returns:
            Tuple[List[str], List[str]]: A tuple of two lists, ORFs before and after the insertion.
        """
        pos -= 1
        num_pre_nts = 20

        seq_t = self.ref_seq[pos - num_pre_nts : pos]
        seq_t_ins = seq_t + ins_seq

        orfs = self.translate_dna_3_frames(seq_t, codon_type="reference")
        m_orfs = self.translate_dna_3_frames(seq_t_ins, codon_type="mutant")

        return orfs, m_orfs

    @staticmethod
    def get_good_orfs(protein_seqs: List[str]) -> List[str]:
        """
        Filter out low-quality ORFs based on certain criteria like length.

        Args:
            protein_seqs (List[str]): List of translated protein sequences (ORFs).

        Returns:
            List[str]: A list of good quality ORFs.
        """
        # Filter based on the length of the ORF and allow some stop codons
        # Only filter out ORFs that are too short or completely invalid
        good_orfs = []
        for orf in protein_seqs:
            # Keep ORFs that have at least 3 amino acids and are not completely invalid
            if (
                len(orf) >= 3 and orf.count("_") < len(orf) * 0.5
            ):  # Allow up to 50% stop codons
                good_orfs.append(orf)
        return good_orfs


class AminoAcidMutationProcessor:
    """
    Class for analyzing possible mutations between codons and amino acids.

    Attributes:
        codon_table (CodonTable): A CodonTable object providing codon-to-amino acid mappings.
    """

    def __init__(self):
        """
        Initialize the CodonMutationAnalyzer with a specific codon table.

        Args:
            codon_table (CodonTable): A CodonTable object for codon-to-amino acid mapping.
        """
        self.codon_table = CodonTable.unambiguous_dna_by_id[1]

    @staticmethod
    def is_mutable_codon(codon_a: str, codon_b: str) -> bool:
        """
        Determine if two codons can mutate into each other by comparing their nucleotides.

        Args:
            codon_a (str): The first codon to compare.
            codon_b (str): The second codon to compare.

        Returns:
            bool: True if the codons can mutate into each other, False otherwise.
        """
        return any(
            nucleotide_a == nucleotide_b
            for nucleotide_a, nucleotide_b in zip(codon_a, codon_b)
        )

    def is_mutable_aa(self, aa1: str, aa2: str) -> bool:
        """
        Determine if there is a possible mutation between two amino acids based on codon mappings.

        Args:
            aa1 (str): The first amino acid to compare.
            aa2 (str): The second amino acid to compare.

        Returns:
            bool: True if a mutation is possible between the two amino acids, False otherwise.
        """
        for codon_a, amino_acid_a in self.codon_table.forward_table.items():
            for codon_b, amino_acid_b in self.codon_table.forward_table.items():
                if (
                    codon_a != codon_b
                    and self.is_mutable_codon(codon_a, codon_b)
                    and amino_acid_a == aa1
                    and amino_acid_b == aa2
                ):
                    return True
        return False

    @staticmethod
    def get_row_aa_mutations(orf: str, morf: str, match_start: int) -> str:
        """
        Get amino acid mutations by comparing original and mutated ORF sequences.

        This method scans through the entire ORF sequence to find all amino acid
        differences between the original and mutated sequences.

        Args:
            orf (str): Original ORF sequence
            morf (str): Mutated ORF sequence
            match_start (int): Starting position of the match in the reference protein

        Returns:
            str: Comma-separated list of mutations in format 'original_aa{position}mutated_aa'
        """
        mutations = []

        # Ensure both sequences have the same length for comparison
        min_length = min(len(orf), len(morf))

        for i in range(min_length):
            if orf[i] != morf[i]:
                # Calculate actual protein position by adding match_start and current position
                actual_protein_pos = match_start + i + 1
                mutation = f"{orf[i]}{actual_protein_pos}{morf[i]}"
                mutations.append(mutation)

        # If there are no mutations found, return empty string
        if not mutations:
            return ""

        # Return comma-separated list of mutations
        return ",".join(mutations)

    @staticmethod
    def get_aa_dels(orf: str, m_orf: str, snp_r: int) -> List[str]:
        """
        Get amino acid deletions from aligned sequences.

        Args:
            orf (str): Original reading frame sequence.
            m_orf (str): Mutated reading frame sequence.
            snp_r (int): SNP position in the reading frame.

        Returns:
            List[str]: List of amino acid deletions in the format '{original_aa}{position}{deletion}'.
        """
        # Check if either sequence is None or empty
        if not orf or not m_orf:
            return []

        try:
            alignments = pairwise2.align.localms(orf, m_orf, 2, -1, -1, -1)
        except (TypeError, ValueError):
            alignments = pairwise2.align.localxx(orf, m_orf)
        del_aas = []

        pre_matched = sum(
            1
            for _ in takewhile(
                lambda x: x[0] == x[1], zip(alignments[0][0], alignments[0][1])
            )
        )
        i_snp_r = snp_r + pre_matched

        m_len = 0
        for a, b in zip(alignments[0][0], alignments[0][1]):
            if a != b:
                m_len += 1
                pos_aa = i_snp_r + m_len
                if b == "-":
                    del_aas.append(a + str(pos_aa) + "del")
                else:
                    del_aas.append(a + str(pos_aa) + b)

        return del_aas

    @staticmethod
    def get_aa_ins(orf: str, m_orf: str, snp_r: int) -> str:
        """
        Get amino acid insertions from aligned sequences.

        Args:
            orf (str): Original reading frame sequence.
            m_orf (str): Mutated reading frame sequence.
            snp_r (int): SNP position in the reading frame.

        Returns:
            str: Amino acid insertions in the format 'ins{position}{insertions}'.
        """
        if orf and m_orf:
            try:
                alignments = pairwise2.align.localms(orf, m_orf, 2, -1, -1, -1)
            except (TypeError, ValueError):
                alignments = pairwise2.align.localxx(orf, m_orf)
            pre_matched = sum(
                a == b for a, b in zip(alignments[0][0], alignments[0][1])
            )
            i_snp_r = snp_r + pre_matched

            ins_aas = "".join(
                b
                for a, b in zip(alignments[0][0], alignments[0][1])
                if a != b and a == "-"
            )
            return f"ins{i_snp_r}{ins_aas}" if ins_aas else ""

        return ""


class AlignmentProcessor:
    """
    Class for processing alignments and identifying mutations in sequences.
    """

    def __init__(
        self, proteome: Proteome, aa_mutation_processor: AminoAcidMutationProcessor
    ):
        self.proteome = proteome
        self.proteome_dir = self.proteome.proteome_dir
        self.aa_mutation_processor = aa_mutation_processor

    @staticmethod
    def align_local(query_sequence: str, found_sequence: str) -> List[str]:
        """
        Perform local alignment between two sequences.

        Args:
            query_sequence (str): The query sequence.
            found_sequence (str): The found sequence.

        Returns:
            List[str]: A list containing the aligned sequences and the score.
        """
        # Check if either sequence is None or empty
        if not query_sequence or not found_sequence:
            return ["", 0]

        # Use a simpler alignment method that's more compatible
        try:
            # Try the most common API
            alignments = pairwise2.align.localms(
                query_sequence, found_sequence, 2, -1, -1, -1
            )
        except (TypeError, ValueError):
            # Fallback to a simpler alignment method
            alignments = pairwise2.align.localxx(query_sequence, found_sequence)
        if not alignments:
            return ["", 0]

        match = [
            "|" if a == b else "*" for a, b in zip(alignments[0][0], alignments[0][1])
        ]
        score = match.count("|")
        aligned_sequences = "\n".join(
            [alignments[0][0], "".join(match), alignments[0][1]]
        )

        return [aligned_sequences, score]

    @staticmethod
    def most_match_orf_align(self, orf: str, ORFs: List[str]) -> str:
        """
        Find the ORF with the most matches to a given ORF.

        Args:
            orf (str): The ORF to compare.
            ORFs (List[str]): A list of ORFs to compare.

        Returns:
            str: The ORF with the most matches.
        """
        max_matches = 0
        best_ORF = None

        for ORF in ORFs:
            alignments = pairwise2.align.localxx(orf, ORF)
            matches = sum(a == b for a, b in zip(alignments[0][0], alignments[0][1]))

            if matches > max_matches:
                max_matches = matches
                best_ORF = ORF

        return best_ORF

    @staticmethod
    def most_match_orf(target_orf: str, candidate_orfs: List[str]) -> str:
        """
        Find the ORF from a list of ORFs that most closely matches the given ORF.

        Args:
            arget_orf: The ORF to match.
            candidate_orfs: A list of candidate ORFs.

        Returns:
            str: The ORF that most closely matches the target ORF, or None if no close match is found.
        """
        max_distance = 5  # detect maximum 5 amino acids deletions
        if not target_orf or not candidate_orfs:
            return None

        distances = [
            min(
                (
                    match.dist
                    for match in find_near_matches(
                        target_orf, orf, max_l_dist=max_distance
                    )
                ),
                default=1000,
            )
            for orf in candidate_orfs
        ]

        min_distance = min(distances)
        if min_distance == 1000:
            return None

        closest_match_index = distances.index(min_distance)

        return candidate_orfs[closest_match_index]

    def align_orfs_to_ref_proteome_snps(
        self, orfs: List[str], m_orfs: List[str], dna_pos: int = None, mutation_info: str = None, ref_seq: str = None
    ) -> List[Dict[str, str]]:
        """
        Align ORFs to the reference proteome and identify SNPs.

        Args:
            orfs (List[str]): Original ORFs.
            m_orfs (List[str]): Mutated ORFs.
            dna_pos (int, optional): DNA position of the mutation for silent mutation detection.
            mutation_info (str, optional): Mutation information in format "A->B".

        Returns:
            List[Dict[str, str]]: List of protein mutations.
        """
        protein_mutations = []
        on_cds_cnt = 0

        # First, try direct DNA-to-protein mapping if DNA position is provided
        if dna_pos is not None and mutation_info is not None:
            direct_mutation = self._get_direct_protein_mutation(dna_pos, mutation_info, ref_seq)
            if direct_mutation:
                protein_mutations.append(direct_mutation)
                on_cds_cnt += 1
                return protein_mutations  # Return immediately if direct mapping succeeds

        # If direct mapping didn't work, try the original ORF alignment approach
        # This should not happen often if coordinates are correct, but serves as fallback
        if on_cds_cnt == 0:
            for orf, m_orf in zip(orfs, m_orfs):
                if not orf or not m_orf:  # Skip empty ORFs
                    continue

                m_orf = self.most_match_orf(orf, m_orfs)
                if not m_orf:  # Skip if no good match found
                    continue

                for record in SeqIO.parse(self.proteome_dir, "fasta"):
                    # Use record.description to get the full header, not record.id which gets truncated at first space
                    header = record.description.split("|")[1]
                    ref_pro_seq = str(record.seq).replace(" ", "")

                    # Try exact match first
                    matches = find_near_matches(orf, ref_pro_seq, max_l_dist=0)

                    # Only allow fuzzy matching for longer proteins and ORFs to prevent false positives
                    # Don't allow fuzzy matching for very short proteins (< 100 aa) or very short ORFs (< 15 aa)
                    if not matches and len(ref_pro_seq) >= 100 and len(orf) >= 15:
                        matches = find_near_matches(orf, ref_pro_seq, max_l_dist=2)

                    if len(matches) == 1:
                        seq_matched = ref_pro_seq[matches[0].start : matches[0].end]
                        aligned_seqs, score = self.align_local(m_orf, seq_matched)

                        if score != 0:
                            # Here, call the function from the GenomeAnalyzer class
                            mutation_aas = self.aa_mutation_processor.get_row_aa_mutations(
                                orf, m_orf, matches[0].start
                            )

                            if mutation_aas:
                                # Regular mutation
                                protein_mutation = {
                                    "protein": header,
                                    "mutation": mutation_aas,
                                }
                                protein_mutations.append(protein_mutation)
                                on_cds_cnt += 1

        if on_cds_cnt == 0:
            protein_mutations.append({"protein": "None-CDS", "mutation": "NA"})

        return protein_mutations

    def _get_direct_protein_mutation(self, dna_pos: int, mutation_info: str, ref_seq: str = None) -> Dict[str, str]:
        """
        Directly map DNA position to protein mutation without relying on ORF alignment.
        Handles both simple ranges and join operations in gene coordinates.

        Args:
            dna_pos (int): DNA position (1-based)
            mutation_info (str): Mutation information in format "A->B"

        Returns:
            Dict[str, str]: Protein mutation dict or None if not found
        """
        # Parse mutation info
        original_nt, mutated_nt = mutation_info.split("->")

        for record in SeqIO.parse(self.proteome_dir, "fasta"):
            # Use record.description to get the full header, not record.id which gets truncated at first space
            parts = record.description.split("|")

            # Handle different FASTA header formats
            if len(parts) >= 3:
                # Format: NP_066243.1|nucleoprotein|470..2689
                header = parts[1]
                coordinates = parts[2]
            elif len(parts) == 2:
                # Format: NP_066243.1|nucleoprotein_470..2689
                header = parts[1]
                # Try to extract coordinates from the header
                if "_" in header:
                    header_parts = header.split("_")
                    if len(header_parts) > 1:
                        coordinates = header_parts[-1]  # Last part should be coordinates
                        header = "_".join(header_parts[:-1])  # Rest is the protein name
                    else:
                        continue  # Skip if we can't parse
                else:
                    continue  # Skip if we can't parse
            else:
                continue  # Skip if format is not recognized

            # Import join operation utilities
            from ..utils.sequence_utils import (
                calculate_amino_acid_position_enhanced,
                create_genome_to_joined_mapping_enhanced,
                parse_gene_coordinates_enhanced,
            )

            try:
                # Parse coordinates (handles simple ranges, join operations, and complement notation)
                coordinate_pairs, is_complement = parse_gene_coordinates_enhanced(coordinates)

                # Check if DNA position is within any of the coordinate pairs
                position_in_gene = False
                for start_pos, end_pos in coordinate_pairs:
                    if start_pos <= dna_pos <= end_pos:
                        position_in_gene = True
                        break

                if position_in_gene:
                    # Calculate amino acid position using enhanced function
                    aa_position = calculate_amino_acid_position_enhanced(dna_pos, coordinates)

                    if aa_position is not None:
                        ref_pro_seq = str(record.seq).replace(" ", "")

                        # Check if amino acid position is valid
                        if 0 <= aa_position - 1 < len(ref_pro_seq):

                            # Get the codon containing this position
                            if ref_seq is not None:
                                # Use enhanced coordinate mapping to get the correct codon
                                genome_to_joined, _ = create_genome_to_joined_mapping_enhanced(coordinates)
                                joined_position = genome_to_joined[dna_pos]

                                # Calculate codon boundaries in the joined sequence
                                codon_number = (joined_position - 1) // 3
                                codon_start_joined = codon_number * 3 + 1
                                codon_end_joined = codon_start_joined + 2

                                # Extract the original gene sequence (properly handling complement)
                                from ..utils.sequence_utils import (
                                    extract_gene_sequence_enhanced,
                                )
                                gene_sequence, was_reverse_complemented = extract_gene_sequence_enhanced(ref_seq, coordinates)

                                # Extract the codon from the gene sequence
                                if codon_start_joined - 1 < len(gene_sequence) and codon_end_joined <= len(gene_sequence):
                                    original_codon = gene_sequence[codon_start_joined - 1:codon_end_joined]

                                    # Calculate which position within the codon is mutated
                                    codon_offset = (joined_position - 1) % 3

                                    # For complement genes, the mutation affects the reverse complement
                                    if is_complement:
                                        # Convert mutation to what it would be on the reverse complement
                                        from Bio.Seq import Seq
                                        complement_map = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
                                        mutated_nt_complement = complement_map.get(mutated_nt, mutated_nt)
                                        # Reverse the codon offset for complement genes
                                        codon_offset = 2 - codon_offset
                                        mutated_codon = original_codon[:codon_offset] + mutated_nt_complement + original_codon[codon_offset + 1:]
                                    else:
                                        mutated_codon = original_codon[:codon_offset] + mutated_nt + original_codon[codon_offset + 1:]

                                    # Translate codons to amino acids
                                    from Bio.Seq import Seq
                                    original_aa_from_codon = str(Seq(original_codon).translate())
                                    mutated_aa_from_codon = str(Seq(mutated_codon).translate())

                                    # Create mutation string
                                    if original_aa_from_codon == mutated_aa_from_codon:
                                        # Silent mutation
                                        mutation_str = f"{original_aa_from_codon}{aa_position}{original_aa_from_codon}"
                                    else:
                                        # Check if this creates a stop codon
                                        if mutated_aa_from_codon == "*":
                                            # Nonsense mutation (stop codon)
                                            mutation_str = f"{original_aa_from_codon}{aa_position}*"
                                        else:
                                            # Regular missense mutation
                                            mutation_str = f"{original_aa_from_codon}{aa_position}{mutated_aa_from_codon}"

                                    return {
                                        "protein": header,
                                        "mutation": mutation_str,
                                    }

            except (ValueError, IndexError):
                # Skip records with invalid coordinate formats
                continue

        return None

    def align_orfs_to_ref_proteome_row(
        self, orfs: List[str], m_orfs: List[str], dna_pos_range: str = None
    ) -> List[Dict[str, str]]:
        """
        Align ORFs to a reference proteome and find mutations.
        Now includes coordinate checking to prevent false positives.

        Args:
            orfs (List[str]): A list of ORFs.
            m_orfs (List[str]): A list of mutated ORFs.
            dna_pos_range (str, optional): DNA position range in format "start:end" for coordinate checking.

        Returns:
            List[Dict[str, str]]: List of protein mutations.
        """
        protein_mutations = []

        # Parse DNA position range if provided
        pos_start = pos_end = None
        if dna_pos_range:
            try:
                pos_parts = dna_pos_range.split(":")
                pos_start = int(pos_parts[0])
                pos_end = int(pos_parts[-1])
            except (ValueError, IndexError):
                # If parsing fails, disable coordinate checking
                pos_start = pos_end = None

        for orf, m_orf in zip(orfs, m_orfs):
            if not orf or not m_orf:  # Skip empty ORFs
                continue

            m_orf = self.most_match_orf(orf, m_orfs)
            if not m_orf:  # Skip if no good match found
                continue

            for record in SeqIO.parse(self.proteome_dir, "fasta"):
                # Use record.description to get the full header, not record.id which gets truncated at first space
                parts = record.description.split("|")
                if len(parts) < 3:
                    continue

                header = parts[1]
                coordinates = parts[2]
                ref_pro_seq = str(record.seq).replace(" ", "")

                # First check coordinate containment if DNA position is provided
                if pos_start is not None and pos_end is not None:
                    from ..utils.sequence_utils import parse_gene_coordinates_enhanced
                    try:
                        coordinate_pairs, is_complement = parse_gene_coordinates_enhanced(coordinates)
                        position_in_gene = False
                        for start_pos, end_pos in coordinate_pairs:
                            if start_pos <= pos_start and pos_end <= end_pos:
                                position_in_gene = True
                                break

                        if not position_in_gene:
                            continue  # Skip this protein if mutation is not within its coordinates

                    except (ValueError, IndexError):
                        continue  # Skip if coordinates can't be parsed

                # Try exact match first
                matches = find_near_matches(orf, ref_pro_seq, max_l_dist=0)

                # Only allow fuzzy matching for longer proteins and ORFs to prevent false positives
                if not matches and len(ref_pro_seq) >= 100 and len(orf) >= 15:
                    matches = find_near_matches(orf, ref_pro_seq, max_l_dist=1)

                    # If still no match, try with more tolerance for longer proteins only
                    if not matches and len(ref_pro_seq) >= 200:
                        matches = find_near_matches(orf, ref_pro_seq, max_l_dist=3)

                if len(matches) == 1:
                    seq_matched = ref_pro_seq[matches[0].start : matches[0].end]
                    aligned_seqs, score = self.align_local(m_orf, seq_matched)

                    if score != 0:
                        # Again, calling the function from the GenomeAnalyzer class
                        mutation_aas = self.aa_mutation_processor.get_row_aa_mutations(
                            orf, m_orf, matches[0].start
                        )

                        # Only add if mutations were found
                        if mutation_aas:
                            protein_mutation = {
                                "protein": header,
                                "mutation": mutation_aas,
                            }
                            protein_mutations.append(protein_mutation)

        return protein_mutations

    def align_orfs_to_ref_proteome_indels(
        self, orfs: List[str], m_orfs: List[str]
    ) -> List[dict]:
        """
        Align open reading frames (ORFs) to a reference proteome.

        Args:
            orfs (List[str]): List of original ORFs.
            m_orfs (List[str]): List of mutated ORFs.
            proteins (str): Name of the reference proteome.

        Returns:
            List[dict]: List of dictionaries with matched ORFs information.
        """
        target_orfs = []

        for orf in orfs:
            for record in SeqIO.parse(self.proteome_dir, "fasta"):
                # Use record.description to get the full header, not record.id which gets truncated at first space
                header = record.description.split("|")[1]
                ref_pro_seq = str(record.seq).replace(" ", "")
                matches = find_near_matches(orf, ref_pro_seq, max_l_dist=0)

                if matches:
                    snp_r = matches[0].start
                    m_orf = self.most_match_orf(orf, m_orfs)
                    target_orfs.append(
                        {"header": header, "SNP_R": snp_r, "ORF": orf, "mORF": m_orf}
                    )
                    break

        return target_orfs


class ProMutationDetector:
    def __init__(
        self,
        orf_processor: ORFProcessor,
        alignment_processor: AlignmentProcessor,
        aa_mutation_processor: AminoAcidMutationProcessor,
        ref_genome_instance: ReferenceGenome,
    ) -> None:
        self.codon_table = CodonTable.unambiguous_dna_by_id[1]
        self.orf_processor = orf_processor
        self.alignment_processor = alignment_processor
        self.aa_mutation_processor = aa_mutation_processor
        self.ref_seq = ref_genome_instance.ref_seq

    def protein_point_snp(self, genome_snp: Dict) -> Dict:
        """
        Get protein mutations from a Single Nucleotide Polymorphism (SNP).

        Args:
            genome_snp (Dict): A SNP dictionary in the genome.

        Returns:
            Dict: Updated SNP dictionary with protein mutation information.
        """
        snp_pos = genome_snp["pos"]
        mutated_to = genome_snp["SNP"][3]

        orfs, m_orfs = self.orf_processor.get_orfs_snps(snp_pos, mutated_to)
        protein_mutation = self.alignment_processor.align_orfs_to_ref_proteome_snps(
            orfs, m_orfs, dna_pos=int(snp_pos), mutation_info=genome_snp["SNP"], ref_seq=self.ref_seq
        )
        # Process protein mutations
        genome_snp.update({"proteinMutation": protein_mutation})
        return genome_snp

    def protein_row_snp(self, genome_snp: Dict) -> Dict:
        """
        Get protein mutations from a Single Nucleotide Polymorphism (SNP) for row-based data.

        Args:
            genome_snp (Dict): A SNP dictionary in the genome.

        Returns:
            Dict: Updated SNP dictionary with protein mutation information.
        """
        pos_range = genome_snp["pos"].split(":")
        pos_a, pos_b = int(pos_range[0]) - 1, int(pos_range[-1])
        mutated_to = genome_snp["SNP"].split("->")[1]

        seq_t = self.ref_seq[pos_a - 2 : pos_b + 25]
        orfs = self.orf_processor.translate_dna_3_frames_row(seq_t)
        seq_m = seq_t[:2] + mutated_to + seq_t[2 + len(mutated_to) :]
        m_orfs = self.orf_processor.translate_dna_3_frames_row(seq_m)

        protein_mutations = self.alignment_processor.align_orfs_to_ref_proteome_row(
            orfs, m_orfs, dna_pos_range=genome_snp["pos"]
        )

        genome_snp["proteinMutation"] = protein_mutations

        return genome_snp

    def protein_point_row_snps(self, genome_snps: List[Dict]) -> List[Dict]:
        """
        Get protein point mutations from a list of Single Nucleotide Polymorphisms (SNPs).

        Args:
            genome_snps (List[Dict]): A list of SNP dictionaries in the genome.

        Returns:
            List[Dict]: List of SNP dictionaries with protein mutation information.
        """
        processed_snps = []
        for genome_snp in genome_snps:
            processed_snp = (
                self.protein_row_snp(genome_snp)
                if genome_snp["type"] == "rowMutation"
                else self.protein_point_snp(genome_snp)
            )
            processed_snps.append(processed_snp)
        # print('=============================================')
        # print('processed_snps', processed_snps)

        return processed_snps

    def protein_dels_snps(self, genome_snps: List[Dict]) -> List[Dict]:
        """
        Get protein deletions and SNPs from genome SNPs.

        Args:
            genome_snps (List[Dict]): List of SNP dictionaries in the genome.

        Returns:
            List[Dict]: List of updated SNP dictionaries with protein mutation information.
        """
        snp_results = []
        for genome_snp in genome_snps:
            del_seq = genome_snp["SNP"]
            frameshift = "yes" if len(del_seq) % 3 != 0 else "no"

            orfs, m_orfs = self.orf_processor.get_orfs_dels(genome_snp["pos"])

            if not orfs or not m_orfs:
                header = (
                    "Invalid protein sequence"
                    if not orfs
                    else "mutated to invalid protein sequence"
                )
                protein_mutations = [
                    {"protein": header, "frameshift": frameshift, "mutation": "NA"}
                ]
            else:
                target_orfs = (
                    self.alignment_processor.align_orfs_to_ref_proteome_indels(
                        orfs, m_orfs
                    )
                )
                protein_mutations = self._process_target_orfs(
                    target_orfs, del_seq, frameshift, mut_type="del"
                )

            genome_snp["proteinMutation"] = protein_mutations
            snp_results.append(genome_snp)

        return snp_results

    def protein_ins_snps(self, genome_snps: List[Dict]) -> List[Dict]:
        """
        Get protein mutations from a list of insertion SNPs.

        Args:
            genome_snps (List[Dict]): List of SNP dictionaries in the genome.

        Returns:
            List[Dict]: List of updated SNP dictionaries with protein mutation information.
        """
        snp_results = []
        for genome_snp in genome_snps:
            ins_seq = genome_snp["SNP"]
            frameshift = "no" if len(ins_seq) % 3 == 0 else "yes"

            orfs, m_orfs = self.orf_processor.get_orfs_ins(
                int(genome_snp["pos"]) + 1, ins_seq
            )
            target_orfs = self.alignment_processor.align_orfs_to_ref_proteome_indels(
                orfs, m_orfs
            )

            protein_mutations = self._process_target_orfs(
                target_orfs, ins_seq, frameshift, mut_type="ins"
            )

            genome_snp["proteinMutation"] = protein_mutations
            snp_results.append(genome_snp)

        return snp_results

    def protein_hot_snps(self, genome_snps: List[Dict]) -> List[Dict]:
        """
        Get protein hot SNPs from genome SNPs using codon-aware mutation mapping.

        Args:
            genome_snps (List[Dict]): List of SNP dictionaries in the genome.

        Returns:
            List[Dict]: List of updated SNP dictionaries with protein mutation information.
        """
        from Bio.Seq import Seq
        snp_results = []

        for genome_snp in genome_snps:
            # Parse mutation information
            pos_range = genome_snp["pos"].split(":")
            pos_start = int(pos_range[0])
            pos_end = int(pos_range[-1])
            orig_bases, mut_bases = genome_snp["SNP"].split("->")

            # Map mutations to their correct positions
            # The mutation notation "start:end,orig->mut" means specific positions mutate
            # e.g., "10447:10449,GC->AA" means 10447G->A and 10449C->A (10448 unchanged)
            mutations = []
            mutation_index = 0  # index for orig_bases/mut_bases
            for pos in range(pos_start, pos_end + 1):
                if mutation_index < len(orig_bases):
                    # Check if this position should have a mutation
                    # For "10447:10449,GC->AA", we want positions 10447 and 10449
                    # So we need to map mutation_index 0 to pos_start, mutation_index 1 to pos_end
                    if mutation_index == 0 and pos == pos_start:
                        orig_base = orig_bases[mutation_index]
                        mut_base = mut_bases[mutation_index]
                        if orig_base != mut_base:
                            mutations.append((pos, orig_base, mut_base))
                        mutation_index += 1
                    elif mutation_index == 1 and pos == pos_end:
                        orig_base = orig_bases[mutation_index]
                        mut_base = mut_bases[mutation_index]
                        if orig_base != mut_base:
                            mutations.append((pos, orig_base, mut_base))
                        mutation_index += 1

            # Process each protein that could be affected
            protein_mutations = []
            for record in SeqIO.parse(self.alignment_processor.proteome.proteome_dir, "fasta"):
                # Use record.description to get the full header, not record.id which gets truncated at first space
                parts = record.description.split("|")

                # Handle FASTA header format: ID|protein_name|coordinates
                if len(parts) >= 3:
                    # Format: YP_009725297.1|leader_nsp1|266..805
                    # Format: NP_066243.1|nucleoprotein|470..2689
                    # Format: YP_009725307.1|RNA-dependent-polymerase|join(13442..13468,13468..16236)
                    header = parts[1]
                    coordinates = parts[2]
                else:
                    continue  # Skip if format is not recognized

                # Import join operation utilities
                from ..utils.sequence_utils import parse_gene_coordinates_enhanced

                try:
                    # Parse coordinates (handles simple ranges, join operations, and complement notation)
                    coordinate_pairs, is_complement = parse_gene_coordinates_enhanced(coordinates)

                    # For hot mutations, check if the mutation range is CONTAINED within the protein range
                    position_in_gene = False
                    for start_pos, end_pos in coordinate_pairs:
                        if start_pos <= pos_start and pos_end <= end_pos:
                            position_in_gene = True
                            protein_start, protein_end = start_pos, end_pos
                            break

                    if not position_in_gene:
                        continue

                except (ValueError, IndexError):
                    # Skip records with invalid coordinate formats
                    continue

                # Redundant check removed - position_in_gene already ensures containment

                # Find affected codons
                affected_codons = {}
                for pos, orig_base, mut_base in mutations:
                    if protein_start <= pos <= protein_end:
                        # Calculate which codon this position belongs to
                        offset_from_start = pos - protein_start
                        codon_index = offset_from_start // 3
                        codon_start = protein_start + (codon_index * 3)

                        if codon_start not in affected_codons:
                            # Get reference codon
                            ref_codon = self.ref_seq[codon_start - 1:codon_start + 2]
                            affected_codons[codon_start] = list(ref_codon)

                        # Apply mutation to codon
                        codon_pos = pos - codon_start
                        affected_codons[codon_start][codon_pos] = mut_base

                # Generate protein mutations for each affected codon
                for codon_start, mutated_codon_list in affected_codons.items():
                    ref_codon = self.ref_seq[codon_start - 1:codon_start + 2]
                    mutated_codon = ''.join(mutated_codon_list)

                    # Calculate protein position
                    protein_pos = ((codon_start - protein_start) // 3) + 1

                    # Translate codons
                    ref_aa = str(Seq(ref_codon).translate())
                    mut_aa = str(Seq(mutated_codon).translate())

                    # Create mutation string
                    if ref_aa == mut_aa:
                        mutation = f"{ref_aa}{protein_pos}{ref_aa}"
                    else:
                        mutation = f"{ref_aa}{protein_pos}{mut_aa}"

                    # Check if this protein is already in mutations
                    existing_mutation = None
                    for pm in protein_mutations:
                        if pm["protein"] == header:
                            existing_mutation = pm
                            break

                    if existing_mutation:
                        # Append to existing mutation string
                        existing_mutation["mutation"] += f",{mutation}"
                    else:
                        # Create new protein mutation
                        protein_mutations.append({
                            "protein": header,
                            "mutation": mutation
                        })

            genome_snp["proteinMutation"] = protein_mutations
            snp_results.append(genome_snp)

        return snp_results

    def _process_target_orfs(self, target_orfs, seq, frameshift, mut_type):
        protein_mutations = []
        for target_orf in target_orfs:
            header = target_orf["header"]
            orf = target_orf["ORF"]
            m_orf = target_orf["mORF"]
            snp_r = target_orf["SNP_R"]
            aligned_seqs, score = self.alignment_processor.align_local(orf, m_orf)

            if mut_type == "del":
                aa_mutations = self.aa_mutation_processor.get_aa_dels(orf, m_orf, snp_r)
            elif mut_type == "ins":
                aa_mutations = self.aa_mutation_processor.get_aa_ins(orf, m_orf, snp_r)

            if (
                len(seq) > 21
            ):  # Handle large insertion/deletion (invalid) >7 amino acids
                protein_mutation = {
                    "protein": "None-CDS",
                    "frameshift": frameshift,
                    "mutation": f"long-{mut_type}",
                }
            else:
                protein_mutation = {
                    "protein": header,
                    "frameshift": frameshift,
                    "mutation": aa_mutations,
                }

            protein_mutations.append(protein_mutation)

        return protein_mutations
