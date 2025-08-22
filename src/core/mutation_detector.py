# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2023-10-31 12:42:13
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-20 09:50:14
Email: wang.rui@nyu.edu
FilePath: /VARIANT/src/core/mutation_detector.py
Description: Mutation detection and analysis module for viral genomes.
'''



from typing import Any, Dict, List, Tuple, Union

import more_itertools as mit
from Bio import AlignIO

from ..utils import sequence_utils as utils
from ..utils.mutation_utils import generate_hash_key


class ReferenceGenome:
    def __init__(self, ref_dir: str):
        fasta_data = utils.read_fasta(ref_dir)[0]
        self.ref_name = fasta_data["seqId"]
        self.ref_seq = fasta_data["seq"]


class MultipleSequenceAlignment:
    def __init__(self, msa_dir: str):
        self.align = AlignIO.read(msa_dir, "clustal")

    def find_ref_msa(self, ref_name: str, virus_config: Dict = None) -> Tuple[int, str]:
        """
        Find the index of the reference genome in the MSA file.

        This method uses a hierarchical approach:
        1. First tries the exact ref_name
        2. Then tries reference names from virus_config
        3. Finally falls back to generic reference patterns

        Args:
            ref_name (str): The name of the reference genome
            virus_config (Dict, optional): Virus configuration containing reference genome info

        Returns:
            Tuple[int, str]: A tuple containing the index of the reference genome and its sequence
        """
        # First try exact match
        ref_matches = [i for i, seq in enumerate(self.align) if ref_name in seq.id]

        if not ref_matches and virus_config:
            # Try reference names from virus config
            ref_names = self._get_reference_names_from_config(virus_config)
            for common_name in ref_names:
                ref_matches = [
                    i for i, seq in enumerate(self.align) if common_name in seq.id
                ]
                if ref_matches:
                    ref_name = common_name
                    break

        if not ref_matches:
            # Try generic reference patterns
            generic_ref_names = ["reference", "ref", "REF", "Reference"]
            for common_name in generic_ref_names:
                ref_matches = [i for i, seq in enumerate(self.align) if common_name in seq.id]
                if ref_matches:
                    ref_name = common_name
                    break

        if not ref_matches:
            # Last resort: use the first sequence as reference
            print(f"Warning: Reference genome '{ref_name}' not found in MSA file. Using first sequence as reference.")
            ref_seq_id = 0
            ref_seq_msa = str(self.align[ref_seq_id].seq)
            return ref_seq_id, ref_seq_msa

        ref_seq_id = ref_matches[0]
        ref_seq_msa = str(self.align[ref_seq_id].seq)
        return ref_seq_id, ref_seq_msa

    def _get_reference_names_from_config(self, virus_config: Dict) -> List[str]:
        """
        Get reference names from virus configuration.
        
        Args:
            virus_config (Dict): Virus configuration dictionary
            
        Returns:
            List[str]: List of reference names to try
        """
        ref_names = []

        # Get reference genome filename from config and extract the reference name
        if 'reference_genome' in virus_config:
            ref_file = virus_config['reference_genome']
            # Extract reference name (e.g., NC_007366.1 from NC_007366.1.fasta)
            if ref_file.endswith('.fasta'):
                ref_name = ref_file[:-6]  # Remove .fasta
            else:
                ref_name = ref_file
            ref_names.append(ref_name)

            # Print helpful message for users
            print(f"Looking for reference genome: {ref_name}")
            print(f"Make sure your MSA file contains a sequence with ID containing '{ref_name}'")

        # Add any additional reference names specified in config
        if 'reference_names' in virus_config:
            ref_names.extend(virus_config['reference_names'])

        return ref_names

    def find_seq_msa(self, genome_id: str) -> Tuple[int, str]:
        """
        This function finds the index of a genome in the MSA file.

        Parameters:
        genome_id (str): The ID of the genome to find in the MSA file.

        Returns:
        tuple: A tuple containing the index of the genome and its sequence in the MSA file: [idx, self.seq_msa]
        """
        for i, seq in enumerate(self.align):
            if genome_id in seq.id:
                return i, str(seq.seq)
        print(f"Error: Genome ID {genome_id} not found in the MSA file.")
        return None, None

    def msa_pos_to_ref_pos(self, pos_msa: int, ref_seq_msa: str) -> int:
        """
        Maps the positions on MSA reference sequence to the position on actual reference genome.

        Parameters:
        pos_msa (int): The position on the MSA reference sequence in the MSA file.
        ref_seq_msa (str): The reference sequence in the MSA file.

        Returns:
        int: The corresponding position on the actual reference genome.
        """
        num_dash = ref_seq_msa[:pos_msa].count("-")
        pos_ref = pos_msa - num_dash
        return pos_ref


class GeneMutationDetector:
    def __init__(
        self,
        ref_genome: ReferenceGenome,
        msa: MultipleSequenceAlignment,
        genome_id: str,
        virus_config: Dict = None,
    ):
        self.genome_id = genome_id
        self.ref_seq = ref_genome
        self.msa = msa
        self.virus_config = virus_config
        self.ref_seq_id, self.ref_seq_msa = self.msa.find_ref_msa(self.ref_seq.ref_name, virus_config)
        self.seq_id, self.seq_msa = self.msa.find_seq_msa(self.genome_id)

    @staticmethod
    def create_mutation_dict(
        mutation_type: str, pos: Union[int, List[int]], snp: Union[str, List[str]]
    ) -> Dict[str, Any]:
        """
        Create a dictionary representing a mutation.

        Args:
            mutation_type: The type of mutation (e.g. 'pointMutation', 'rowMutation', etc.).
            pos: The position(s) of the mutation.
            snp: The nucleotide change(s) of the mutation.

        Returns:
            A dictionary representing the mutation.

        Raises:
            TypeError: If `pos` or `snp` is not a list or a single value.

        """
        if isinstance(pos, list):
            pos_str = f"{pos[0]}:{pos[-1]}"
            snp_from = "".join(s.split("->")[0] for s in snp)
            snp_to = "".join(s.split("->")[1] for s in snp)
            snp_str = f"{snp_from}->{snp_to}"
        else:
            pos_str = str(pos)
            snp_str = snp

        mut_data = f"{pos_str}{snp_str}"
        return {
            "key": generate_hash_key(mut_data),
            "type": mutation_type,
            "pos": pos_str,
            "SNP": snp_str,
        }

    def insertions_snps(self) -> List[Tuple[str, int]]:
        """
        Extracts all insertion sequences and their indices in the reference genome from two aligned sequences.

        Returns:
        list: A list of tuples containing an insertion sequence and its index in the reference genome.
        """
        pos_start = self.ref_seq_msa.find("-")
        pos_end = self.ref_seq_msa.rfind("-")

        # Extract the insertion sequences and their indices
        ins_seqs = []
        ins_idxs = []
        ins_seq = ""
        ins_idx = None
        for i in range(pos_start, pos_end + 1):
            if self.ref_seq_msa[i] == "-" and self.seq_msa[i] != "-":
                ins_seq += self.seq_msa[i]
                if ins_idx is None:
                    ins_idx = i
            else:
                if ins_seq:
                    ins_seqs.append(ins_seq)
                    ins_idxs.append(ins_idx)
                    ins_seq = ""
                    ins_idx = None
        if ins_seq:
            ins_seqs.append(ins_seq)
            ins_idxs.append(ins_idx)
        return list(zip(ins_seqs, ins_idxs))

    def get_insertions(self) -> List[Dict[str, Any]]:
        """
        Find insertions in the query sequence relative to the reference sequence.

        Returns:
            A list of dictionaries representing the insertions in the query sequence.
        """
        insert_mutations = []
        for ins_seq, pos_insert in self.insertions_snps():
            pos_ref = self.msa.msa_pos_to_ref_pos(pos_insert, self.ref_seq_msa)
            insert_mutations.append(
                self.create_mutation_dict("insertion", pos_ref, ins_seq)
            )
        return insert_mutations

    def get_deletions(self) -> List[Dict[str, Any]]:
        """
        Find row mutations in the query sequence relative to the reference sequence.

        Returns:
            A tuple of two lists of dictionaries representing the row mutations in the query sequence.
            The first list contains all row mutations, and the second list contains only consecutive row mutations.
        """
        mutations = []
        del_seq = ""
        prev_base = ""
        for i, (ref_base, seq_base) in enumerate(zip(self.ref_seq_msa, self.seq_msa)):
            if ref_base != "-" and seq_base == "-":
                del_seq += ref_base
            else:
                if del_seq:
                    pos_msa = i - len(del_seq)
                    pos_ref = self.msa.msa_pos_to_ref_pos(pos_msa, self.ref_seq_msa)
                    end_ref = pos_ref + len(del_seq)
                    if prev_base != "N" and self.seq_msa[i : i + 1] != "N":
                        mutations.append(
                            self.create_mutation_dict(
                                "deletion", f"{pos_ref + 1}:{end_ref}", del_seq
                            )
                        )
                    del_seq = ""
            prev_base = seq_base
        return mutations

    def get_others_mutations(
        self,
    ) -> Tuple[
        List[Dict[str, Union[str, int]]],
        List[Dict[str, Union[str, int]]],
        List[Dict[str, Union[str, int]]],
    ]:
        """
        Find row, hot, and point mutations in the query sequence relative to the reference sequence.

        Returns:
            A list of dictionaries representing the point mutations in the query sequence.
        """
        row_mutations = []
        hot_mutations = []
        point_mutations = []

        not_row_mutations_dict = {}  # not row but point ref: C--ATCTA seq: GTTCTCTA CT->GC is not a row mutation but point mutation
        other_mutations_dict = {}  # (hot mutations + point mutations + row mutations) - not row mutations
        for pos_msa, (ref_base, seq_base) in enumerate(
            zip(self.ref_seq_msa, self.seq_msa)
        ):
            if ref_base in "ATGC" and seq_base in "ATGC" and ref_base != seq_base:
                pos_ref = self.msa.msa_pos_to_ref_pos(pos_msa, self.ref_seq_msa)
                pos_genome = pos_ref + 1
                mut_nt = f"{ref_base}->{seq_base}"
                # Check bounds before accessing adjacent positions
                ref_prev = self.ref_seq_msa[pos_msa - 1] if pos_msa > 0 else "-"
                ref_next = self.ref_seq_msa[pos_msa + 1] if pos_msa < len(self.ref_seq_msa) - 1 else "-"
                seq_prev = self.seq_msa[pos_msa - 1] if pos_msa > 0 else "-"
                seq_next = self.seq_msa[pos_msa + 1] if pos_msa < len(self.seq_msa) - 1 else "-"

                if (
                    ref_prev != "-"
                    and ref_next != "-"
                    and seq_prev != "-"
                    and seq_next != "-"
                ):
                    other_mutations_dict[pos_genome] = mut_nt
                else:
                    not_row_mutations_dict[pos_genome] = mut_nt
        point_mutations_dict = {**other_mutations_dict, **not_row_mutations_dict}

        # ======================================= Row Mutations ======================================
        # Group consecutive SNPs into row mutations
        cons = [
            list(group)
            for group in mit.consecutive_groups(sorted(other_mutations_dict.keys()))
        ]
        row_poss = [x for x in cons if len(x) > 1]
        row_snps = [[other_mutations_dict[y] for y in x] for x in row_poss]
        for idx, _ in enumerate(row_poss):
            row_mutations.append(
                self.create_mutation_dict("rowMutation", row_poss[idx], row_snps[idx])
            )

        # ======================================= Hot Mutations ======================================
        # Find mutations that are 2 positions apart (flanking 1 conserved base)
        # But exclude positions that are already part of row mutations
        hot_poss = []
        row_positions = set()
        for row_pos in row_poss:
            row_positions.update(row_pos)

        for key in sorted(other_mutations_dict.keys()):
            # Check if there's a mutation 2 positions away
            if key + 2 in other_mutations_dict:
                # Only create hot mutation if neither position is part of a row mutation
                if key not in row_positions and key + 2 not in row_positions:
                    hot_poss.append([key, key + 2])

        # Remove duplicates and ensure we don't double-count
        unique_hot_poss = []
        for pos_pair in hot_poss:
            if pos_pair not in unique_hot_poss:
                unique_hot_poss.append(pos_pair)

        hot_snps = [
            (other_mutations_dict[key1], other_mutations_dict[key2])
            for key1, key2 in unique_hot_poss
        ]
        for idx, _ in enumerate(unique_hot_poss):
            hot_mutations.append(
                self.create_mutation_dict("hotMutation", unique_hot_poss[idx], hot_snps[idx])
            )

        # ======================================= Point Mutations ======================================
        point_poss = [
            key
            for key in point_mutations_dict
            if key not in [x for y in row_poss for x in y]
            and key not in [x for y in hot_poss for x in y]
        ]
        point_snps = [point_mutations_dict[key] for key in point_poss]
        for idx, _ in enumerate(point_poss):
            point_mutations.append(
                self.create_mutation_dict(
                    "pointMutation", point_poss[idx], point_snps[idx]
                )
            )

        return row_mutations, hot_mutations, point_mutations


# if __name__ == "__main__":
#     genome_id = "EPI_ISL_16327572"
#     ref_genome = ReferenceGenome("../data/refs/NC_045512.fasta")
#     seq_msa = MultipleSequenceAlignment(
#         "../data/clustalW/gisaid_hcov-19_20221201_20221230_China_0_msa.txt"
#     )
#     mutation_finder = GeneMutationDetector(ref_genome, seq_msa, genome_id)

#     print("==============================Insertions==============================")
#     insertions = mutation_finder.get_insertions()
#     print(insertions)

#     print("==============================Deletions==============================")
#     deletions = mutation_finder.get_deletions()
#     print(deletions)

#     print("==============================Row Mutations==============================")
#     (
#         row_mutations,
#         hot_mutations,
#         point_mutations,
#     ) = mutation_finder.get_others_mutations()
#     print(row_mutations)

#     print("==============================Hot Mutations==============================")
#     print(hot_mutations)

#     print("==============================Point Mutations==============================")
#     print(len(point_mutations))
#     # print(point_mutations)
