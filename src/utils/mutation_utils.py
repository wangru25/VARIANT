# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-08-19 17:30:00
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-20 09:47:20
Email: wang.rui@nyu.edu
FilePath: /VARIANT/src/utils/mutation_utils.py
Description: Mutation analysis utilities and helper functions.
'''

import hashlib
from typing import Any, Dict, List, Tuple


def generate_hash_key(mut_data: str) -> int:
    '''
    Generate a hash key for mutation data.

    Args:
        mut_data (str): Mutation data string.

    Returns:
        int: Hash key.
    '''
    return int(hashlib.md5(mut_data.encode()).hexdigest(), 16)


def create_mutation_dict(pos_list: List[int], snp_list: List[str], mutation_type: str, hash_function: callable = generate_hash_key) -> Dict[str, Any]:
    '''
    Create a standardized mutation dictionary.

    Args:
        pos_list: List of positions affected by the mutation
        snp_list: List of SNP changes (format: "from->to")
        mutation_type: Type of mutation (e.g., "insertion", "deletion")
        hash_function: Function to generate hash keys

    Returns:
        Dictionary containing mutation information

    Example:
        >>> create_mutation_dict([100, 101], ["A->G", "T->C"], "point_mutation")
        {'key': 12345678, 'type': 'point_mutation', 'pos': '100:101', 'SNP': 'AT->GC'}
    '''
    # Create position string
    if len(pos_list) > 1:
        pos_str = f"{pos_list[0]}:{pos_list[-1]}"
    else:
        pos_str = str(pos_list[0])

    # Create SNP string
    snp_from = "".join(snp.split("->")[0] for snp in snp_list)
    snp_to = "".join(snp.split("->")[1] for snp in snp_list)
    snp_str = f"{snp_from}->{snp_to}"

    # Create mutation data string for hashing
    mutation_data = f"{pos_str}{snp_str}"

    return {
        "key": hash_function(mutation_data),
        "type": mutation_type,
        "pos": pos_str,
        "SNP": snp_str,
    }


def is_subset(snp_keys_current: List[int], snp_keys_next: List[int]) -> bool:
    '''
    Check if the first list is a subset of the second list.

    Args:
        snp_keys_current: List of current SNP keys
        snp_keys_next: List of next SNP keys

    Returns:
        True if current keys are a subset of next keys

    Example:
        >>> is_subset([10, 5], [9, 4, 5, 8, 10])
        True
    '''
    return set(snp_keys_current).issubset(set(snp_keys_next))


def find_common_elements(lists: List[List[int]]) -> List[int]:
    '''
    Find common elements across multiple lists.

    Args:
        lists (List[List[int]]): List of lists to find common elements from.

    Returns:
        List[int]: List of common elements.
    '''
    if not lists:
        return []

    common = set(lists[0])
    for lst in lists[1:]:
        common &= set(lst)

    return list(common)


def sort_dict_by_consecutive_keys(dictionary: Dict[int, str]) -> Dict[int, str]:
    '''
    Sort dictionary by consecutive keys.

    Args:
        dictionary (Dict[int, str]): Dictionary to sort.

    Returns:
        Dict[int, str]: Sorted dictionary.
    '''
    sorted_keys = sorted(dictionary.keys())
    return {key: dictionary[key] for key in sorted_keys}


def is_keys_subset(sub_list: List[Any], test_list: List[Any]) -> Tuple[bool, List[Any]]:
    '''
    Check if the first list is a proper subset of the second list.

    Args:
        sub_list: List to check if it's a subset
        test_list: List to check against

    Returns:
        Tuple of (is_subset, difference_list)

    Example:
        >>> is_keys_subset([1, 2], [1, 2, 3, 4])
        (True, [3, 4])
    '''
    is_subset_result = set(sub_list).issubset(set(test_list)) and len(test_list) > len(
        sub_list
    )

    if is_subset_result:
        difference = [x for x in test_list if x not in sub_list]
    else:
        difference = []

    return is_subset_result, difference


def convert_key_to_int(key: str) -> int:
    '''
    Convert a string key to an integer.
    Handles position ranges like "11288:11296" by using the start position.

    Args:
        key (str): String key.

    Returns:
        int: Integer key.
    '''
    try:
        # Handle position ranges like "11288:11296"
        if ':' in key:
            start_pos = key.split(':')[0]
            return int(start_pos)
        else:
            return int(key)
    except ValueError:
        return 0


def is_next_genome(keys0: List[int], keys1: List[int]) -> bool:
    '''
    Check if keys1 is a subset of keys0.

    Args:
        keys0 (List[int]): First list of keys.
        keys1 (List[int]): Second list of keys.

    Returns:
        bool: True if keys1 is a subset of keys0.
    '''
    return all(key in keys0 for key in keys1)


def contains_element(test_string: str, test_list: List[str]) -> bool:
    '''
    Check if a string contains any element from a list.

    Args:
        test_string: String to search in
        test_list: List of strings to search for

    Returns:
        True if any element from test_list is found in test_string
    '''
    return any(element in test_string for element in test_list)


def validate_mutation_data(mutation_dict: Dict[str, Any]) -> Tuple[bool, List[str]]:
    '''
    Validate mutation dictionary data.

    Args:
        mutation_dict: Mutation dictionary to validate

    Returns:
        Tuple of (is_valid, list_of_errors)
    '''
    errors = []
    required_keys = ["key", "type", "pos", "SNP"]

    # Check required keys
    for key in required_keys:
        if key not in mutation_dict:
            errors.append(f"Missing required key: {key}")

    # Validate mutation type
    valid_types = {
        "insertion",
        "deletion",
        "point_mutation",
        "row_mutation",
        "hot_mutation",
    }
    if "type" in mutation_dict and mutation_dict["type"] not in valid_types:
        errors.append(f"Invalid mutation type: {mutation_dict['type']}")

    # Validate SNP format
    if "SNP" in mutation_dict:
        snp = mutation_dict["SNP"]
        if "->" not in snp:
            errors.append(f"Invalid SNP format: {snp}")

    return len(errors) == 0, errors


def classify_mutation_type(protein_mutation: List[Dict]) -> str:
    '''
    Classify mutation by biological type based on protein mutation information.

    Args:
        protein_mutation (List[Dict]): List of protein mutation dictionaries.

    Returns:
        str: Biological classification (silent, missense, nonsense, deletion, insertion, frameshift)
    '''
    if not protein_mutation:
        return "unknown"

    # Check for frameshift mutations first
    for mutation in protein_mutation:
        if isinstance(mutation, dict) and mutation.get('frameshift') == 'yes':
            return "frameshift"

    # Check for deletions
    for mutation in protein_mutation:
        if isinstance(mutation, dict) and mutation.get('mutation'):
            mut_str = mutation.get('mutation')
            if isinstance(mut_str, list):
                # Check if any mutation in the list is a deletion
                if any('del' in str(m) for m in mut_str):
                    return "deletion"
            elif isinstance(mut_str, str) and 'del' in mut_str:
                return "deletion"

    # Check for insertions
    for mutation in protein_mutation:
        if isinstance(mutation, dict) and mutation.get('mutation'):
            mut_str = mutation.get('mutation')
            if isinstance(mut_str, list):
                # Check if any mutation in the list is an insertion
                if any('ins' in str(m) for m in mut_str):
                    return "insertion"
            elif isinstance(mut_str, str) and 'ins' in mut_str:
                return "insertion"

    # Check for nonsense mutations (stop codons)
    for mutation in protein_mutation:
        if isinstance(mutation, dict) and mutation.get('mutation'):
            mut_str = mutation.get('mutation')
            if isinstance(mut_str, list):
                # Check if any mutation in the list creates a stop codon
                if any('*' in str(m) or 'X' in str(m) for m in mut_str):
                    return "nonsense"
            elif isinstance(mut_str, str) and ('*' in mut_str or 'X' in mut_str):
                return "nonsense"

    # Check for silent vs missense mutations
    for mutation in protein_mutation:
        if isinstance(mutation, dict) and mutation.get('mutation'):
            mut_str = mutation.get('mutation')
            if isinstance(mut_str, list):
                # Check if all mutations are silent (same amino acid)
                all_silent = True
                for m in mut_str:
                    if isinstance(m, str) and len(m) >= 3:
                        # Check if it's a silent mutation (e.g., "A123A")
                        parts = m.split('del')[0] if 'del' in m else m
                        if len(parts) >= 3 and parts[0] == parts[-1] and parts[1:-1].isdigit():
                            continue
                        else:
                            all_silent = False
                            break
                if all_silent:
                    return "silent"
            elif isinstance(mut_str, str) and len(mut_str) >= 3:
                # Check if it's a silent mutation (e.g., "A123A")
                parts = mut_str.split('del')[0] if 'del' in mut_str else mut_str
                if len(parts) >= 3 and parts[0] == parts[-1] and parts[1:-1].isdigit():
                    return "silent"

    # If not silent, it's missense
    return "missense"


def detect_hot_mutation(original_seq: str, mutated_seq: str) -> bool:
    '''
    Detect if a mutation is a hot mutation (exactly 2 substitutions flanking exactly 1 conserved base).

    Args:
        original_seq (str): Original nucleotide sequence.
        mutated_seq (str): Mutated nucleotide sequence.

    Returns:
        bool: True if it's a hot mutation, False otherwise.
    '''
    if len(original_seq) != len(mutated_seq):
        return False

    # Find all positions where nucleotides differ
    substitutions = []
    for i, (orig, mut) in enumerate(zip(original_seq, mutated_seq)):
        if orig != mut:
            substitutions.append(i)

    # Hot mutation requires exactly 2 substitutions
    if len(substitutions) != 2:
        return False

    # Check if there's exactly 1 conserved base between the two substitutions
    pos1, pos2 = substitutions[0], substitutions[1]

    # For hot mutation, the two substitutions should be at positions 0 and 2 (flanking position 1)
    # This means the sequence should be 3 nucleotides long with the middle one conserved
    if len(original_seq) == 3 and pos1 == 0 and pos2 == 2:
        # Check if the middle position (position 1) is conserved
        return original_seq[1] == mutated_seq[1]

    # For other cases, count conserved bases between the two substitution positions
    conserved_count = 0
    for i in range(pos1 + 1, pos2):
        if original_seq[i] == mutated_seq[i]:
            conserved_count += 1

    # Hot mutation has exactly 1 conserved base between substitutions
    return conserved_count == 1


def split_multi_protein_mutations(genome_snpp: Dict) -> List[Dict]:
    '''
    Split mutations affecting multiple proteins into separate mutation records.
    Note: Hot mutations and row mutations are not split by protein as they represent single genomic events.

    Args:
        genome_snpp (Dict): Original mutation record.

    Returns:
        List[Dict]: List of separated mutation records.
    '''
    protein_mutation = genome_snpp.get('proteinMutation', [])
    mutation_type = genome_snpp.get('type', '')

    # For hot mutations and row mutations, don't split by protein - treat as single genomic events
    if mutation_type in ['hotMutation', 'rowMutation']:
        return [genome_snpp]

    # If there's only one protein affected, return as is
    if len(protein_mutation) <= 1:
        return [genome_snpp]

    # Split mutations affecting multiple proteins
    separated_mutations = []
    for protein_info in protein_mutation:
        if isinstance(protein_info, dict):
            separated_mutation = genome_snpp.copy()
            separated_mutation['proteinMutation'] = [protein_info]
            separated_mutations.append(separated_mutation)

    return separated_mutations if separated_mutations else [genome_snpp]


def format_amino_acid_change_for_csv(amino_acid_change) -> str:
    '''
    Format amino acid change for CSV output as a quoted list string, even for single mutations.

    Args:
        amino_acid_change: Amino acid change (can be string, list, or other format)

    Returns:
        str: Formatted amino acid change for CSV (always a quoted list string like "['S375F', 'T376A']")
    '''
    if isinstance(amino_acid_change, list):
        return f'"{str(amino_acid_change)}"'
    elif isinstance(amino_acid_change, str):
        # If it's a comma-separated string, convert to list
        if ',' in amino_acid_change:
            items = [item.strip() for item in amino_acid_change.split(',')]
            return f'"{str(items)}"'
        # If it's a single mutation string, wrap in a list
        return f'"{str([amino_acid_change])}"'
    else:
        return f'"{str([str(amino_acid_change)])}"'
