# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-06-23 14:27:07
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-19 17:42:40
Email: wang.rui@nyu.edu
FilePath: /VARIANT/src/utils/__init__.py
Description: Python module for __init__.
'''

from .list_utils import (
    contains_element,
    delete_char_from_string,
    delete_element_from_list,
    is_keys_subset,
)
from .mutation_utils import (
    create_mutation_dict,
    find_common_elements,
    generate_hash_key,
    is_subset,
    sort_dict_by_consecutive_keys,
)
from .sequence_utils import clustal_genomes, has_only_valid_nts, read_fasta, write_fasta

__all__ = [
    # Sequence utilities
    "read_fasta",
    "write_fasta",
    "clustal_genomes",
    "has_only_valid_nts",
    # Mutation utilities
    "create_mutation_dict",
    "generate_hash_key",
    "is_subset",
    "find_common_elements",
    "sort_dict_by_consecutive_keys",
    # List utilities
    "delete_element_from_list",
    "delete_char_from_string",
    "is_keys_subset",
    "contains_element",

]
