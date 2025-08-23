# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-08-19 17:30:00
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-23 11:05:49
Email: rw3594@nyu.edu
FilePath: /VARIANT/src/__init__.py
Description: Python module for __init__.
'''

__version__ = "1.0.0"
__author__ = "Rui Wang"
__email__ = "rw3594@nyu.edu"

# Import main functionality
from .core import (
    AlignmentProcessor,
    AminoAcidMutationProcessor,
    GeneMutationDetector,
    GenomeSNPProcessor,
    MultipleSequenceAlignment,
    ORFProcessor,
    ProMutationDetector,
    Proteome,
    ReferenceGenome,
)

# Import utilities
from .utils import list_utils, mutation_utils, sequence_utils

# Import PRF analyzer


__all__ = [
    # Core classes
    "ReferenceGenome",
    "GeneMutationDetector",
    "MultipleSequenceAlignment",
    "Proteome",
    "ProMutationDetector",
    "ORFProcessor",
    "AminoAcidMutationProcessor",
    "AlignmentProcessor",
    "GenomeSNPProcessor",
    # Utility modules
    "mutation_utils",
    "sequence_utils",
    "list_utils",
    # PRF analyzer

]
