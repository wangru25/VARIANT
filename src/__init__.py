# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-08-19 17:30:00
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-19 17:30:00
Email: wang.rui@nyu.edu
FilePath: /VARIANT/./src/__init__.py
Description: Python module for __init__.
'''


"""
MutParser: A Comprehensive Framework for Multi-Scale Viral Mutation Analysis.

This package provides tools for analyzing viral mutations at multiple scales,
from nucleotide changes to protein impacts, including programmed ribosomal
frameshifting analysis for SARS-CoV-2 and other viruses.
"""

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
