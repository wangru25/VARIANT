"""
Author: Rui Wang
Date: 2025-06-23 14:33:45
LastModifiedBy: Rui Wang
LastEditTime: 2025-06-23 14:36:23
Email: wang.rui@nyu.edu
FilePath: /7_MutParser/src/core/__init__.py
Description:
"""

"""
Core functionality for MutParser.

This package contains the main classes and functions for virus mutation analysis.
"""

from .genome_processor import GenomeSNPProcessor
from .mutation_detector import GeneMutationDetector, MultipleSequenceAlignment
from .protein_analyzer import (
    AlignmentProcessor,
    AminoAcidMutationProcessor,
    ORFProcessor,
    ProMutationDetector,
    Proteome,
)
from .reference_genome import ReferenceGenome
from .mutation_processor import MutationProcessor

__all__ = [
    "ReferenceGenome",
    "GeneMutationDetector",
    "MultipleSequenceAlignment",
    "Proteome",
    "ProMutationDetector",
    "ORFProcessor",
    "AminoAcidMutationProcessor",
    "AlignmentProcessor",
    "GenomeSNPProcessor",
    "MutationProcessor",
]
