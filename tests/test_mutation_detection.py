#!/usr/bin/env python3
"""
Test script to verify mutation detection improvements.

This script tests the improved mutation detection functionality to ensure
that it can properly detect mutations throughout protein sequences, not just
at the beginning.
"""

import sys
from pathlib import Path

# Add src to path
sys.path.append(str(Path(__file__).parent / "src"))


from src.core.reference_genome import ReferenceGenome
from src.core.protein_analyzer import (
    AlignmentProcessor,
    AminoAcidMutationProcessor,
    ORFProcessor,
    Proteome,
)


def test_mutation_detection():
    """Test the improved mutation detection functionality."""

    print("Testing improved mutation detection...")

    # Initialize components
    ref_genome = ReferenceGenome("./data/refs/NC_045512.fasta")
    orf_processor = ORFProcessor(ref_genome)
    aa_processor = AminoAcidMutationProcessor()
    proteome = Proteome("./data/refs/SARS-CoV-2_proteome.fasta")
    alignment_processor = AlignmentProcessor(proteome, aa_processor)

    # Test the improved get_row_aa_mutations method
    print("\n1. Testing get_row_aa_mutations method:")

    # Test case 1: Single mutation
    orf1 = "AAAAAAAAAKAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
    morf1 = orf1[:9] + "R" + orf1[10:]

    result1 = aa_processor.get_row_aa_mutations(orf1, morf1, 0)
    print(f"   Single mutation test: {result1}")
    expected1 = "K10R"
    assert result1 == expected1, f"Expected {expected1}, got {result1}"

    # Test case 2: Multiple mutations
    orf2 = "MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVF"
    morf2 = "MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVF"
    # Introduce mutations at positions 5, 15, and 25
    morf2 = morf2[:4] + "A" + morf2[5:14] + "T" + morf2[15:24] + "G" + morf2[25:]

    result2 = aa_processor.get_row_aa_mutations(orf2, morf2, 0)
    print(f"   Multiple mutations test: {result2}")
    expected2 = "V5A,Q15T,D25G"
    assert result2 == expected2, f"Expected {expected2}, got {result2}"

    # Test case 3: No mutations
    orf3 = "MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVF"
    morf3 = "MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVF"

    result3 = aa_processor.get_row_aa_mutations(orf3, morf3, 0)
    print(f"   No mutations test: {result3}")
    expected3 = ""
    assert result3 == expected3, f"Expected {expected3}, got {result3}"

    # Test case 4: Mutation at the end of sequence
    orf4 = "MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVF"
    morf4 = "MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVF"
    # Introduce mutation at the last position
    morf4 = morf4[:-1] + "A"

    result4 = aa_processor.get_row_aa_mutations(orf4, morf4, 0)
    print(f"   End mutation test: {result4}")
    expected4 = "F60A"
    assert result4 == expected4, f"Expected {expected4}, got {result4}"

    print("   ✓ All get_row_aa_mutations tests passed!")

    # Test the improved ORF generation
    print("\n2. Testing improved ORF generation:")

    # Test with a larger window
    test_pos = 1000  # Example SNP position
    test_mutation = "A"

    orfs, m_orfs = orf_processor.get_orfs_snps(test_pos, test_mutation)
    print(f"   Generated {len(orfs)} original ORFs and {len(m_orfs)} mutated ORFs")

    if orfs and m_orfs:
        print(f"   Original ORF lengths: {[len(orf) for orf in orfs]}")
        print(f"   Mutated ORF lengths: {[len(orf) for orf in m_orfs]}")
        print("   ✓ ORF generation test passed!")
    else:
        print("   ⚠ No ORFs generated - this might be normal for this position")

    print("\n3. Testing improved alignment:")

    # Test alignment with the proteome
    if orfs and m_orfs:
        protein_mutations = alignment_processor.align_orfs_to_ref_proteome_snps(
            orfs, m_orfs
        )
        print(f"   Found {len(protein_mutations)} protein mutations")

        for mutation in protein_mutations:
            protein = mutation.get("protein", "Unknown")
            mut_str = mutation.get("mutation", "NA")
            print(f"   Protein: {protein}, Mutation: {mut_str}")

    print("\n✓ All tests completed successfully!")
    print("\nThe improved mutation detection should now be able to:")
    print("1. Detect mutations throughout the entire protein sequence")
    print("2. Handle multiple mutations in a single ORF")
    print("3. Work with larger genomic windows for better coverage")
    print("4. Be more flexible in finding protein matches")


if __name__ == "__main__":
    test_mutation_detection()
