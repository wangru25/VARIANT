"""
Simple comprehensive test for mutation_detector module using real data.
"""

import pytest
import os
from src.core.mutation_detector import (
    ReferenceGenome,
    MultipleSequenceAlignment,
    GeneMutationDetector,
)


def test_mutation_detector_with_real_data():
    """Test the complete mutation detection workflow with real SARS-CoV-2 data."""
    
    # Check if data files exist
    ref_file = "data/SARS-CoV-2/refs/NC_045512.fasta"
    msa_file = "data/SARS-CoV-2/clustalW/test_msa_1.txt"
    
    # Get absolute paths
    current_dir = os.path.dirname(os.path.abspath(__file__))
    project_root = os.path.dirname(current_dir)
    ref_file = os.path.join(project_root, ref_file)
    msa_file = os.path.join(project_root, msa_file)
    
    if not (os.path.exists(ref_file) and os.path.exists(msa_file)):
        pytest.skip(f"Required data files not found. Looking for:\n  {ref_file}\n  {msa_file}")
    
    print("Testing mutation detector with real SARS-CoV-2 data...")
    
    # 1. Test ReferenceGenome
    print("  - Loading reference genome...")
    ref_genome = ReferenceGenome(ref_file)
    assert ref_genome.ref_name == "NC_045512"
    assert len(ref_genome.ref_seq) == 29903  # SARS-CoV-2 genome length
    print(f"    ✓ Reference genome loaded: {len(ref_genome.ref_seq)} bp")
    
    # 2. Test MultipleSequenceAlignment
    print("  - Loading MSA file...")
    msa = MultipleSequenceAlignment(msa_file)
    assert len(msa.align) == 2
    print(f"    ✓ MSA loaded: {len(msa.align)} sequences")
    
    # 3. Test finding sequences
    print("  - Finding reference and query sequences...")
    ref_idx, ref_seq = msa.find_ref_msa("NC_045512")
    query_idx, query_seq = msa.find_seq_msa("EPI_ISL_16327572")
    assert ref_idx == 0
    assert query_idx == 1
    print(f"    ✓ Reference sequence found at index {ref_idx}")
    print(f"    ✓ Query sequence found at index {query_idx}")
    
    # 4. Test GeneMutationDetector
    print("  - Creating mutation detector...")
    detector = GeneMutationDetector(ref_genome, msa, "EPI_ISL_16327572")
    assert detector.genome_id == "EPI_ISL_16327572"
    print(f"    ✓ Mutation detector created for {detector.genome_id}")
    
    # 5. Test mutation detection
    print("  - Detecting mutations...")
    print("==============================Insertions==============================")
    # Insertions
    insertions = detector.get_insertions()
    print(f"    ✓ Insertions found: {len(insertions)}")
    print(insertions)
    print("==============================Deletions==============================")
    # Deletions
    deletions = detector.get_deletions()
    print(f"    ✓ Deletions found: {len(deletions)}")
    print(deletions)
    print("==============================Row Mutations==============================")
    # Other mutations (row, hot, point)
    row_mutations, hot_mutations, point_mutations = detector.get_others_mutations()
    print(row_mutations)
    print("==============================Hot Mutations==============================")
    print(hot_mutations)
    print("==============================Point Mutations==============================")
    print(point_mutations)
    print(f"    ✓ Point mutations found: {len(point_mutations)}")

    # 6. Validate results
    total_mutations = len(insertions) + len(deletions) + len(row_mutations) + len(hot_mutations) + len(point_mutations)
    print(f"    ✓ Total mutations detected: {total_mutations}")
    
    # Verify we found some mutations (real data should have differences)
    assert total_mutations > 0, "Should detect mutations in real data"
    
    # 7. Validate mutation structure
    all_mutations = insertions + deletions + row_mutations + hot_mutations + point_mutations
    
    for i, mutation in enumerate(all_mutations):
        assert isinstance(mutation, dict), f"Mutation {i} should be a dictionary"
        assert "key" in mutation, f"Mutation {i} missing 'key' field"
        assert "type" in mutation, f"Mutation {i} missing 'type' field"
        assert "pos" in mutation, f"Mutation {i} missing 'pos' field"
        assert "SNP" in mutation, f"Mutation {i} missing 'SNP' field"
    
    print("    ✓ All mutations have correct structure")
    
    # 8. Show some example mutations
    if total_mutations > 0:
        print("\n  Example mutations found:")
        for i, mutation in enumerate(all_mutations[:3]):  # Show first 3
            print(f"    {i+1}. {mutation['type']}: {mutation['pos']} {mutation['SNP']}")
        if total_mutations > 3:
            print(f"    ... and {total_mutations - 3} more")
    
    print("\n✅ All tests passed! Mutation detector is working correctly.")


def test_mutation_dict_creation():
    """Test the create_mutation_dict static method."""
    
    # Test single position mutation
    result = GeneMutationDetector.create_mutation_dict("pointMutation", 241, "C->T")
    assert result["type"] == "pointMutation"
    assert result["pos"] == "241"
    assert result["SNP"] == "C->T"
    assert "key" in result
    
    # Test multiple position mutation
    result = GeneMutationDetector.create_mutation_dict(
        "rowMutation", [241, 242, 243], ["C->T", "T->G", "A->C"]
    )
    assert result["type"] == "rowMutation"
    assert result["pos"] == "241:243"
    assert result["SNP"] == "CTA->TGC"
    assert "key" in result
    
    print("✅ Mutation dictionary creation tests passed!")


if __name__ == "__main__":

    test_mutation_detector_with_real_data()

    

