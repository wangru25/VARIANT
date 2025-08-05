#!/usr/bin/env python3
"""
Debug script to understand why ORF generation can produce None values
"""

import sys
import os

# Set up path
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_dir)

def debug_orf_none_generation():
    """Debug why ORFs become None"""
    
    print("=== ORF NONE GENERATION DEBUG ===")
    
    from src.core.reference_genome import ReferenceGenome
    from src.core.protein_analyzer import ORFProcessor
    
    # Test both viruses to compare
    configs = {
        "SARS-CoV-2": "data/SARS-CoV-2/refs/NC_045512.fasta",
        "ZaireEbola": "data/ZaireEbola/refs/NC_002549.1.fasta"
    }
    
    for virus_name, ref_path in configs.items():
        print(f"\n--- {virus_name} ORF Analysis ---")
        
        try:
            ref_genome = ReferenceGenome(ref_path)
            orf_processor = ORFProcessor(ref_genome)
            
            print(f"Reference sequence length: {len(ref_genome.ref_seq)}")
            
            # Test various deletion scenarios that might produce None ORFs
            test_positions = [
                "1:5",      # Very start of genome
                "10:15",    # Early in genome
                "100:105",  # Middle
                f"{len(ref_genome.ref_seq)-10}:{len(ref_genome.ref_seq)-5}",  # Near end
                f"{len(ref_genome.ref_seq)-5}:{len(ref_genome.ref_seq)}",     # Very end
            ]
            
            for pos in test_positions:
                print(f"\n  Testing deletion at position {pos}")
                
                try:
                    orfs_result = orf_processor.get_orfs_dels(pos)
                    original_orfs, mutated_orfs = orfs_result
                    
                    print(f"    Original ORFs: {len(original_orfs)}")
                    print(f"    Mutated ORFs: {len(mutated_orfs)}")
                    
                    # Check for None values
                    orig_none_count = sum(1 for orf in original_orfs if orf is None)
                    mut_none_count = sum(1 for orf in mutated_orfs if orf is None)
                    orig_empty_count = sum(1 for orf in original_orfs if orf == "")
                    mut_empty_count = sum(1 for orf in mutated_orfs if orf == "")
                    
                    print(f"    Original: {orig_none_count} None, {orig_empty_count} empty")
                    print(f"    Mutated: {mut_none_count} None, {mut_empty_count} empty")
                    
                    if orig_none_count > 0 or mut_none_count > 0:
                        print(f"    ⚠️ Found None ORFs!")
                        
                        # Show details
                        for i, (orig, mut) in enumerate(zip(original_orfs, mutated_orfs)):
                            if orig is None or mut is None:
                                print(f"      ORF {i}: orig={orig}, mut={mut}")
                    
                    # Check for very short sequences that might cause issues
                    for i, (orig, mut) in enumerate(zip(original_orfs, mutated_orfs)):
                        if orig is not None and len(orig) < 3:
                            print(f"    ⚠️ Very short original ORF {i}: '{orig}' (len={len(orig)})")
                        if mut is not None and len(mut) < 3:
                            print(f"    ⚠️ Very short mutated ORF {i}: '{mut}' (len={len(mut)})")
                    
                except Exception as e:
                    print(f"    ❌ Error: {e}")
            
            # Test extreme cases that are more likely to produce None
            print(f"\n  Testing extreme cases...")
            
            # Test deletions that might remove too much sequence
            extreme_cases = [
                "1:50",     # Large deletion at start
                "50:100",   # Large deletion in middle
                f"{len(ref_genome.ref_seq)-50}:{len(ref_genome.ref_seq)}",  # Large deletion at end
            ]
            
            for pos in extreme_cases:
                print(f"\n    Extreme case: deletion at {pos}")
                
                try:
                    # Manually check what get_orfs_dels does step by step
                    pos_a, pos_b = map(int, pos.split(":"))
                    pos_a -= 1
                    
                    num_pre_nts = 15
                    num_post_nts = 15
                    
                    # Check bounds
                    start_pos = pos_a - num_pre_nts
                    end_pos = pos_b + num_post_nts
                    
                    print(f"      Extracting sequence from {start_pos} to {end_pos}")
                    print(f"      Genome length: {len(ref_genome.ref_seq)}")
                    
                    if start_pos < 0:
                        print(f"      ⚠️ Start position {start_pos} is negative!")
                    if end_pos > len(ref_genome.ref_seq):
                        print(f"      ⚠️ End position {end_pos} exceeds genome length!")
                    
                    # Extract sequences
                    seq_t = ref_genome.ref_seq[max(0, start_pos):min(len(ref_genome.ref_seq), end_pos)]
                    seq_t_del = (
                        ref_genome.ref_seq[max(0, start_pos):pos_a] + 
                        ref_genome.ref_seq[pos_b:min(len(ref_genome.ref_seq), end_pos)]
                    )
                    
                    print(f"      Original sequence length: {len(seq_t)}")
                    print(f"      Deleted sequence length: {len(seq_t_del)}")
                    
                    if len(seq_t) < 3:
                        print(f"      ⚠️ Original sequence too short for translation!")
                    if len(seq_t_del) < 3:
                        print(f"      ⚠️ Deleted sequence too short for translation!")
                    
                    # Try translation
                    if len(seq_t) >= 3:
                        original_orfs = orf_processor.translate_dna_3_frames(seq_t, codon_type="reference")
                        print(f"      Original ORFs: {len(original_orfs)} - {[len(orf) if orf else 'None' for orf in original_orfs]}")
                    
                    if len(seq_t_del) >= 3:
                        mutated_orfs = orf_processor.translate_dna_3_frames(seq_t_del, codon_type="mutant")
                        print(f"      Mutated ORFs: {len(mutated_orfs)} - {[len(orf) if orf else 'None' for orf in mutated_orfs]}")
                    
                except Exception as e:
                    print(f"      ❌ Extreme case error: {e}")
                    import traceback
                    traceback.print_exc()
            
        except Exception as e:
            print(f"❌ Error with {virus_name}: {e}")

def test_translate_dna_3_frames():
    """Test the translate_dna_3_frames method directly with edge cases"""
    
    print(f"\n=== TRANSLATE DNA 3 FRAMES TEST ===")
    
    from src.core.reference_genome import ReferenceGenome
    from src.core.protein_analyzer import ORFProcessor
    
    # Use ZaireEbola since it has more issues
    ref_genome = ReferenceGenome("data/ZaireEbola/refs/NC_002549.1.fasta")
    orf_processor = ORFProcessor(ref_genome)
    
    # Test edge cases that might produce None or empty ORFs
    test_sequences = [
        ("Empty", ""),
        ("Too short", "AT"),
        ("3 nucleotides", "ATG"),
        ("6 nucleotides", "ATGCCC"),
        ("Stop codons", "TAATAGTAG"),
        ("Invalid codons", "NNNNNNNNN"),
        ("Mixed valid/invalid", "ATGNNNTAG"),
        ("Very short after filtering", "TAGTAGTAG"),  # All stop codons
    ]
    
    for name, seq in test_sequences:
        print(f"\n--- Testing: {name} ('{seq}') ---")
        
        try:
            # Test both reference and mutant modes
            for codon_type in ["reference", "mutant"]:
                print(f"  {codon_type} mode:")
                
                if len(seq) < 3:
                    print(f"    Sequence too short for translation")
                    continue
                
                result = orf_processor.translate_dna_3_frames(seq, codon_type=codon_type)
                print(f"    Result: {result}")
                print(f"    Types: {[type(orf) for orf in result]}")
                print(f"    Lengths: {[len(orf) if orf else 'None' for orf in result]}")
                
                # Check for None values
                none_count = sum(1 for orf in result if orf is None)
                empty_count = sum(1 for orf in result if orf == "")
                
                if none_count > 0:
                    print(f"    ⚠️ Found {none_count} None ORFs")
                if empty_count > 0:
                    print(f"    ⚠️ Found {empty_count} empty ORFs")
        
        except Exception as e:
            print(f"    ❌ Error: {e}")

if __name__ == "__main__":
    debug_orf_none_generation()
    test_translate_dna_3_frames()