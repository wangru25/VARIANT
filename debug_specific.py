#!/usr/bin/env python3
"""
Debug script to trace the specific path where query_sequence becomes None
"""

import sys
import os

# Set up path for imports
current_dir = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, current_dir)

def debug_deletion_processing():
    """Debug the specific deletion processing that causes the error"""
    
    print("=== DELETION PROCESSING DEBUG ===")
    
    from src.core.reference_genome import ReferenceGenome
    from src.core.mutation_detector import MultipleSequenceAlignment, GeneMutationDetector
    from src.core.protein_analyzer import ProMutationDetector, ORFProcessor, AlignmentProcessor, AminoAcidMutationProcessor
    from src.core.reference_genome import Proteome
    
    # ZaireEbola configuration (the one that fails)
    ref_genome_path = "data/ZaireEbola/refs/NC_002549.1.fasta"
    proteome_path = "data/ZaireEbola/refs/ZaireEbola_proteome.fasta"
    msa_path = "data/ZaireEbola/clustalW/zaire_ebola_mas.txt"
    genome_id = "AB050936v1"
    
    print(f"Analyzing ZaireEbola deletion processing...")
    
    # Initialize components
    ref_genome = ReferenceGenome(ref_genome_path)
    proteome = Proteome(proteome_path)
    msa = MultipleSequenceAlignment(msa_path)
    gene_mut_detector = GeneMutationDetector(ref_genome, msa, genome_id)
    
    print(f"Reference genome length: {len(ref_genome.ref_seq)}")
    
    # Get deletions - this should work without the protein analysis
    print("\n--- Getting deletions ---")
    try:
        deletions = gene_mut_detector.get_deletions()
        print(f"✅ Found {len(deletions)} deletions")
        
        if deletions:
            del_mutation = deletions[0]
            print(f"First deletion: {del_mutation}")
            
            # Now trace what happens in protein_dels_snps
            print(f"\n--- Tracing protein_dels_snps ---")
            
            # Create the protein mutation detector
            orf_processor = ORFProcessor(ref_genome)
            aa_mutation_processor = AminoAcidMutationProcessor()
            alignment_processor = AlignmentProcessor(proteome, aa_mutation_processor)
            pro_mut_detector = ProMutationDetector(
                orf_processor, alignment_processor, aa_mutation_processor, ref_genome
            )
            
            # This is where the error occurs - let's trace step by step
            print(f"Calling protein_dels_snps...")
            
            try:
                # This should trigger the error
                result = pro_mut_detector.protein_dels_snps([del_mutation])
                print(f"✅ protein_dels_snps succeeded: {result}")
            except Exception as e:
                print(f"❌ protein_dels_snps failed: {e}")
                
                # Let's manually trace what happens inside protein_dels_snps
                print(f"\n--- Manual tracing of protein_dels_snps ---")
                
                # Get the deletion info
                pos_start = del_mutation['pos_start']
                pos_end = del_mutation['pos_end']
                del_seq = del_mutation['seq']
                
                print(f"Deletion: pos {pos_start}-{pos_end}, seq '{del_seq}'")
                
                # Get original ORFs
                print(f"Getting original ORFs...")
                try:
                    original_orfs = orf_processor.get_orfs_dels(f"{pos_start}:{pos_end}")
                    print(f"Original ORFs: {len(original_orfs[0])} ORFs")
                    print(f"First few original ORFs: {[orf[:20] if orf else 'None' for orf in original_orfs[0][:3]]}")
                except Exception as orf_e:
                    print(f"❌ Error getting original ORFs: {orf_e}")
                    return
                
                # Create mutated sequence
                print(f"Creating mutated sequence...")
                ref_seq = ref_genome.ref_seq
                mutated_seq = ref_seq[:pos_start-1] + ref_seq[pos_end:]
                print(f"Original length: {len(ref_seq)}, Mutated length: {len(mutated_seq)}")
                
                # Create temporary reference genome for mutated sequence
                print(f"Creating mutated reference genome...")
                # This is likely where the issue occurs
                
                # Let's see what happens when we create a new ReferenceGenome with mutated sequence
                # We need to create a temporary file for this
                import tempfile
                with tempfile.NamedTemporaryFile(mode='w', suffix='.fasta', delete=False) as temp_file:
                    temp_file.write(f">temp_mutated\n{mutated_seq}\n")
                    temp_fasta_path = temp_file.name
                
                try:
                    mutated_ref_genome = ReferenceGenome(temp_fasta_path)
                    print(f"✅ Created mutated reference genome: {len(mutated_ref_genome.ref_seq)}")
                    
                    # Create mutated ORF processor
                    mutated_orf_processor = ORFProcessor(mutated_ref_genome)
                    print(f"✅ Created mutated ORF processor")
                    
                    # Get mutated ORFs
                    print(f"Getting mutated ORFs...")
                    try:
                        mutated_orfs = mutated_orf_processor.get_orfs_dels(f"{pos_start}:{pos_end}")
                        print(f"Mutated ORFs: {len(mutated_orfs[0])} ORFs")
                        print(f"First few mutated ORFs: {[orf[:20] if orf else 'None' for orf in mutated_orfs[0][:3]]}")
                        
                        # Now check what happens in _process_target_orfs
                        print(f"\n--- Checking _process_target_orfs ---")
                        
                        # Find the first non-None ORF pair to test
                        for i, (orf, m_orf) in enumerate(zip(original_orfs[0], mutated_orfs[0])):
                            print(f"ORF pair {i}: orig={type(orf)} {len(orf) if orf else 'None'}, mut={type(m_orf)} {len(m_orf) if m_orf else 'None'}")
                            
                            if orf is None or m_orf is None:
                                print(f"  ❌ Found None ORF at index {i}")
                                if i >= 5:  # Stop after checking first few
                                    break
                            else:
                                print(f"  ✅ Valid ORF pair at index {i}")
                                
                                # This is where align_local gets called with None values
                                print(f"    Calling align_local with orf='{orf[:20]}...', m_orf='{m_orf[:20]}...'")
                                break
                    
                    except Exception as mut_orf_e:
                        print(f"❌ Error getting mutated ORFs: {mut_orf_e}")
                    
                except Exception as ref_e:
                    print(f"❌ Error creating mutated reference genome: {ref_e}")
                finally:
                    # Clean up temp file
                    try:
                        os.unlink(temp_fasta_path)
                    except:
                        pass
        
    except Exception as del_e:
        print(f"❌ Error getting deletions: {del_e}")
        import traceback
        traceback.print_exc()

if __name__ == "__main__":
    debug_deletion_processing()