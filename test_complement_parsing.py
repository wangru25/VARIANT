#!/usr/bin/env python3
"""
Test script to verify complement coordinate parsing in mutation_summary.py
"""

import sys
from pathlib import Path

# Add src to path
sys.path.append(str(Path(__file__).parent / "src"))

from src.utils.mutation_summary import _validate_protein_assignment

def test_complement_parsing():
    """Test complement coordinate parsing with HIV-1 proteome"""
    proteome_file = "data/HIV-1/refs/HIV-1_proteome.fasta"
    virus_name = "HIV-1"
    
    # Test positions that should fall within the Asp protein complement(6919..7488)
    test_cases = [
        (6919, "Asp"),  # Start position
        (7200, "Asp"),  # Middle position  
        (7488, "Asp"),  # End position
        (6918, "Envelope surface glycoprotein gp160"),  # Just before start (still in gp160)
        (7489, "Envelope surface glycoprotein gp160"),  # Just after end (still in gp160)
        (5000, "Vif"),  # Should fall in Vif protein (4587..5165)
        (1000, "Pr55(Gag)"),  # Should fall in Pr55(Gag) (336..1838)
        (7950, "Tat"),  # Should fall in Tat's second range (7925..7970)
        (5500, "Tat"),  # Should fall in Tat's first range (5377..5591)
    ]
    
    print("Testing complement coordinate parsing...")
    print(f"Using proteome file: {proteome_file}")
    print("=" * 60)
    
    for position, expected_protein in test_cases:
        result = _validate_protein_assignment(virus_name, position, "test_protein", proteome_file)
        status = "✓" if result == expected_protein else "✗"
        print(f"{status} Position {position:4d}: Expected '{expected_protein}', Got '{result}'")
    
    print("=" * 60)
    print("Test completed!")

if __name__ == "__main__":
    test_complement_parsing()
