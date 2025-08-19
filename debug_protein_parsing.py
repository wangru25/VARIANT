#!/usr/bin/env python3
"""
Debug script to understand protein boundary parsing
"""

import sys
import os
from pathlib import Path
from Bio import SeqIO
import re

# Add src to path
sys.path.append(str(Path(__file__).parent / "src"))

def debug_protein_boundaries():
    """Debug protein boundary parsing"""
    proteome_file = "data/HIV-1/refs/HIV-1_proteome.fasta"
    
    if not os.path.exists(proteome_file):
        print(f"Error: {proteome_file} not found!")
        return
    
    protein_boundaries = {}
    
    print("Parsing HIV-1 proteome file...")
    print("=" * 80)
    
    for record in SeqIO.parse(proteome_file, "fasta"):
        header = record.description
        print(f"Header: {header}")
        
        # Handle join() format first
        join_match = re.search(r'\|join\((\d+)\.\.(\d+),(\d+)\.\.(\d+)\)$', header)
        if join_match:
            start1, end1, start2, end2 = map(int, join_match.groups())
            start, end = start1, end2  # Full range coverage
            protein_name_match = re.search(r'\|([^|]+)\|join\(', header)
            if protein_name_match:
                protein_name = protein_name_match.group(1)
                protein_boundaries[protein_name] = (start, end)
                print(f"  → JOIN: {protein_name} = ({start}, {end}) [covers {start1}-{end1}, {start2}-{end2}]")
        else:
            # Handle complement format
            complement_match = re.search(r'\|complement\((\d+)\.\.(\d+)\)$', header)
            if complement_match:
                start, end = int(complement_match.group(1)), int(complement_match.group(2))
                protein_name_match = re.search(r'\|([^|]+)\|complement\(', header)
                if protein_name_match:
                    protein_name = protein_name_match.group(1)
                    protein_boundaries[protein_name] = (start, end)
                    print(f"  → COMPLEMENT: {protein_name} = ({start}, {end})")
            else:
                # Handle simple format
                location_match = re.search(r'\|(\d+)\.\.(\d+)$', header)
                if location_match:
                    start, end = int(location_match.group(1)), int(location_match.group(2))
                    protein_name_match = re.search(r'\|([^|]+)\|\d+\.\.\d+$', header)
                    if protein_name_match:
                        protein_name = protein_name_match.group(1)
                        protein_boundaries[protein_name] = (start, end)
                        print(f"  → SIMPLE: {protein_name} = ({start}, {end})")
        print()
    
    print("=" * 80)
    print("Final protein boundaries:")
    for protein, (start, end) in sorted(protein_boundaries.items(), key=lambda x: x[1][0]):
        print(f"  {protein}: {start}-{end}")
    
    print("\n" + "=" * 80)
    print("Testing specific positions:")
    test_positions = [6919, 7200, 7488, 7925, 7970]
    
    for pos in test_positions:
        matches = []
        for protein, (start, end) in protein_boundaries.items():
            if start <= pos <= end:
                matches.append(protein)
        print(f"Position {pos}: {matches}")

if __name__ == "__main__":
    debug_protein_boundaries()
