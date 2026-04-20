# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2026-04-15 12:22:57
LastModifiedBy: Rui Wang
LastEditTime: 2026-04-15 16:40:36
Email: wang.rui@nyu.edu
FilePath: /Test-EDS/convert_dbn_to_dssr.py
Description: 
'''
#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Convert .dbn file to X3DNA-DSSR format
"""

import sys
import os
from pathlib import Path

def infer_pdb_id_from_path(dbn_file):
    stem = Path(dbn_file).stem
    base = stem.split("-")[0].split("_")[0].strip()
    return base.upper() if base else "INPUT"


def convert_dbn_to_dssr(dbn_file, output_file, pdb_id=None):
    """
    Convert .dbn file to X3DNA-DSSR format
    
    .dbn files typically have format:
    >header
    sequence
    dot-bracket
    """
    if pdb_id is None:
        pdb_id = infer_pdb_id_from_path(dbn_file)
    with open(dbn_file, 'r') as f:
        lines = [line.strip() for line in f.readlines() if line.strip()]
    
    print(f"Reading {dbn_file}...")
    print(f"Found {len(lines)} non-empty lines")
    
    # Parse .dbn format
    sequence = None
    dot_bracket = None
    header = None
    
    i = 0
    while i < len(lines):
        line = lines[i]
        
        # Check for header (starts with >)
        if line.startswith('>'):
            header = line[1:].strip()
            print(f"Found header: {header}")
            i += 1
            continue
        
        # Next non-header line is usually sequence
        if sequence is None and not line.startswith('>'):
            # Check if it looks like a sequence (only A, U, G, C, T, N, etc.)
            if all(c.upper() in 'AUCGTN-&' for c in line):
                sequence = line
                print(f"Found sequence: {sequence[:50]}...")
                i += 1
                continue
        
        # Next line after sequence is usually dot-bracket
        if sequence is not None and dot_bracket is None:
            # Check if it looks like dot-bracket (has parentheses, brackets, dots)
            if any(c in '().[]{}&' for c in line):
                dot_bracket = line
                print(f"Found dot-bracket: {dot_bracket[:50]}...")
                i += 1
                continue
        
        i += 1
    
    if not sequence or not dot_bracket:
        print("ERROR: Could not find sequence and/or dot-bracket notation")
        print("File contents:")
        for i, line in enumerate(lines[:10], 1):
            print(f"  {i}: {line}")
        return False
    
    # Parse multiple strands if present
    strands = []
    current_strand = None
    
    with open(dbn_file, 'r') as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith('>'):
                if current_strand:
                    strands.append(current_strand)
                current_strand = {'header': line[1:], 'sequence': None, 'dot_bracket': None}
            elif current_strand:
                if current_strand['sequence'] is None:
                    # Check if it's sequence (mostly A, U, G, C)
                    if all(c.upper() in 'AUCGTN-&' for c in line):
                        current_strand['sequence'] = line
                elif current_strand['dot_bracket'] is None:
                    # Next line is dot-bracket
                    if any(c in '().[]{}&' for c in line):
                        current_strand['dot_bracket'] = line
    
    if current_strand:
        strands.append(current_strand)
    
    # Create X3DNA-DSSR format
    output_lines = []
    output_lines.append("Secondary structures in dot-bracket notation (dbn) as a whole and per chain\n")
    
    # Whole structure (combine all strands)
    if len(strands) > 0:
        whole_seq = '&'.join([s['sequence'] for s in strands if s['sequence']])
        whole_db = '&'.join([s['dot_bracket'] for s in strands if s['dot_bracket']])
        output_lines.append(f">{pdb_id}\n")
        output_lines.append(f"{whole_seq}\n")
        output_lines.append(f"{whole_db}\n")
        
        # Individual chains
        for idx, strand in enumerate(strands, 1):
            chain_id = chr(64 + idx)  # A, B, C, etc.
            if strand['sequence'] and strand['dot_bracket']:
                output_lines.append(f">Chain {chain_id}\n")
                output_lines.append(f"{strand['sequence']}\n")
                output_lines.append(f"{strand['dot_bracket']}\n")
    else:
        # Fallback to original method
        output_lines.append(f">{pdb_id}\n")
        output_lines.append(f"{sequence}\n")
        output_lines.append(f"{dot_bracket}\n")
        output_lines.append(">Chain A\n")
        output_lines.append(f"{sequence}\n")
        output_lines.append(f"{dot_bracket}\n")
    
    # Write output
    output_dir = os.path.dirname(output_file)
    if output_dir:
        os.makedirs(output_dir, exist_ok=True)
    with open(output_file, 'w') as f:
        f.writelines([line + '\n' for line in output_lines])
    
    print("\n✅ Successfully converted!")
    print(f"Output file: {output_file}")
    print("\nFile preview:")
    with open(output_file, 'r') as f:
        preview = f.read()
        print(preview[:500])
        if len(preview) > 500:
            print("...")
    
    return True

if __name__ == "__main__":
    if len(sys.argv) < 2:
        print("Usage: python convert_dbn_to_dssr.py <input.dbn> [output.out] [pdb_id]")
        print("\nExample:")
        print("  python convert_dbn_to_dssr.py ~/Downloads/ABCD-2D.dbn")
        sys.exit(1)
    
    dbn_file = sys.argv[1]
    inferred_id = infer_pdb_id_from_path(dbn_file)
    output_file = sys.argv[2] if len(sys.argv) > 2 else os.path.join("PDB_DSSR", f"{inferred_id}.out")
    pdb_id = sys.argv[3] if len(sys.argv) > 3 else inferred_id
    
    if not os.path.isfile(dbn_file):
        print(f"ERROR: File not found: {dbn_file}")
        sys.exit(1)
    
    success = convert_dbn_to_dssr(dbn_file, output_file, pdb_id)
    
    if success:
        print("\n" + "="*60)
        print("Next steps:")
        print("="*60)
        print("1. Verify the output file looks correct")
        print("2. Run: python3 PDBto2D.py")
        print("3. Run: python3 Dual_Library.py")
    else:
        sys.exit(1)
