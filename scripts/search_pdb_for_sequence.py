#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-01-15 14:30:00
LastModifiedBy: Rui Wang
LastEditTime: 2025-01-15 14:30:00
Email: rw3594@nyu.edu
FilePath: /VARIANT/scripts/search_pdb_for_sequence.py
Description: Script to search PDB database for structures matching a given RNA sequence.
'''

import sys
from pathlib import Path
from Bio import SeqIO
from Bio.Seq import Seq
import urllib.request
import urllib.parse

def rna_to_dna(seq_rna: str) -> str:
    '''Convert RNA sequence to DNA sequence.'''
    return seq_rna.replace('U', 'T').replace('u', 't')

def search_pdb_blast(sequence: str, sequence_type: str = 'rna') -> None:
    '''
    Search PDB using BLAST for matching structures.
    
    Args:
        sequence: RNA or DNA sequence
        sequence_type: 'rna' or 'dna'
    '''
    print(f"Searching PDB for sequence: {sequence[:50]}...")
    print(f"Sequence length: {len(sequence)} nucleotides")
    
    # Convert RNA to DNA for BLAST search
    if sequence_type == 'rna':
        dna_seq = rna_to_dna(sequence)
    else:
        dna_seq = sequence
    
    # PDB BLAST URL
    blast_url = "https://www.rcsb.org/search/sequence"
    
    print(f"\nTo search PDB for this sequence:")
    print(f"1. Go to: {blast_url}")
    print(f"2. Paste your sequence (DNA format):")
    print(f"   {dna_seq}")
    print(f"\nOr use the RCSB PDB API:")
    print(f"   https://search.rcsb.org/#sequence")
    
    # Try to use BioPython to search (if internet available)
    try:
        print(f"\nAttempting to search PDB via API...")
        # This would require the RCSB PDB API client
        print("Note: Direct API search requires additional setup.")
        print("Recommended: Use the web interface at https://www.rcsb.org/search/sequence")
    except Exception as e:
        print(f"API search not available: {e}")

def extract_sequence_from_fasta(fasta_file: Path) -> list:
    '''Extract sequences from FASTA file.'''
    sequences = []
    try:
        for record in SeqIO.parse(fasta_file, "fasta"):
            seq_str = str(record.seq)
            sequences.append({
                'id': record.id,
                'description': record.description,
                'sequence': seq_str,
                'length': len(seq_str)
            })
    except Exception as e:
        print(f"Error reading FASTA file: {e}")
    return sequences

def main():
    '''Main function to search PDB for sequences.'''
    if len(sys.argv) < 2:
        print("Usage: python search_pdb_for_sequence.py <fasta_file> [sequence_id]")
        print("\nExample:")
        print("  python search_pdb_for_sequence.py data/RSV-1/fasta/sequences_w_reference.fasta")
        print("  python search_pdb_for_sequence.py data/RSV-1/fasta/sequences_w_reference.fasta 'RSV-1|Reference genome'")
        sys.exit(1)
    
    fasta_file = Path(sys.argv[1])
    if not fasta_file.exists():
        print(f"Error: FASTA file not found: {fasta_file}")
        sys.exit(1)
    
    sequences = extract_sequence_from_fasta(fasta_file)
    
    if not sequences:
        print("No sequences found in FASTA file")
        sys.exit(1)
    
    # If specific sequence ID provided
    if len(sys.argv) > 2:
        seq_id = sys.argv[2]
        matching_seqs = [s for s in sequences if seq_id in s['id'] or seq_id in s['description']]
        if matching_seqs:
            sequences = matching_seqs
        else:
            print(f"Warning: Sequence ID '{seq_id}' not found, showing all sequences")
    
    print(f"Found {len(sequences)} sequence(s) in {fasta_file}\n")
    
    for seq_info in sequences:
        print("=" * 60)
        print(f"Sequence ID: {seq_info['id']}")
        print(f"Description: {seq_info['description']}")
        print(f"Length: {seq_info['length']} nucleotides")
        print(f"\nRNA Sequence:")
        print(f"{seq_info['sequence']}")
        print(f"\nDNA Sequence (for PDB search):")
        print(f"{rna_to_dna(seq_info['sequence'])}")
        print("\n" + "-" * 60)
        
        # Search PDB
        search_pdb_blast(seq_info['sequence'], sequence_type='rna')
        print()

if __name__ == "__main__":
    main()