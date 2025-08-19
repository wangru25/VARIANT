"""
Mutation Summary Utilities

This module provides functionality to extract and format mutation information
from analysis output files into structured CSV format.
"""

import os
import re
import csv
from typing import Dict, List, Tuple, Optional
from datetime import datetime


def extract_mutation_summary_to_csv(virus_name: str, genome_id: str, segment: Optional[str] = None, proteome_file: Optional[str] = None) -> None:
    """
    Extract mutation information from .txt output files and generate a comprehensive CSV file.
    
    Parameters:
    -----
    virus_name : str
        Name of the virus (e.g., 'SARS-CoV-2')
    genome_id : str
        Sample identifier (e.g., 'EPI_ISL_16327572')
    segment : Optional[str]
        Segment name for multi-segment viruses
    proteome_file : Optional[str]
        Path to the proteome FASTA file for protein boundary validation
    """
    print(f"Extracting mutation summary for {genome_id}...")
    
    # File paths
    from datetime import date
    today = date.today().strftime("%Y%m%d")
    today = 20250807
    
    if segment:
        txt_file = f"result/{virus_name}/{segment}/{genome_id}_{today}.txt"
        csv_file = f"result/{virus_name}/{segment}/{genome_id}_mutation_summary_{today}.csv"
    else:
        txt_file = f"result/{virus_name}/{genome_id}_{today}.txt"
        csv_file = f"result/{virus_name}/{genome_id}_mutation_summary_{today}.csv"
    
    if not os.path.exists(txt_file):
        print(f"Error: {txt_file} not found!")
        return
    
    mutations_data = []
    
    with open(txt_file, 'r') as f:
        content = f.read()
    
    # Parse different mutation types
    mutations_data.extend(_parse_point_mutations(content, virus_name, proteome_file))
    mutations_data.extend(_parse_deletions(content, virus_name, proteome_file))
    mutations_data.extend(_parse_insertions(content, virus_name, proteome_file))
    mutations_data.extend(_parse_frameshifts(content, virus_name, proteome_file))
    
    # Write to CSV
    _write_mutations_to_csv(mutations_data, csv_file)
    
    print(f"Mutation summary saved to: {csv_file}")
    print(f"Total mutations processed: {len(mutations_data)}")


def _parse_point_mutations(content: str, virus_name: str, proteome_file: Optional[str] = None) -> List[Dict]:
    """Parse point mutations from the content."""
    mutations = []
    pattern = r'(missense|silent)\s+(\d+(?::\d+)?)\s+([ATCG]+)->([ATCG]+)\s+\[(.*?)\]'
    
    for match in re.finditer(pattern, content, re.MULTILINE):
        mutation_type = match.group(1)
        nt_position = match.group(2)
        nt_ref = match.group(3)
        nt_alt = match.group(4)
        protein_mutations_str = match.group(5)
        
        # Parse protein mutations
        protein_mutations = _parse_protein_mutations(protein_mutations_str)
        
        for protein_name, mutation_desc in protein_mutations:
            # Validate protein assignment
            # Extract position for validation
            pos_match = re.search(r'(\d+)', nt_position)
            if pos_match:
                pos = int(pos_match.group(1))
                corrected_protein = _validate_protein_assignment(virus_name, pos, protein_name, proteome_file)
            else:
                corrected_protein = protein_name
            if corrected_protein == 'None-CDS':
                continue  # Skip if position is in non-coding region
            
            # Determine mutation type
            actual_type = _determine_mutation_type(mutation_desc, mutation_type)
            
            # Format NT mutation
            if ':' in nt_position:
                # Handle range positions
                nt_desc = f"{nt_position}{nt_ref}->{nt_alt}"
            else:
                nt_desc = f"{nt_ref}{nt_position}{nt_alt}"
            
            mutations.append({
                'protein_name': corrected_protein,
                'protein_mutation': mutation_desc,
                'aa_mutation_type': actual_type,
                'nt_mutation': nt_desc,
                'caused_by_nt_mutation_type': 'NT changes'
            })
    
    return mutations


def _parse_deletions(content: str, virus_name: str, proteome_file: Optional[str] = None) -> List[Dict]:
    """Parse deletions from the content."""
    mutations = []
    pattern = r'deletion\s+(\d+):(\d+)\s+([ATCG]+)\s+\[(.*)\]'
    
    for match in re.finditer(pattern, content, re.MULTILINE):
        start_pos = int(match.group(1))
        end_pos = int(match.group(2))
        deleted_seq = match.group(3)
        protein_mutations_str = match.group(4)
        
        # Parse protein mutations
        protein_mutations = _parse_protein_mutations(protein_mutations_str)
        
        for protein_name, mutation_desc in protein_mutations:
            # Validate protein assignment
            corrected_protein = _validate_protein_assignment(virus_name, start_pos, protein_name, proteome_file)
            if corrected_protein == 'None-CDS':
                continue
            
            # Determine mutation type
            actual_type = _determine_mutation_type(mutation_desc, 'deletion')
            
            # Format NT mutation
            nt_desc = f"{start_pos}:{end_pos}{deleted_seq}"
            
            mutations.append({
                'protein_name': corrected_protein,
                'protein_mutation': mutation_desc,
                'aa_mutation_type': actual_type,
                'nt_mutation': nt_desc,
                'caused_by_nt_mutation_type': 'deletion'
            })
    
    return mutations


def _parse_frameshifts(content: str, virus_name: str, proteome_file: Optional[str] = None) -> List[Dict]:
    """Parse frameshifts from the content."""
    mutations = []
    pattern = r'frameshift\s+(\d+)\s+([ATCG]+)\s+\[(.*)\]'
    
    for match in re.finditer(pattern, content, re.MULTILINE):
        nt_position = match.group(1)
        nt_seq = match.group(2)
        protein_mutations_str = match.group(3)
        
        # Parse protein mutations
        protein_mutations = _parse_protein_mutations(protein_mutations_str)
        
        for protein_name, mutation_desc in protein_mutations:
            # Validate protein assignment
            pos = int(nt_position)
            corrected_protein = _validate_protein_assignment(virus_name, pos, protein_name, proteome_file)
            if corrected_protein == 'None-CDS':
                continue
            
            # Determine mutation type
            actual_type = _determine_mutation_type(mutation_desc, 'frameshift')
            
            # Format NT mutation
            nt_desc = f"{nt_position}{nt_seq}"
            
            mutations.append({
                'protein_name': corrected_protein,
                'protein_mutation': mutation_desc,
                'aa_mutation_type': actual_type,
                'nt_mutation': nt_desc,
                'caused_by_nt_mutation_type': 'frameshift'
            })
    
    return mutations


def _parse_insertions(content: str, virus_name: str, proteome_file: Optional[str] = None) -> List[Dict]:
    """Parse insertions from the content."""
    mutations = []
    pattern = r'insertion\s+(\d+):(\d+)\s+->([ATCG]+)\s+\[(.*?)\]'
    
    for match in re.finditer(pattern, content, re.MULTILINE):
        start_pos = int(match.group(1))
        end_pos = int(match.group(2))
        inserted_seq = match.group(3)
        protein_mutations_str = match.group(4)
        
        # Parse protein mutations
        protein_mutations = _parse_protein_mutations(protein_mutations_str)
        
        for protein_name, mutation_desc in protein_mutations:
            # Validate protein assignment
            corrected_protein = _validate_protein_assignment(virus_name, start_pos, protein_name, proteome_file)
            if corrected_protein == 'None-CDS':
                continue
            
            # Determine mutation type
            actual_type = _determine_mutation_type(mutation_desc, 'insertion')
            
            # Format NT mutation
            nt_desc = f"{start_pos}:{end_pos}->{inserted_seq}"
            
            mutations.append({
                'protein_name': corrected_protein,
                'protein_mutation': mutation_desc,
                'aa_mutation_type': actual_type,
                'nt_mutation': nt_desc,
                'caused_by_nt_mutation_type': 'insertion'
            })
    
    return mutations


def _parse_protein_mutations(protein_mutations_str: str) -> List[Tuple[str, str]]:
    """Parse protein mutations from the string format."""
    mutations = []
    
    # Handle list format: [{'protein': 'name', 'mutation': 'desc'}]
    if protein_mutations_str.startswith('[') and protein_mutations_str.endswith(']'):
        # Extract individual mutation dictionaries
        pattern = r"\{'protein':\s*'([^']+)',\s*'mutation':\s*'([^']+)'\}"
        for match in re.finditer(pattern, protein_mutations_str):
            protein_name = match.group(1)
            mutation_desc = match.group(2)
            
            # Handle comma-separated mutations
            if ',' in mutation_desc and not mutation_desc.endswith('del'):
                for mut in mutation_desc.split(','):
                    mutations.append((protein_name, mut.strip()))
            else:
                mutations.append((protein_name, mutation_desc))
    else:
        # Handle single dictionary format: {'protein': 'name', 'mutation': 'desc'}
        # Also handle list format: {'protein': 'name', 'mutation': ['item1', 'item2']}
        
        # First try list format (with optional fields between protein and mutation)
        list_pattern = r"'protein':\s*'([^']+)',\s*.*?'mutation':\s*\[(.*?)\]"
        list_match = re.search(list_pattern, protein_mutations_str)
        if list_match:
            protein_name = list_match.group(1)
            list_content = list_match.group(2)
            # Parse list items
            list_items = [item.strip().strip("'\"") for item in list_content.split(',')]
            for item in list_items:
                mutations.append((protein_name, item))
        else:
            # Try single mutation format
            pattern = r"'protein':\s*'([^']+)',\s*'mutation':\s*'([^']*)'"
            match = re.search(pattern, protein_mutations_str)
            if match:
                protein_name = match.group(1)
                mutation_desc = match.group(2)
                
                # Skip empty mutations
                if mutation_desc:
                    # Handle comma-separated mutations
                    if ',' in mutation_desc and not mutation_desc.endswith('del'):
                        for mut in mutation_desc.split(','):
                            mutations.append((protein_name, mut.strip()))
                    else:
                        mutations.append((protein_name, mutation_desc))
    
    return mutations


def _determine_mutation_type(mutation_desc: str, original_type: str) -> str:
    """Determine the actual mutation type based on the mutation description."""
    if mutation_desc.endswith('del'):
        return 'deletion'
    elif mutation_desc.endswith('_'):  # Stop codon mutation (e.g., K2_)
        return 'stop'
    elif len(mutation_desc) >= 3 and mutation_desc[0].isalpha() and mutation_desc[-1].isalpha() and mutation_desc[1:-1].isdigit():
        # Check if it's a silent mutation (same amino acid)
        if mutation_desc[0] == mutation_desc[-1]:
            return 'silent'
        else:
            return 'missense'
    elif '->' in mutation_desc:
        return 'missense'
    else:
        return original_type


def _validate_protein_assignment(virus_name: str, nt_position: int, assigned_protein: str, proteome_file: Optional[str] = None) -> str:
    """Validate if a given nucleotide position falls within the genomic coordinates of an assigned protein."""
    if proteome_file and os.path.exists(proteome_file):
        try:
            from Bio import SeqIO
            protein_boundaries = {}
            
            # Parse proteome FASTA file to extract protein boundaries
            for record in SeqIO.parse(proteome_file, "fasta"):
                # Extract coordinates from FASTA header
                # Format varies, but typically contains location info
                header = record.description
                
                # Look for location patterns in the header
                # Format: ">YP_009725297.1|leader_nsp1|266..805" or ">YP_009725307.1|RNA-dependent-polymerase|join(13442..13468,13468..16236)"
                import re
                
                # Handle join() format first
                join_match = re.search(r'\|join\((\d+)\.\.(\d+),(\d+)\.\.(\d+)\)$', header)
                if join_match:
                    # For join format, use the full range from first start to second end
                    start1, end1, start2, end2 = map(int, join_match.groups())
                    start, end = start1, end2  # Full range coverage
                    # Extract protein name
                    protein_name_match = re.search(r'\|([^|]+)\|join\(', header)
                    if protein_name_match:
                        protein_name = protein_name_match.group(1)
                        protein_boundaries[protein_name] = (start, end)
                else:
                    # Handle simple format
                    location_match = re.search(r'\|(\d+)\.\.(\d+)$', header)
                    if location_match:
                        start, end = int(location_match.group(1)), int(location_match.group(2))
                        # Extract protein name from header (e.g., "leader_nsp1" from "YP_009725297.1|leader_nsp1|266..805")
                        protein_name_match = re.search(r'\|([^|]+)\|\d+\.\.\d+$', header)
                        if protein_name_match:
                            protein_name = protein_name_match.group(1)
                            protein_boundaries[protein_name] = (start, end)
            
            # Check if the nucleotide position falls within any protein's boundaries
            for protein, (start, end) in protein_boundaries.items():
                if start <= nt_position <= end:
                    return protein
            
            # If position is not within any protein, return None-CDS
            return 'None-CDS'
            
        except Exception as e:
            print(f"Warning: Could not parse proteome file {proteome_file}: {e}")
            return assigned_protein
    
    # If no proteome file or parsing failed, return the assigned protein as-is
    return assigned_protein


def _write_mutations_to_csv(mutations_data: List[Dict], csv_file: str) -> None:
    """Write mutations data to CSV file."""
    with open(csv_file, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(['Protein Name', 'Protein Mutations', 'AA Mutation Type', 'NT Mutation', 'Caused by NT Mutation Type'])
        
        for mutation in mutations_data:
            writer.writerow([
                mutation['protein_name'],
                mutation['protein_mutation'],
                mutation['aa_mutation_type'],
                mutation['nt_mutation'],
                mutation['caused_by_nt_mutation_type']
            ])
