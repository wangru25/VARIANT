# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-08-28 15:00:00
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-29 10:01:23
Email: rw3594@nyu.edu
FilePath: /VARIANT/src/visualization/figure1_mutation_analysis.py
Description: Simple combined genome organization and protein mutation analysis visualization
'''

import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
import numpy as np
from Bio import SeqIO
import re
import json
import pandas as pd
from typing import Dict, List, Tuple, Optional, Union
import os
import argparse
import glob

# Set default template
pio.templates.default = "simple_white"

class SimpleCombinedAnalyzer:
    """Simple combined genome organization and protein mutation analysis visualizations."""
    
    def __init__(self, virus_name: str, sample_id: str = None):
        """
        Initialize the simple combined analyzer.

        Args:
            virus_name (str): Name of the virus (e.g., 'SARS-CoV-2', 'HIV-1')
            sample_id (str, optional): Specific sample ID to analyze
        """
        self.virus_name = virus_name
        self.sample_id = sample_id
        
        # Auto-detect file paths
        self.proteome_path = self._auto_detect_proteome_path()
        self.mutation_csv_path = self._auto_detect_mutation_csv_path()
        self.reference_genome_path = self._auto_detect_reference_genome_path()
        
        # Load data
        self.proteins = self._load_proteins()
        self.mutations = self._parse_mutations_from_csv()
        self.genome_length = self._get_genome_length()
        self.genome_organization = self._create_genome_organization()
        
    def _auto_detect_proteome_path(self) -> str:
        """Auto-detect proteome file path."""
        # First check if virus has segment-based structure
        segment_dirs = glob.glob(f"data/{self.virus_name}/segment_*")
        
        if segment_dirs:
            # Segment-based virus - use first segment
            segment_name = os.path.basename(sorted(segment_dirs)[0])
            possible_paths = [
                f"data/{self.virus_name}/{segment_name}/refs/{segment_name}_proteome.fasta",
                f"data/{self.virus_name}/{segment_name}/refs/*proteome*.fasta"
            ]
        else:
            # Regular virus structure
            possible_paths = [
                f"data/{self.virus_name}/refs/{self.virus_name}_proteome.fasta",
                f"data/{self.virus_name}/refs/{self.virus_name}_proteome_all_possible.fasta",
                f"data/{self.virus_name}/refs/*proteome*.fasta"
            ]
        
        for path in possible_paths:
            if '*' in path:
                matches = glob.glob(path)
                if matches:
                    return matches[0]
            elif os.path.exists(path):
                return path
        
        raise FileNotFoundError(f"Proteome file not found for {self.virus_name}")
    
    def _auto_detect_mutation_csv_path(self) -> str:
        """Auto-detect mutation CSV file path."""
        # First check if virus has segment-based structure
        segment_dirs = glob.glob(f"result/{self.virus_name}/segment_*")
        
        if segment_dirs:
            # Segment-based virus - use first segment for panel 1 (genome organization)
            segment_name = os.path.basename(sorted(segment_dirs)[0])
            result_dir = f"result/{self.virus_name}/{segment_name}/"
        else:
            # Regular virus structure
            result_dir = f"result/{self.virus_name}/"
        
        if not os.path.exists(result_dir):
            raise FileNotFoundError(f"Result directory not found: {result_dir}")
        
        if self.sample_id:
            specific_path = os.path.join(result_dir, f"{self.sample_id}_mutation_summary.csv")
            if os.path.exists(specific_path):
                return specific_path
        
        csv_files = glob.glob(os.path.join(result_dir, "*_mutation_summary.csv"))
        if not csv_files:
            raise FileNotFoundError(f"No mutation summary CSV files found in {result_dir}")
        
        return csv_files[0]
    
    def _get_all_segment_data(self) -> Dict[str, Dict]:
        """Get mutation data from all segments for multi-segment viruses."""
        segment_dirs = glob.glob(f"result/{self.virus_name}/segment_*")
        if not segment_dirs:
            return {}
        
        all_segment_data = {}
        
        for segment_dir in sorted(segment_dirs):
            segment_name = os.path.basename(segment_dir)
            csv_files = glob.glob(os.path.join(segment_dir, "*_mutation_summary.csv"))
            
            if csv_files:
                csv_path = csv_files[0]
                try:
                    # Parse mutations from this segment
                    mutations = self._parse_mutations_from_csv_path(csv_path)
                    
                    # Load proteome for this segment
                    proteome_path = f"data/{self.virus_name}/{segment_name}/refs/{segment_name}_proteome.fasta"
                    if os.path.exists(proteome_path):
                        proteins = self._load_proteins_from_path(proteome_path)
                        
                        # Calculate conservation data for this segment
                        segment_data = self._calculate_protein_conservation_for_segment(mutations, proteins)
                        all_segment_data[segment_name] = segment_data
                        
                except Exception as e:
                    print(f"Warning: Could not process segment {segment_name}: {e}")
        
        return all_segment_data
    
    def _load_proteins_from_path(self, proteome_path: str) -> Dict[str, Dict]:
        """Load protein sequences from a specific FASTA file path."""
        proteins = {}
        
        with open(proteome_path, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                header_parts = record.description.split('|')
                if len(header_parts) >= 2:
                    protein_name = header_parts[1]
                    sequence = str(record.seq)
                    
                    position_range = None
                    if len(header_parts) >= 3:
                        pos_str = header_parts[2]
                        if '..' in pos_str and not pos_str.startswith('join'):
                            try:
                                parts = pos_str.split('..')
                                if len(parts) == 2:
                                    position_range = (int(parts[0]), int(parts[1]))
                            except (ValueError, IndexError):
                                pass
                    
                    proteins[protein_name] = {
                        'sequence': sequence,
                        'length': len(sequence),
                        'position_range': position_range,
                        'id': header_parts[0]
                    }
        
        return proteins
    
    def _calculate_protein_conservation_for_segment(self, mutations: List[Dict], proteins: Dict[str, Dict]) -> Dict[str, Dict]:
        """Calculate conservation metrics for a specific segment."""
        conservation_data = {}
        
        protein_mutations = {}
        for mutation in mutations:
            protein_name = mutation['protein']
            if protein_name not in protein_mutations:
                protein_mutations[protein_name] = []
            protein_mutations[protein_name].append(mutation)
        
        for protein_name, protein_data in proteins.items():
            if protein_name not in protein_mutations:
                protein_mutations[protein_name] = []
            
            mutations = protein_mutations[protein_name]
            protein_length = protein_data['length']
            genomic_length = protein_data['position_range'][1] - protein_data['position_range'][0] if protein_data['position_range'] else 0
            
            # Filter out silent mutations for protein-level mutation rate
            protein_level_mutations = [mut for mut in mutations if mut['type'] != 'silent']
            protein_mutation_count = len(protein_level_mutations)
            
            # Total mutation count (including silent) for other metrics
            total_mutation_count = len(mutations)
            
            # Protein-level mutation rate (excluding silent mutations)
            mutation_rate = protein_mutation_count / protein_length if protein_length > 0 else 0
            mutation_density = total_mutation_count / genomic_length if genomic_length > 0 else 0
            conservation_ratio = (protein_length - protein_mutation_count) / protein_length if protein_length > 0 else 1.0
            
            mutation_details = {}
            for mutation in mutations:
                mut_type = mutation['type']
                if mut_type not in mutation_details:
                    mutation_details[mut_type] = []
                mutation_details[mut_type].append(mutation['mutation'])

            conservation_data[protein_name] = {
                'mutation_count': protein_mutation_count,
                'total_mutation_count': total_mutation_count,
                'mutation_rate': mutation_rate,
                'mutation_density': mutation_density,
                'conservation_ratio': conservation_ratio,
                'protein_length': protein_length,
                'genomic_length': genomic_length,
                'mutation_details': mutation_details
            }
        
        return conservation_data
    
    def _auto_detect_reference_genome_path(self) -> Optional[str]:
        """Auto-detect reference genome file path."""
        # First check if virus has segment-based structure
        segment_dirs = glob.glob(f"data/{self.virus_name}/segment_*")
        
        if segment_dirs:
            # Segment-based virus - use first segment
            segment_name = os.path.basename(sorted(segment_dirs)[0])
            possible_paths = [
                f"data/{self.virus_name}/{segment_name}/refs/NC_*.fasta",
                f"data/{self.virus_name}/{segment_name}/refs/*.fasta"
            ]
        else:
            # Regular virus structure
            possible_paths = [
                f"data/{self.virus_name}/refs/NC_*.fasta",
                f"data/{self.virus_name}/refs/*.fasta"
            ]
        
        for path in possible_paths:
            matches = glob.glob(path)
            if matches:
                return matches[0]
        
        return None
        
    def _load_proteins(self) -> Dict[str, Dict]:
        """Load protein sequences from FASTA file."""
        proteins = {}
        
        with open(self.proteome_path, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                header_parts = record.description.split('|')
                if len(header_parts) >= 2:
                    protein_name = header_parts[1]
                    sequence = str(record.seq)
                    
                    position_range = None
                    if len(header_parts) >= 3:
                        pos_str = header_parts[2]
                        if '..' in pos_str and not pos_str.startswith('join'):
                            try:
                                parts = pos_str.split('..')
                                if len(parts) == 2:
                                    position_range = (int(parts[0]), int(parts[1]))
                            except (ValueError, IndexError):
                                pass
                        elif pos_str.startswith('join'):
                            try:
                                start_bracket = pos_str.find('(')
                                end_bracket = pos_str.find(')')
                                if start_bracket != -1 and end_bracket != -1:
                                    join_content = pos_str[start_bracket+1:end_bracket]
                                    first_range = join_content.split(',')[0]
                                    if '..' in first_range:
                                        parts = first_range.split('..')
                                        if len(parts) == 2:
                                            start_pos = int(parts[0])
                                            last_range = join_content.split(',')[-1]
                                            if '..' in last_range:
                                                end_parts = last_range.split('..')
                                                if len(end_parts) == 2:
                                                    end_pos = int(end_parts[1])
                                                    position_range = (start_pos, end_pos)
                            except (ValueError, IndexError):
                                pass
                    
                    proteins[protein_name] = {
                        'sequence': sequence,
                        'length': len(sequence),
                        'position_range': position_range,
                        'id': header_parts[0]
                    }
        
        return proteins
    
    def _parse_mutations_from_csv(self) -> List[Dict]:
        """Parse mutation data from CSV file."""
        return self._parse_mutations_from_csv_path(self.mutation_csv_path)
    
    def _parse_mutations_from_csv_path(self, csv_path: str) -> List[Dict]:
        """Parse mutation data from a specific CSV file path."""
        mutations = []
        
        try:
            df = pd.read_csv(csv_path)
            
            for _, row in df.iterrows():
                protein_name = row['Protein Name']
                protein_mutation = row['Protein Mutations']
                aa_mutation_type = row['AA Mutation Type']
                nt_mutation = row['NT Mutation']
                
                if protein_name in ["Invalid protein sequence", "None-CDS"]:
                    continue
                
                position = self._parse_nt_position(nt_mutation)
                if position is None:
                    continue
                
                mutations.append({
                    'type': aa_mutation_type,
                    'position': position,
                    'protein': protein_name,
                    'mutation': protein_mutation,
                    'nt_mutation': nt_mutation,
                    'raw_line': f"{aa_mutation_type} {protein_mutation} in {protein_name}"
                })
                
        except Exception as e:
            print(f"Error reading CSV file {csv_path}: {e}")
            return []
        
        return mutations
    
    def _parse_nt_position(self, nt_mutation: str) -> Optional[Union[int, Tuple[int, int]]]:
        """Parse nucleotide position from NT Mutation column."""
        try:
            if ':' in nt_mutation:
                parts = nt_mutation.split(':')
                if len(parts) >= 2:
                    start_pos = int(parts[0])
                    end_part = parts[1]
                    end_match = re.search(r'^(\d+)', end_part)
                    if end_match:
                        end_pos = int(end_match.group(1))
                        return (start_pos, end_pos)
                    else:
                        return start_pos
            
            match = re.search(r'[ACGT](\d+)[ACGT]', nt_mutation)
            if match:
                return int(match.group(1))
            
            numbers = re.findall(r'\d+', nt_mutation)
            if numbers:
                return int(numbers[0])
                
        except (ValueError, IndexError) as e:
            print(f"Warning: Could not parse position from {nt_mutation}: {e}")
        
        return None
    
    def _get_genome_length(self) -> int:
        """Get genome length from reference genome or estimate from protein positions."""
        if self.reference_genome_path and os.path.exists(self.reference_genome_path):
            try:
                with open(self.reference_genome_path, 'r') as f:
                    for record in SeqIO.parse(f, 'fasta'):
                        return len(record.seq)
            except Exception as e:
                print(f"Warning: Could not read reference genome: {e}")
        
        max_position = 0
        for protein_data in self.proteins.values():
            if protein_data['position_range']:
                max_position = max(max_position, protein_data['position_range'][1])
        
        if max_position == 0:
            max_position = 30000
        
        return max_position
    
    def _create_genome_organization(self) -> Dict:
        """Create generic genome organization data structure for any virus."""
        # Check if this is a multi-segment virus
        segment_dirs = glob.glob(f"data/{self.virus_name}/segment_*")
        
        if segment_dirs:
            # Multi-segment virus - combine all segments
            genome_org = {
                'genes': [],
                'regulatory_regions': []
            }
            
            current_position = 1  # Start at position 1
            for segment_dir in sorted(segment_dirs):
                segment_name = os.path.basename(segment_dir)
                proteome_path = f"data/{self.virus_name}/{segment_name}/refs/{segment_name}_proteome.fasta"
                
                if os.path.exists(proteome_path):
                    segment_proteins = self._load_proteins_from_path(proteome_path)
                    
                    # Get segment mutations
                    csv_files = glob.glob(f"result/{self.virus_name}/{segment_name}/*_mutation_summary.csv")
                    if csv_files:
                        segment_mutations = self._parse_mutations_from_csv_path(csv_files[0])
                        
                        # Add genes from this segment
                        sorted_proteins = []
                        for protein_name, protein_data in segment_proteins.items():
                            if protein_data['position_range']:
                                sorted_proteins.append((protein_name, protein_data))
                        
                        sorted_proteins.sort(key=lambda x: x[1]['position_range'][0])
                        
                        for protein_name, protein_data in sorted_proteins:
                            original_start, original_end = protein_data['position_range']
                            protein_length = original_end - original_start
                            
                            gene_info = {
                                'name': f"{segment_name}_{protein_name}",
                                'start': current_position,  # Sequential positioning
                                'end': current_position + protein_length,
                                'length': protein_length,
                                'protein_length': protein_data['length'],
                                'color': self._get_protein_color(f"{segment_name}_{protein_name}"),
                                'segment': segment_name,
                                'original_start': original_start,
                                'original_end': original_end
                            }
                            
                            genome_org['genes'].append(gene_info)
                            
                            # Update position for next protein: previous end + 1
                            current_position = gene_info['end'] + 1
        else:
            # Single-segment virus - use existing logic
            genome_org = {
                'genes': [],
                'regulatory_regions': []
            }
            
            sorted_proteins = []
            for protein_name, protein_data in self.proteins.items():
                if protein_data['position_range']:
                    sorted_proteins.append((protein_name, protein_data))
            
            sorted_proteins.sort(key=lambda x: x[1]['position_range'][0])
            
            for protein_name, protein_data in sorted_proteins:
                start_pos, end_pos = protein_data['position_range']
                
                gene_info = {
                    'name': protein_name,
                    'start': start_pos,
                    'end': end_pos,
                    'length': end_pos - start_pos,
                    'protein_length': protein_data['length'],
                    'color': self._get_protein_color(protein_name)
                }
                
                genome_org['genes'].append(gene_info)
            
            if self.genome_length > 0:
                utr_size = max(100, int(self.genome_length * 0.05))
                genome_org['regulatory_regions'] = [
                    {'name': '5\' UTR', 'start': 1, 'end': utr_size, 'type': 'UTR'},
                    {'name': '3\' UTR', 'start': self.genome_length - utr_size, 'end': self.genome_length, 'type': 'UTR'}
                ]
        
        return genome_org
    
    def _get_protein_color(self, protein_name: str) -> str:
        """Get color for protein using Plotly continuous color palette."""
        import plotly.colors as pc
        
        # Get the number of proteins to determine how many colors we need
        num_proteins = len(self.genome_organization['genes']) if hasattr(self, 'genome_organization') else 10
        
        # Sample colors based on the number of proteins needed
        # Use a larger scale to get better color distribution, then sample the middle portion
        scale_size = max(50, num_proteins * 2)  # Ensure we have enough colors to sample from
        start_offset = int(scale_size * 0.15)  # Skip the lightest 15%
        end_offset = int(scale_size * 0.15)    # Skip the darkest 15%
        
        colors = pc.sample_colorscale('twilight', scale_size)[start_offset:scale_size-end_offset]
        
        # If we still have more colors than needed, sample evenly
        if len(colors) > num_proteins:
            step = len(colors) / num_proteins
            colors = [colors[int(i * step)] for i in range(num_proteins)]
        
        total = sum(ord(c) for c in protein_name)
        color_index = total % len(colors)
        return colors[color_index]
    
    def _get_mutation_color(self, mutation_type: str) -> str:
        """Get color for mutation type."""
        color_map = {
            'missense': '#FF0000',      # Bright red
            'silent': '#778DA9',        # Grey
            'deletion': '#14746F',      # Dark green
            'insertion': '#8A2BE2',     # Blue violet
            'nonsense': '#4361EE',      # Deep blue
            'frameshift': '#8B4513',    # Saddle brown
        }
        return color_map.get(mutation_type, '#1f77b4')
    
    def _add_mutation_legend(self, fig):
        """Add mutation type legend in top left corner."""
        mutation_types_in_data = set()
        for mutation in self.mutations:
            mutation_types_in_data.add(mutation['type'])
        
        mutation_type_info = {
            'missense': ('Missense', '#FF0000', 'Amino acid change'),
            'silent': ('Silent', '#778DA9', 'No amino acid change'),
            'deletion': ('Deletion', '#14746F', 'Nucleotide deletion'),
            'insertion': ('Insertion', '#8A2BE2', 'Nucleotide insertion'),
            'nonsense': ('Nonsense', '#4361EE', 'Stop codon introduction'),
            'frameshift': ('Frameshift', '#8B4513', 'Reading frame shift'),
        }
        
        legend_text = "<b>Mutation Types:</b><br>"
        main_types = ['missense', 'deletion', 'insertion', 'nonsense']

        for mut_type in main_types:
            if mut_type in mutation_type_info:
                display_name, color, description = mutation_type_info[mut_type]
                legend_text += f"<span style='color:{color};'>|</span> <span style='color:{color};'>{display_name}</span>: <span style='color:{color};'>{description}</span><br>"
        
        for mut_type in sorted(mutation_types_in_data):
            if mut_type not in main_types and mut_type in mutation_type_info:
                display_name, color, description = mutation_type_info[mut_type]
                legend_text += f"<span style='color:{color};'>|</span> <span style='color:{color};'>{display_name}</span>: <span style='color:{color};'>{description}</span><br>"
        
        fig.add_annotation(
            x=0.02,
            y=0.98,
            xref='paper',
            yref='paper',
            text=legend_text,
            showarrow=False,
            bordercolor='black',
            borderwidth=1,
            font=dict(size=10),
            align='left'
        )
    
    def _format_mutation_breakdown(self, mutation_details):
        """Format mutation details for hover display."""
        if not mutation_details:
            return "No mutations"
        
        breakdown = []
        for mutation_type, mutations in mutation_details.items():
            if mutations:
                # Show all mutations with compact formatting
                mutations_display = ', '.join(mutations)
                # Split into multiple lines if too long
                if len(mutations_display) > 60:
                    parts = mutations_display.split(', ')
                    lines = []
                    current_line = ""
                    for part in parts:
                        if len(current_line + part) > 60:
                            if current_line:
                                lines.append(current_line.rstrip())
                            current_line = part + ", "
                        else:
                            current_line += part + ", "
                    if current_line:
                        lines.append(current_line.rstrip(', '))
                    mutations_display = '<br>'.join(lines)
                
                breakdown.append(f"<b>{mutation_type}</b>: {mutations_display}")
        
        return '<br>'.join(breakdown) if breakdown else "No mutations"
    
    def _calculate_protein_conservation(self) -> Dict[str, Dict]:
        """Calculate conservation metrics for each protein."""
        # Check if this is a multi-segment virus
        segment_dirs = glob.glob(f"result/{self.virus_name}/segment_*")
        
        if segment_dirs:
            # Multi-segment virus - combine data from all segments
            all_segment_data = self._get_all_segment_data()
            conservation_data = {}
            
            for segment_name, segment_data in all_segment_data.items():
                for protein_name, protein_data in segment_data.items():
                    # Add segment prefix to protein name to avoid conflicts
                    full_protein_name = f"{segment_name}_{protein_name}"
                    conservation_data[full_protein_name] = protein_data
            
            return conservation_data
        else:
            # Regular virus - use existing logic
            conservation_data = {}
            
            protein_mutations = {}
            for mutation in self.mutations:
                protein_name = mutation['protein']
                if protein_name not in protein_mutations:
                    protein_mutations[protein_name] = []
                protein_mutations[protein_name].append(mutation)
            
            for protein_name, protein_data in self.proteins.items():
                if protein_name not in protein_mutations:
                    protein_mutations[protein_name] = []
                
                mutations = protein_mutations[protein_name]
                protein_length = protein_data['length']
                genomic_length = protein_data['position_range'][1] - protein_data['position_range'][0] if protein_data['position_range'] else 0
                
                # Filter out silent mutations for protein-level mutation rate
                protein_level_mutations = [mut for mut in mutations if mut['type'] != 'silent']
                protein_mutation_count = len(protein_level_mutations)
                
                # Total mutation count (including silent) for other metrics
                total_mutation_count = len(mutations)
                
                # Protein-level mutation rate (excluding silent mutations)
                mutation_rate = protein_mutation_count / protein_length if protein_length > 0 else 0
                mutation_density = total_mutation_count / genomic_length if genomic_length > 0 else 0
                conservation_ratio = (protein_length - protein_mutation_count) / protein_length if protein_length > 0 else 1.0
                
                mutation_details = {}
                for mutation in mutations:
                    mut_type = mutation['type']
                    if mut_type not in mutation_details:
                        mutation_details[mut_type] = []
                    mutation_details[mut_type].append(mutation['mutation'])

                conservation_data[protein_name] = {
                    'mutation_count': protein_mutation_count,  # Protein-level mutations (excluding silent)
                    'total_mutation_count': total_mutation_count,  # All mutations including silent
                    'mutation_rate': mutation_rate,
                    'mutation_density': mutation_density,
                    'conservation_ratio': conservation_ratio,
                    'protein_length': protein_length,
                    'genomic_length': genomic_length,
                    'mutation_details': mutation_details
                }
            
            return conservation_data
    
    def create_combined_chart(self, output_path: str = None):
        """Create combined genome organization and protein mutation analysis visualization."""
        
        if output_path is None:
            genome_id = os.path.basename(self.mutation_csv_path).replace('_mutation_summary.csv', '')
            output_dir = f"imgs/visualizations/{self.virus_name}/"
            os.makedirs(output_dir, exist_ok=True)
            output_path = f"{output_dir}{genome_id}_combined_analysis.html"
        
        # Calculate dynamic panel heights based on number of proteins
        num_proteins = len(self.proteins)
        
        # Adjust panel heights based on protein count
        if num_proteins <= 10:
            # For viruses with few proteins (like ZaireEbola), give more space to panels 2 and 3
            row_heights = [0.3, 0.35, 0.35]
        elif num_proteins <= 20:
            # For medium viruses
            row_heights = [0.45, 0.275, 0.275]
        else:
            # For viruses with many proteins (like SARS-CoV-2), give more space to panel 1
            row_heights = [0.5, 0.25, 0.25]
        
        # Create subplots with dynamic heights
        fig = make_subplots(
            rows=3, cols=1,
            subplot_titles=('', '', ''),
            vertical_spacing=0.2,
            row_heights=row_heights,
            specs=[[{"secondary_y": False}], [{"secondary_y": False}], [{"secondary_y": False}]]
        )
        
        # Add mutation legend
        self._add_mutation_legend(fig)
        
        # Calculate conservation data
        conservation_data = self._calculate_protein_conservation()
        
        # For multi-segment viruses, ensure we use the combined data
        segment_dirs = glob.glob(f"result/{self.virus_name}/segment_*")
        if segment_dirs:
            print(f"Multi-segment virus detected: {len(segment_dirs)} segments")
            # Get combined data from all segments
            all_segment_data = self._get_all_segment_data()
            if all_segment_data:
                print(f"Successfully collected data from {len(all_segment_data)} segments")
                conservation_data = {}
                for segment_name, segment_data in all_segment_data.items():
                    for protein_name, protein_data in segment_data.items():
                        # Add segment prefix to protein name to avoid conflicts
                        full_protein_name = f"{segment_name}_{protein_name}"
                        conservation_data[full_protein_name] = protein_data
                print(f"Combined data contains {len(conservation_data)} proteins from all segments")
            else:
                print("Warning: Could not collect data from segments")
        else:
            print("Single-segment virus detected")
        
        # === PANEL 1: GENOME ORGANIZATION (from Figure 1) ===
        all_genes = []
        
        for i, gene in enumerate(self.genome_organization['genes']):
            bar_height = 0.6
            gap_size = 0.6
            y_pos = i * (bar_height + gap_size)
            
            fig.add_trace(go.Bar(
                x=[gene['end'] - gene['start']],
                y=[y_pos],
                orientation='h',
                name=gene['name'],
                marker_color=gene['color'],
                base=gene['start'],
                width=0.6,
                hovertemplate=f"<b>{gene['name']}</b><br>"
                            f"Position: {gene['start']}-{gene['end']}<br>"
                            f"Length: {gene['length']} bp<br>"
                            f"Protein: {gene['protein_length']} aa<extra></extra>",
                showlegend=False,
                opacity=0.4
            ), row=1, col=1)
            all_genes.append(gene)
        
        # Add regulatory regions
        for region in self.genome_organization['regulatory_regions']:
            mid_pos = (region['start'] + region['end']) / 2
            fig.add_annotation(
                x=mid_pos,
                y=1.04,
                xref='x',
                yref='paper',
                text=region['name'],
                showarrow=False,
                font=dict(size=10, color='gray'),
                textangle=-45
            )
        
        # Add mutations (from Figure 1)
        mutation_data = []
        
        # Check if this is a multi-segment virus
        segment_dirs = glob.glob(f"result/{self.virus_name}/segment_*")
        
        if segment_dirs:
            # Multi-segment virus - add mutations from all segments
            # Use the same offset calculation as in genome organization
            segment_offset = 0
            for segment_dir in sorted(segment_dirs):
                segment_name = os.path.basename(segment_dir)
                csv_files = glob.glob(os.path.join(segment_dir, "*_mutation_summary.csv"))
                
                if csv_files:
                    segment_mutations = self._parse_mutations_from_csv_path(csv_files[0])
                    
                    for mutation in segment_mutations:
                        protein_name = f"{segment_name}_{mutation['protein']}"
                        gene = None
                        for g in self.genome_organization['genes']:
                            if g['name'] == protein_name:
                                gene = g
                                break
                        
                        if gene is not None:
                            if isinstance(mutation['position'], tuple):
                                mut_pos = (mutation['position'][0] + mutation['position'][1]) / 2
                            else:
                                mut_pos = mutation['position']
                            
                            # Calculate relative position within the protein and map to new sequential position
                            relative_pos = mut_pos - gene['original_start']
                            gene_pos = gene['start'] + relative_pos
                            
                            gene_index = None
                            for i, g in enumerate(all_genes):
                                if g['name'] == gene['name']:
                                    gene_index = i
                                    break
                            
                            if gene_index is not None:
                                bar_height = 0.6
                                gap_size = 0.6
                                y_pos = gene_index * (bar_height + gap_size)
                                
                                mutation_data.append({
                                    'x': gene_pos,
                                    'y': y_pos,
                                    'type': mutation['type'],
                                    'mutation': mutation['mutation'],
                                    'nt_mutation': mutation['nt_mutation'],
                                    'position': mutation['position'],
                                    'gene': gene['name']
                                })
                    
                    # Update offset for next segment (same calculation as genome organization)
                    segment_proteins = self._load_proteins_from_path(f"data/{self.virus_name}/{segment_name}/refs/{segment_name}_proteome.fasta")
                    segment_length = max([p['position_range'][1] for p in segment_proteins.values() if p['position_range']])
                    segment_offset += segment_length + 1000  # Add gap between segments
        else:
            # Single-segment virus - use existing logic
            for mutation in self.mutations:
                protein_name = mutation['protein']
                gene = None
                for g in self.genome_organization['genes']:
                    if g['name'] == protein_name:
                        gene = g
                        break
                
                if gene is not None:
                    if isinstance(mutation['position'], tuple):
                        mut_pos = (mutation['position'][0] + mutation['position'][1]) / 2
                    else:
                        mut_pos = mutation['position']
                    
                    if gene['start'] <= mut_pos <= gene['end']:
                        relative_pos = mut_pos - gene['start']
                        gene_pos = gene['start'] + relative_pos
                        
                        gene_index = None
                        for i, g in enumerate(all_genes):
                            if g['name'] == gene['name']:
                                gene_index = i
                                break
                        
                        if gene_index is not None:
                            bar_height = 0.6
                            gap_size = 0.6
                            y_pos = gene_index * (bar_height + gap_size)
                            
                            mutation_data.append({
                                'x': gene_pos,
                                'y': y_pos,
                                'type': mutation['type'],
                                'mutation': mutation['mutation'],
                                'nt_mutation': mutation['nt_mutation'],
                                'position': mutation['position'],
                                'gene': gene['name']
                            })
        
        # Plot mutations as lines (from Figure 1)
        mutation_types = set(mut['type'] for mut in mutation_data)
        for mut_type in mutation_types:
            type_mutations = [mut for mut in mutation_data if mut['type'] == mut_type]
            
            # Add lines
            for mut in type_mutations:
                fig.add_shape(
                    type="line",
                    x0=mut['x'],
                    y0=mut['y'] - 0.3,
                    x1=mut['x'],
                    y1=mut['y'] + 0.3,
                    line=dict(
                        color=self._get_mutation_color(mut_type),
                        width=1.5
                    ),
                    opacity=1.0,
                    layer='above'
                )
            
            # Add hover
            x_positions = [mut['x'] for mut in type_mutations]
            y_positions = [mut['y'] for mut in type_mutations]
            
            colored_hover_texts = []
            mutation_color = self._get_mutation_color(mut_type)
            for mut in type_mutations:
                hover_text = f"<b>Mutation: {mut['type']}</b><br>Position: {mut['position']}<br>Nucleotide Change: {mut['nt_mutation']}<br>Protein Change: {mut['mutation']}<br>Gene: {mut['gene']}"
                colored_text = f'<span style="background-color:{mutation_color}; color:white; padding:5px; border-radius:3px;">{hover_text}</span>'
                colored_hover_texts.append(colored_text)
            
            fig.add_trace(go.Scatter(
                x=x_positions,
                y=y_positions,
                mode='markers',
                marker=dict(
                    size=1,
                    opacity=0,
                    color=mutation_color
                ),
                hovertemplate='%{text}<extra></extra>',
                text=colored_hover_texts,
                showlegend=False
            ), row=1, col=1)
        
        # Update y-axis for panel 1
        bar_height = 0.6
        gap_size = 0.6
        max_y_pos = (len(all_genes) - 1) * (bar_height + gap_size) + bar_height
        y_range = [-0.2, max_y_pos + 0.2]
        tick_positions = [i * (bar_height + gap_size) for i in range(len(all_genes))]
        tick_text_with_colors = []
        for g in all_genes:
            tick_text_with_colors.append(f'<span style="color:{g["color"]}">{g["name"]}</span>')
        
        fig.update_yaxes(
            tickmode='array',
            ticktext=tick_text_with_colors,
            tickvals=tick_positions,
            tickangle=0,
            tickfont=dict(size=11, family='Arial'),
            linecolor='black',
            linewidth=1,
            range=y_range,
            dtick=1,
            row=1, col=1
        )
        
        # === PANEL 2: GENE-LEVEL MUTATION ANALYSIS ===
        gene_proteins = []
        gene_rates = []
        gene_colors = []
        gene_mutation_counts = []
        
        for gene in self.genome_organization['genes']:
            conservation_info = conservation_data.get(gene['name'], {})
            gene_proteins.append(gene['name'])
            # Gene-level mutation rate includes all mutations (including silent)
            gene_rate = conservation_info.get('total_mutation_count', 0) / conservation_info.get('genomic_length', 1) if conservation_info.get('genomic_length', 0) > 0 else 0
            gene_rates.append(gene_rate)
            gene_mutation_counts.append(conservation_info.get('total_mutation_count', 0))
            gene_colors.append(gene['color'])
        
        # Add gene-level mutation bars (panel 2)
        fig.add_trace(go.Bar(
            x=gene_proteins,
            y=gene_rates,
            name='Gene Mutation Rate',
            marker_color=gene_colors,
            width=0.3,
            hovertemplate='<b>%{x}</b><br>' +
                         'Gene Mutation Rate: %{y:.3f}<br>' +
                         'Genomic Length: %{customdata[0]} bp<br>' +
                         'Total Mutations: %{customdata[1]}<extra></extra>',
            customdata=[[
                conservation_data.get(protein, {}).get('genomic_length', 0),
                conservation_data.get(protein, {}).get('total_mutation_count', 0)
            ] for protein in gene_proteins],
            opacity=0.3
        ), row=2, col=1)
        
        # Update x-axis for panel 2 (colored protein names)
        fig.update_xaxes(
            tickangle=45,
            tickmode='array',
            ticktext=[f'<span style="color:{color}">{protein}</span>' for protein, color in zip(gene_proteins, gene_colors)],
            tickvals=list(range(len(gene_proteins))),
            tickfont=dict(size=11, family='Arial'),
            linecolor='black',
            linewidth=1,
            title_text="",
            row=2, col=1
        )
        
        # === PANEL 3: PROTEIN MUTATION ANALYSIS (from Figure 2) ===
        conservation_proteins = []
        conservation_rates = []
        conservation_colors = []
        mutation_counts = []
        
        for gene in self.genome_organization['genes']:
            conservation_info = conservation_data.get(gene['name'], {})
            conservation_proteins.append(gene['name'])
            conservation_rates.append(conservation_info.get('mutation_rate', 0))
            mutation_counts.append(conservation_info.get('mutation_count', 0))
            conservation_colors.append(gene['color'])
        
        # Add protein mutation bars (panel 3)
        fig.add_trace(go.Bar(
            x=conservation_proteins,
            y=conservation_rates,
            name='Mutation Rate',
            marker_color=conservation_colors,
            width=0.3,
            hovertemplate='<b>%{x}</b><br>' +
                         'Protein Mutation Rate: %{y:.3f}<br>' +
                         'Protein Length: %{customdata[0]} aa<br>' +
                         'Protein Mutations (without silent mutations): %{customdata[1]}<br>' +
                         'Total Mutations (including silent mutations): %{customdata[2]}<br><br>' +
                         '<b>Mutations:</b><br>' +
                         '%{customdata[3]}<extra></extra>',
            customdata=[[
                conservation_data.get(protein, {}).get('protein_length', 0),
                conservation_data.get(protein, {}).get('mutation_count', 0),
                conservation_data.get(protein, {}).get('total_mutation_count', 0),
                self._format_mutation_breakdown(conservation_data.get(protein, {}).get('mutation_details', {}))
            ] for protein in conservation_proteins],
            opacity=0.4
        ), row=3, col=1)
        
        # Calculate simple height based on content
        num_proteins = len(all_genes)
        
        # Panel 1 height based on number of proteins
        baseline_y_range = 36.8
        current_y_range = (len(all_genes) - 1) * (0.6 + 0.6) + 0.6 + 0.4
        panel1_height = int(600 * (current_y_range / baseline_y_range))
        
        # Fixed heights for panels 2 and 3
        panel2_height = 300  # Gene-level analysis
        panel3_height = 300  # Protein analysis
        
        # Total height with spacing
        total_height = panel1_height + panel2_height + panel3_height + 400  # Extra space for margins and spacing
        
        # Update layout
        genome_id = os.path.basename(self.mutation_csv_path).replace('_mutation_summary.csv', '')
        fig.update_layout(
            title=dict(
                text=f"{self.virus_name} ({genome_id}) Genome and Protein Mutation Analysis",
                font=dict(size=16, family='Arial', color='black'),
                x=0.5,
                xanchor='center'
            ),
            width=1200,
            height=total_height,
            margin=dict(l=60, r=60, t=100, b=80),
            showlegend=False,
            plot_bgcolor='white',
            paper_bgcolor='white'
        )
        
        # Update axes for panel 1
        # Calculate the maximum x position from all genes
        max_x_pos = max([gene['end'] for gene in all_genes]) if all_genes else self.genome_length
        
        fig.update_xaxes(
            range=[0, max_x_pos],
            showgrid=False,
            linecolor='black',
            linewidth=1,
            tickfont=dict(size=11, family='Arial', color='black'),
            title_text="Genome Position (base pairs)",
            row=1, col=1
        )
        
        # Update y-axis for panel 2
        fig.update_yaxes(
            linecolor='black',
            linewidth=1,
            tickfont=dict(size=11, family='Arial', color='black'),
            title_text="Gene Mutation Rate",
            row=2, col=1
        )
        
        # Update x-axis for panel 3 (colored protein names)
        fig.update_xaxes(
            tickangle=45,
            tickmode='array',
            ticktext=[f'<span style="color:{color}">{protein}</span>' for protein, color in zip(gene_proteins, gene_colors)],
            tickvals=list(range(len(gene_proteins))),
            tickfont=dict(size=11, family='Arial'),
            linecolor='black',
            linewidth=1,
            title_text="",
            row=3, col=1
        )
        
        # Update x-axis for panel 3 (colored protein names)
        fig.update_xaxes(
            tickangle=45,
            tickmode='array',
            ticktext=[f'<span style="color:{color}">{protein}</span>' for protein, color in zip(conservation_proteins, conservation_colors)],
            tickvals=list(range(len(conservation_proteins))),
            tickfont=dict(size=11, family='Arial'),
            linecolor='black',
            linewidth=1,
            title_text="",
            row=3, col=1
        )
        
        # Update y-axis for panel 3
        fig.update_yaxes(
            linecolor='black',
            linewidth=1,
            tickfont=dict(size=11, family='Arial', color='black'),
            title_text="Protein Mutation Rate",
            row=3, col=1
        )

        # Save
        fig.write_html(output_path)
        fig.write_image(output_path.replace('.html', '.pdf'), width=1200, height=total_height)
        
        print(f"Combined analysis chart saved to: {output_path}")
        
        # Print summary
        print(f"\nCombined Analysis Summary:")
        print(f"  Virus: {self.virus_name}")
        print(f"  Genome ID: {genome_id}")
        
        # For multi-segment viruses, show combined protein count
        if segment_dirs:
            total_proteins = len(conservation_data)
            total_mutations = sum(data.get('total_mutation_count', 0) for data in conservation_data.values())
            print(f"  Proteins (all segments): {total_proteins}")
            print(f"  Mutations (all segments): {total_mutations}")
        else:
            print(f"  Proteins: {len(self.proteins)}")
            print(f"  Mutations: {len(self.mutations)}")
        
        print(f"  Genome length: {self.genome_length:,} bp")
        
        print(f"\nProtein Mutation Analysis Summary:")
        if segment_dirs:
            # For multi-segment viruses, show all proteins from all segments
            for protein_name, rate in zip(conservation_proteins, conservation_rates):
                protein_count = conservation_data.get(protein_name, {}).get('mutation_count', 0)
                total_count = conservation_data.get(protein_name, {}).get('total_mutation_count', 0)
                print(f"  {protein_name}: {rate:.3f} ({protein_count} protein mutations, {total_count} total mutations)")
        else:
            # For single-segment viruses
            for protein_name, rate in zip(conservation_proteins, conservation_rates):
                protein_count = conservation_data.get(protein_name, {}).get('mutation_count', 0)
                total_count = conservation_data.get(protein_name, {}).get('total_mutation_count', 0)
                print(f"  {protein_name}: {rate:.3f} ({protein_count} protein mutations, {total_count} total mutations)")
        
        print(f"\nGene-Level Mutation Analysis Summary:")
        if segment_dirs:
            # For multi-segment viruses, show all genes from all segments
            for protein_name, rate in zip(gene_proteins, gene_rates):
                total_count = conservation_data.get(protein_name, {}).get('total_mutation_count', 0)
                genomic_length = conservation_data.get(protein_name, {}).get('genomic_length', 0)
                print(f"  {protein_name}: {rate:.3f} ({total_count} mutations in {genomic_length} bp)")
        else:
            # For single-segment viruses
            for protein_name, rate in zip(gene_proteins, gene_rates):
                total_count = conservation_data.get(protein_name, {}).get('total_mutation_count', 0)
                genomic_length = conservation_data.get(protein_name, {}).get('genomic_length', 0)
                print(f"  {protein_name}: {rate:.3f} ({total_count} mutations in {genomic_length} bp)")
        
        return fig

def process_multiple_viruses():
    """Process multiple virus datasets with auto-detected paths."""
    data_dir = "data"
    if not os.path.exists(data_dir):
        print(f"Data directory not found: {data_dir}")
        return
    
    virus_dirs = [d for d in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir, d)) and not d.startswith('.')]
    
    print(f"Found virus directories: {virus_dirs}")
    
    for virus_name in virus_dirs:
        print(f"\n{'='*50}")
        print(f"Processing {virus_name}")
        print(f"{'='*50}")
        
        try:
            analyzer = SimpleCombinedAnalyzer(virus_name)
            analyzer.create_combined_chart()
            
            print(f"Successfully processed {virus_name}")
            
            # Check if multi-segment virus
            segment_dirs = glob.glob(f"data/{virus_name}/segment_*")
            if segment_dirs:
                # Count proteins from all segments
                total_proteins = len(analyzer.genome_organization['genes'])
                print(f"  Proteins (all segments): {total_proteins}")
            else:
                print(f"  Proteins: {len(analyzer.proteins)}")
            
            print(f"  Mutations: {len(analyzer.mutations)}")
            print(f"  Genome length: {analyzer.genome_length:,} bp")
            
        except Exception as e:
            print(f"Error processing {virus_name}: {e}")
            continue

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Generate simple combined genome organization and protein mutation analysis visualization',
        formatter_class=argparse.RawDescriptionHelpFormatter)
    
    parser.add_argument(
        '--virus', 
        type=str, 
        help='Virus name (e.g., SARS-CoV-2, HIV-1, Chikungunya, ZaireEbola)'
    )
    
    parser.add_argument(
        '--genome-id', 
        type=str, 
        help='Specific genome ID to analyze (optional, will auto-detect if not provided)'
    )
    
    parser.add_argument(
        '--output', 
        type=str, 
        help='Output file path (optional, will auto-generate if not provided)'
    )
    
    parser.add_argument(
        '--list-viruses', 
        action='store_true',
        help='List available viruses in the data directory'
    )
    
    parser.add_argument(
        '--process-all', 
        action='store_true',
        help='Process all available viruses'
    )
    
    return parser.parse_args()

def main():
    """Main function to run the simple combined analysis."""
    args = parse_arguments()
    
    if args.list_viruses:
        data_dir = "data"
        if os.path.exists(data_dir):
            virus_dirs = [d for d in os.listdir(data_dir) if os.path.isdir(os.path.join(data_dir, d)) and not d.startswith('.')]
            print("Available viruses:")
            for virus in virus_dirs:
                print(f"  - {virus}")
        else:
            print("Data directory not found")
        return
    
    if args.process_all:
        process_multiple_viruses()
        return
    
    if not args.virus:
        print("Error: Please specify a virus name or use --list-viruses to see available options")
        print("Example: python figure1_2_combined_simple.py --virus SARS-CoV-2")
        return
    
    try:
        analyzer = SimpleCombinedAnalyzer(args.virus, args.genome_id)
        analyzer.create_combined_chart(args.output)
        
        print(f"\n✅ Successfully generated simple combined analysis plot:")
        print(f"   Virus: {args.virus}")
        print(f"   Genome ID: {args.genome_id if args.genome_id else 'Auto-detected'}")
        print(f"   Output: {args.output if args.output else 'Auto-generated'}")
        
        # Check if multi-segment virus
        segment_dirs = glob.glob(f"data/{args.virus}/segment_*")
        if segment_dirs:
            # Count proteins from all segments
            total_proteins = len(analyzer.genome_organization['genes'])
            print(f"   Proteins (all segments): {total_proteins}")
        else:
            print(f"   Proteins: {len(analyzer.proteins)}")
        
        print(f"   Mutations: {len(analyzer.mutations)}")
        print(f"   Genome length: {analyzer.genome_length:,} bp")
        
    except FileNotFoundError as e:
        print(f"❌ Error: {e}")
        print("\nTroubleshooting:")
        print("1. Check that the virus name is correct (use --list-viruses to see available viruses)")
        print("2. Ensure that mutation analysis has been run for this virus")
        print("3. Verify that proteome and reference genome files exist")
    except Exception as e:
        print(f"❌ Error: {e}")

if __name__ == "__main__":
    main()
