# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-08-28 16:30:00
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-29 10:01:20
Email: rw3594@nyu.edu
FilePath: /VARIANT/src/visualization/figure2_row_hot_mutations.py
Description: Panel 1 genome organization visualization using row_hot_mutations.csv data
'''

import os
import glob
import csv
import re
import argparse
import pandas as pd
from typing import Dict, List, Optional, Union, Tuple
from Bio import SeqIO
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

# Set default template
pio.templates.default = "simple_white"

class RowHotMutationVisualizer:
    def __init__(self, virus_name: str, sample_id: str = None):
        self.virus_name = virus_name
        self.sample_id = sample_id
        
        # Auto-detect paths
        self.proteome_path = self._auto_detect_proteome_path()
        self.row_hot_csv_path = self._auto_detect_row_hot_csv_path()
        self.reference_genome_path = self._auto_detect_reference_genome_path()
        
        # Load data
        self.proteins = self._load_proteins()
        self.row_hot_mutations = self._parse_row_hot_mutations_from_csv()
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
    
    def _auto_detect_row_hot_csv_path(self) -> str:
        """Auto-detect row_hot_mutations CSV file path."""
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
            specific_path = os.path.join(result_dir, f"{self.sample_id}_row_hot_mutations.csv")
            if os.path.exists(specific_path):
                return specific_path
        
        csv_files = glob.glob(os.path.join(result_dir, "*_row_hot_mutations.csv"))
        if not csv_files:
            raise FileNotFoundError(f"No row_hot_mutations CSV files found in {result_dir}")

        if self.sample_id:
            sample_id_lower = self.sample_id.lower()
            contains_matches = [
                path for path in csv_files
                if sample_id_lower in os.path.basename(path).lower()
            ]
            if contains_matches:
                return sorted(contains_matches)[0]

            normalized_sample = re.sub(r'[^a-z0-9]+', '', sample_id_lower)
            normalized_matches = []
            for path in csv_files:
                stem = os.path.basename(path).replace('_row_hot_mutations.csv', '')
                normalized_stem = re.sub(r'[^a-z0-9]+', '', stem.lower())
                if normalized_sample and normalized_sample in normalized_stem:
                    normalized_matches.append(path)
            if normalized_matches:
                return sorted(normalized_matches)[0]

            raise FileNotFoundError(
                f"No row/hot mutation CSV matched genome ID '{self.sample_id}' in {result_dir}"
            )

        return sorted(csv_files)[0]
    
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
    
    def _parse_row_hot_mutations_from_csv(self) -> List[Dict]:
        """Parse row_hot_mutations data from CSV file."""
        return self._parse_row_hot_mutations_from_csv_path(self.row_hot_csv_path)
    
    def _parse_row_hot_mutations_from_csv_path(self, csv_path: str) -> List[Dict]:
        """Parse row_hot_mutations data from a specific CSV file path."""
        mutations = []
        
        try:
            df = pd.read_csv(csv_path)
            
            for _, row in df.iterrows():
                mutation_type = row['mutation_type']  # 'row' or 'hot'
                position_str = row['position']
                nucleotide_change = row['nucleotide_change']
                protein_affected = row['protein_affected']
                amino_acid_change = row['amino_acid_change']
                biological_classification = row['biological_classification']
                
                # Parse position
                position = self._parse_nt_position(position_str)
                if position is None:
                    continue
                
                mutations.append({
                    'mutation_type': mutation_type,  # 'row' or 'hot'
                    'position': position,
                    'nucleotide_change': nucleotide_change,
                    'protein_affected': protein_affected,
                    'amino_acid_change': amino_acid_change,
                    'biological_classification': biological_classification
                })
                
        except Exception as e:
            print(f"Error reading row_hot_mutations CSV file {csv_path}: {e}")
            return []
        
        return mutations
    
    def _parse_nt_position(self, nt_mutation: str) -> Optional[Union[int, Tuple[int, int]]]:
        """Parse nucleotide position from position string."""
        try:
            if ':' in nt_mutation:
                parts = nt_mutation.split(':')
                if len(parts) >= 2:
                    start_pos = int(parts[0])
                    end_pos = int(parts[1])
                    return (start_pos, end_pos)
            
            # Try to extract single position
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
        """Get color for mutation type (row or hot)."""
        # Use fully opaque RGBA colors to avoid any translucency
        color_map = {
            'row': 'rgba(106, 0, 255, 1.0)',   # Deep violet
            'hot': 'rgba(194, 24, 7, 1.0)',    # Deep crimson
        }
        return color_map.get(mutation_type, '#1f77b4')

    def _add_row_hot_legend(self, fig, mutation_types: set) -> None:
        """Add mutation type legend (same style as figure1) in the top-left."""
        # Keep order consistent: row first, then hot, then any others if ever present
        ordered_types = [t for t in ['row', 'hot'] if t in mutation_types]
        ordered_types += [t for t in sorted(mutation_types) if t not in ('row', 'hot')]
        
        legend_text = "<b>Mutation Types:</b><br>"
        name_map = {
            'row': 'Row',
            'hot': 'Hot'
        }
        for t in ordered_types:
            color = self._get_mutation_color(t)
            display = name_map.get(t, t.title())
            # Use a colored vertical bar '|' plus colored label text (matches figure1 approach)
            legend_text += f"<span style='color:{color};'>|</span> <span style='color:{color};'>{display} Mutations</span><br>"
        
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
    
    def create_row_hot_mutation_chart(self, output_path: str = None) -> go.Figure:
        """Create genome organization chart with row and hot mutations."""
        
        # Auto-detect genome ID from CSV filename
        genome_id = os.path.basename(self.row_hot_csv_path).replace('_row_hot_mutations.csv', '')
        
        # Auto-generate output path if not provided
        if output_path is None:
            os.makedirs('imgs/visualizations', exist_ok=True)
            virus_dir = f'imgs/visualizations/{self.virus_name}'
            os.makedirs(virus_dir, exist_ok=True)
            output_path = f'{virus_dir}/{genome_id}_row_hot_mutations.html'
        
        # Create figure
        fig = go.Figure()
        
        # === PANEL 1: GENOME ORGANIZATION ===
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
                opacity=0.3
            ))
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
        
        # Add row and hot mutations
        mutation_data = []
        
        # Check if this is a multi-segment virus
        segment_dirs = glob.glob(f"result/{self.virus_name}/segment_*")
        
        if segment_dirs:
            # Multi-segment virus - add mutations from all segments
            segment_offset = 0
            for segment_dir in sorted(segment_dirs):
                segment_name = os.path.basename(segment_dir)
                csv_files = glob.glob(os.path.join(segment_dir, "*_row_hot_mutations.csv"))
                
                if csv_files:
                    segment_mutations = self._parse_row_hot_mutations_from_csv_path(csv_files[0])
                    
                    for mutation in segment_mutations:
                        protein_name = f"{segment_name}_{mutation['protein_affected']}"
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
                                    'type': mutation['mutation_type'],
                                    'nucleotide_change': mutation['nucleotide_change'],
                                    'amino_acid_change': mutation['amino_acid_change'],
                                    'biological_classification': mutation['biological_classification'],
                                    'position': mutation['position'],
                                    'gene': gene['name']
                                })
                    
                    # Update offset for next segment
                    segment_proteins = self._load_proteins_from_path(f"data/{self.virus_name}/{segment_name}/refs/{segment_name}_proteome.fasta")
                    segment_length = max([p['position_range'][1] for p in segment_proteins.values() if p['position_range']])
                    segment_offset += segment_length + 1000  # Add gap between segments
        else:
            # Single-segment virus - use existing logic
            for mutation in self.row_hot_mutations:
                protein_name = mutation['protein_affected']
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
                                'x': mut_pos,
                                'y': y_pos,
                                'type': mutation['mutation_type'],
                                'nucleotide_change': mutation['nucleotide_change'],
                                'amino_acid_change': mutation['amino_acid_change'],
                                'biological_classification': mutation['biological_classification'],
                                'position': mutation['position'],
                                'gene': gene['name']
                            })
        
        # Plot mutations as lines (exact same approach as figure1_mutation_analysis.py)
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
                        width=2
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
                hover_text = f"<b>{mut_type.title()} Mutation</b><br>Position: {mut['position']}<br>Gene: {mut['gene']}<br>Nucleotide Change: {mut['nucleotide_change']}<br>Amino Acid Change: {mut['amino_acid_change']}<br>Classification: {mut['biological_classification']}"
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
            ))
        
        # Figure 1-style annotation legend (colored '|' lines with labels)
        self._add_row_hot_legend(fig, mutation_types)
        
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
            title_text="Genes"
        )
        
        # Calculate the maximum x position from all genes
        max_x_pos = max([gene['end'] for gene in all_genes]) if all_genes else self.genome_length
        
        # Update x-axis
        fig.update_xaxes(
            range=[0, max_x_pos],
            showgrid=False,
            linecolor='black',
            linewidth=1,
            tickfont=dict(size=11, family='Arial', color='black'),
            title_text="Genome Position (base pairs)"
        )
        
        # Calculate dynamic height based on number of proteins (same as figure1)
        num_proteins = len(all_genes)
        
        # Panel height based on number of proteins
        baseline_y_range = 36.8
        current_y_range = (len(all_genes) - 1) * (0.6 + 0.6) + 0.6 + 0.4
        panel_height = int(600 * (current_y_range / baseline_y_range))
        
        # Total height with margins
        total_height = panel_height + 200  # Extra space for margins and title
        
        # Update layout
        fig.update_layout(
            title=f"Row and Hot Mutations - {self.virus_name} ({genome_id})",
            width=1200,
            height=total_height,
            margin=dict(l=50, r=50, t=100, b=50),
            showlegend=False
        )
        
        # Save
        fig.write_html(output_path)
        fig.write_image(output_path.replace('.html', '.pdf'), width=1200, height=total_height)
        
        print(f"Row/Hot mutation chart saved to: {output_path}")
        
        # Print summary
        print(f"\nRow/Hot Mutation Analysis Summary:")
        print(f"  Virus: {self.virus_name}")
        print(f"  Genome ID: {genome_id}")
        print(f"  Proteins: {len(all_genes)}")
        
        # Count mutations by type
        row_count = len([m for m in mutation_data if m['type'] == 'row'])
        hot_count = len([m for m in mutation_data if m['type'] == 'hot'])
        
        print(f"  Row Mutations: {row_count}")
        print(f"  Hot Mutations: {hot_count}")
        print(f"  Total Mutations: {len(mutation_data)}")
        print(f"  Genome length: {self.genome_length:,} bp")
        
        return fig

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Generate genome organization visualization with row and hot mutations',
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
    
    return parser.parse_args()

def main():
    """Main function to run the row/hot mutation analysis."""
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
    
    if not args.virus:
        print("Error: Please specify a virus name or use --list-viruses to see available options")
        print("Example: python figure2_row_hot_mutations.py --virus SARS-CoV-2")
        return
    
    try:
        visualizer = RowHotMutationVisualizer(args.virus, args.genome_id)
        visualizer.create_row_hot_mutation_chart(args.output)
        
        print(f"\n✅ Successfully generated row/hot mutation analysis plot:")
        print(f"   Virus: {args.virus}")
        print(f"   Genome ID: {args.genome_id if args.genome_id else 'Auto-detected'}")
        print(f"   Output: {args.output if args.output else 'Auto-generated'}")
        
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
