# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-08-19 16:07:24
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-19 17:37:22
Email: wang.rui@nyu.edu
FilePath: /VARIANT/plot/genome_organization_visualization.py
Description: 
'''
# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-08-19 17:30:00
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-19 17:30:00
Email: wang.rui@nyu.edu
FilePath: /VARIANT/plot/genome_organization_visualization.py
Description: Genome organization visualization for viral genomes with mutation analysis.
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

# Set default template
pio.templates.default = "simple_white"

class GenomeOrganizationVisualizer:
    """Create genome organization visualizations showing genes across reading frames."""
    
    def __init__(self, proteome_path: str, mutation_csv_path: str, reference_genome_path: str = None):
        """
        Initialize the genome organization visualizer.

        Args:
            proteome_path (str): Path to the proteome FASTA file
            mutation_csv_path (str): Path to mutation summary CSV file
            reference_genome_path (str, optional): Path to reference genome FASTA
        """
        self.proteome_path = proteome_path
        self.mutation_csv_path = mutation_csv_path
        self.reference_genome_path = reference_genome_path
        self.proteins = self._load_proteins()
        self.mutations = self._parse_mutations_from_csv()
        self.genome_length = self._get_genome_length()
        self.genome_organization = self._create_genome_organization()
        
    def _load_proteins(self) -> Dict[str, Dict]:
        """Load protein sequences from FASTA file."""
        proteins = {}
        
        with open(self.proteome_path, 'r') as f:
            for record in SeqIO.parse(f, 'fasta'):
                # Parse header: YP_009725297.1|leader_nsp1|266..805
                header_parts = record.description.split('|')
                if len(header_parts) >= 2:
                    protein_name = header_parts[1]
                    sequence = str(record.seq)
                    
                    # Extract position range if available - handle complex cases
                    position_range = None
                    if len(header_parts) >= 3:
                        pos_str = header_parts[2]
                        # Handle simple ranges like "266..805"
                        if '..' in pos_str and not pos_str.startswith('join'):
                            try:
                                parts = pos_str.split('..')
                                if len(parts) == 2:
                                    position_range = (int(parts[0]), int(parts[1]))
                            except (ValueError, IndexError):
                                pass
                        # Handle complex cases like "join(13442..13468,13468..16236)"
                        elif pos_str.startswith('join'):
                            try:
                                # Extract the first range from join statement
                                # join(13442..13468,13468..16236) -> get 13442..13468
                                start_bracket = pos_str.find('(')
                                end_bracket = pos_str.find(')')
                                if start_bracket != -1 and end_bracket != -1:
                                    join_content = pos_str[start_bracket+1:end_bracket]
                                    # Get the first range (before the comma)
                                    first_range = join_content.split(',')[0]
                                    if '..' in first_range:
                                        parts = first_range.split('..')
                                        if len(parts) == 2:
                                            start_pos = int(parts[0])
                                            # For the end, we'll use the last range's end position
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
        mutations = []
        
        try:
            # Read CSV file
            df = pd.read_csv(self.mutation_csv_path)
            
            for _, row in df.iterrows():
                protein_name = row['Protein Name']
                protein_mutation = row['Protein Mutations']
                aa_mutation_type = row['AA Mutation Type']
                nt_mutation = row['NT Mutation']
                
                # Skip invalid proteins
                if protein_name in ["Invalid protein sequence", "None-CDS"]:
                    continue
                
                # Parse nucleotide position from NT Mutation column
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
            print(f"Error reading CSV file: {e}")
            return []
        
        return mutations
    
    def _parse_nt_position(self, nt_mutation: str) -> Optional[Union[int, Tuple[int, int]]]:
        """Parse nucleotide position from NT Mutation column."""
        try:
            # Handle range mutations like "10447:10449GC->AA" or "21633:21641TACCCCCTG"
            if ':' in nt_mutation:
                # Extract position range before the nucleotide change
                parts = nt_mutation.split(':')
                if len(parts) >= 2:
                    start_pos = int(parts[0])
                    # Extract end position - look for the number before any letters
                    end_part = parts[1]
                    # Find the first sequence of digits
                    end_match = re.search(r'^(\d+)', end_part)
                    if end_match:
                        end_pos = int(end_match.group(1))
                        return (start_pos, end_pos)
                    else:
                        # Fallback to single position
                        return start_pos
            
            # Handle single position mutations like "T670G" or "C378A"
            # Extract position from the mutation string
            # Look for patterns like "C378A" where 378 is the position
            match = re.search(r'[ACGT](\d+)[ACGT]', nt_mutation)
            if match:
                return int(match.group(1))
            
            # If no pattern found, try to extract any number
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
        
        # Estimate from protein position ranges
        max_position = 0
        for protein_data in self.proteins.values():
            if protein_data['position_range']:
                max_position = max(max_position, protein_data['position_range'][1])
        
        # If no position ranges, use a typical SARS-CoV-2 length
        if max_position == 0:
            max_position = 30000  # Typical SARS-CoV-2 genome length
        
        return max_position
    
    def _create_genome_organization(self) -> Dict:
        """Create generic genome organization data structure for any virus."""
        
        # Generic genome organization - works for any virus
        genome_org = {
            'genes': [],
            'regulatory_regions': []
        }
        
        # Sort proteins by position
        sorted_proteins = []
        for protein_name, protein_data in self.proteins.items():
            if protein_data['position_range']:
                sorted_proteins.append((protein_name, protein_data))
        
        sorted_proteins.sort(key=lambda x: x[1]['position_range'][0])
        
        # Add all proteins as individual genes
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
        
        # Add basic regulatory regions if genome length is available
        if self.genome_length > 0:
            # Estimate UTR regions (first 5% and last 5% of genome)
            utr_size = max(100, int(self.genome_length * 0.05))
            genome_org['regulatory_regions'] = [
                {'name': '5\' UTR', 'start': 1, 'end': utr_size, 'type': 'UTR'},
                {'name': '3\' UTR', 'start': self.genome_length - utr_size, 'end': self.genome_length, 'type': 'UTR'}
            ]
        
        return genome_org
    
    def _get_protein_color(self, protein_name: str) -> str:
        """Get color for protein using Plotly continuous color palette."""
        import plotly.colors as pc
        
        # Use Plotly's Viridis continuous color palette
        # This creates a smooth, professional color gradient
        colors = pc.cyclical.Twilight_r
        # colors = pc.qualitative.Set3
        # colors = pc.diverging.Spectral

        # Use deterministic color assignment based on protein name
        total = sum(ord(c) for c in protein_name)
        color_index = total % len(colors)
        return colors[color_index]
    
    def _get_mutation_color(self, mutation_type: str) -> str:
        """Get color for mutation type."""
        color_map = {
            'missense': '#d62728',      # Red for missense
            'silent': '#2ca02c',        # Green for silent
            'deletion': '#ff7f0e',      # Orange for deletion
            'insertion': '#9467bd',     # Purple for insertion
            'nonsense': '#e377c2',      # Pink for nonsense
            'frameshift': '#8c564b',    # Brown for frameshift
        }
        return color_map.get(mutation_type, '#1f77b4')  # Default blue for any other types
    
    def _format_mutation_breakdown(self, mutation_details: dict) -> str:
        """Format mutation details breakdown for hover display."""
        if not mutation_details:
            return "No mutations"
        
        breakdown_lines = []
        for mut_type, mutations in sorted(mutation_details.items()):
            # Capitalize first letter and show mutations
            formatted_type = mut_type.capitalize()
            count = len(mutations)
            
            if count <= 5:
                # Show all mutations if 5 or fewer
                mutations_str = ", ".join(mutations)
                breakdown_lines.append(f"{formatted_type}: {mutations_str}")
            else:
                # Show all mutations, 5 per line
                breakdown_lines.append(f"{formatted_type} ({count}):")
                for i in range(0, len(mutations), 5):
                    chunk = mutations[i:i+5]
                    mutations_str = ", ".join(chunk)
                    breakdown_lines.append(f"  {mutations_str}")
        
        return "<br>".join(breakdown_lines)
    
    def _add_mutation_legend(self, fig):
        """Add mutation type legend in top left corner."""
        # Get unique mutation types from the data
        mutation_types_in_data = set()
        for mutation in self.mutations:
            mutation_types_in_data.add(mutation['type'])
        
        # Define mutation types and their colors - prioritize the main four types
        mutation_type_info = {
            'missense': ('Missense', '#d62728', 'Amino acid change'),
            'deletion': ('Deletion', '#ff7f0e', 'Nucleotide deletion'),
            'insertion': ('Insertion', '#9467bd', 'Nucleotide insertion'),
            'nonsense': ('Nonsense', '#e377c2', 'Stop codon introduction'),
            'silent': ('Silent', '#2ca02c', 'No amino acid change'),
            'frameshift': ('Frameshift', '#8c564b', 'Reading frame shift')
        }
        
        # Create legend text - show main four types first, then others present in data
        legend_text = "<b>Mutation Types:</b><br>"
        
        # Define the main four mutation types in order
        main_types = ['missense', 'deletion', 'insertion', 'nonsense']
        
        # Show main types first (all as solid dots since they're valid mutation types)
        for mut_type in main_types:
            if mut_type in mutation_type_info:
                display_name, color, description = mutation_type_info[mut_type]
                legend_text += f"<span style='color:{color};'>●</span> {display_name}: {description}<br>"
        
        # Show other types present in data
        for mut_type in sorted(mutation_types_in_data):
            if mut_type not in main_types and mut_type in mutation_type_info:
                display_name, color, description = mutation_type_info[mut_type]
                legend_text += f"<span style='color:{color};'>●</span> {display_name}: {description}<br>"
        
        # Add annotation for legend
        fig.add_annotation(
            x=0.02,  # Left position
            y=0.98,  # Top position
            xref='paper',
            yref='paper',
            text=legend_text,
            showarrow=False,
            bgcolor='rgba(255, 255, 255, 0.9)',
            bordercolor='black',
            borderwidth=1,
            font=dict(size=10),
            align='left'
        )
    
    def _calculate_protein_conservation(self) -> Dict[str, Dict]:
        """Calculate conservation metrics for each protein."""
        conservation_data = {}
        
        # Count mutations per protein
        protein_mutations = {}
        for mutation in self.mutations:
            protein_name = mutation['protein']
            if protein_name not in protein_mutations:
                protein_mutations[protein_name] = []
            protein_mutations[protein_name].append(mutation)
        
        # Calculate metrics for each protein
        for protein_name, protein_data in self.proteins.items():
            if protein_name not in protein_mutations:
                protein_mutations[protein_name] = []
            
            mutations = protein_mutations[protein_name]
            protein_length = protein_data['length']
            genomic_length = protein_data['position_range'][1] - protein_data['position_range'][0] if protein_data['position_range'] else 0
            
            # Calculate conservation metrics
            mutation_count = len(mutations)
            mutation_rate = mutation_count / protein_length if protein_length > 0 else 0
            mutation_density = mutation_count / genomic_length if genomic_length > 0 else 0
            conservation_ratio = (protein_length - mutation_count) / protein_length if protein_length > 0 else 1.0
            
            # Group mutations by type with their details
            mutation_details = {}
            for mutation in mutations:
                mut_type = mutation['type']
                if mut_type not in mutation_details:
                    mutation_details[mut_type] = []
                mutation_details[mut_type].append(mutation['mutation'])

            conservation_data[protein_name] = {
                'mutation_count': mutation_count,
                'mutation_rate': mutation_rate,
                'mutation_density': mutation_density,
                'conservation_ratio': conservation_ratio,
                'protein_length': protein_length,
                'genomic_length': genomic_length,
                'mutation_details': mutation_details
            }
        
        return conservation_data
    
    def create_genome_organization_chart(self, output_path: str = "genome_organization.html"):
        """Create genome organization visualization with protein conservation analysis."""
        
        # Create subplots: proteins and conservation
        fig = make_subplots(
            rows=2, cols=1,
            subplot_titles=('All Proteins', 'Protein Conservation'),
            vertical_spacing=0.08,
            row_heights=[0.7, 0.3],  # Adjust heights for two panels
            specs=[[{"secondary_y": False}], [{"secondary_y": False}]]
        )
        
        # Add mutation type legend in top left corner
        self._add_mutation_legend(fig)
        
        # Calculate protein conservation metrics
        conservation_data = self._calculate_protein_conservation()
        
        # Track all genes for mutation overlay
        all_genes = []
        
        # Add all proteins to the first panel with unique Y positions
        print(f"Adding proteins to visualization:")
        y_offset = 0
        
        for gene in self.genome_organization['genes']:
            print(f"  Adding: {gene['name']} ({gene['start']}-{gene['end']})")
            
            # Use simple numeric Y position for better axis handling
            y_pos = y_offset
            y_offset += 1
            
            # Get conservation data for this protein
            conservation_info = conservation_data.get(gene['name'], {})
            mutation_count = conservation_info.get('mutation_count', 0)
            mutation_rate = conservation_info.get('mutation_rate', 0)
            
            fig.add_trace(go.Bar(
                x=[gene['end'] - gene['start']],
                y=[y_pos],
                orientation='h',
                name=gene['name'],
                marker_color=gene['color'],  # Use original colorful protein colors
                base=gene['start'],
                hovertemplate=f"<b>{gene['name']}</b><br>"
                            f"Position: {gene['start']}-{gene['end']}<br>"
                            f"Length: {gene['length']} bp<br>"
                            f"Protein: {gene['protein_length']} aa<br>"
                            f"Mutations: {mutation_count}<br>"
                            f"Mutation Rate: {mutation_rate:.3f}<extra></extra>",
                showlegend=False
            ), row=1, col=1)
            all_genes.append(gene)
        print(f"Total genes added to visualization: {len(all_genes)}")
        
        # 2. Protein conservation panel
        print("Adding protein conservation panel...")
        conservation_proteins = []
        conservation_rates = []
        conservation_colors = []
        
        for gene in self.genome_organization['genes']:
            conservation_info = conservation_data.get(gene['name'], {})
            conservation_proteins.append(gene['name'])
            conservation_rates.append(conservation_info.get('mutation_rate', 0))
            # Use Nature Reviews color palette instead of conservation colors
            conservation_colors.append(gene['color'])
        
        # Add conservation bar chart with enhanced hover data
        fig.add_trace(go.Bar(
            x=conservation_proteins,
            y=conservation_rates,
            name='Mutation Rate',
            marker_color=conservation_colors,
            hovertemplate='<b>%{x}</b><br>' +
                         'Mutation Rate: %{y:.3f}<br>' +
                         'Protein Length: %{customdata[0]} aa<br>' +
                         'Total Mutations: %{customdata[1]}<br>' +
                         '%{customdata[2]}<extra></extra>',
            customdata=[[
                conservation_data.get(protein, {}).get('protein_length', 0),
                conservation_data.get(protein, {}).get('mutation_count', 0),
                self._format_mutation_breakdown(conservation_data.get(protein, {}).get('mutation_details', {}))
            ] for protein in conservation_proteins]
        ), row=2, col=1)
        
        # Print conservation analysis summary
        print(f"\nConservation Analysis Summary:")
        for protein_name, rate in zip(conservation_proteins, conservation_rates):
            print(f"  {protein_name}: {rate:.3f}")
        
        # Add regulatory regions as annotations
        for region in self.genome_organization['regulatory_regions']:
            mid_pos = (region['start'] + region['end']) / 2
            fig.add_annotation(
                x=mid_pos,
                y=1.02,
                xref='x',
                yref='paper',
                text=region['name'],
                showarrow=False,
                font=dict(size=10, color='gray'),
                textangle=-45
            )
        
        # Add mutation markers on genes
        for mutation in self.mutations:
            protein_name = mutation['protein']
            # Find the gene in our list
            gene = None
            for g in self.genome_organization['genes']:
                if g['name'] == protein_name:
                    gene = g
                    break
            
            if gene is not None:
                
                # Calculate mutation position within gene
                if isinstance(mutation['position'], tuple):
                    mut_pos = (mutation['position'][0] + mutation['position'][1]) / 2
                else:
                    mut_pos = mutation['position']
                
                # Map to gene position
                if gene['start'] <= mut_pos <= gene['end']:
                    relative_pos = mut_pos - gene['start']
                    gene_pos = gene['start'] + relative_pos
                    
                    # Find the corresponding Y position for this gene
                    gene_index = None
                    for i, g in enumerate(all_genes):
                        if g['name'] == gene['name']:
                            gene_index = i
                            break
                    
                    if gene_index is not None:
                        y_pos = gene_index
                        
                        fig.add_trace(go.Scatter(
                            x=[gene_pos],
                            y=[y_pos],
                            mode='markers',
                            marker=dict(
                                size=8,
                                color=self._get_mutation_color(mutation['type']),
                                symbol='circle',
                                line=dict(color='white', width=1)
                            ),
                            name=f"{mutation['type']} - {mutation['mutation']}",
                            hovertemplate=f"<b>Mutation: {mutation['type']}</b><br>"
                                        f"Position: {mutation['position']}<br>"
                                        f"Change: {mutation['mutation']}<br>"
                                        f"Gene: {gene['name']}<extra></extra>",
                            showlegend=False
                        ), row=1, col=1)
        
        # Update layout - Nature Reviews style
        fig.update_layout(
            title=dict(
                text="Viral Genome Organization with Mutation Analysis",
                font=dict(size=14, color='black'),
                x=0.5,
                xanchor='center'
            ),
            xaxis_title="Genome Position (base pairs)",
            width=1200,
            height=1200,
            margin=dict(l=60, r=60, t=80, b=60),  # More margin for clarity
            showlegend=False,
            plot_bgcolor='white',  # Clean white background
            paper_bgcolor='white'
        )
        
        # Update axes for both panels - Nature Reviews style
        for i in range(1, 3):
            if i == 1:
                # Protein panel - genome position on X-axis
                fig.update_xaxes(
                    title_text="Genome Position (base pairs)", 
                    row=i, col=1, 
                    range=[0, self.genome_length],
                    showgrid=True,
                    gridcolor='#f0f0f0',
                    linecolor='black',
                    tickfont=dict(size=10, color='black')
                )
                fig.update_yaxes(
                    title_text="Proteins", 
                    row=i, col=1,
                    tickmode='array',
                    ticktext=[g['name'] for g in all_genes],
                    tickvals=list(range(len(all_genes))),
                    tickangle=0,
                    tickfont=dict(size=9, color='black'),
                    linecolor='black'
                )
            elif i == 2:
                # Conservation panel - protein names on X-axis
                fig.update_xaxes(
                    title_text="Proteins", 
                    row=i, col=1, 
                    tickangle=45,
                    tickfont=dict(size=8, color='black'),
                    linecolor='black'
                )
                fig.update_yaxes(
                    title_text="Mutation Rate", 
                    row=i, col=1,
                    linecolor='black',
                    tickfont=dict(size=10, color='black')
                )
            else:
                # Mutation density panel - genome position on X-axis
                fig.update_xaxes(
                    title_text="Genome Position (base pairs)", 
                    row=i, col=1, 
                    range=[0, self.genome_length],
                    linecolor='black',
                    tickfont=dict(size=10, color='black')
                )
                fig.update_yaxes(
                    title_text="Mutation Count", 
                    row=i, col=1,
                    linecolor='black',
                    tickfont=dict(size=10, color='black')
                )

        # Save the figure
        fig.write_html(output_path)
        fig.write_image(output_path.replace('.html', '.pdf'), width=1200, height=1200)
        
        print(f"Genome organization chart saved to: {output_path}")
        return fig

def process_multiple_viruses():
    """
    Process multiple virus datasets and create genome organization visualizations.
    
    Automatically processes CSV files for SARS-CoV-2, HIV-1, Chikungunya, and ZaireEbola.
    Creates HTML and PDF visualizations for each virus dataset.
    """
    
    # Define virus configurations
    virus_configs = {
        'SARS-CoV-2': {
            'proteome_path': 'data/SARS-CoV-2/refs/SARS-CoV-2_proteome.fasta',
            'reference_genome_path': 'data/SARS-CoV-2/refs/NC_045512.fasta',
            'result_dir': 'result/SARS-CoV-2/',
            'output_dir': 'plot/'
        },
        'HIV-1': {
            'proteome_path': 'data/HIV-1/refs/HIV-1_proteome.fasta',
            'reference_genome_path': 'data/HIV-1/refs/NC_001802.1.fasta',
            'result_dir': 'result/HIV-1/',
            'output_dir': 'plot/'
        },
        'Chikungunya': {
            'proteome_path': 'data/Chikungunya/refs/Chikungunya_proteome.fasta',
            'reference_genome_path': 'data/Chikungunya/refs/NC_004162.fasta',
            'result_dir': 'result/Chikungunya/',
            'output_dir': 'plot/'
        },
        'ZaireEbola': {
            'proteome_path': 'data/ZaireEbola/refs/ZaireEbola_proteome.fasta',
            'reference_genome_path': 'data/ZaireEbola/refs/NC_002549.1.fasta',
            'result_dir': 'result/ZaireEbola/',
            'output_dir': 'plot/'
        }
    }
    
    for virus_name, config in virus_configs.items():
        print(f"\n{'='*50}")
        print(f"Processing {virus_name}")
        print(f"{'='*50}")
        
        # Check if proteome file exists
        if not os.path.exists(config['proteome_path']):
            print(f"Warning: Proteome file not found: {config['proteome_path']}")
            continue
        
        # Find CSV files in result directory
        csv_files = []
        if os.path.exists(config['result_dir']):
            for file in os.listdir(config['result_dir']):
                if file.endswith('_mutation_summary.csv'):
                    csv_files.append(os.path.join(config['result_dir'], file))
        
        if not csv_files:
            print(f"No mutation summary CSV files found in {config['result_dir']}")
            continue
        
        print(f"Found {len(csv_files)} CSV files")
        
        # Process each CSV file
        for csv_file in csv_files[:3]:  # Limit to first 3 files per virus
            try:
                sample_name = os.path.basename(csv_file).replace('_mutation_summary.csv', '')
                output_file = f"{config['output_dir']}{virus_name}_{sample_name}_genome_organization.html"
                
                print(f"Processing {sample_name}...")
                
                # Create visualizer
                visualizer = GenomeOrganizationVisualizer(
                    config['proteome_path'], 
                    csv_file, 
                    config['reference_genome_path']
                )
                
                # Create visualization
                visualizer.create_genome_organization_chart(output_file)
                
                print(f"Created: {output_file}")
                print(f"  Proteins: {len(visualizer.proteins)}")
                print(f"  Mutations: {len(visualizer.mutations)}")
                print(f"  Genome length: {visualizer.genome_length:,} bp")
                
            except Exception as e:
                print(f"Error processing {csv_file}: {e}")
                continue

def main():
    """
    Main function to test the genome organization visualization.
    
    Creates genome organization visualizations for SARS-CoV-2 data.
    """
    
    # File paths - using CSV files instead of .txt files
    proteome_path = "data/SARS-CoV-2/refs/SARS-CoV-2_proteome.fasta"
    mutation_csv_path = "result/SARS-CoV-2/EPI_ISL_16127650_mutation_summary.csv"
    reference_genome_path = "data/SARS-CoV-2/refs/NC_045512.fasta"
    
    # Check if files exist
    if not os.path.exists(proteome_path):
        print(f"Error: Proteome file not found: {proteome_path}")
        return
    
    if not os.path.exists(mutation_csv_path):
        print(f"Error: Mutation CSV file not found: {mutation_csv_path}")
        return
    
    # Create visualizer
    visualizer = GenomeOrganizationVisualizer(proteome_path, mutation_csv_path, reference_genome_path)
    
    # Create genome organization visualizations
    print("Creating genome organization chart...")
    visualizer.create_genome_organization_chart("plot/genome_organization_analysis.html")
    
    # Print summary statistics
    print(f"\nGenome Organization Analysis Summary:")
    print(f"Proteins loaded: {len(visualizer.proteins)}")
    print(f"Mutations parsed: {len(visualizer.mutations)}")
    print(f"Genome length: {visualizer.genome_length:,} base pairs")
    
    # Show gene distribution
    print(f"Total genes: {len(visualizer.genome_organization['genes'])}")
    for gene in visualizer.genome_organization['genes']:
        print(f"  {gene['name']}: {gene['start']}-{gene['end']} ({gene['length']} bp)")

if __name__ == "__main__":
    # Process multiple viruses
    process_multiple_viruses()
    
    # Also run the original single virus example
    print(f"\n{'='*50}")
    print("Running single virus example (SARS-CoV-2)")
    print(f"{'='*50}")
    main()
