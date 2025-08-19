#!/usr/bin/env python3
"""
Figure 3: Genome-Wide Mutation Landscape
=======================================

This script creates a comprehensive genome-wide mutation landscape showing:
- Mutation frequency across the genome
- Protein region annotations
- Known variant positions
- Mutation type distribution
Generic version that works with any virus data.
"""

import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
import pandas as pd
import numpy as np
import os
import glob
from collections import defaultdict
import re

# Set default template
pio.templates.default = "simple_white"

def get_genome_regions(virus_name):
    """
    Get genome regions for a specific virus. This is a placeholder function.
    In a real implementation, you would load this from a configuration file
    or database based on the virus name.
    
    Args:
        virus_name (str): Name of the virus
        
    Returns:
        dict: Genome regions with coordinates
    """
    # Default regions - in practice, this would be loaded from virus-specific config
    default_regions = {
        'gene_1': (1, 1000),
        'gene_2': (1001, 2000),
        'gene_3': (2001, 3000),
        'gene_4': (3001, 4000),
        'gene_5': (4001, 5000)
    }
    
    # Virus-specific regions (add more as needed)
    virus_regions = {
        'SARS-CoV-2': {
            'leader_nsp1': (266, 805),
            'nsp2': (806, 2719),
            'nsp3': (2720, 8554),
            'nsp4': (8555, 10054),
            '3C-likease': (10055, 10972),
            'nsp6': (10973, 11842),
            'nsp7': (11843, 12091),
            'nsp8': (12092, 12685),
            'nsp9': (12686, 13024),
            'nsp10': (13025, 13441),
            'RNA-dependent-polymerase': (13442, 16236),
            'helicase': (16237, 18039),
            "3'-to-5'_exonuclease": (18040, 19620),
            'endoRNAse': (19621, 20658),
            'nsp16': (20659, 21552),
            'spike_surface_glycoprotein': (21563, 25384),
            'ORF3a': (25393, 26220),
            'envelope': (26245, 26472),
            'membrane_glycoprotein': (26523, 27191),
            'ORF6': (27202, 27387),
            'ORF7a': (27394, 27759),
            'ORF7b': (27756, 27887),
            'ORF8': (27894, 28259),
            'nucleocapsid_phosphoprotein': (28274, 29533),
            'ORF10': (29558, 29674)
        }
    }
    
    return virus_regions.get(virus_name, default_regions)

def get_known_variants(virus_name):
    """
    Get known variants for a specific virus. This is a placeholder function.
    
    Args:
        virus_name (str): Name of the virus
        
    Returns:
        dict: Known variants with positions and descriptions
    """
    # Default variants - in practice, this would be loaded from virus-specific config
    default_variants = {
        100: ('Variant1', 'Gene1', 'Common variant'),
        500: ('Variant2', 'Gene2', 'Important variant'),
        1000: ('Variant3', 'Gene3', 'Major variant')
    }
    
    # Virus-specific variants (add more as needed)
    virus_variants = {
        'SARS-CoV-2': {
            23403: ('D614G', 'Spike', 'Major variant'),
            23063: ('N501Y', 'Spike', 'Alpha/Beta/Delta'),
            22813: ('K417N', 'Spike', 'Beta/Delta'),
            22917: ('L452R', 'Spike', 'Delta'),
            22995: ('T478K', 'Spike', 'Delta'),
            23013: ('E484A', 'Spike', 'Beta/Delta'),
            23055: ('Q498R', 'Spike', 'Omicron'),
            23075: ('Y505H', 'Spike', 'Omicron'),
            14408: ('P323L', 'RdRP', 'Major variant'),
            3037: ('C241T', 'NSP3', 'Silent'),
            27964: ('Q57H', 'ORF3a', 'Common variant')
        }
    }
    
    return virus_variants.get(virus_name, default_variants)

def parse_mutation_positions(file_path):
    """
    Parse mutation positions from a mutation file.
    
    Args:
        file_path (str): Path to the mutation file
        
    Returns:
        list: List of mutation positions and types
    """
    mutations = []
    
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    parts = line.split()
                    if len(parts) >= 2:
                        mutation_type = parts[0]
                        
                        # Extract position
                        pos_str = parts[1]
                        if ':' in pos_str:
                            # Range mutation (e.g., 10447:10449)
                            start_pos = int(pos_str.split(':')[0])
                            end_pos = int(pos_str.split(':')[1])
                            for pos in range(start_pos, end_pos + 1):
                                mutations.append((pos, mutation_type))
                        else:
                            # Single position
                            try:
                                pos = int(pos_str)
                                mutations.append((pos, mutation_type))
                            except ValueError:
                                continue
    except FileNotFoundError:
        print(f"Warning: File not found {file_path}")
    
    return mutations

def get_mutation_data(virus_name):
    """
    Collect mutation data from virus result files.
    
    Args:
        virus_name (str): Name of the virus to analyze
        
    Returns:
        tuple: (mutation_data, total_files)
    """
    result_dir = f"../result/{virus_name}"
    
    if not os.path.exists(result_dir):
        print(f"❌ Result directory not found: {result_dir}")
        return {}, 0
    
    mutation_files = glob.glob(os.path.join(result_dir, '*_20250807.txt'))
    
    # Filter out summary and PRF files
    mutation_files = [f for f in mutation_files if 'mutation_summary' not in f and 'potential_PRF' not in f]
    
    # Collect all mutations
    all_mutations = []
    file_count = 0
    
    for file_path in mutation_files:
        mutations = parse_mutation_positions(file_path)
        all_mutations.extend(mutations)
        file_count += 1
    
    # Count mutation frequency at each position
    position_counts = defaultdict(lambda: defaultdict(int))
    total_files = file_count
    
    for pos, mutation_type in all_mutations:
        position_counts[pos][mutation_type] += 1
    
    return dict(position_counts), total_files

def create_genome_landscape_figure(mutation_data, total_files, virus_name, genome_length=30000):
    """
    Create the genome-wide mutation landscape figure.
    
    Args:
        mutation_data (dict): Mutation frequency data
        total_files (int): Total number of files analyzed
        virus_name (str): Name of the virus
        genome_length (int): Length of the genome
        
    Returns:
        plotly.graph_objects.Figure: The created figure
    """
    
    # Get virus-specific data
    genome_regions = get_genome_regions(virus_name)
    known_variants = get_known_variants(virus_name)
    
    # Create subplots
    fig = make_subplots(
        rows=3, cols=1,
        subplot_titles=(
            f"Mutation Frequency Across {virus_name} Genome",
            "Mutation Types by Position",
            "Protein Regions and Known Variants"
        ),
        specs=[[{"secondary_y": False}],
               [{"secondary_y": False}],
               [{"secondary_y": False}]],
        vertical_spacing=0.1,
        row_heights=[0.4, 0.3, 0.3]
    )
    
    # Genome positions
    genome_positions = list(range(1, genome_length + 1))
    
    # Plot 1: Overall mutation frequency
    mutation_frequency = []
    for pos in genome_positions:
        if pos in mutation_data:
            total_mutations = sum(mutation_data[pos].values())
            frequency = (total_mutations / total_files) * 100
            mutation_frequency.append(frequency)
        else:
            mutation_frequency.append(0)
    
    fig.add_trace(
        go.Scatter(
            x=genome_positions,
            y=mutation_frequency,
            mode='lines',
            name='Mutation Frequency (%)',
            line=dict(color='#1f77b4', width=2),
            fill='tonexty',
            fillcolor='rgba(31, 119, 180, 0.3)'
        ),
        row=1, col=1
    )
    
    # Plot 2: Mutation types
    mutation_types = ['missense', 'silent', 'deletion', 'frameshift', 'insertion']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    
    for i, mutation_type in enumerate(mutation_types):
        type_frequency = []
        for pos in genome_positions:
            if pos in mutation_data:
                count = mutation_data[pos].get(mutation_type, 0)
                frequency = (count / total_files) * 100
                type_frequency.append(frequency)
            else:
                type_frequency.append(0)
        
        fig.add_trace(
            go.Scatter(
                x=genome_positions,
                y=type_frequency,
                mode='lines',
                name=f'{mutation_type.capitalize()} (%)',
                line=dict(color=colors[i], width=1.5),
                opacity=0.8
            ),
            row=2, col=1
        )
    
    # Plot 3: Protein regions
    region_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
                     '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    
    for i, (region_name, (start, end)) in enumerate(genome_regions.items()):
        fig.add_shape(
            type="rect",
            x0=start,
            y0=0,
            x1=end,
            y1=1,
            fillcolor=region_colors[i % len(region_colors)],
            opacity=0.3,
            line=dict(color="black", width=1),
            row=3, col=1
        )
        
        # Add region label
        fig.add_annotation(
            x=(start + end) / 2,
            y=0.5,
            text=region_name.replace('_', '<br>'),
            showarrow=False,
            font=dict(size=8, color="black"),
            xanchor="center",
            yanchor="middle",
            row=3, col=1
        )
    
    # Add known variants
    for pos, (variant_name, protein, description) in known_variants.items():
        fig.add_trace(
            go.Scatter(
                x=[pos],
                y=[0.8],
                mode='markers+text',
                name=variant_name,
                text=[variant_name],
                textposition='top center',
                marker=dict(
                    symbol='diamond',
                    size=12,
                    color='red',
                    line=dict(color='black', width=2)
                ),
                showlegend=False
            ),
            row=3, col=1
        )
        
        fig.add_annotation(
            x=pos,
            y=0.6,
            text=f"{variant_name}<br>{protein}<br>{description}",
            showarrow=True,
            arrowhead=2,
            arrowsize=1,
            arrowwidth=2,
            arrowcolor="red",
            font=dict(size=8, color="red"),
            xanchor="center",
            yanchor="top",
            row=3, col=1
        )
    
    # Update layout
    fig.update_layout(
        title={
            'text': f"{virus_name} Genome-Wide Mutation Landscape",
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 20, 'family': 'Arial Black'}
        },
        width=1000,
        height=1000,
        margin=dict(l=50, r=50, t=100, b=50),
        showlegend=True,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        )
    )
    
    # Update axes
    fig.update_xaxes(title_text="Genome Position", row=3, col=1)
    fig.update_yaxes(title_text="Mutation Frequency (%)", row=1, col=1)
    fig.update_yaxes(title_text="Mutation Type Frequency (%)", row=2, col=1)
    fig.update_yaxes(title_text="Protein Regions", row=3, col=1, range=[0, 1])
    
    # Add summary statistics
    total_mutations = sum(sum(counts.values()) for counts in mutation_data.values())
    unique_positions = len(mutation_data)
    
    fig.add_annotation(
        x=0.02,
        y=0.98,
        xref="paper",
        yref="paper",
        text=f"Total Files: {total_files}<br>Total Mutations: {total_mutations}<br>Unique Positions: {unique_positions}",
        showarrow=False,
        font=dict(size=12, color="black"),
        xanchor="left",
        yanchor="top",
        bgcolor="white",
        bordercolor="black",
        borderwidth=1
    )
    
    return fig

def create_mutation_hotspot_figure(mutation_data, total_files, virus_name):
    """
    Create a focused view of mutation hotspots.
    
    Args:
        mutation_data (dict): Mutation frequency data
        total_files (int): Total number of files analyzed
        virus_name (str): Name of the virus
        
    Returns:
        plotly.graph_objects.Figure: The created figure
    """
    
    # Find hotspots (positions with high mutation frequency)
    hotspot_threshold = 0.1  # 10% of files
    hotspots = []
    
    for pos, counts in mutation_data.items():
        total_mutations = sum(counts.values())
        frequency = total_mutations / total_files
        if frequency >= hotspot_threshold:
            hotspots.append((pos, frequency, counts))
    
    # Sort by frequency
    hotspots.sort(key=lambda x: x[1], reverse=True)
    
    # Create figure
    fig = go.Figure()
    
    positions = [h[0] for h in hotspots[:20]]  # Top 20 hotspots
    frequencies = [h[1] * 100 for h in hotspots[:20]]
    mutation_types = [list(h[2].keys()) for h in hotspots[:20]]
    
    # Create color-coded bars
    colors = []
    for mut_types in mutation_types:
        if 'missense' in mut_types:
            colors.append('#1f77b4')
        elif 'silent' in mut_types:
            colors.append('#ff7f0e')
        elif 'deletion' in mut_types:
            colors.append('#2ca02c')
        elif 'frameshift' in mut_types:
            colors.append('#d62728')
        else:
            colors.append('#9467bd')
    
    fig.add_trace(go.Bar(
        x=positions,
        y=frequencies,
        marker_color=colors,
        text=[f"{freq:.1f}%" for freq in frequencies],
        textposition='auto',
        hovertemplate='Position: %{x}<br>Frequency: %{y:.1f}%<br>Types: %{customdata}<extra></extra>',
        customdata=mutation_types
    ))
    
    fig.update_layout(
        title={
            'text': f"Top 20 Mutation Hotspots in {virus_name}",
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 18, 'family': 'Arial Black'}
        },
        xaxis_title="Genome Position",
        yaxis_title="Mutation Frequency (%)",
        width=800,
        height=800,
        margin=dict(l=50, r=50, t=100, b=50)
    )
    
    return fig

def main():
    """Main function to create and save the genome landscape figures."""
    print("Creating Figure 3: Genome-Wide Mutation Landscape...")
    
    # Get available viruses
    result_dir = "../result"
    if not os.path.exists(result_dir):
        print(f"❌ Result directory not found: {result_dir}")
        return
    
    virus_dirs = [d for d in os.listdir(result_dir) 
                  if os.path.isdir(os.path.join(result_dir, d)) and not d.startswith('.')]
    
    if not virus_dirs:
        print("❌ No virus data found in result directory")
        return
    
    print(f"Available viruses: {virus_dirs}")
    
    # For now, use the first available virus (you can modify this to analyze all or specific ones)
    virus_name = virus_dirs[0]
    print(f"Analyzing {virus_name}...")
    
    # Prioritize SARS-CoV-2 if available
    if 'SARS-CoV-2' in virus_dirs:
        virus_name = 'SARS-CoV-2'
        print(f"Prioritizing {virus_name} due to comprehensive data...")
    
    # Get mutation data
    mutation_data, total_files = get_mutation_data(virus_name)
    
    if not mutation_data:
        print(f"❌ No mutation data found for {virus_name}")
        return
    
    print(f"Analyzed {total_files} {virus_name} samples")
    
    # Estimate genome length from mutation data
    genome_length = max(mutation_data.keys()) if mutation_data else 30000
    
    # Create main landscape figure
    fig1 = create_genome_landscape_figure(mutation_data, total_files, virus_name, genome_length)
    fig1.write_html(f"figure3a_{virus_name}_genome_landscape.html")
    fig1.write_image(f"figure3a_{virus_name}_genome_landscape.pdf", width=1000, height=1000)
    
    # Create hotspot figure
    fig2 = create_mutation_hotspot_figure(mutation_data, total_files, virus_name)
    fig2.write_html(f"figure3b_{virus_name}_mutation_hotspots.html")
    fig2.write_image(f"figure3b_{virus_name}_mutation_hotspots.pdf", width=800, height=800)
    
    print(f"✅ Figure 3 saved as:")
    print(f"   - figure3a_{virus_name}_genome_landscape.html")
    print(f"   - figure3a_{virus_name}_genome_landscape.pdf")
    print(f"   - figure3b_{virus_name}_mutation_hotspots.html")
    print(f"   - figure3b_{virus_name}_mutation_hotspots.pdf")
    
    # Print hotspot summary
    print(f"\n🔥 Top 10 Mutation Hotspots in {virus_name}:")
    hotspot_threshold = 0.1
    hotspots = []
    for pos, counts in mutation_data.items():
        total_mutations = sum(counts.values())
        frequency = total_mutations / total_files
        if frequency >= hotspot_threshold:
            hotspots.append((pos, frequency, counts))
    
    hotspots.sort(key=lambda x: x[1], reverse=True)
    for i, (pos, freq, counts) in enumerate(hotspots[:10]):
        print(f"  {i+1}. Position {pos}: {freq*100:.1f}% ({dict(counts)})")
    
    return fig1, fig2

if __name__ == "__main__":
    main()
