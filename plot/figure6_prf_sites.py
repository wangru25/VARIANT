#!/usr/bin/env python3
"""
Figure 6: PRF Site Distribution
==============================

This script creates comprehensive analysis of PRF (Programmed Ribosomal Frameshifting) sites showing:
- PRF site distribution across the genome
- Sequence context and types (-1, +1 PRF)
- Overlapping gene regions
- Biological significance
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

def parse_prf_data(file_path):
    """
    Parse PRF data from the potential_PRF CSV file.
    
    Args:
        file_path (str): Path to the PRF file
        
    Returns:
        list: List of PRF site data
    """
    prf_sites = []
    
    try:
        df = pd.read_csv(file_path)
        for _, row in df.iterrows():
            prf_data = {
                'position': row['position'],
                'end_position': row['end_position'],
                'sequence': row['sequence'],
                'type': row['type']
            }
            prf_sites.append(prf_data)
    except FileNotFoundError:
        print(f"Warning: File not found {file_path}")
    except Exception as e:
        print(f"Warning: Error parsing {file_path}: {e}")
    
    return prf_sites

def get_prf_data(virus_name):
    """
    Collect PRF data from virus result files.
    
    Args:
        virus_name (str): Name of the virus to analyze
        
    Returns:
        list: PRF site data
    """
    result_dir = f"../result/{virus_name}"
    
    if not os.path.exists(result_dir):
        print(f"❌ Result directory not found: {result_dir}")
        return []
    
    prf_files = glob.glob(os.path.join(result_dir, '*potential_PRF_20250807.csv'))
    
    all_prf_sites = []
    
    for file_path in prf_files:
        prf_sites = parse_prf_data(file_path)
        all_prf_sites.extend(prf_sites)
    
    return all_prf_sites

def create_prf_genome_map(prf_sites, virus_name):
    """
    Create a genome map showing PRF site distribution.
    
    Args:
        prf_sites (list): PRF site data
        virus_name (str): Name of the virus
        
    Returns:
        plotly.graph_objects.Figure: The created figure
    """
    
    # Get virus-specific genome regions
    genome_regions = get_genome_regions(virus_name)
    
    # Create subplots
    fig = make_subplots(
        rows=2, cols=1,
        subplot_titles=(
            f"PRF Site Distribution Across {virus_name} Genome",
            "PRF Sites by Type and Protein Region"
        ),
        specs=[[{"secondary_y": False}],
               [{"secondary_y": False}]],
        vertical_spacing=0.15,
        row_heights=[0.6, 0.4]
    )
    
    # Plot 1: Genome-wide PRF distribution
    positions = []
    end_positions = []
    sequences = []
    types = []
    colors = []
    
    for site in prf_sites:
        positions.append(site['position'])
        end_positions.append(site['end_position'])
        sequences.append(site['sequence'])
        types.append(site['type'])
        
        # Color coding by PRF type
        if '-1 PRF' in site['type']:
            colors.append('#d62728')  # Red for -1 PRF
        elif '+1 PRF' in site['type']:
            colors.append('#1f77b4')  # Blue for +1 PRF
        else:
            colors.append('#ff7f0e')  # Orange for others
    
    # Create scatter plot for PRF sites
    fig.add_trace(
        go.Scatter(
            x=positions,
            y=[1] * len(positions),
            mode='markers',
            marker=dict(
                size=15,
                color=colors,
                symbol='diamond',
                line=dict(color='black', width=2)
            ),
            text=[f"Pos: {pos}<br>Seq: {seq}<br>Type: {type_}" 
                  for pos, seq, type_ in zip(positions, sequences, types)],
            hovertemplate='%{text}<extra></extra>',
            name='PRF Sites'
        ),
        row=1, col=1
    )
    
    # Add protein region annotations
    region_colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', 
                     '#8c564b', '#e377c2', '#7f7f7f', '#bcbd22', '#17becf']
    
    for i, (region_name, (start, end)) in enumerate(genome_regions.items()):
        fig.add_shape(
            type="rect",
            x0=start,
            y0=0.5,
            x1=end,
            y1=1.5,
            fillcolor=region_colors[i % len(region_colors)],
            opacity=0.3,
            line=dict(color="black", width=1),
            row=1, col=1
        )
        
        # Add region label
        fig.add_annotation(
            x=(start + end) / 2,
            y=0.75,
            text=region_name.replace('_', '<br>'),
            showarrow=False,
            font=dict(size=8, color="black"),
            xanchor="center",
            yanchor="middle",
            row=1, col=1
        )
    
    # Plot 2: PRF sites by type and region
    prf_by_region = defaultdict(lambda: defaultdict(int))
    
    for site in prf_sites:
        pos = site['position']
        prf_type = site['type']
        
        # Find which region this PRF site is in
        region_found = False
        for region_name, (start, end) in genome_regions.items():
            if start <= pos <= end:
                prf_by_region[region_name][prf_type] += 1
                region_found = True
                break
        
        if not region_found:
            prf_by_region['Intergenic'][prf_type] += 1
    
    # Create stacked bar chart
    regions = list(prf_by_region.keys())
    prf_types = ['-1 PRF', '+1 PRF']
    
    for prf_type in prf_types:
        values = [prf_by_region[region].get(prf_type, 0) for region in regions]
        color = '#d62728' if prf_type == '-1 PRF' else '#1f77b4'
        
        fig.add_trace(
            go.Bar(
                name=prf_type,
                x=regions,
                y=values,
                marker_color=color,
                opacity=0.8
            ),
            row=2, col=1
        )
    
    # Update layout
    fig.update_layout(
        title={
            'text': f"PRF Site Distribution and Analysis - {virus_name}",
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
    fig.update_xaxes(title_text="Genome Position", row=1, col=1)
    fig.update_yaxes(title_text="PRF Sites", row=1, col=1, range=[0, 2])
    fig.update_xaxes(title_text="Protein Region", row=2, col=1)
    fig.update_yaxes(title_text="Number of PRF Sites", row=2, col=1)
    
    return fig

def create_prf_sequence_analysis(prf_sites, virus_name):
    """
    Create analysis of PRF sequence patterns.
    
    Args:
        prf_sites (list): PRF site data
        virus_name (str): Name of the virus
        
    Returns:
        plotly.graph_objects.Figure: The created figure
    """
    
    # Analyze sequence patterns
    sequence_patterns = defaultdict(int)
    type_sequences = defaultdict(list)
    
    for site in prf_sites:
        seq = site['sequence']
        prf_type = site['type']
        sequence_patterns[seq] += 1
        type_sequences[prf_type].append(seq)
    
    # Create subplots
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            "PRF Sequence Frequency",
            "PRF Types Distribution",
            "Sequence Length Analysis",
            "PRF Site Position Distribution"
        ),
        specs=[[{"type": "bar"}, {"type": "pie"}],
               [{"type": "histogram"}, {"type": "scatter"}]]
    )
    
    # Plot 1: Sequence frequency
    sequences = list(sequence_patterns.keys())
    frequencies = [sequence_patterns[seq] for seq in sequences]
    
    fig.add_trace(
        go.Bar(
            x=sequences,
            y=frequencies,
            marker_color='#1f77b4',
            text=[f"{val}" for val in frequencies],
            textposition='auto',
            hovertemplate='Sequence: %{x}<br>Frequency: %{y}<extra></extra>'
        ),
        row=1, col=1
    )
    
    # Plot 2: PRF type distribution
    type_counts = defaultdict(int)
    for site in prf_sites:
        type_counts[site['type']] += 1
    
    fig.add_trace(
        go.Pie(
            labels=list(type_counts.keys()),
            values=list(type_counts.values()),
            marker_colors=['#d62728', '#1f77b4'],
            hole=0.3
        ),
        row=1, col=2
    )
    
    # Plot 3: Sequence length analysis
    seq_lengths = [len(site['sequence']) for site in prf_sites]
    
    fig.add_trace(
        go.Histogram(
            x=seq_lengths,
            nbinsx=10,
            marker_color='#2ca02c',
            name='Sequence Length'
        ),
        row=2, col=1
    )
    
    # Plot 4: Position distribution
    positions = [site['position'] for site in prf_sites]
    types = [site['type'] for site in prf_sites]
    
    colors = ['#d62728' if '-1 PRF' in t else '#1f77b4' for t in types]
    
    fig.add_trace(
        go.Scatter(
            x=positions,
            y=[1] * len(positions),
            mode='markers',
            marker=dict(
                size=10,
                color=colors,
                opacity=0.7
            ),
            text=[f"Pos: {pos}<br>Type: {type_}" for pos, type_ in zip(positions, types)],
            hovertemplate='%{text}<extra></extra>',
            name='PRF Positions'
        ),
        row=2, col=2
    )
    
    fig.update_layout(
        title={
            'text': f"PRF Sequence and Pattern Analysis - {virus_name}",
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 18, 'family': 'Arial Black'}
        },
        width=1000,
        height=800,
        margin=dict(l=50, r=50, t=100, b=50),
        showlegend=False
    )
    
    return fig

def create_prf_biological_significance(prf_sites, virus_name):
    """
    Create analysis of biological significance of PRF sites.
    
    Args:
        prf_sites (list): PRF site data
        virus_name (str): Name of the virus
        
    Returns:
        plotly.graph_objects.Figure: The created figure
    """
    
    # Get virus-specific genome regions
    genome_regions = get_genome_regions(virus_name)
    
    # Analyze PRF sites by protein region and biological function
    region_analysis = defaultdict(lambda: defaultdict(int))
    functional_impact = defaultdict(int)
    
    for site in prf_sites:
        pos = site['position']
        prf_type = site['type']
        
        # Find protein region
        for region_name, (start, end) in genome_regions.items():
            if start <= pos <= end:
                region_analysis[region_name][prf_type] += 1
                
                # Assess functional impact
                region_lower = region_name.lower()
                if 'spike' in region_lower:
                    functional_impact['Spike Protein'] += 1
                elif 'polymerase' in region_lower or 'rdrp' in region_lower:
                    functional_impact['RNA Polymerase'] += 1
                elif 'nsp' in region_lower:
                    functional_impact['Non-structural Proteins'] += 1
                elif 'orf' in region_lower:
                    functional_impact['Accessory Proteins'] += 1
                else:
                    functional_impact['Other Proteins'] += 1
                break
        else:
            region_analysis['Intergenic'][prf_type] += 1
            functional_impact['Intergenic'] += 1
    
    # Create figure
    fig = go.Figure()
    
    # Create stacked bar chart
    regions = list(region_analysis.keys())
    prf_types = ['-1 PRF', '+1 PRF']
    
    for prf_type in prf_types:
        values = [region_analysis[region].get(prf_type, 0) for region in regions]
        color = '#d62728' if prf_type == '-1 PRF' else '#1f77b4'
        
        fig.add_trace(go.Bar(
            name=prf_type,
            x=regions,
            y=values,
            marker_color=color,
            opacity=0.8,
            hovertemplate='Region: %{x}<br>Type: %{name}<br>Count: %{y}<extra></extra>'
        ))
    
    fig.update_layout(
        title={
            'text': f"Biological Significance of PRF Sites by Protein Region - {virus_name}",
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 18, 'family': 'Arial Black'}
        },
        xaxis_title="Protein Region",
        yaxis_title="Number of PRF Sites",
        barmode='stack',
        width=800,
        height=800,
        margin=dict(l=50, r=50, t=100, b=50),
        showlegend=True,
        xaxis_tickangle=-45
    )
    
    return fig

def main():
    """Main function to create and save the PRF figures."""
    print("Creating Figure 6: PRF Site Distribution...")
    
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
    
    # Get PRF data
    prf_sites = get_prf_data(virus_name)
    
    if not prf_sites:
        print(f"❌ No PRF data found for {virus_name}")
        return
    
    print(f"Found {len(prf_sites)} PRF sites")
    
    # Create genome map
    fig1 = create_prf_genome_map(prf_sites, virus_name)
    fig1.write_html(f"figure6a_{virus_name}_prf_genome_map.html")
    fig1.write_image(f"figure6a_{virus_name}_prf_genome_map.pdf", width=1000, height=1000)
    
    # Create sequence analysis
    fig2 = create_prf_sequence_analysis(prf_sites, virus_name)
    fig2.write_html(f"figure6b_{virus_name}_prf_sequence_analysis.html")
    fig2.write_image(f"figure6b_{virus_name}_prf_sequence_analysis.pdf", width=1000, height=800)
    
    # Create biological significance analysis
    fig3 = create_prf_biological_significance(prf_sites, virus_name)
    fig3.write_html(f"figure6c_{virus_name}_prf_biological_significance.html")
    fig3.write_image(f"figure6c_{virus_name}_prf_biological_significance.pdf", width=800, height=800)
    
    print(f"✅ Figure 6 saved as:")
    print(f"   - figure6a_{virus_name}_prf_genome_map.html")
    print(f"   - figure6a_{virus_name}_prf_genome_map.pdf")
    print(f"   - figure6b_{virus_name}_prf_sequence_analysis.html")
    print(f"   - figure6b_{virus_name}_prf_sequence_analysis.pdf")
    print(f"   - figure6c_{virus_name}_prf_biological_significance.html")
    print(f"   - figure6c_{virus_name}_prf_biological_significance.pdf")
    
    # Print PRF summary
    print(f"\n🔄 PRF Site Summary for {virus_name}:")
    
    # Count by type
    type_counts = defaultdict(int)
    for site in prf_sites:
        type_counts[site['type']] += 1
    
    for prf_type, count in type_counts.items():
        print(f"  {prf_type}: {count} sites")
    
    # Count by sequence
    sequence_counts = defaultdict(int)
    for site in prf_sites:
        sequence_counts[site['sequence']] += 1
    
    print(f"\n🧬 Top PRF Sequences in {virus_name}:")
    sorted_sequences = sorted(sequence_counts.items(), key=lambda x: x[1], reverse=True)
    for i, (seq, count) in enumerate(sorted_sequences[:5]):
        print(f"  {i+1}. {seq}: {count} occurrences")
    
    # Analyze by protein region
    genome_regions = get_genome_regions(virus_name)
    print(f"\n📊 PRF Sites by Protein Region in {virus_name}:")
    region_counts = defaultdict(int)
    for site in prf_sites:
        pos = site['position']
        for region_name, (start, end) in genome_regions.items():
            if start <= pos <= end:
                region_counts[region_name] += 1
                break
        else:
            region_counts['Intergenic'] += 1
    
    sorted_regions = sorted(region_counts.items(), key=lambda x: x[1], reverse=True)
    for region, count in sorted_regions:
        print(f"  {region}: {count} PRF sites")
    
    return fig1, fig2, fig3

if __name__ == "__main__":
    main()
