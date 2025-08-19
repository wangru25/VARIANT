#!/usr/bin/env python3
"""
Figure 5: Hot Mutation Analysis
==============================

This script creates comprehensive analysis of hot mutations showing:
- Mutation frequency across genome positions
- Biological classification of mutations
- Hot mutation patterns and significance
- Sample distribution of hot mutations
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

def parse_hot_mutations(file_path):
    """
    Parse hot mutations from a row_hot_mutations CSV file.
    
    Args:
        file_path (str): Path to the hot mutations file
        
    Returns:
        list: List of hot mutation data
    """
    hot_mutations = []
    
    try:
        df = pd.read_csv(file_path)
        for _, row in df.iterrows():
            mutation_data = {
                'genome_id': row['genome_id'],
                'mutation_type': row['mutation_type'],
                'position': row['position'],
                'nucleotide_change': row['nucleotide_change'],
                'protein_affected': row['protein_affected'],
                'amino_acid_change': row['amino_acid_change'],
                'biological_classification': row['biological_classification']
            }
            hot_mutations.append(mutation_data)
    except FileNotFoundError:
        print(f"Warning: File not found {file_path}")
    except Exception as e:
        print(f"Warning: Error parsing {file_path}: {e}")
    
    return hot_mutations

def get_hot_mutation_data(virus_name):
    """
    Collect hot mutation data from virus result files.
    
    Args:
        virus_name (str): Name of the virus to analyze
        
    Returns:
        tuple: (hot_mutations, total_files)
    """
    result_dir = f"../result/{virus_name}"
    
    if not os.path.exists(result_dir):
        print(f"❌ Result directory not found: {result_dir}")
        return [], 0
    
    hot_mutation_files = glob.glob(os.path.join(result_dir, '*_row_hot_mutations_20250807.csv'))
    
    all_hot_mutations = []
    file_count = 0
    
    for file_path in hot_mutation_files:
        hot_mutations = parse_hot_mutations(file_path)
        all_hot_mutations.extend(hot_mutations)
        file_count += 1
    
    return all_hot_mutations, file_count

def create_hot_mutation_scatter(hot_mutations, total_files, virus_name):
    """
    Create a scatter plot showing hot mutation frequency and distribution.
    
    Args:
        hot_mutations (list): Hot mutation data
        total_files (int): Total number of files analyzed
        virus_name (str): Name of the virus
        
    Returns:
        plotly.graph_objects.Figure: The created figure
    """
    
    # Count mutations by position
    position_counts = defaultdict(int)
    position_types = defaultdict(list)
    position_proteins = defaultdict(list)
    
    for mutation in hot_mutations:
        pos_str = mutation['position']
        if ':' in pos_str:
            # Range mutation, use start position
            pos = int(pos_str.split(':')[0])
        else:
            pos = int(pos_str)
        
        position_counts[pos] += 1
        position_types[pos].append(mutation['mutation_type'])
        position_proteins[pos].append(mutation['protein_affected'])
    
    # Prepare data for plotting
    positions = list(position_counts.keys())
    frequencies = [position_counts[pos] for pos in positions]
    
    # Create color coding based on mutation type
    colors = []
    for pos in positions:
        types = position_types[pos]
        if 'hot' in types:
            colors.append('#d62728')  # Red for hot mutations
        elif 'row' in types:
            colors.append('#1f77b4')  # Blue for row mutations
        else:
            colors.append('#ff7f0e')  # Orange for others
    
    # Create size based on frequency
    sizes = [max(10, freq * 5) for freq in frequencies]
    
    # Create hover text
    hover_text = []
    for pos in positions:
        types = position_types[pos]
        proteins = list(set(position_proteins[pos]))
        text = f"Position: {pos}<br>Frequency: {position_counts[pos]}<br>Types: {', '.join(types)}<br>Proteins: {', '.join(proteins)}"
        hover_text.append(text)
    
    fig = go.Figure()
    
    fig.add_trace(go.Scatter(
        x=positions,
        y=frequencies,
        mode='markers',
        marker=dict(
            size=sizes,
            color=colors,
            opacity=0.7,
            line=dict(color='black', width=1)
        ),
        text=hover_text,
        hovertemplate='%{text}<extra></extra>',
        name='Hot Mutations'
    ))
    
    fig.update_layout(
        title={
            'text': f"Hot Mutation Frequency Across {virus_name} Genome",
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 18, 'family': 'Arial Black'}
        },
        xaxis_title="Genome Position",
        yaxis_title="Number of Samples with Mutation",
        width=800,
        height=800,
        margin=dict(l=50, r=50, t=100, b=50),
        showlegend=False
    )
    
    return fig

def create_mutation_type_distribution(hot_mutations, virus_name):
    """
    Create a chart showing distribution of hot mutation types.
    
    Args:
        hot_mutations (list): Hot mutation data
        virus_name (str): Name of the virus
        
    Returns:
        plotly.graph_objects.Figure: The created figure
    """
    
    # Count by mutation type
    type_counts = defaultdict(int)
    type_positions = defaultdict(list)
    
    for mutation in hot_mutations:
        mut_type = mutation['mutation_type']
        type_counts[mut_type] += 1
        type_positions[mut_type].append(mutation['position'])
    
    # Create bar chart
    types = list(type_counts.keys())
    counts = [type_counts[t] for t in types]
    
    colors = ['#d62728' if t == 'hot' else '#1f77b4' if t == 'row' else '#ff7f0e' for t in types]
    
    fig = go.Figure(data=go.Bar(
        x=types,
        y=counts,
        marker_color=colors,
        text=[f"{val}" for val in counts],
        textposition='auto',
        hovertemplate='Type: %{x}<br>Count: %{y}<extra></extra>'
    ))
    
    fig.update_layout(
        title={
            'text': f"Distribution of Hot Mutation Types - {virus_name}",
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 18, 'family': 'Arial Black'}
        },
        xaxis_title="Mutation Type",
        yaxis_title="Number of Mutations",
        width=800,
        height=800,
        margin=dict(l=50, r=50, t=100, b=50)
    )
    
    return fig

def create_protein_impact_analysis(hot_mutations, virus_name):
    """
    Create analysis of hot mutations by protein.
    
    Args:
        hot_mutations (list): Hot mutation data
        virus_name (str): Name of the virus
        
    Returns:
        plotly.graph_objects.Figure: The created figure
    """
    
    # Count by protein
    protein_counts = defaultdict(int)
    protein_mutations = defaultdict(list)
    
    for mutation in hot_mutations:
        protein = mutation['protein_affected']
        protein_counts[protein] += 1
        protein_mutations[protein].append(mutation)
    
    # Sort by count
    sorted_proteins = sorted(protein_counts.items(), key=lambda x: x[1], reverse=True)
    
    proteins = [p[0] for p in sorted_proteins[:15]]  # Top 15 proteins
    counts = [p[1] for p in sorted_proteins[:15]]
    
    # Create color coding
    colors = []
    for protein in proteins:
        protein_lower = protein.lower()
        if 'spike' in protein_lower:
            colors.append('#d62728')  # Red for spike
        elif 'polymerase' in protein_lower or 'rdrp' in protein_lower:
            colors.append('#1f77b4')  # Blue for polymerase
        elif 'nsp' in protein_lower:
            colors.append('#ff7f0e')  # Orange for NSPs
        elif 'orf' in protein_lower:
            colors.append('#2ca02c')  # Green for ORFs
        else:
            colors.append('#9467bd')  # Purple for others
    
    fig = go.Figure(data=go.Bar(
        x=proteins,
        y=counts,
        marker_color=colors,
        text=[f"{val}" for val in counts],
        textposition='auto',
        hovertemplate='Protein: %{x}<br>Hot Mutations: %{y}<extra></extra>'
    ))
    
    fig.update_layout(
        title={
            'text': f"Hot Mutations by Protein (Top 15) - {virus_name}",
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 18, 'family': 'Arial Black'}
        },
        xaxis_title="Protein",
        yaxis_title="Number of Hot Mutations",
        width=800,
        height=800,
        margin=dict(l=50, r=50, t=100, b=50),
        xaxis_tickangle=-45
    )
    
    return fig

def create_biological_classification_analysis(hot_mutations, virus_name):
    """
    Create analysis of biological classification of hot mutations.
    
    Args:
        hot_mutations (list): Hot mutation data
        virus_name (str): Name of the virus
        
    Returns:
        plotly.graph_objects.Figure: The created figure
    """
    
    # Count by biological classification
    bio_counts = defaultdict(int)
    bio_details = defaultdict(list)
    
    for mutation in hot_mutations:
        bio_class = mutation['biological_classification']
        bio_counts[bio_class] += 1
        bio_details[bio_class].append(mutation)
    
    # Create pie chart
    labels = list(bio_counts.keys())
    values = [bio_counts[label] for label in labels]
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    
    fig = go.Figure(data=go.Pie(
        labels=labels,
        values=values,
        marker_colors=colors[:len(labels)],
        hole=0.3,
        textinfo='label+percent+value',
        textposition='inside'
    ))
    
    fig.update_layout(
        title={
            'text': f"Biological Classification of Hot Mutations - {virus_name}",
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 18, 'family': 'Arial Black'}
        },
        width=800,
        height=800,
        margin=dict(l=50, r=50, t=100, b=50),
        showlegend=True
    )
    
    return fig

def create_sample_distribution_heatmap(hot_mutations, total_files, virus_name):
    """
    Create a heatmap showing sample distribution of hot mutations.
    
    Args:
        hot_mutations (list): Hot mutation data
        total_files (int): Total number of files analyzed
        virus_name (str): Name of the virus
        
    Returns:
        plotly.graph_objects.Figure: The created figure
    """
    
    # Group mutations by position and count samples
    position_samples = defaultdict(set)
    position_types = defaultdict(set)
    
    for mutation in hot_mutations:
        pos_str = mutation['position']
        if ':' in pos_str:
            pos = int(pos_str.split(':')[0])
        else:
            pos = int(pos_str)
        
        position_samples[pos].add(mutation['genome_id'])
        position_types[pos].add(mutation['mutation_type'])
    
    # Get top positions by sample count
    position_counts = [(pos, len(samples)) for pos, samples in position_samples.items()]
    position_counts.sort(key=lambda x: x[1], reverse=True)
    
    top_positions = [p[0] for p in position_counts[:20]]  # Top 20 positions
    
    # Create sample presence matrix
    sample_ids = set()
    for mutation in hot_mutations:
        sample_ids.add(mutation['genome_id'])
    
    sample_ids = sorted(list(sample_ids))
    
    # Create presence matrix
    presence_matrix = []
    for pos in top_positions:
        row = []
        for sample_id in sample_ids:
            if sample_id in position_samples[pos]:
                row.append(1)
            else:
                row.append(0)
        presence_matrix.append(row)
    
    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=presence_matrix,
        x=sample_ids,
        y=[f"Pos {pos}" for pos in top_positions],
        colorscale='Blues',
        showscale=True,
        hovertemplate='Position: %{y}<br>Sample: %{x}<br>Present: %{z}<extra></extra>'
    ))
    
    fig.update_layout(
        title={
            'text': f"Sample Distribution of Top 20 Hot Mutation Positions - {virus_name}",
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 18, 'family': 'Arial Black'}
        },
        xaxis_title="Sample ID",
        yaxis_title="Genome Position",
        width=1000,
        height=800,
        margin=dict(l=50, r=50, t=100, b=50)
    )
    
    return fig

def main():
    """Main function to create and save the hot mutation figures."""
    print("Creating Figure 5: Hot Mutation Analysis...")
    
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
    
    # Get hot mutation data
    hot_mutations, total_files = get_hot_mutation_data(virus_name)
    
    if not hot_mutations:
        print(f"❌ No hot mutation data found for {virus_name}")
        return
    
    print(f"Analyzed {total_files} {virus_name} samples")
    print(f"Found {len(hot_mutations)} hot mutations")
    
    # Create scatter plot
    fig1 = create_hot_mutation_scatter(hot_mutations, total_files, virus_name)
    fig1.write_html(f"figure5a_{virus_name}_hot_mutation_scatter.html")
    fig1.write_image(f"figure5a_{virus_name}_hot_mutation_scatter.pdf", width=800, height=800)
    
    # Create mutation type distribution
    fig2 = create_mutation_type_distribution(hot_mutations, virus_name)
    fig2.write_html(f"figure5b_{virus_name}_mutation_type_distribution.html")
    fig2.write_image(f"figure5b_{virus_name}_mutation_type_distribution.pdf", width=800, height=800)
    
    # Create protein impact analysis
    fig3 = create_protein_impact_analysis(hot_mutations, virus_name)
    fig3.write_html(f"figure5c_{virus_name}_protein_impact_analysis.html")
    fig3.write_image(f"figure5c_{virus_name}_protein_impact_analysis.pdf", width=800, height=800)
    
    # Create biological classification analysis
    fig4 = create_biological_classification_analysis(hot_mutations, virus_name)
    fig4.write_html(f"figure5d_{virus_name}_biological_classification.html")
    fig4.write_image(f"figure5d_{virus_name}_biological_classification.pdf", width=800, height=800)
    
    # Create sample distribution heatmap
    fig5 = create_sample_distribution_heatmap(hot_mutations, total_files, virus_name)
    fig5.write_html(f"figure5e_{virus_name}_sample_distribution_heatmap.html")
    fig5.write_image(f"figure5e_{virus_name}_sample_distribution_heatmap.pdf", width=1000, height=800)
    
    print(f"✅ Figure 5 saved as:")
    print(f"   - figure5a_{virus_name}_hot_mutation_scatter.html")
    print(f"   - figure5a_{virus_name}_hot_mutation_scatter.pdf")
    print(f"   - figure5b_{virus_name}_mutation_type_distribution.html")
    print(f"   - figure5b_{virus_name}_mutation_type_distribution.pdf")
    print(f"   - figure5c_{virus_name}_protein_impact_analysis.html")
    print(f"   - figure5c_{virus_name}_protein_impact_analysis.pdf")
    print(f"   - figure5d_{virus_name}_biological_classification.html")
    print(f"   - figure5d_{virus_name}_biological_classification.pdf")
    print(f"   - figure5e_{virus_name}_sample_distribution_heatmap.html")
    print(f"   - figure5e_{virus_name}_sample_distribution_heatmap.pdf")
    
    # Print summary statistics
    print(f"\n🔥 Hot Mutation Summary for {virus_name}:")
    
    # Count by type
    type_counts = defaultdict(int)
    for mutation in hot_mutations:
        type_counts[mutation['mutation_type']] += 1
    
    for mut_type, count in type_counts.items():
        print(f"  {mut_type.capitalize()}: {count} mutations")
    
    # Count by protein
    protein_counts = defaultdict(int)
    for mutation in hot_mutations:
        protein_counts[mutation['protein_affected']] += 1
    
    print(f"\n📊 Top 5 Proteins with Hot Mutations in {virus_name}:")
    sorted_proteins = sorted(protein_counts.items(), key=lambda x: x[1], reverse=True)
    for i, (protein, count) in enumerate(sorted_proteins[:5]):
        print(f"  {i+1}. {protein}: {count} hot mutations")
    
    return fig1, fig2, fig3, fig4, fig5

if __name__ == "__main__":
    main()
