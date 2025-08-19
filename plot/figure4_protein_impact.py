#!/usr/bin/env python3
"""
Figure 4: Protein-Level Mutation Impact
=====================================

This script creates a comprehensive analysis of protein-level mutation impact showing:
- Which proteins are most affected by mutations
- Distribution of mutation types across proteins
- Functional impact assessment
- Conservation analysis
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
import csv

# Set default template
pio.templates.default = "simple_white"

def parse_protein_mutations(file_path):
    """
    Parse protein-level mutations from a mutation summary CSV file.
    
    Args:
        file_path (str): Path to the mutation summary CSV file
        
    Returns:
        dict: Protein mutation data
    """
    protein_mutations = defaultdict(lambda: defaultdict(int))
    
    try:
        with open(file_path, 'r') as f:
            reader = csv.DictReader(f)
            for row in reader:
                protein_name = row['Protein Name']
                mutation_type = row['AA Mutation Type']
                
                if protein_name and mutation_type:
                    protein_mutations[protein_name][mutation_type] += 1
    except FileNotFoundError:
        print(f"Warning: File not found {file_path}")
    except Exception as e:
        print(f"Warning: Error parsing {file_path}: {e}")
    
    return dict(protein_mutations)

def get_protein_mutation_data(virus_name):
    """
    Collect protein-level mutation data from virus result files.
    
    Args:
        virus_name (str): Name of the virus to analyze
        
    Returns:
        tuple: (protein_data, total_files)
    """
    result_dir = f"../result/{virus_name}"
    
    if not os.path.exists(result_dir):
        print(f"❌ Result directory not found: {result_dir}")
        return {}, 0
    
    # Look for mutation summary CSV files
    mutation_summary_files = glob.glob(os.path.join(result_dir, '*mutation_summary_20250807.csv'))
    
    # If no mutation summary files, try to generate them
    if not mutation_summary_files:
        print(f"No mutation summary files found for {virus_name}, attempting to generate...")
        try:
            # Import the mutation summary function
            import sys
            sys.path.append('../src')
            from utils.mutation_summary import extract_mutation_summary_to_csv
            
            # Find available txt files
            txt_files = glob.glob(os.path.join(result_dir, '*_20250807.txt'))
            for txt_file in txt_files[:5]:  # Process first 5 files
                genome_id = os.path.basename(txt_file).replace('_20250807.txt', '')
                try:
                    extract_mutation_summary_to_csv(virus_name, genome_id)
                except Exception as e:
                    print(f"Warning: Could not generate summary for {genome_id}: {e}")
            
            # Check again for generated files
            mutation_summary_files = glob.glob(os.path.join(result_dir, '*mutation_summary_20250807.csv'))
        except Exception as e:
            print(f"Could not generate mutation summaries: {e}")
            return {}, 0
    
    # Collect all protein mutations
    all_protein_mutations = defaultdict(lambda: defaultdict(int))
    file_count = 0
    
    for file_path in mutation_summary_files:
        protein_mutations = parse_protein_mutations(file_path)
        for protein, mutations in protein_mutations.items():
            for mutation_type, count in mutations.items():
                all_protein_mutations[protein][mutation_type] += count
        file_count += 1
    
    return dict(all_protein_mutations), file_count

def create_protein_impact_heatmap(protein_data, total_files, virus_name):
    """
    Create a heatmap showing mutation impact across proteins.
    
    Args:
        protein_data (dict): Protein mutation data
        total_files (int): Total number of files analyzed
        virus_name (str): Name of the virus
        
    Returns:
        plotly.graph_objects.Figure: The created figure
    """
    
    # Define mutation types and proteins
    mutation_types = ['missense', 'silent', 'deletion', 'frameshift', 'insertion']
    
    # Filter proteins with significant mutations
    significant_proteins = []
    for protein, mutations in protein_data.items():
        total_mutations = sum(mutations.values())
        if total_mutations >= 5:  # At least 5 mutations
            significant_proteins.append(protein)
    
    # Sort proteins by total mutation count
    significant_proteins.sort(key=lambda p: sum(protein_data[p].values()), reverse=True)
    
    # Prepare data for heatmap
    z_data = []
    for protein in significant_proteins:
        row = []
        for mutation_type in mutation_types:
            count = protein_data[protein].get(mutation_type, 0)
            # Normalize by total files
            frequency = (count / total_files) * 100
            row.append(frequency)
        z_data.append(row)
    
    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=z_data,
        x=mutation_types,
        y=significant_proteins,
        colorscale='Viridis',
        text=[[f"{val:.1f}%" for val in row] for row in z_data],
        texttemplate="%{text}",
        textfont={"size": 10},
        hoverongaps=False,
        hovertemplate='Protein: %{y}<br>Mutation Type: %{x}<br>Frequency: %{z:.2f}%<extra></extra>'
    ))
    
    fig.update_layout(
        title={
            'text': f"Protein-Level Mutation Impact Heatmap - {virus_name}",
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 18, 'family': 'Arial Black'}
        },
        xaxis_title="Mutation Type",
        yaxis_title="Protein",
        width=800,
        height=800,
        margin=dict(l=50, r=50, t=100, b=50)
    )
    
    return fig

def create_protein_mutation_summary(protein_data, total_files, virus_name):
    """
    Create a summary chart showing total mutations per protein.
    
    Args:
        protein_data (dict): Protein mutation data
        total_files (int): Total number of files analyzed
        virus_name (str): Name of the virus
        
    Returns:
        plotly.graph_objects.Figure: The created figure
    """
    
    # Calculate total mutations per protein
    protein_totals = []
    for protein, mutations in protein_data.items():
        total = sum(mutations.values())
        if total > 0:
            protein_totals.append((protein, total))
    
    # Sort by total mutations
    protein_totals.sort(key=lambda x: x[1], reverse=True)
    
    proteins = [p[0] for p in protein_totals[:20]]  # Top 20 proteins
    totals = [p[1] for p in protein_totals[:20]]
    
    # Create color coding based on protein function
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
        y=totals,
        marker_color=colors,
        text=[f"{val}" for val in totals],
        textposition='auto',
        hovertemplate='Protein: %{x}<br>Total Mutations: %{y}<extra></extra>'
    ))
    
    fig.update_layout(
        title={
            'text': f"Total Mutations by Protein (Top 20) - {virus_name}",
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 18, 'family': 'Arial Black'}
        },
        xaxis_title="Protein",
        yaxis_title="Total Mutations",
        width=800,
        height=800,
        margin=dict(l=50, r=50, t=100, b=50),
        xaxis_tickangle=-45
    )
    
    return fig

def create_mutation_type_distribution(protein_data, total_files, virus_name):
    """
    Create a stacked bar chart showing mutation type distribution across proteins.
    
    Args:
        protein_data (dict): Protein mutation data
        total_files (int): Total number of files analyzed
        virus_name (str): Name of the virus
        
    Returns:
        plotly.graph_objects.Figure: The created figure
    """
    
    # Get top 10 proteins by total mutations
    protein_totals = []
    for protein, mutations in protein_data.items():
        total = sum(mutations.values())
        if total > 0:
            protein_totals.append((protein, total))
    
    protein_totals.sort(key=lambda x: x[1], reverse=True)
    top_proteins = [p[0] for p in protein_totals[:10]]
    
    # Define mutation types and colors
    mutation_types = ['missense', 'silent', 'deletion', 'frameshift', 'insertion']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    
    # Create traces for each mutation type
    traces = []
    
    for i, mutation_type in enumerate(mutation_types):
        values = []
        for protein in top_proteins:
            count = protein_data[protein].get(mutation_type, 0)
            values.append(count)
        
        traces.append(go.Bar(
            name=mutation_type.capitalize(),
            x=top_proteins,
            y=values,
            marker_color=colors[i],
            opacity=0.8,
            hovertemplate='<b>%{x}</b><br>' +
                         f'{mutation_type.capitalize()}: %{{y}}<br>' +
                         '<extra></extra>'
        ))
    
    fig = go.Figure(data=traces)
    
    fig.update_layout(
        title={
            'text': f"Mutation Type Distribution Across Top 10 Proteins - {virus_name}",
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 18, 'family': 'Arial Black'}
        },
        xaxis_title="Protein",
        yaxis_title="Number of Mutations",
        barmode='stack',
        width=800,
        height=800,
        margin=dict(l=50, r=50, t=100, b=50),
        showlegend=True,
        legend=dict(
            orientation="h",
            yanchor="bottom",
            y=1.02,
            xanchor="right",
            x=1
        ),
        xaxis_tickangle=-45
    )
    
    return fig

def create_functional_impact_radar(protein_data, total_files, virus_name):
    """
    Create a radar chart showing functional impact scores.
    
    Args:
        protein_data (dict): Protein mutation data
        total_files (int): Total number of files analyzed
        virus_name (str): Name of the virus
        
    Returns:
        plotly.graph_objects.Figure: The created figure
    """
    
    # Define impact categories and their weights
    impact_categories = {
        'Mutation Frequency': 0.3,
        'Functional Impact': 0.4,
        'Conservation': 0.2,
        'Structural Impact': 0.1
    }
    
    # Calculate impact scores for top proteins
    protein_totals = []
    for protein, mutations in protein_data.items():
        total = sum(mutations.values())
        if total > 0:
            protein_totals.append((protein, total))
    
    protein_totals.sort(key=lambda x: x[1], reverse=True)
    top_proteins = [p[0] for p in protein_totals[:8]]  # Top 8 for radar chart
    
    # Calculate impact scores
    impact_scores = {}
    max_mutations = max([p[1] for p in protein_totals]) if protein_totals else 1
    
    for protein in top_proteins:
        mutations = protein_data[protein]
        total = sum(mutations.values())
        
        # Mutation frequency score (normalized)
        freq_score = (total / max_mutations) * 100
        
        # Functional impact score (weighted by mutation types)
        functional_score = 0
        if total > 0:
            missense_ratio = mutations.get('missense', 0) / total
            frameshift_ratio = mutations.get('frameshift', 0) / total
            deletion_ratio = mutations.get('deletion', 0) / total
            functional_score = (missense_ratio * 0.6 + frameshift_ratio * 0.8 + deletion_ratio * 0.7) * 100
        
        # Conservation score (inverse of silent mutations)
        silent_ratio = mutations.get('silent', 0) / total if total > 0 else 0
        conservation_score = (1 - silent_ratio) * 100
        
        # Structural impact score (based on frameshifts and deletions)
        structural_score = ((mutations.get('frameshift', 0) + mutations.get('deletion', 0)) / total) * 100 if total > 0 else 0
        
        impact_scores[protein] = {
            'Mutation Frequency': freq_score,
            'Functional Impact': functional_score,
            'Conservation': conservation_score,
            'Structural Impact': structural_score
        }
    
    # Create radar chart
    fig = go.Figure()
    
    categories = list(impact_categories.keys())
    
    for protein in top_proteins:
        scores = [impact_scores[protein][cat] for cat in categories]
        scores.append(scores[0])  # Close the polygon
        
        fig.add_trace(go.Scatterpolar(
            r=scores,
            theta=categories + [categories[0]],
            fill='toself',
            name=protein,
            line=dict(width=2)
        ))
    
    fig.update_layout(
        polar=dict(
            radialaxis=dict(
                visible=True,
                range=[0, 100]
            )),
        title={
            'text': f"Functional Impact Assessment of Top Proteins - {virus_name}",
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

def main():
    """Main function to create and save the protein impact figures."""
    print("Creating Figure 4: Protein-Level Mutation Impact...")
    
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
    
    # Get protein mutation data
    protein_data, total_files = get_protein_mutation_data(virus_name)
    
    if not protein_data:
        print(f"❌ No protein mutation data found for {virus_name}")
        # Try SARS-CoV-2 if available
        if 'SARS-CoV-2' in virus_dirs and virus_name != 'SARS-CoV-2':
            virus_name = 'SARS-CoV-2'
            print(f"Trying {virus_name}...")
            protein_data, total_files = get_protein_mutation_data(virus_name)
            if not protein_data:
                print(f"❌ No protein mutation data found for {virus_name}")
                return
        else:
            return
    
    print(f"Analyzed {total_files} {virus_name} samples")
    print(f"Found mutations in {len(protein_data)} proteins")
    
    # Create heatmap
    fig1 = create_protein_impact_heatmap(protein_data, total_files, virus_name)
    fig1.write_html(f"figure4a_{virus_name}_protein_impact_heatmap.html")
    fig1.write_image(f"figure4a_{virus_name}_protein_impact_heatmap.pdf", width=800, height=800)
    
    # Create summary chart
    fig2 = create_protein_mutation_summary(protein_data, total_files, virus_name)
    fig2.write_html(f"figure4b_{virus_name}_protein_mutation_summary.html")
    fig2.write_image(f"figure4b_{virus_name}_protein_mutation_summary.pdf", width=800, height=800)
    
    # Create mutation type distribution
    fig3 = create_mutation_type_distribution(protein_data, total_files, virus_name)
    fig3.write_html(f"figure4c_{virus_name}_mutation_type_distribution.html")
    fig3.write_image(f"figure4c_{virus_name}_mutation_type_distribution.pdf", width=800, height=800)
    
    # Create functional impact radar
    fig4 = create_functional_impact_radar(protein_data, total_files, virus_name)
    fig4.write_html(f"figure4d_{virus_name}_functional_impact_radar.html")
    fig4.write_image(f"figure4d_{virus_name}_functional_impact_radar.pdf", width=800, height=800)
    
    print(f"✅ Figure 4 saved as:")
    print(f"   - figure4a_{virus_name}_protein_impact_heatmap.html")
    print(f"   - figure4a_{virus_name}_protein_impact_heatmap.pdf")
    print(f"   - figure4b_{virus_name}_protein_mutation_summary.html")
    print(f"   - figure4b_{virus_name}_protein_mutation_summary.pdf")
    print(f"   - figure4c_{virus_name}_mutation_type_distribution.html")
    print(f"   - figure4c_{virus_name}_mutation_type_distribution.pdf")
    print(f"   - figure4d_{virus_name}_functional_impact_radar.html")
    print(f"   - figure4d_{virus_name}_functional_impact_radar.pdf")
    
    # Print protein summary
    print(f"\n📊 Top 10 Most Mutated Proteins in {virus_name}:")
    protein_totals = []
    for protein, mutations in protein_data.items():
        total = sum(mutations.values())
        if total > 0:
            protein_totals.append((protein, total, mutations))
    
    protein_totals.sort(key=lambda x: x[1], reverse=True)
    for i, (protein, total, mutations) in enumerate(protein_totals[:10]):
        print(f"  {i+1}. {protein}: {total} mutations ({dict(mutations)})")
    
    return fig1, fig2, fig3, fig4

if __name__ == "__main__":
    main()