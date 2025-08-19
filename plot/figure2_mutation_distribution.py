#!/usr/bin/env python3
"""
Figure 2: Mutation Type Distribution Across Viruses
==================================================

This script creates a stacked bar chart comparing mutation types (missense, silent, 
deletion, frameshift) across different viral families. Generic version that works with any virus data.
"""

import plotly.graph_objects as go
import plotly.io as pio
import pandas as pd
import numpy as np
import os
import glob
from collections import defaultdict
from plotly.subplots import make_subplots

# Set default template
pio.templates.default = "simple_white"

def parse_mutation_file(file_path):
    """
    Parse a mutation file and count mutation types.
    
    Args:
        file_path (str): Path to the mutation file
        
    Returns:
        dict: Count of each mutation type
    """
    mutation_counts = defaultdict(int)
    
    try:
        with open(file_path, 'r') as f:
            for line in f:
                line = line.strip()
                if line:
                    # Extract mutation type from the beginning of the line
                    parts = line.split()
                    if parts:
                        mutation_type = parts[0]
                        mutation_counts[mutation_type] += 1
    except FileNotFoundError:
        print(f"Warning: File not found {file_path}")
        return {}
    
    return dict(mutation_counts)

def get_virus_data():
    """
    Collect mutation data from all virus result directories.
    Automatically detects available virus data.
    
    Returns:
        dict: Virus data with mutation counts
    """
    virus_data = {}
    result_dir = "../result"
    
    if not os.path.exists(result_dir):
        print(f"❌ Result directory not found: {result_dir}")
        return {}
    
    # Automatically detect virus directories
    virus_dirs = [d for d in os.listdir(result_dir) 
                  if os.path.isdir(os.path.join(result_dir, d)) and not d.startswith('.')]
    
    print(f"Found virus directories: {virus_dirs}")
    
    for virus_dir in virus_dirs:
        virus_path = os.path.join(result_dir, virus_dir)
        
        # Find all mutation files (excluding summary and PRF files)
        mutation_files = []
        
        # Check for segment subdirectories (like H3N2)
        segment_dirs = [d for d in os.listdir(virus_path) 
                       if os.path.isdir(os.path.join(virus_path, d)) and 'segment' in d.lower()]
        
        if segment_dirs:
            # Handle segmented viruses like H3N2
            for segment_dir in segment_dirs:
                segment_path = os.path.join(virus_path, segment_dir)
                segment_files = glob.glob(os.path.join(segment_path, '*_20250807.txt'))
                mutation_files.extend([f for f in segment_files 
                                     if 'mutation_summary' not in f and 'potential_PRF' not in f])
        else:
            # Handle non-segmented viruses
            all_files = glob.glob(os.path.join(virus_path, '*_20250807.txt'))
            mutation_files = [f for f in all_files 
                            if 'mutation_summary' not in f and 'potential_PRF' not in f]
        
        # Aggregate mutation counts from all files
        total_counts = defaultdict(int)
        file_count = 0
        
        for file_path in mutation_files:
            counts = parse_mutation_file(file_path)
            for mutation_type, count in counts.items():
                total_counts[mutation_type] += count
            file_count += 1
        
        if total_counts:
            virus_data[virus_dir] = {
                'counts': dict(total_counts),
                'file_count': file_count
            }
    
    return virus_data

def create_mutation_distribution_figure(virus_data):
    """
    Create a stacked bar chart showing mutation type distribution across viruses.
    
    Args:
        virus_data (dict): Virus data with mutation counts
        
    Returns:
        plotly.graph_objects.Figure: The created figure
    """
    
    # Define mutation types and colors
    mutation_types = ['missense', 'silent', 'deletion', 'frameshift', 'insertion']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    
    # Prepare data for plotting
    viruses = list(virus_data.keys())
    
    if not viruses:
        print("❌ No virus data available for plotting")
        return None
    
    # Create traces for each mutation type
    traces = []
    
    for i, mutation_type in enumerate(mutation_types):
        values = []
        for virus in viruses:
            if virus in virus_data:
                counts = virus_data[virus]['counts']
                values.append(counts.get(mutation_type, 0))
            else:
                values.append(0)
        
        traces.append(go.Bar(
            name=mutation_type.capitalize(),
            x=viruses,
            y=values,
            marker_color=colors[i],
            opacity=0.8,
            hovertemplate='<b>%{x}</b><br>' +
                         f'{mutation_type.capitalize()}: %{{y}}<br>' +
                         '<extra></extra>'
        ))
    
    # Create figure
    fig = go.Figure(data=traces)
    
    # Update layout
    fig.update_layout(
        title={
            'text': "Mutation Type Distribution Across Viral Families",
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 18, 'family': 'Arial Black'}
        },
        xaxis_title="Viral Family",
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
        )
    )
    
    # Add sample size annotations
    for i, virus in enumerate(viruses):
        if virus in virus_data:
            file_count = virus_data[virus]['file_count']
            total_mutations = sum(virus_data[virus]['counts'].values())
            
            fig.add_annotation(
                x=i,
                y=total_mutations + max([virus_data[v]['counts'].get(mt, 0) for v in viruses for mt in mutation_types]) * 0.05,
                text=f"n={file_count}<br>Total: {total_mutations}",
                showarrow=False,
                font=dict(size=10, color="gray"),
                xanchor="center",
                yanchor="bottom"
            )
    
    return fig

def create_mutation_type_pie_charts(virus_data):
    """
    Create individual pie charts for each virus showing mutation type proportions.
    
    Args:
        virus_data (dict): Virus data with mutation counts
        
    Returns:
        plotly.graph_objects.Figure: The created figure
    """
    
    mutation_types = ['missense', 'silent', 'deletion', 'frameshift', 'insertion']
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd']
    
    # Calculate number of rows and columns for subplots
    n_viruses = len(virus_data)
    if n_viruses == 0:
        print("❌ No virus data available for pie charts")
        return None
    
    n_cols = min(3, n_viruses)
    n_rows = (n_viruses + n_cols - 1) // n_cols
    
    fig = make_subplots(
        rows=n_rows, cols=n_cols,
        subplot_titles=list(virus_data.keys()),
        specs=[[{"type": "pie"}] * n_cols] * n_rows
    )
    
    virus_list = list(virus_data.keys())
    
    for i, virus in enumerate(virus_list):
        row = (i // n_cols) + 1
        col = (i % n_cols) + 1
        
        counts = virus_data[virus]['counts']
        values = [counts.get(mt, 0) for mt in mutation_types]
        labels = [f"{mt.capitalize()}<br>({val})" for mt, val in zip(mutation_types, values)]
        
        # Only show labels for non-zero values
        non_zero_indices = [j for j, val in enumerate(values) if val > 0]
        non_zero_values = [values[j] for j in non_zero_indices]
        non_zero_labels = [labels[j] for j in non_zero_indices]
        non_zero_colors = [colors[j] for j in non_zero_indices]
        
        fig.add_trace(
            go.Pie(
                labels=non_zero_labels,
                values=non_zero_values,
                marker_colors=non_zero_colors,
                hole=0.3,
                textinfo='label+percent',
                textposition='inside'
            ),
            row=row, col=col
        )
    
    fig.update_layout(
        title={
            'text': "Mutation Type Proportions by Viral Family",
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 18, 'family': 'Arial Black'}
        },
        width=800,
        height=800,
        margin=dict(l=50, r=50, t=100, b=50),
        showlegend=False
    )
    
    return fig

def main():
    """Main function to create and save the mutation distribution figures."""
    print("Creating Figure 2: Mutation Type Distribution Across Viruses...")
    
    # Get virus data
    virus_data = get_virus_data()
    
    if not virus_data:
        print("❌ No virus data found. Please ensure result files exist.")
        return
    
    print(f"Found data for {len(virus_data)} viruses: {list(virus_data.keys())}")
    
    # Create stacked bar chart
    fig1 = create_mutation_distribution_figure(virus_data)
    if fig1:
        fig1.write_html("figure2a_mutation_distribution_stacked.html")
        fig1.write_image("figure2a_mutation_distribution_stacked.pdf", width=800, height=800)
    
    # Create pie charts
    fig2 = create_mutation_type_pie_charts(virus_data)
    if fig2:
        fig2.write_html("figure2b_mutation_distribution_pie.html")
        fig2.write_image("figure2b_mutation_distribution_pie.pdf", width=800, height=800)
    
    print("✅ Figure 2 saved as:")
    print("   - figure2a_mutation_distribution_stacked.html")
    print("   - figure2a_mutation_distribution_stacked.pdf")
    print("   - figure2b_mutation_distribution_pie.html")
    print("   - figure2b_mutation_distribution_pie.pdf")
    
    # Print summary statistics
    print("\n📊 Summary Statistics:")
    for virus, data in virus_data.items():
        total = sum(data['counts'].values())
        print(f"  {virus}: {data['file_count']} samples, {total} total mutations")
        for mt, count in data['counts'].items():
            if count > 0:
                percentage = (count / total) * 100
                print(f"    {mt.capitalize()}: {count} ({percentage:.1f}%)")
    
    return fig1, fig2

if __name__ == "__main__":
    main()
