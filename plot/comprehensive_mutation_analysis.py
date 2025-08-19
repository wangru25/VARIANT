#!/usr/bin/env python3
"""
Comprehensive Mutation Analysis Visualization
Generates multiple plots analyzing mutation patterns across all viruses in the result folder.
"""

import os
import glob
import pandas as pd
import numpy as np
import plotly.graph_objects as go
import plotly.express as px
import plotly.io as pio
from plotly.subplots import make_subplots
from collections import defaultdict, Counter
import re
from pathlib import Path

# Set default template
pio.templates.default = "simple_white"


def load_all_mutation_data():
    """Load all mutation summary CSV files from the result directory."""
    mutation_files = glob.glob("result/**/*_mutation_summary.csv", recursive=True)
    
    all_data = []
    virus_info = []
    
    for file_path in mutation_files:
        try:
            # Extract virus name and genome ID from file path
            parts = file_path.split('/')
            virus_name = parts[1]  # result/virus_name/...
            genome_id = os.path.basename(file_path).replace('_mutation_summary.csv', '')
            
            # Handle multi-segment viruses (H3N2)
            segment = None
            if len(parts) > 3 and parts[2].startswith('segment_'):
                segment = parts[2]
                virus_name = f"{virus_name}_{segment}"
            
            # Load CSV data
            df = pd.read_csv(file_path)
            df['virus'] = virus_name
            df['genome_id'] = genome_id
            df['segment'] = segment
            
            all_data.append(df)
            virus_info.append({
                'virus': virus_name,
                'genome_id': genome_id,
                'segment': segment,
                'file_path': file_path,
                'mutation_count': len(df)
            })
            
        except Exception as e:
            print(f"Error loading {file_path}: {e}")
    
    if all_data:
        combined_df = pd.concat(all_data, ignore_index=True)
        virus_summary = pd.DataFrame(virus_info)
        return combined_df, virus_summary
    else:
        return pd.DataFrame(), pd.DataFrame()


def create_mutation_type_distribution_plot(data_df):
    """Create a stacked bar plot showing mutation type distribution across viruses."""
    
    # Count mutation types by virus
    mutation_counts = data_df.groupby(['virus', 'AA Mutation Type']).size().unstack(fill_value=0)
    
    # Calculate percentages
    mutation_percentages = mutation_counts.div(mutation_counts.sum(axis=1), axis=0) * 100
    
    # Create stacked bar plot
    fig = go.Figure()
    
    colors = ['#1f77b4', '#ff7f0e', '#2ca02c', '#d62728', '#9467bd', '#8c564b']
    
    for i, mutation_type in enumerate(mutation_percentages.columns):
        fig.add_trace(go.Bar(
            name=mutation_type,
            x=mutation_percentages.index,
            y=mutation_percentages[mutation_type],
            marker_color=colors[i % len(colors)],
            opacity=0.8
        ))
    
    fig.update_layout(
        title="Mutation Type Distribution Across Viruses",
        xaxis_title="Virus",
        yaxis_title="Percentage (%)",
        barmode='stack',
        width=1000,
        height=800,
        margin=dict(l=50, r=50, t=50, b=50),
        showlegend=True
    )
    
    # Save plots
    fig.write_html("plot/mutation_type_distribution.html")
    fig.write_image("plot/mutation_type_distribution.pdf", width=1000, height=800)
    
    return fig


def create_protein_mutation_heatmap(data_df):
    """Create a heatmap showing mutation frequency by protein across viruses."""
    
    # Count mutations by protein and virus
    protein_counts = data_df.groupby(['virus', 'Protein Name']).size().unstack(fill_value=0)
    
    # Fill NaN values with 0
    protein_counts = protein_counts.fillna(0)
    
    # Create heatmap
    fig = go.Figure(data=go.Heatmap(
        z=protein_counts.values,
        x=protein_counts.columns,
        y=protein_counts.index,
        colorscale='Viridis',
        showscale=True,
        colorbar=dict(title="Mutation Count")
    ))
    
    fig.update_layout(
        title="Protein Mutation Frequency Heatmap",
        xaxis_title="Protein",
        yaxis_title="Virus",
        width=1200,
        height=800,
        margin=dict(l=50, r=50, t=50, b=50)
    )
    
    # Save plots
    fig.write_html("plot/protein_mutation_heatmap.html")
    fig.write_image("plot/protein_mutation_heatmap.pdf", width=1200, height=800)
    
    return fig


def create_mutation_count_comparison(virus_summary):
    """Create a bar plot comparing total mutation counts across viruses."""
    
    # Group by virus (combining segments for H3N2)
    virus_totals = virus_summary.groupby('virus')['mutation_count'].sum().sort_values(ascending=False)
    
    fig = go.Figure()
    
    fig.add_trace(go.Bar(
        x=virus_totals.index,
        y=virus_totals.values,
        marker_color='#1f77b4',
        opacity=0.8
    ))
    
    fig.update_layout(
        title="Total Mutation Count by Virus",
        xaxis_title="Virus",
        yaxis_title="Total Mutations",
        width=800,
        height=800,
        margin=dict(l=50, r=50, t=50, b=50)
    )
    
    # Save plots
    fig.write_html("plot/mutation_count_comparison.html")
    fig.write_image("plot/mutation_count_comparison.pdf", width=800, height=800)
    
    return fig


def create_nt_mutation_analysis(data_df):
    """Create plots analyzing nucleotide mutation patterns."""
    
    # Extract nucleotide changes
    data_df['nt_change'] = data_df['NT Mutation'].str.extract(r'([A-Z]+->[A-Z]+)')
    
    # Count nucleotide changes
    nt_counts = data_df['nt_change'].value_counts().head(10)
    
    fig = go.Figure()
    
    fig.add_trace(go.Bar(
        x=nt_counts.index,
        y=nt_counts.values,
        marker_color='#2ca02c',
        opacity=0.8
    ))
    
    fig.update_layout(
        title="Top 10 Nucleotide Changes",
        xaxis_title="Nucleotide Change",
        yaxis_title="Count",
        width=800,
        height=800,
        margin=dict(l=50, r=50, t=50, b=50)
    )
    
    # Save plots
    fig.write_html("plot/nt_mutation_analysis.html")
    fig.write_image("plot/nt_mutation_analysis.pdf", width=800, height=800)
    
    return fig


def create_virus_comparison_dashboard(data_df, virus_summary):
    """Create a comprehensive dashboard comparing viruses."""
    
    # Create subplots
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            'Mutation Count by Virus',
            'Mutation Type Distribution',
            'Top Proteins by Mutation Count',
            'Nucleotide Change Patterns'
        ),
        specs=[[{"secondary_y": False}, {"secondary_y": False}],
               [{"secondary_y": False}, {"secondary_y": False}]]
    )
    
    # Plot 1: Mutation count by virus
    virus_totals = virus_summary.groupby('virus')['mutation_count'].sum().sort_values(ascending=False)
    fig.add_trace(
        go.Bar(x=virus_totals.index, y=virus_totals.values, name="Total Mutations", marker_color='#1f77b4'),
        row=1, col=1
    )
    
    # Plot 2: Mutation type distribution
    mutation_counts = data_df['AA Mutation Type'].value_counts()
    fig.add_trace(
        go.Pie(labels=mutation_counts.index, values=mutation_counts.values, name="Mutation Types"),
        row=1, col=2
    )
    
    # Plot 3: Top proteins
    protein_counts = data_df['Protein Name'].value_counts().head(10)
    fig.add_trace(
        go.Bar(x=protein_counts.index, y=protein_counts.values, name="Protein Mutations", marker_color='#ff7f0e'),
        row=2, col=1
    )
    
    # Plot 4: Nucleotide changes
    nt_counts = data_df['NT Mutation'].str.extract(r'([A-Z]+->[A-Z]+)')[0].value_counts().head(8)
    fig.add_trace(
        go.Bar(x=nt_counts.index, y=nt_counts.values, name="NT Changes", marker_color='#2ca02c'),
        row=2, col=2
    )
    
    fig.update_layout(
        title="Comprehensive Virus Mutation Analysis Dashboard",
        width=1200,
        height=1200,
        margin=dict(l=50, r=50, t=50, b=50),
        showlegend=False
    )
    
    # Save plots
    fig.write_html("plot/virus_comparison_dashboard.html")
    fig.write_image("plot/virus_comparison_dashboard.pdf", width=1200, height=1200)
    
    return fig


def create_sars_cov2_detailed_analysis(data_df):
    """Create detailed analysis plots for SARS-CoV-2."""
    
    # Filter SARS-CoV-2 data
    sars_data = data_df[data_df['virus'].str.contains('SARS-CoV-2', na=False)]
    
    if len(sars_data) == 0:
        print("No SARS-CoV-2 data found")
        return None
    
    # Create subplots for SARS-CoV-2 analysis
    fig = make_subplots(
        rows=2, cols=2,
        subplot_titles=(
            'SARS-CoV-2 Protein Mutation Distribution',
            'SARS-CoV-2 Mutation Types',
            'SARS-CoV-2 Top Mutated Proteins',
            'SARS-CoV-2 Nucleotide Changes'
        ),
        specs=[[{"secondary_y": False}, {"secondary_y": False}],
               [{"secondary_y": False}, {"secondary_y": False}]]
    )
    
    # Plot 1: Protein distribution
    protein_counts = sars_data['Protein Name'].value_counts().head(15)
    fig.add_trace(
        go.Bar(x=protein_counts.index, y=protein_counts.values, marker_color='#1f77b4'),
        row=1, col=1
    )
    
    # Plot 2: Mutation types
    mutation_counts = sars_data['AA Mutation Type'].value_counts()
    fig.add_trace(
        go.Pie(labels=mutation_counts.index, values=mutation_counts.values),
        row=1, col=2
    )
    
    # Plot 3: Top mutated proteins
    top_proteins = sars_data['Protein Name'].value_counts().head(10)
    fig.add_trace(
        go.Bar(x=top_proteins.index, y=top_proteins.values, marker_color='#ff7f0e'),
        row=2, col=1
    )
    
    # Plot 4: Nucleotide changes
    nt_counts = sars_data['NT Mutation'].str.extract(r'([A-Z]+->[A-Z]+)')[0].value_counts().head(10)
    fig.add_trace(
        go.Bar(x=nt_counts.index, y=nt_counts.values, marker_color='#2ca02c'),
        row=2, col=2
    )
    
    fig.update_layout(
        title="SARS-CoV-2 Detailed Mutation Analysis",
        width=1200,
        height=1200,
        margin=dict(l=50, r=50, t=50, b=50),
        showlegend=False
    )
    
    # Save plots
    fig.write_html("plot/sars_cov2_detailed_analysis.html")
    fig.write_image("plot/sars_cov2_detailed_analysis.pdf", width=1200, height=1200)
    
    return fig


def main():
    """Main function to generate all plots."""
    print("🔬 Loading mutation data...")
    
    # Load all mutation data
    data_df, virus_summary = load_all_mutation_data()
    
    if len(data_df) == 0:
        print("❌ No mutation data found!")
        return
    
    print(f"✅ Loaded {len(data_df)} mutations from {len(virus_summary)} genomes")
    print(f"📊 Viruses found: {', '.join(virus_summary['virus'].unique())}")
    
    # Create plots directory if it doesn't exist
    os.makedirs("plot", exist_ok=True)
    
    print("\n📈 Generating plots...")
    
    # Generate all plots
    plots = {}
    
    print("  📊 Creating mutation type distribution plot...")
    plots['mutation_type_dist'] = create_mutation_type_distribution_plot(data_df)
    
    print("  🔥 Creating protein mutation heatmap...")
    plots['protein_heatmap'] = create_protein_mutation_heatmap(data_df)
    
    print("  📊 Creating mutation count comparison...")
    plots['mutation_count'] = create_mutation_count_comparison(virus_summary)
    
    print("  🧬 Creating nucleotide mutation analysis...")
    plots['nt_analysis'] = create_nt_mutation_analysis(data_df)
    
    print("  📋 Creating virus comparison dashboard...")
    plots['dashboard'] = create_virus_comparison_dashboard(data_df, virus_summary)
    
    print("  🦠 Creating SARS-CoV-2 detailed analysis...")
    plots['sars_detailed'] = create_sars_cov2_detailed_analysis(data_df)
    
    print("\n✅ All plots generated successfully!")
    print("\n📁 Generated files:")
    print("  - mutation_type_distribution.html/pdf")
    print("  - protein_mutation_heatmap.html/pdf")
    print("  - mutation_count_comparison.html/pdf")
    print("  - nt_mutation_analysis.html/pdf")
    print("  - virus_comparison_dashboard.html/pdf")
    print("  - sars_cov2_detailed_analysis.html/pdf")
    
    # Print summary statistics
    print(f"\n📊 Summary Statistics:")
    print(f"  Total mutations: {len(data_df)}")
    print(f"  Total genomes: {len(virus_summary)}")
    print(f"  Viruses analyzed: {len(virus_summary['virus'].unique())}")
    print(f"  Unique proteins: {data_df['Protein Name'].nunique()}")
    print(f"  Mutation types: {', '.join(data_df['AA Mutation Type'].unique())}")


if __name__ == "__main__":
    main()
