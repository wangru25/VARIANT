# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-08-19 17:30:00
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-19 17:30:00
Email: wang.rui@nyu.edu
FilePath: /VARIANT/plot/tmp/figure1_workflow_overview.py
Description: Figure 1: VARIANT Workflow Overview - creates a comprehensive workflow diagram showing the complete VARIANT pipeline from input data to final output, highlighting all major components and data flows.
'''

import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots
import numpy as np

# Set default template
pio.templates.default = "simple_white"

def create_workflow_overview():
    """
    Create a comprehensive workflow diagram for the VARIANT pipeline.
    """
    
    # Define workflow stages and their components
    stages = [
        {
            'name': 'Input Data',
            'components': ['FASTA Files', 'MSA Files', 'Reference Genomes'],
            'color': '#1f77b4',
            'x': 0.1
        },
        {
            'name': 'Data Processing',
            'components': ['Sequence Alignment', 'Quality Control', 'Format Validation'],
            'color': '#ff7f0e',
            'x': 0.25
        },
        {
            'name': 'Core Analysis',
            'components': ['Mutation Detection', 'Protein Translation', 'ORF Analysis'],
            'color': '#2ca02c',
            'x': 0.4
        },
        {
            'name': 'Specialized Analysis',
            'components': ['PRF Detection', 'Hot Mutation Analysis', 'Row Mutation Analysis'],
            'color': '#d62728',
            'x': 0.55
        },
        {
            'name': 'Output Generation',
            'components': ['Mutation Reports', 'CSV Files', 'Summary Statistics'],
            'color': '#9467bd',
            'x': 0.7
        },
        {
            'name': 'Visualization',
            'components': ['Genome Maps', 'Mutation Plots', 'Statistical Charts'],
            'color': '#8c564b',
            'x': 0.85
        }
    ]
    
    # Create figure
    fig = go.Figure()
    
    # Add workflow boxes
    for i, stage in enumerate(stages):
        # Main stage box
        fig.add_shape(
            type="rect",
            x0=stage['x'] - 0.08,
            y0=0.6,
            x1=stage['x'] + 0.08,
            y1=0.9,
            fillcolor=stage['color'],
            opacity=0.8,
            line=dict(color="black", width=2)
        )
        
        # Stage title
        fig.add_annotation(
            x=stage['x'],
            y=0.85,
            text=stage['name'],
            showarrow=False,
            font=dict(size=14, color="white", family="Arial Black"),
            xanchor="center",
            yanchor="middle"
        )
        
        # Component boxes
        for j, component in enumerate(stage['components']):
            y_pos = 0.5 - j * 0.08
            fig.add_shape(
                type="rect",
                x0=stage['x'] - 0.06,
                y0=y_pos - 0.02,
                x1=stage['x'] + 0.06,
                y1=y_pos + 0.02,
                fillcolor=stage['color'],
                opacity=0.6,
                line=dict(color="black", width=1)
            )
            
            fig.add_annotation(
                x=stage['x'],
                y=y_pos,
                text=component,
                showarrow=False,
                font=dict(size=10, color="black"),
                xanchor="center",
                yanchor="middle"
            )
    
    # Add connecting arrows
    for i in range(len(stages) - 1):
        fig.add_annotation(
            x=stages[i]['x'] + 0.08,
            y=0.75,
            xref="x",
            yref="y",
            axref="x",
            ayref="y",
            ax=stages[i+1]['x'] - 0.08,
            ay=0.75,
            arrowhead=2,
            arrowsize=1,
            arrowwidth=2,
            arrowcolor="black"
        )
    
    # Add data flow indicators
    data_flows = [
        {'from': (0.1, 0.3), 'to': (0.25, 0.3), 'label': 'Sequences'},
        {'from': (0.25, 0.2), 'to': (0.4, 0.2), 'label': 'Aligned Data'},
        {'from': (0.4, 0.1), 'to': (0.55, 0.1), 'label': 'Mutations'},
        {'from': (0.55, 0.0), 'to': (0.7, 0.0), 'label': 'Analysis Results'},
        {'from': (0.7, -0.1), 'to': (0.85, -0.1), 'label': 'Reports'}
    ]
    
    for flow in data_flows:
        fig.add_annotation(
            x=flow['from'][0],
            y=flow['from'][1],
            xref="x",
            yref="y",
            axref="x",
            ayref="y",
            ax=flow['to'][0],
            ay=flow['to'][1],
            arrowhead=1,
            arrowsize=0.8,
            arrowwidth=1,
            arrowcolor="gray",
            opacity=0.7
        )
        
        fig.add_annotation(
            x=(flow['from'][0] + flow['to'][0]) / 2,
            y=flow['from'][1] - 0.02,
            text=flow['label'],
            showarrow=False,
            font=dict(size=9, color="gray"),
            xanchor="center",
            yanchor="top"
        )
    
    # Add CLI interface indicator
    fig.add_shape(
        type="rect",
        x0=0.02,
        y0=0.02,
        x1=0.98,
        y1=0.15,
        fillcolor="lightblue",
        opacity=0.3,
        line=dict(color="blue", width=2, dash="dash")
    )
    
    fig.add_annotation(
        x=0.5,
        y=0.12,
        text="Command Line Interface",
        showarrow=False,
        font=dict(size=12, color="blue", family="Arial Black"),
        xanchor="center",
        yanchor="middle"
    )
    
    cli_commands = [
        "variant setup --virus <VIRUS_NAME>",
        "variant analyze --virus <VIRUS_NAME> --msa data/clustalW/msa.txt",
        "variant prf --genome data/refs/reference.fasta"
    ]
    
    for i, cmd in enumerate(cli_commands):
        fig.add_annotation(
            x=0.5,
            y=0.08 - i * 0.02,
            text=f"$ {cmd}",
            showarrow=False,
            font=dict(size=10, color="blue", family="Courier New"),
            xanchor="center",
            yanchor="middle"
        )
    
    # Update layout
    fig.update_layout(
        title={
            'text': "VARIANT: Comprehensive Viral Mutation Analysis Pipeline",
            'x': 0.5,
            'xanchor': 'center',
            'font': {'size': 20, 'family': 'Arial Black'}
        },
        xaxis=dict(
            showgrid=False,
            showticklabels=False,
            range=[0, 1],
            zeroline=False
        ),
        yaxis=dict(
            showgrid=False,
            showticklabels=False,
            range=[-0.2, 1],
            zeroline=False
        ),
        width=1000,
        height=800,
        margin=dict(l=50, r=50, t=100, b=50),
        showlegend=False,
        plot_bgcolor='white'
    )
    
    return fig

def main():
    """Main function to create and save the workflow figure."""
    print("Creating Figure 1: Workflow Overview...")
    
    # Create the workflow diagram
    fig = create_workflow_overview()
    
    # Save the figure
    fig.write_html("figure1_workflow_overview.html")
    fig.write_image("figure1_workflow_overview.pdf", width=1000, height=800)
    
    print("✅ Figure 1 saved as:")
    print("   - figure1_workflow_overview.html")
    print("   - figure1_workflow_overview.pdf")
    
    return fig

if __name__ == "__main__":
    main()
