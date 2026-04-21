# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-08-29 08:51:53
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-29 10:00:26
Email: rw3594@nyu.edu
FilePath: /VARIANT/src/visualization/figure3_PRF.py
Description: Figure 3 - Plot PRF regions on genome organization bars, using PRF records as input.
'''

import os
import glob
import csv
import json
import re
import argparse
import pandas as pd
import yaml
from typing import Dict, List, Optional, Union, Tuple
from Bio import SeqIO
import plotly.graph_objects as go
import plotly.io as pio
from plotly.subplots import make_subplots

# Set default template
pio.templates.default = "simple_white"

class PRFVisualizer:
    def __init__(self, virus_name: str, sample_id: str = None):
        self.virus_name = virus_name
        self.sample_id = sample_id
        
        # Auto-detect paths
        self.proteome_path = self._auto_detect_proteome_path()
        self.prf_path = self._auto_detect_prf_path()
        self.reference_genome_path = self._auto_detect_reference_genome_path()
        
        # Load data
        self.proteins = self._load_proteins()
        self.prf_records = self._parse_prf_records()
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
    
    def _auto_detect_prf_path(self) -> str:
        """Auto-detect PRF records file path."""
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
            specific_paths = [
                os.path.join(result_dir, f"{self.sample_id}_prf*.csv"),
                os.path.join(result_dir, f"{self.sample_id}_prf*.json")
            ]
            for pattern in specific_paths:
                matches = glob.glob(pattern)
                if matches:
                    return matches[0]
        
        # Look for any PRF files (including custom-virus outputs like
        # "prf_analysis.prf_candidates.csv" and legacy "potential_PRF.csv").
        prf_files = (
            glob.glob(os.path.join(result_dir, "*prf*.csv")) +
            glob.glob(os.path.join(result_dir, "*prf*.json")) +
            glob.glob(os.path.join(result_dir, "potential_PRF.csv"))
        )
        if not prf_files:
            raise FileNotFoundError(f"No potential PRF record files found in {result_dir}")
        
        return prf_files[0]
    
    def _auto_detect_reference_genome_path(self) -> Optional[str]:
        """Auto-detect reference genome file path."""
        # Prefer the explicit reference genome from virus_config.yaml.
        try:
            with open("virus_config.yaml", "r") as f:
                config = yaml.safe_load(f) or {}
            virus_cfg = (config.get("viruses") or {}).get(self.virus_name, {})
            configured_ref = virus_cfg.get("reference_genome")
            if configured_ref:
                configured_path = f"data/{self.virus_name}/refs/{configured_ref}"
                if os.path.exists(configured_path):
                    return configured_path
        except Exception:
            pass

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
    
    def _parse_prf_records(self) -> List[Dict]:
        """Parse PRF records from file."""
        return self._parse_prf_records_from_path(self.prf_path)
    
    def _parse_prf_records_from_path(self, prf_path: str) -> List[Dict]:
        """Parse PRF records from a specific file path."""
        records = []
        
        try:
            if prf_path.endswith('.json'):
                with open(prf_path, 'r') as f:
                    data = json.load(f)
                items = data.get('sites', data if isinstance(data, list) else [])
                for item in items:
                    records.append(self._normalize_prf_record(item))
            else:
                with open(prf_path, 'r') as f:
                    reader = csv.DictReader(f)
                    for row in reader:
                        records.append(self._normalize_prf_record(row))
        except Exception as e:
            print(f"Error reading PRF file {prf_path}: {e}")
            return []
        
        return [r for r in records if r.get('start') is not None and r.get('end') is not None]
    
    def _normalize_prf_record(self, item: Dict) -> Dict:
        """Normalize PRF record to standard format."""
        def to_int(val) -> Optional[int]:
            if val is None:
                return None
            if isinstance(val, int):
                return val
            if isinstance(val, str):
                m = re.findall(r'\d+', val)
                if m:
                    return int(m[0])
            return None
        
        motif = (item.get('sequence') or item.get('slippery_motif') or item.get('motif') or '').strip()
        start = to_int(
            item.get('start') or
            item.get('position') or
            item.get('site_start_1based')
        )
        end = to_int(
            item.get('end') or
            item.get('end_position') or
            item.get('site_end_1based')
        )
        if end is None and start is not None and motif:
            end = start + len(motif) - 1
        if end is None and start is not None:
            end = start
        
        prf_type = (item.get('type') or item.get('site_type') or '').strip()
        normalized_type = prf_type
        if prf_type in ['-1', '-1PRF', '-1 prf']:
            normalized_type = '-1 PRF'
        elif prf_type in ['+1', '+1PRF', '+1 prf']:
            normalized_type = '+1 PRF'
        elif not prf_type:
            normalized_type = '-1 PRF'

        return {
            'start': start,
            'end': end,
            'sequence': motif,
            'type': normalized_type,
            'protein': item.get('protein') or item.get('gene')
        }
    
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
    
    def _get_prf_color(self, prf_type: str):
        """Get color for PRF type."""
        mapping = {
            '-1 PRF': ('rgba(255, 0, 0, 1.0)', 0.8),      # Bright red
            '+1 PRF': ('rgba(0, 0, 255, 1.0)', 0.8),      # Bright blue
            'potential -1 PRF': ('rgba(255, 0, 0, 1.0)', 0.8),      # Bright red
            'potential +1 PRF': ('rgba(0, 0, 255, 1.0)', 0.8),      # Bright blue
        }
        return mapping.get(prf_type, ('rgba(255, 165, 0, 1.0)', 0.8))  # Orange for others

    def _add_prf_legend(self, fig, prf_types: set) -> None:
        """Add PRF type legend in the top-left."""
        # Keep order consistent
        ordered_types = ['-1 PRF', '+1 PRF']
        ordered_types = [t for t in ordered_types if t in prf_types]
        ordered_types += [t for t in sorted(prf_types) if t not in ordered_types]
        
        legend_text = "<b>Potential PRF Types:</b><br>"
        for t in ordered_types:
            color, _ = self._get_prf_color(t)
            # Use a colored rectangle '█' plus colored label text
            legend_text += f"<span style='color:{color};'>█</span> <span style='color:{color};'>{t}</span><br>"
        
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
    
    def create_prf_chart(self, output_path: str = None) -> go.Figure:
        """Create genome organization chart with PRF regions."""
        
        # Auto-detect genome ID from PRF filename
        genome_id = os.path.basename(self.prf_path)
        genome_id = re.sub(r'_prf.*', '', genome_id)
        genome_id = genome_id.replace('.csv', '').replace('.json', '')
        
        # Auto-generate output path if not provided
        if output_path is None:
            os.makedirs('imgs/visualizations', exist_ok=True)
            virus_dir = f'imgs/visualizations/{self.virus_name}'
            os.makedirs(virus_dir, exist_ok=True)
            output_path = f'{virus_dir}/{genome_id}_prf_regions.html'
        
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
        
        # Add PRF regions
        prf_data = []
        
        # Check if this is a multi-segment virus
        segment_dirs = glob.glob(f"result/{self.virus_name}/segment_*")
        
        if segment_dirs:
            # Multi-segment virus - add PRF regions from all segments
            for segment_dir in sorted(segment_dirs):
                segment_name = os.path.basename(segment_dir)
                prf_files = glob.glob(os.path.join(segment_dir, "*_prf*.csv")) + glob.glob(os.path.join(segment_dir, "*_prf*.json"))
                
                if prf_files:
                    segment_prf_records = self._parse_prf_records_from_path(prf_files[0])
                    
                    for prf_record in segment_prf_records:
                        protein_name = f"{segment_name}_{prf_record['protein']}" if prf_record.get('protein') else None
                        gene = None
                        
                        # Try to find matching gene
                        if protein_name:
                            for g in self.genome_organization['genes']:
                                if g['name'] == protein_name:
                                    gene = g
                                    break
                        
                        if gene is None:
                            # Try to find gene by position
                            for g in self.genome_organization['genes']:
                                if 'original_start' in g and g['original_start'] <= prf_record['start'] <= g['original_end']:
                                    gene = g
                                    break
                            if gene is None:
                                for g in self.genome_organization['genes']:
                                    if g['start'] <= prf_record['start'] <= g['end']:
                                        gene = g
                                        break
                        
                        if gene is not None:
                            # Calculate position on the adjusted bar
                            if 'original_start' in gene:
                                x0 = gene['start'] + max(0, prf_record['start'] - gene['original_start'])
                                x1 = gene['start'] + max(0, prf_record['end'] - gene['original_start'])
                            else:
                                x0 = max(gene['start'], prf_record['start'])
                                x1 = min(gene['end'], prf_record['end'])
                            
                            if x1 < x0:
                                x0, x1 = x1, x0
                            
                            gene_index = None
                            for i, g in enumerate(all_genes):
                                if g['name'] == gene['name']:
                                    gene_index = i
                                    break
                            
                            if gene_index is not None:
                                bar_height = 0.6
                                gap_size = 0.6
                                y_center = gene_index * (bar_height + gap_size)
                                
                                prf_data.append({
                                    'x0': x0,
                                    'x1': x1,
                                    'y_center': y_center,
                                    'type': prf_record['type'],
                                    'sequence': prf_record['sequence'],
                                    'start': prf_record['start'],
                                    'end': prf_record['end'],
                                    'gene': gene['name']
                                })
        else:
            # Single-segment virus - use existing logic
            for prf_record in self.prf_records:
                protein_name = prf_record.get('protein')
                gene = None
                
                # Try to find matching gene
                if protein_name:
                    for g in self.genome_organization['genes']:
                        if g['name'] == protein_name:
                            gene = g
                            break
                
                if gene is None:
                    # Try to find gene by position
                    for g in self.genome_organization['genes']:
                        if g['start'] <= prf_record['start'] <= g['end']:
                            gene = g
                            break
                
                if gene is not None:
                    x0 = max(gene['start'], prf_record['start'])
                    x1 = min(gene['end'], prf_record['end'])
                    
                    if x1 < x0:
                        x0, x1 = x1, x0
                    
                    gene_index = None
                    for i, g in enumerate(all_genes):
                        if g['name'] == gene['name']:
                            gene_index = i
                            break
                    
                    if gene_index is not None:
                        bar_height = 0.6
                        gap_size = 0.6
                        y_center = gene_index * (bar_height + gap_size)
                        
                        prf_data.append({
                            'x0': x0,
                            'x1': x1,
                            'y_center': y_center,
                            'type': prf_record['type'],
                            'sequence': prf_record['sequence'],
                            'start': prf_record['start'],
                            'end': prf_record['end'],
                            'gene': gene['name']
                        })
        
        # Plot PRF regions as rectangles
        for prf in prf_data:
            color, _ = self._get_prf_color(prf['type'])
            
            # Add rectangle for each PRF region
            fig.add_shape(
                type="rect",
                x0=prf['start'],                    # Start position
                x1=prf['end'],                      # End position  
                y0=prf['y_center'] - 0.2,           # Bottom of bar
                y1=prf['y_center'] + 0.2,           # Top of bar
                line=dict(color=color, width=2),     # Thicker border
                fillcolor=color,
                opacity=0.9,                        # Higher opacity
                layer='above'
            )
            
            # Add hover functionality with matching color
            hover_text = f"<b>{prf['type']}</b><br>Region: {prf['start']}-{prf['end']}<br>Gene: {prf['gene']}<br>Sequence: {prf['sequence']}"
            colored_hover_text = f'<span style="background-color:{color}; color:white; padding:5px; border-radius:3px;">{hover_text}</span>'
            
            fig.add_trace(go.Scatter(
                x=[(prf['start'] + prf['end']) / 2],  # Center of region
                y=[prf['y_center']],                  # Same y as protein bar
                mode='markers',
                marker=dict(size=1, opacity=0, color=color),       # Invisible marker
                hovertemplate='%{text}<extra></extra>',
                text=[colored_hover_text],
                showlegend=False
            ))
        
        # PRF legend
        prf_types = set(prf['type'] for prf in prf_data)
        self._add_prf_legend(fig, prf_types)
        
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
            title=f"Potential PRF Regions - {self.virus_name} ({genome_id})",
            width=1200,
            height=total_height,
            margin=dict(l=50, r=50, t=100, b=50),
            showlegend=False
        )
        
        # Save
        fig.write_html(output_path)
        fig.write_image(output_path.replace('.html', '.pdf'), width=1200, height=total_height)
        
        print(f"PRF chart saved to: {output_path}")
        
        # Print summary
        print(f"\nPotential PRF Analysis Summary:")
        print(f"  Virus: {self.virus_name}")
        print(f"  Genome ID: {genome_id}")
        print(f"  Proteins: {len(all_genes)}")
        
        # Count PRF regions by type
        prf_counts = {}
        for prf in prf_data:
            prf_type = prf['type']
            prf_counts[prf_type] = prf_counts.get(prf_type, 0) + 1
        
        for prf_type, count in prf_counts.items():
            print(f"  {prf_type}: {count}")
        print(f"  Total Potential PRF Regions: {len(prf_data)}")
        print(f"  Genome length: {self.genome_length:,} bp")
        
        return fig

def parse_arguments():
    """Parse command line arguments."""
    parser = argparse.ArgumentParser(
        description='Generate genome organization visualization with PRF regions',
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
    """Main function to run the PRF analysis."""
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
        print("Example: python figure3_PRF.py --virus SARS-CoV-2")
        return
    
    try:
        visualizer = PRFVisualizer(args.virus, args.genome_id)
        visualizer.create_prf_chart(args.output)
        
        print(f"\n✅ Successfully generated PRF analysis plot:")
        print(f"   Virus: {args.virus}")
        print(f"   Genome ID: {args.genome_id if args.genome_id else 'Auto-detected'}")
        print(f"   Output: {args.output if args.output else 'Auto-generated'}")
        
    except FileNotFoundError as e:
        print(f"❌ Error: {e}")
        print("\nTroubleshooting:")
        print("1. Check that the virus name is correct (use --list-viruses to see available viruses)")
        print("2. Ensure that PRF analysis has been run for this virus")
        print("3. Verify that proteome and reference genome files exist")
    except Exception as e:
        print(f"❌ Error: {e}")

if __name__ == "__main__":
    main()
