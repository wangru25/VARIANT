# VARIANT Visualization Suite

This folder contains all visualization scripts and outputs for the VARIANT project.

## Quick Start

```bash
# Run from project root
python plot/genome_organization_visualization.py

# Or navigate to plot folder
cd plot/
python genome_organization_visualization.py
```

## Available Visualizations

1. **Basic Protein Alignment** - `test_protein_visualization.py`
2. **Enhanced Protein Analysis** - `enhanced_protein_visualization.py`
3. **Genome-Based Protein Visualization** - `enhanced_genome_protein_visualization.py`
4. **MSA Conservation Analysis** - `msa_conservation_visualization.py`
5. **Genome Organization** - `genome_organization_visualization.py` ⭐

## Documentation

- `COMPREHENSIVE_VISUALIZATION_SUMMARY.md` - Complete documentation
- `PROTEIN_VISUALIZATION_SUMMARY.md` - Protein visualization overview

## Output Files

All visualizations generate both HTML (interactive) and PDF (static) files.

## Data Sources

- Proteome: `../data/SARS-CoV-2/refs/SARS-CoV-2_proteome.fasta`
- Mutations: `../result/SARS-CoV-2/EPI_ISL_16127660_20250807.txt`
- Reference: `../data/SARS-CoV-2/refs/NC_045512.fasta`
