# VARIANT Visualization Module

This module contains visualization tools for generating publication-ready figures from VARIANT analysis results.

## Files

- `figure1_genome_analysis.py` - Genome organization visualization with mutation overlays
- `__init__.py` - Module initialization and exports

## Usage

### Command Line Interface

```bash
# From project root
python generate_genome_plot.py --virus SARS-CoV-2 --genome-id EPI_ISL_16127650
python generate_genome_plot.py --virus HIV-1
python generate_genome_plot.py --list-viruses
python generate_genome_plot.py --process-all
```

### Python API

```python
from src.visualization import GenomeOrganizationVisualizer

# Create visualization for specific virus and genome
visualizer = GenomeOrganizationVisualizer('SARS-CoV-2', 'EPI_ISL_16127650')
visualizer.create_genome_organization_chart('my_plot.html')
```

## Features

- **Auto-detection**: Automatically finds proteome and mutation CSV files
- **Multi-virus support**: Works with SARS-CoV-2, HIV-1, Chikungunya, ZaireEbola
- **Publication-ready**: High-resolution HTML and PDF output
- **Interactive elements**: Hover details and zoom capabilities
- **Consistent styling**: Professional appearance across all viruses

## Output

Generates both interactive HTML files and static PDF files suitable for:
- Scientific publications
- Web applications
- Presentations
- Further analysis
