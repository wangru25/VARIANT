# Protein Visualization with Plotly

This document summarizes the protein alignment visualization capabilities created for the VARIANT project using Plotly. The visualizations leverage existing mutation data from `.txt` output files and reference proteome FASTA files to create interactive protein-level analysis charts.

## Overview

The visualization system provides multiple views of protein-level mutations:

1. **Protein Alignment Charts** - Horizontal bar charts showing protein lengths with mutation markers
2. **Protein Sequence Heatmaps** - Detailed sequence-level visualization with mutation annotations
3. **Mutation Density Analysis** - Statistical analysis of mutation distribution across proteins
4. **Comprehensive Protein Comparison** - Multi-panel analysis of protein characteristics

## Files Created

### Test Visualizations (`test_protein_visualization.py`)
- `test_protein_alignment.html/pdf` - Basic protein alignment with mutation markers
- `test_mutation_summary.html/pdf` - Stacked bar chart of mutation types by protein

### Enhanced Visualizations (`enhanced_protein_visualization.py`)
- `enhanced_protein_heatmap.html/pdf` - Protein sequence heatmaps with mutation annotations
- `enhanced_mutation_density.html/pdf` - Mutation density analysis with dual-axis charts
- `enhanced_protein_comparison.html/pdf` - Comprehensive 4-panel protein analysis

## Data Sources

### Input Files
- **Proteome FASTA**: `data/SARS-CoV-2/refs/SARS-CoV-2_proteome.fasta`
  - Contains reference protein sequences
  - Includes protein names, IDs, and genome position ranges
  - 31 proteins loaded successfully

- **Mutation Data**: `result/SARS-CoV-2/EPI_ISL_16127660_20250807.txt`
  - Parsed from existing VARIANT output
  - Contains genome-level mutations with protein-level annotations
  - 83 mutations parsed successfully

### Data Structure
```python
# Protein data structure
proteins = {
    'protein_name': {
        'sequence': 'MESLVPGFNEKTHVQLSLPVLQVRDVLVRGFGDSVEEVLSEARQHLKDGTCGLVEVEKGVLPQLEQPYVF...',
        'length': 266,
        'position_range': (266, 805),  # Genome coordinates
        'id': 'YP_009725297.1'
    }
}

# Mutation data structure
mutations = [
    {
        'type': 'missense',
        'position': 670,  # or (start, end) for range mutations
        'protein': 'leader_nsp1',
        'mutation': 'S135R',
        'raw_line': 'missense 670 T->G [{\'protein\': \'leader_nsp1\', \'mutation\': \'S135R\'}]'
    }
]
```

## Visualization Features

### 1. Protein Alignment Chart
- **Purpose**: Overview of all proteins with mutation locations
- **Features**:
  - Horizontal bars representing protein lengths
  - Color-coded by protein function
  - Diamond markers showing mutation positions
  - Interactive hover information with protein details and mutations
  - Mutation type color coding (missense=red, silent=green, deletion=orange)

### 2. Protein Sequence Heatmap
- **Purpose**: Detailed sequence-level analysis
- **Features**:
  - Amino acid sequence visualization using color-coded heatmap
  - Mutation markers overlaid on sequence positions
  - One subplot per protein with mutations
  - Genome-to-protein position mapping
  - Detailed hover information for each mutation

### 3. Mutation Density Analysis
- **Purpose**: Statistical analysis of mutation distribution
- **Features**:
  - Mutation density (mutations per amino acid) by protein
  - Dual-axis chart showing density and absolute counts
  - Sorted by mutation density
  - Color-coded by protein function
  - Comprehensive hover information

### 4. Comprehensive Protein Comparison
- **Purpose**: Multi-dimensional protein analysis
- **Features**:
  - 4-panel layout:
    1. Protein lengths
    2. Mutation counts
    3. Mutation density
    4. Length vs mutations scatter plot
  - Consistent color coding across panels
  - Interactive hover information

## Key Findings from Sample Data

### Mutation Distribution
- **Total mutations**: 83
- **Proteins affected**: 18 out of 31
- **Mutation types**: missense, silent, deletion, unknown

### Top Mutated Proteins
1. **spike_surface_glycoprotein**: 34 mutations (density: 0.0267 mutations/aa)
2. **nsp3**: 9 mutations (density: 0.0046 mutations/aa)
3. **nucleocapsid_phosphoprotein**: 5 mutations (density: 0.0119 mutations/aa)
4. **ORF3a**: 4 mutations (density: 0.0145 mutations/aa)
5. **RNA-dependent-polymerase**: 4 mutations (density: 0.0043 mutations/aa)

### Mutation Density Analysis
- **Highest density**: envelope and spike_surface_glycoprotein (0.0267 mutations/aa)
- **Lowest density**: helicase (0.0017 mutations/aa)
- **Average density**: ~0.01 mutations/aa across affected proteins

## Technical Implementation

### Color Schemes
```python
# Protein colors by function
protein_colors = {
    'leader_nsp1': '#1f77b4',      # Blue
    'spike_surface_glycoprotein': '#c5b0d5',  # Purple
    'nsp3': '#2ca02c',             # Green
    # ... etc
}

# Mutation type colors
mutation_colors = {
    'missense': '#d62728',         # Red
    'silent': '#2ca02c',           # Green
    'deletion': '#ff7f0e',         # Orange
    'unknown': '#7f7f7f'           # Gray
}
```

### Position Mapping
- **Genome positions** from mutation data are mapped to **protein positions**
- Uses protein genome position ranges when available
- Fallback to proportional mapping for complex cases (e.g., join regions)

### Interactive Features
- **Hover information**: Detailed mutation and protein information
- **Zoom and pan**: Full Plotly interactivity
- **Export options**: HTML and PDF formats
- **Responsive design**: Adapts to different screen sizes

## Usage Instructions

### Running the Visualizations
```bash
# Basic visualizations
python test_protein_visualization.py

# Enhanced visualizations
python enhanced_protein_visualization.py
```

### Customizing for Different Viruses
1. Update file paths in the script:
   ```python
   proteome_path = "data/[VIRUS]/refs/[VIRUS]_proteome.fasta"
   mutation_data_path = "result/[VIRUS]/[SAMPLE]_[DATE].txt"
   ```

2. Adjust color schemes for virus-specific proteins
3. Modify position mapping if genome structure differs

### Integration with VARIANT Pipeline
The visualization scripts can be integrated into the VARIANT pipeline by:
1. Adding visualization calls after mutation analysis
2. Creating virus-specific visualization functions
3. Adding CLI options for visualization generation

## Future Enhancements

### Potential Improvements
1. **Multiple sample comparison**: Overlay mutations from multiple samples
2. **Time series analysis**: Track mutations over time
3. **Structural annotations**: Add protein structure information
4. **Conservation analysis**: Compare with reference sequences
5. **Functional impact**: Annotate mutations with functional predictions

### Advanced Features
1. **3D protein visualization**: Integrate with molecular visualization tools
2. **Network analysis**: Show protein interaction networks with mutations
3. **Machine learning integration**: Predict mutation effects
4. **Real-time updates**: Live visualization during analysis

## Conclusion

The protein visualization system provides comprehensive, interactive analysis of protein-level mutations using existing VARIANT output data. The visualizations help researchers understand:

- **Mutation distribution** across viral proteins
- **Sequence-level changes** and their positions
- **Statistical patterns** in mutation accumulation
- **Protein-specific mutation rates** and densities

This system leverages the existing data infrastructure while providing new insights into protein-level mutation analysis, making it a valuable addition to the VARIANT toolkit.
