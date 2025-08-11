# Comprehensive Visualization Suite for VARIANT

This document provides a complete overview of all visualization capabilities created for the VARIANT project, leveraging existing mutation data and reference proteome files to create sophisticated bioinformatics visualizations.

**📁 All visualization files are located in the `plot/` folder to maintain a clean root directory structure.**

## Overview

The VARIANT visualization suite provides multiple complementary views of viral mutation analysis:

1. **Basic Protein Alignment** - Horizontal bar charts with mutation markers
2. **Enhanced Protein Analysis** - Sequence heatmaps and statistical analysis
3. **Genome-Based Protein Visualization** - Genome-positioned proteins with mutation highlighting
4. **MSA Conservation Analysis** - Multiple sequence alignment with conservation scores
5. **Genome Organization** - Reading frame-based genome organization (HIV-style)

## File Organization

```
VARIANT/
├── plot/                          # 📁 All visualization files
│   ├── test_protein_visualization.py
│   ├── enhanced_protein_visualization.py
│   ├── enhanced_genome_protein_visualization.py
│   ├── msa_conservation_visualization.py
│   ├── genome_organization_visualization.py
│   ├── PROTEIN_VISUALIZATION_SUMMARY.md
│   ├── COMPREHENSIVE_VISUALIZATION_SUMMARY.md
│   └── *.html, *.pdf              # Output files
├── data/                          # Input data
├── result/                        # Mutation analysis results
└── ...                           # Other project files
```

## Visualization Types

### 1. Basic Protein Alignment (`plot/test_protein_visualization.py`)

**Purpose**: Overview of protein mutations with basic alignment visualization

**Features**:
- Horizontal protein bars with mutation markers
- Color-coded by protein function
- Interactive hover information
- Mutation type color coding

**Files Created**:
- `plot/test_protein_alignment.html/pdf` - Basic protein alignment
- `plot/test_mutation_summary.html/pdf` - Stacked bar chart of mutation types

**Key Statistics**:
- 31 proteins loaded
- 83 mutations parsed
- 18 proteins with mutations
- Top mutated: spike_surface_glycoprotein (34 mutations)

### 2. Enhanced Protein Analysis (`plot/enhanced_protein_visualization.py`)

**Purpose**: Detailed protein-level analysis with sequence visualization

**Features**:
- Protein sequence heatmaps with mutation annotations
- Mutation density analysis with dual-axis charts
- Comprehensive 4-panel protein comparison
- Statistical analysis of mutation distribution

**Files Created**:
- `plot/enhanced_protein_heatmap.html/pdf` - Sequence-level mutation visualization
- `plot/enhanced_mutation_density.html/pdf` - Mutation density analysis
- `plot/enhanced_protein_comparison.html/pdf` - Multi-panel protein analysis

**Key Insights**:
- Mutation density ranges from 0.0017 to 0.0267 mutations/aa
- Highest density: envelope and spike_surface_glycoprotein
- Statistical correlation between protein length and mutation count

### 3. Genome-Based Protein Visualization (`plot/enhanced_genome_protein_visualization.py`)

**Purpose**: Genome-positioned protein analysis with reference overlay

**Features**:
- X-axis based on genome sequence length
- Protein regions mapped to genome positions
- Reference wild-type proteins overlaid on top
- Mutated sequences with darker colors for mutations
- Detailed mutation analysis showing sequence changes

**Files Created**:
- `plot/enhanced_genome_protein_alignment.html/pdf` - Genome-based protein alignment
- `plot/enhanced_detailed_mutation_analysis.html/pdf` - Sequence-level mutation comparison

**Key Features**:
- Genome length: 29,903 base pairs
- Proteins positioned by genome coordinates
- Reference vs. mutated sequence comparison
- Mutation position mapping from genome to protein

### 4. MSA Conservation Analysis (`plot/msa_conservation_visualization.py`)

**Purpose**: Multiple sequence alignment with conservation analysis

**Features**:
- Conservation bar chart with gap frequency
- Multiple sequence alignment visualization
- Amino acid color coding by conservation
- Mutation annotations overlaid on sequences
- Conservation-mutation correlation analysis

**Files Created**:
- `plot/msa_conservation_analysis.html/pdf` - MSA with conservation scores
- `plot/conservation_mutation_correlation.html/pdf` - Correlation analysis

**Key Features**:
- Synthetic MSA data for demonstration
- Conservation scores calculated per position
- Gap frequency analysis
- Mutation density correlation with conservation

### 5. Genome Organization (`plot/genome_organization_visualization.py`)

**Purpose**: HIV-style genome organization across reading frames

**Features**:
- Genes organized by reading frames
- Regulatory region annotations
- Mutation density histogram
- Gene overlap analysis
- Reading frame-specific visualization

**Files Created**:
- `plot/genome_organization_analysis.html/pdf` - Reading frame organization
- `plot/gene_overlap_analysis.html/pdf` - Gene overlap statistics

**Key Statistics**:
- **Reading Frame 1 (Replicase)**: 12 genes
- **Reading Frame 2 (Structural)**: 3 genes  
- **Reading Frame 3 (Accessory)**: 15 genes
- **Regulatory Regions**: 5' UTR, 3' UTR

## Data Sources

### Input Files
- **Proteome FASTA**: `data/SARS-CoV-2/refs/SARS-CoV-2_proteome.fasta`
  - 31 proteins with sequences and position ranges
  - Includes protein IDs and genome coordinates

- **Mutation Data**: `result/SARS-CoV-2/EPI_ISL_16127660_20250807.txt`
  - 83 parsed mutations
  - Genome-level positions with protein-level annotations
  - Mutation types: missense, silent, deletion, unknown

- **Reference Genome**: `data/SARS-CoV-2/refs/NC_045512.fasta`
  - 29,903 base pairs
  - Used for genome length and positioning

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

### Interactive Features
- **Hover information**: Detailed mutation and protein information
- **Zoom and pan**: Full Plotly interactivity
- **Export options**: HTML (interactive) and PDF (static) formats
- **Responsive design**: Adapts to different screen sizes

### Position Mapping
- **Genome positions** mapped to **protein positions**
- Uses protein genome position ranges when available
- Fallback to proportional mapping for complex cases
- Handles range mutations and single position mutations

## Key Findings

### Mutation Distribution
- **Total mutations**: 83
- **Proteins affected**: 18 out of 31
- **Mutation types**: missense (majority), silent, deletion, unknown

### Top Mutated Proteins
1. **spike_surface_glycoprotein**: 34 mutations (density: 0.0267 mutations/aa)
2. **nsp3**: 9 mutations (density: 0.0046 mutations/aa)
3. **nucleocapsid_phosphoprotein**: 5 mutations (density: 0.0119 mutations/aa)
4. **ORF3a**: 4 mutations (density: 0.0145 mutations/aa)
5. **RNA-dependent-polymerase**: 4 mutations (density: 0.0043 mutations/aa)

### Genome Organization
- **Reading Frame 1**: Replicase proteins (nsp1-nsp11, helicase)
- **Reading Frame 2**: Structural proteins (spike, envelope, membrane)
- **Reading Frame 3**: Accessory proteins (ORF3a-ORF10, nucleocapsid)

### Mutation Patterns
- **Highest density**: envelope and spike_surface_glycoprotein (0.0267 mutations/aa)
- **Lowest density**: helicase (0.0017 mutations/aa)
- **Average density**: ~0.01 mutations/aa across affected proteins

## Usage Instructions

### Running All Visualizations
```bash
# Navigate to plot directory
cd plot/

# Basic visualizations
python test_protein_visualization.py

# Enhanced protein analysis
python enhanced_protein_visualization.py

# Genome-based protein visualization
python enhanced_genome_protein_visualization.py

# MSA conservation analysis
python msa_conservation_visualization.py

# Genome organization
python genome_organization_visualization.py
```

### Running from Root Directory
```bash
# Run from project root (recommended)
python plot/test_protein_visualization.py
python plot/enhanced_protein_visualization.py
python plot/enhanced_genome_protein_visualization.py
python plot/msa_conservation_visualization.py
python plot/genome_organization_visualization.py
```

### Customizing for Different Viruses
1. Update file paths in each script:
   ```python
   proteome_path = "../data/[VIRUS]/refs/[VIRUS]_proteome.fasta"
   mutation_data_path = "../result/[VIRUS]/[SAMPLE]_[DATE].txt"
   reference_genome_path = "../data/[VIRUS]/refs/[REFERENCE].fasta"
   ```

2. Adjust color schemes for virus-specific proteins
3. Modify position mapping if genome structure differs
4. Update reading frame assignments based on virus organization

### Integration with VARIANT Pipeline
The visualization scripts can be integrated into the VARIANT pipeline by:
1. Adding visualization calls after mutation analysis
2. Creating virus-specific visualization functions
3. Adding CLI options for visualization generation
4. Automating visualization creation for batch analysis

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
5. **Comparative genomics**: Compare mutations across viral strains

## File Summary

### Scripts Created (in `plot/` folder)
- `test_protein_visualization.py` - Basic protein alignment (16KB)
- `enhanced_protein_visualization.py` - Enhanced protein analysis (21KB)
- `enhanced_genome_protein_visualization.py` - Genome-based visualization (23KB)
- `msa_conservation_visualization.py` - MSA conservation analysis (23KB)
- `genome_organization_visualization.py` - Genome organization (20KB)

### Output Files (in `plot/` folder)
- **HTML files**: Interactive visualizations (4.4-4.5MB each)
- **PDF files**: Static visualizations for publication (34-192KB each)
- **Documentation**: Comprehensive summaries and usage guides

### Documentation Files
- `plot/PROTEIN_VISUALIZATION_SUMMARY.md` - Protein visualization overview
- `plot/COMPREHENSIVE_VISUALIZATION_SUMMARY.md` - Complete documentation (this file)

## Conclusion

The VARIANT visualization suite provides comprehensive, interactive analysis of viral mutations across multiple dimensions:

- **Protein-level analysis** with sequence visualization
- **Genome-level organization** with reading frame analysis
- **Conservation analysis** with MSA visualization
- **Statistical analysis** of mutation patterns and distributions

This system leverages existing data infrastructure while providing new insights into viral mutation analysis, making it a valuable addition to the VARIANT toolkit for bioinformatics research and viral evolution studies.

The visualizations are designed to be:
- **Interactive**: Full Plotly functionality with hover information
- **Comprehensive**: Multiple complementary views of the same data
- **Extensible**: Easy to adapt for different viruses and analysis types
- **Publication-ready**: High-quality PDF exports for manuscripts
- **User-friendly**: Clear documentation and usage instructions
- **Well-organized**: All files contained in the `plot/` folder for clean project structure
