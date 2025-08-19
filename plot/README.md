# VARIANT Visualization Suite

This folder contains all visualization scripts and outputs for the VARIANT project. All scripts are now **generic and virus-agnostic**, automatically detecting available virus data and creating appropriate visualizations.

## 🚀 Quick Start

### Generate All Figures
```bash
# Run from project root
python plot/create_all_figures.py

# Or navigate to plot folder
cd plot/
python create_all_figures.py
```

### Generate Individual Figures
```bash
# Run from project root
python plot/figure1_workflow_overview.py
python plot/figure2_mutation_distribution.py
python plot/figure3_genome_landscape.py
python plot/figure4_protein_impact.py
python plot/figure5_hot_mutations.py
python plot/figure6_prf_sites.py

# Or navigate to plot folder
cd plot/
python figure1_workflow_overview.py
# ... etc
```

## 📊 Available Visualizations

### **Figure 1: Workflow Overview** - `figure1_workflow_overview.py`
- **Purpose**: Comprehensive workflow diagram showing the complete VARIANT pipeline
- **Content**: Input data → Processing → Analysis → Output → Visualization
- **Output**: `figure1_workflow_overview.html/pdf`

### **Figure 2: Mutation Distribution** - `figure2_mutation_distribution.py`
- **Purpose**: Compare mutation types across different viral families
- **Content**: Stacked bar charts and pie charts showing mutation type distribution
- **Output**: `figure2a_mutation_distribution_stacked.html/pdf`, `figure2b_mutation_distribution_pie.html/pdf`

### **Figure 3: Genome Landscape** - `figure3_genome_landscape.py`
- **Purpose**: Genome-wide mutation frequency and hotspot analysis
- **Content**: Linear genome plots, mutation frequency, protein regions, known variants
- **Output**: `figure3a_{virus}_genome_landscape.html/pdf`, `figure3b_{virus}_mutation_hotspots.html/pdf`

### **Figure 4: Protein Impact** - `figure4_protein_impact.py`
- **Purpose**: Protein-level mutation impact analysis
- **Content**: Heatmaps, mutation summaries, functional impact assessment
- **Output**: `figure4a_{virus}_protein_impact_heatmap.html/pdf`, `figure4b_{virus}_protein_mutation_summary.html/pdf`, etc.

### **Figure 5: Hot Mutations** - `figure5_hot_mutations.py`
- **Purpose**: Analysis of mutation hotspots and biological significance
- **Content**: Scatter plots, protein impact, biological classification, sample distribution
- **Output**: `figure5a_{virus}_hot_mutation_scatter.html/pdf`, `figure5b_{virus}_mutation_type_distribution.html/pdf`, etc.

### **Figure 6: PRF Sites** - `figure6_prf_sites.py`
- **Purpose**: Programmed Ribosomal Frameshifting site analysis
- **Content**: Genome maps, sequence analysis, biological significance
- **Output**: `figure6a_{virus}_prf_genome_map.html/pdf`, `figure6b_{virus}_prf_sequence_analysis.html/pdf`, etc.

## 🔧 Generic Features

### **Automatic Virus Detection**
- Scripts automatically detect available virus data in `../result/`
- Works with any virus: SARS-CoV-2, HIV-1, H3N2, Chikungunya, ZaireEbola, etc.
- Handles both segmented (H3N2) and non-segmented viruses

### **Flexible Configuration**
- Virus-specific genome regions and known variants can be added
- Default configurations work for any virus
- Easy to extend for new viruses

### **Robust Error Handling**
- Graceful handling of missing data
- Informative error messages
- Continues processing even if some data is unavailable

## 📁 Output Files

### **File Naming Convention**
- Generic figures: `figure1_workflow_overview.html/pdf`
- Virus-specific figures: `figure3a_SARS-CoV-2_genome_landscape.html/pdf`

### **File Types**
- **HTML**: Interactive visualizations (open in web browser)
- **PDF**: Static figures ready for publication

### **Generated Files Example**
```
figure1_workflow_overview.html
figure1_workflow_overview.pdf
figure2a_mutation_distribution_stacked.html
figure2a_mutation_distribution_stacked.pdf
figure2b_mutation_distribution_pie.html
figure2b_mutation_distribution_pie.pdf
figure3a_SARS-CoV-2_genome_landscape.html
figure3a_SARS-CoV-2_genome_landscape.pdf
figure3b_SARS-CoV-2_mutation_hotspots.html
figure3b_SARS-CoV-2_mutation_hotspots.pdf
...
```

## 📋 Requirements

### **Python Packages**
```bash
pip install plotly pandas numpy
```

### **Data Structure**
```
../result/
├── SARS-CoV-2/
│   ├── *_20250807.txt          # Mutation files
│   ├── *_row_hot_mutations_20250807.csv  # Hot mutation files
│   └── *potential_PRF_20250807.csv       # PRF files
├── HIV-1/
│   └── ...
├── H3N2/
│   ├── segment_1/
│   ├── segment_2/
│   └── ...
└── ...
```

## 🎯 Usage Examples

### **Generate All Figures for Available Viruses**
```bash
python plot/create_all_figures.py
```

### **Generate Specific Figure for All Viruses**
```bash
python plot/figure2_mutation_distribution.py
```

### **Customize for Specific Virus**
Edit the scripts to analyze specific viruses by modifying the virus selection logic in the `main()` functions.

## 📚 Documentation

- `COMPREHENSIVE_VISUALIZATION_SUMMARY.md` - Complete documentation
- `PROTEIN_VISUALIZATION_SUMMARY.md` - Protein visualization overview

## 🔄 Legacy Scripts

The following legacy scripts are still available but are virus-specific:
- `genome_organization_visualization.py` - Original genome visualization
- `genome_organization_analysis.html/pdf` - Legacy output files

## 🆕 New Generic Scripts

All new scripts (figure1-6) are:
- ✅ **Virus-agnostic**
- ✅ **Automatically detect data**
- ✅ **Handle multiple viruses**
- ✅ **Robust error handling**
- ✅ **Publication-ready output**

## 🎨 Customization

### **Adding New Viruses**
1. Add virus-specific genome regions to `get_genome_regions()` function
2. Add known variants to `get_known_variants()` function
3. Place virus data in `../result/{virus_name}/`

### **Modifying Visualizations**
- Edit individual figure scripts for custom plots
- Modify color schemes, layouts, and annotations
- Add new plot types and analyses

### **Batch Processing**
- Use `create_all_figures.py` for comprehensive analysis
- Modify to process specific viruses or figure types
- Integrate with automated workflows
