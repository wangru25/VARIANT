# VARIANT: Viral mutAtion trackeR aImed At GeNome and proTein-level

A comprehensive Python framework for analyzing viral mutations, supporting both single-segment and multi-segment viruses. The framework provides detailed analysis of nucleotide changes, protein impacts, and mutation classifications including missense, nonsense, and frameshift mutations.

## Features

- **🦠 Multi-virus support**: SARS-CoV-2, ZaireEbola, Chikungunya, HIV-1, H3N2, and easily extensible
- **🧬 Multi-segment virus support** with automatic structure detection (e.g., H3N2 influenza)
- **🔬 Comprehensive mutation analysis**: Point mutations, insertions, deletions, row and hot mutations
- **⚗️ Advanced mutation classification**: Missense, nonsense, silent, and frameshift mutations
- **🧪 Protein-level impact analysis**: Amino acid changes with biological significance
- **🤖 Automatic mutation summary generation**: No manual flags required
- **📊 Complex protein coordinate handling**: Supports `join()` coordinates for viral polyproteins
- **🔄 Programmed Ribosomal Frameshifting (PRF) detection** (+1/-1 frameshifting)
- **📈 Structured output**: Text and CSV formats for analysis and visualization
- **⚙️ Configurable per-virus settings**: Easy setup for new viruses
- **🌐 Web Application**: Modern NYU-themed web interface with bulk ZIP downloads
- **🗜️ Bulk File Management**: Organized ZIP downloads for complete datasets

## Web Application

VARIANT now includes a modern web-based interface with professional design and streamlined file management:

### **🎨 Professional Interface**
- **Clean Design**: Modern blue color scheme
- **Responsive Layout**: Works on desktop, tablet, and mobile devices
- **Accessible Design**: High contrast and readable typography

### **📁 Streamlined File Management**
- **Bulk ZIP Downloads**: Download complete datasets as organized ZIP files
- **Clean File Listings**: Simplified file display without redundant information
- **Upload Progress Tracking**: Visual indicators for file upload status
- **Custom Virus Support**: Add new virus types through the web interface

### **🔬 Analysis Features**
- **Real-time Monitoring**: Live status updates during analysis
- **Background Processing**: Non-blocking analysis with job queue
- **Job History**: Track and manage previous analysis jobs
- **Multi-virus Support**: Analyze all supported viruses through the web interface

### **Quick Start**
```bash
# Deploy with Railway (Recommended)
git clone https://github.com/wangru25/VARIANT.git
cd VARIANT

# Deploy to Railway
# 1. Go to railway.app
# 2. Connect your GitHub repository
# 3. Railway will auto-deploy your app

# Access the application
# Web Interface: https://your-app-name.up.railway.app
# API Documentation: https://your-app-name.up.railway.app/docs
```

📖 **[Complete Web Application Guide](docs/README_WEB.md)** | 🚀 **[Railway Deployment Guide](docs/RAILWAY_DEPLOYMENT.md)**

## Supported Viruses

VARIANT currently supports comprehensive analysis for the following viruses:

| **Virus** | **Type** | **Key Proteins** | **Special Features** |
|-----------|----------|------------------|----------------------|
| **SARS-CoV-2** | Single-segment | RNA-dependent-polymerase, spike, nucleocapsid | Complex `join()` coordinates, 130+ genomes |
| **ZaireEbola** | Single-segment | Nucleoprotein, polymerase, spike glycoprotein | Large genome analysis |
| **Chikungunya** | Single-segment | nsp1-4, capsid, E1/E2/E3, 6K proteins | Multi-protein analysis |
| **HIV-1** | Single-segment | Pr55(Gag), reverse transcriptase, integrase, envelope | Polyprotein processing |
| **H3N2** | Multi-segment (8) | PB1/PB2/PA, hemagglutinin, neuraminidase, matrix | Influenza segment analysis |

### Analysis Capabilities per Virus
- **Mutation detection**: Point mutations, insertions, deletions
- **Protein impact**: Amino acid changes, functional predictions
- **Evolutionary analysis**: Row and hot mutation patterns
- **Summary generation**: Automated CSV reports
- **Frameshift detection**: PRF analysis for all supported viruses

## Installation

VARIANT can be installed and used in two ways:

| Feature | CLI Installation | Conda Environment |
|---------|------------------|-------------------|
| **Installation** | `pip install -e .` | `conda env create -f environment.yml` |
| **Command Style** | `variant prf --genome file.fasta` | `python main.py --virus SARS-CoV-2 --detect-frameshifts` |
| **Dependencies** | Auto-installed via pip | Managed by conda |
| **Best For** | Modern CLI users, CI/CD | Bioinformatics workflows, isolated environments |
| **Entry Point** | `variant` command | `python main.py` |

### Option 1: CLI Installation (Recommended)

**Prerequisites**: Python 3.8 or higher

```bash
# Clone the repository
git clone https://github.com/wangru25/VARIANT.git
cd VARIANT

# Install the package and dependencies
pip install -e .

# Verify installation
variant --help
```

### Option 2: Conda Environment

**Prerequisites**: Conda or Miniconda installed

```bash
# Clone the repository
git clone https://github.com/wangru25/VARIANT.git
cd VARIANT

# Create and activate conda environment
conda env create -f environment.yaml
conda activate variant

# Make scripts executable (if needed)
chmod +x main.py

# Verify installation
python main.py --help
```

### Dependencies
- BioPython >= 1.79
- NumPy >= 1.21.0  
- Pandas >= 1.3.0
- PyYAML >= 6.0
- fuzzysearch == 0.7.3
- more-itertools >= 8.12.0
- ViennaRNA >= 2.4.0 (for PRF scanning and RNA structure prediction)

## Usage

Choose your preferred method based on your installation:

### Option 1: CLI Usage (After `pip install -e .`)

The CLI provides modern, user-friendly commands:

```bash
# Show all available commands
variant --help

# PRF Analysis
variant prf --genome data/SARS-CoV-2/refs/NC_045512.fasta --output results.csv

# Mutation Analysis
variant analyze --virus SARS-CoV-2 --msa data/SARS-CoV-2/clustalW/test_msa_2.txt

# Dataset Setup
variant setup --virus SARS-CoV-2 --config virus_config.yaml
variant setup --virus HIV-1
```

### Option 2: Conda Usage (After `conda activate variant`)

Use the traditional Python scripts:

```bash
# Show available options
python main.py --help

# Analyze specific genome
python main.py --virus SARS-CoV-2 --genome-id EPI_ISL_16327572

# Process all genomes in MSA (automatically generates mutation summaries)
python main.py --virus SARS-CoV-2 --process-all

# PRF analysis
python main.py --virus SARS-CoV-2 --detect-frameshifts

# Multi-segment virus analysis (processes all segments)
python main.py --virus H3N2 --process-all

# Specific segment analysis
python main.py --virus H3N2 --segment segment_1 --process-all

# List available viruses
python main.py --list-viruses
```

### Alternative Methods (Both Options)

```bash
# Python module approach
python -m src.cli.commands --help

# Direct script execution (conda users)
python main.py --virus SARS-CoV-2 --genome-id EPI_ISL_16327572
```

## Quick Start Examples

```bash
# Process all SARS-CoV-2 genomes (most common use case)
python main.py --virus SARS-CoV-2 --process-all

# Process all H3N2 influenza segments
python main.py --virus H3N2 --process-all

# Analyze a specific genome
python main.py --virus SARS-CoV-2 --genome-id EPI_ISL_16327572

# Detect frameshift sites
python main.py --virus SARS-CoV-2 --detect-frameshifts

# Advanced PRF scanning with RNA structure prediction
python src/core/prf_scanner.py --fasta data/SARS-CoV-2/refs/NC_045512.fasta --out results/sars_cov2_prf --use-rnafold

# See what viruses are available
python main.py --list-viruses
```

### Command Line Arguments

- `--virus`: Name of the virus to analyze
- `--genome-id`: Specific genome ID to process
- `--process-all`: Process all genomes in the MSA file (automatically generates mutation summaries)
- `--segment`: For multi-segment viruses, specify segment to analyze
- `--detect-frameshifts`: Detect potential frameshifting sites (PRF analysis only)
- `--msa-file`: Custom MSA file path
- `--config`: Custom configuration file path

## Data Organization

### Directory Structure
```
data/
├── SingleSegmentVirus/
│   ├── clustalW/     # MSA files
│   ├── refs/         # Reference genomes and proteomes
│   └── fasta/        # FASTA sequence files
└── MultiSegmentVirus/
    ├── segment_1/
    │   ├── clustalW/
    │   ├── refs/
    │   └── fasta/
    ├── segment_2/
    └── ...

result/
├── SingleSegmentVirus/
│   └── analysis_results
└── MultiSegmentVirus/
    ├── segment_1/
    ├── segment_2/
    └── ...
```

### Configuration (virus_config.yaml)

```yaml
viruses:
  SingleSegmentVirus:
    reference_genome: "reference.fasta"
    proteome_file: "proteome.fasta"
    codon_table_id: 1
    description: "Single segment virus"
    default_msa_file: "msa.txt"
    
  MultiSegmentVirus:
    codon_table_id: 1
    description: "Multi-segment virus"
    segments:
      segment_1:
        reference_genome: "segment1_ref.fasta"
        proteome_file: "segment1_proteome.fasta"
        default_msa_file: "segment1_msa.txt"
      segment_2:
        reference_genome: "segment2_ref.fasta"
        proteome_file: "segment2_proteome.fasta"
        default_msa_file: "segment2_msa.txt"
```

## Output Files

### Text Output (genome_id_date.txt)
Contains detailed mutation analysis including:
- Point mutations
- Insertions/Deletions
- Row mutations (consecutive changes)
- Hot mutations (non-consecutive changes)
- Protein impact classification:
  - Silent mutations
  - Missense mutations
  - Nonsense mutations (stop codons)
  - Frameshift mutations

### Process All Output
When using `--process-all`:
- Analyzes all genomes in the MSA file
- Automatically generates mutation summary CSV files for all processed genomes
- Prints processing progress and summary statistics
- Creates comprehensive mutation analysis for each genome

### PRF Analysis Output
When using `--detect-frameshifts`:

**Potential PRF Analysis (potential_PRF_date.csv):**
- CSV format with columns: `position,end_position,sequence,type`
- -1 PRF sites (slippery sequences)
- +1 PRF sites (shifty stop codons)
- Stem-loop structures (stimulatory elements)

## PRF Scanner: Advanced Programmed Ribosomal Frameshifting Detection

VARIANT includes a comprehensive PRF scanner (`src/core/prf_scanner.py`) that detects candidate programmed ribosomal frameshifting sites with advanced RNA secondary structure prediction and tRNA interaction validation.

### What is Programmed Ribosomal Frameshifting (PRF)?

Programmed ribosomal frameshifting is a biological mechanism where ribosomes "slip" or "shift" their reading frame during translation, allowing a single mRNA to code for multiple proteins. This is particularly common in viruses as a way to maximize their coding capacity and control protein stoichiometry.

### Types of Frameshifting Detected

#### 1. -1 PRF (Backward Frameshifting)
- **Mechanism**: Ribosome shifts backward by one nucleotide
- **Trigger**: Slippery sequences following the pattern `XXXYYYZ`
- **Examples**: 
  - `TTTTTTA` (HIV gag-pol frameshift)
  - `TTTAAAC` (SARS-CoV-2 frameshift)
  - `GGGTTTA` (Coronavirus frameshift)

#### 2. +1 PRF (Forward Frameshifting)
- **Mechanism**: Ribosome shifts forward by one nucleotide
- **Trigger**: Shifty stop codons with specific context
- **Examples**: 
  - `TAA` in `TCCT` context
  - `TAG` in `TCCG` context

### Biological Requirements for PRF

Our detection algorithm validates frameshift sites based on established biological criteria:

1. **Slippery Sequence Pattern**: Must follow `XXXYYYZ` pattern
2. **tRNA Pairing Compatibility**: P-site and A-site tRNAs must maintain pairing after shift
3. **Optimal Spacing**: Downstream stem-loop should be 5-9 nucleotides away
4. **Structural Context**: Stem-loop structures that pause ribosome translocation
5. **Reading Frame Context**: Shift must place ribosome in viable alternative frame

### tRNA Pairing Validation

The algorithm validates that tRNA molecules can maintain proper pairing before and after frameshifting:

- **P-site tRNA**: Must maintain pairing with first two positions
- **A-site tRNA**: Must maintain pairing with first two positions and allow wobble pairing at third position
- **Wobble Compatibility**: Supports non-canonical base pairs observed in viral PRF (G-U, A-C, etc.)

### Example Output Interpretation

```
position,end_position,sequence,type
1630,1637,TTTTTTA,-1 PRF
13461,13468,TTTAAAC,-1 PRF
```

**Breakdown:**
- **Position**: 1630-1637 (nucleotides in genome)
- **Sequence**: `TTTTTTA` (slippery sequence)
- **Type**: `-1 PRF` (backward frameshifting)

### Known Viral Frameshift Sites

Our algorithm successfully detects well-characterized viral frameshift sites:

- **HIV-1**: `TTTTTTA` at position ~1637 (gag-pol frameshift)
- **SARS-CoV-2**: `TTTAAAC` at position ~13462 (ORF1a-ORF1ab frameshift)
- **Coronaviruses**: Various slippery sequences in replicase genes
- **Retroviruses**: Multiple frameshift sites for polyprotein production

### PRF Scanner Usage

The PRF scanner can be used as a standalone tool or integrated with the main VARIANT pipeline:

#### Standalone PRF Scanning

```bash
# Basic PRF scanning with RNA structure prediction
python src/core/prf_scanner.py --fasta data/SARS-CoV-2/refs/NC_045512.fasta --out results/sars_cov2_prf --use-rnafold

# Advanced scanning with custom parameters
python src/core/prf_scanner.py \
    --fasta data/SARS-CoV-2/refs/NC_045512.fasta \
    --out results/sars_cov2_prf \
    --spacer-min 5 \
    --spacer-max 9 \
    --window 120 \
    --use-rnafold \
    --organism sars_cov2

# With GFF annotation for frame context
python src/core/prf_scanner.py \
    --fasta data/SARS-CoV-2/refs/NC_045512.fasta \
    --out results/sars_cov2_prf \
    --use-rnafold \
    --gff data/SARS-CoV-2/annotations.gff

# With tRNA abundance data
python src/core/prf_scanner.py \
    --fasta data/SARS-CoV-2/refs/NC_045512.fasta \
    --out results/sars_cov2_prf \
    --use-rnafold \
    --trna data/tRNA_abundance.csv \
    --organism sars_cov2
```

#### Integration with Main Pipeline

```bash
# PRF detection through main pipeline (now uses advanced PRF scanner)
python main.py --virus SARS-CoV-2 --detect-frameshifts

# Process all genomes with PRF detection
python main.py --virus SARS-CoV-2 --process-all --detect-frameshifts

# PRF detection for other viruses
python main.py --virus HIV-1 --detect-frameshifts
python main.py --virus Chikungunya --detect-frameshifts
python main.py --virus ZaireEbola --detect-frameshifts
```

**Note**: The `--detect-frameshifts` flag now uses the advanced PRF scanner with RNA structure prediction, tRNA validation, and comprehensive analysis. This replaces the previous simple frameshift detection.

### PRF Scanner Parameters

| Parameter | Description | Default | Example |
|-----------|-------------|---------|---------|
| `--fasta` | Input genome FASTA file | Required | `data/SARS-CoV-2/refs/NC_045512.fasta` |
| `--out` | Output file prefix | Required | `results/prf_analysis` |
| `--spacer-min` | Minimum spacer length (nt) | 5 | `--spacer-min 3` |
| `--spacer-max` | Maximum spacer length (nt) | 9 | `--spacer-max 12` |
| `--window` | Downstream fold window size (nt) | 120 | `--window 150` |
| `--use-rnafold` | Enable RNA structure prediction | False | `--use-rnafold` |
| `--gff` | GFF annotation file for frame context | None | `--gff annotations.gff` |
| `--organism` | Organism for tRNA data | human | `--organism sars_cov2` |
| `--trna` | tRNA abundance CSV file | None | `--trna trna_data.csv` |

### PRF Scanner Output

The PRF scanner generates two output files:

#### 1. CSV Output (`*.prf_candidates.csv`)

Comprehensive candidate data with the following columns:

| Column | Description | Example |
|--------|-------------|---------|
| `seqid` | Sequence identifier | `NC_045512` |
| `site_start_1based` | PRF site start position (1-based) | `13462` |
| `slippery_motif` | Slippery sequence motif | `TTTAAAC` |
| `type` | PRF type | `-1` |
| `spacer_nt` | Spacer length (nt) | `7` |
| `fold_window_start` | Structure window start | `13470` |
| `fold_window_end` | Structure window end | `13590` |
| `fold_window_seq` | Downstream sequence | `GUGGCUGUCACUCGGC...` |
| `rnafold_structure` | RNA secondary structure | `((((.((.((((.(((.....))).)))))).))))` |
| `rnafold_mfe` | Minimum free energy (kcal/mol) | `-30.4` |
| `rnafold_type` | Structure classification | `stem_loop` |
| `probknot_structure` | ProbKnot structure prediction | `(((((.........[[[[[.........))))).]]]]]` |
| `probknot_type` | ProbKnot classification | `pseudoknot` |
| `pknots_structure` | PKNOTS structure prediction | `(((((.........[[[[[.........))))).]]]]]` |
| `pknots_type` | PKNOTS classification | `pseudoknot` |
| `frame_context` | Reading frame context | `CDS_frame0` |
| `codon1` | P-site codon | `TTT` |
| `codon2` | A-site codon | `AAA` |
| `trna1_abundance` | P-site tRNA abundance | `0.8` |
| `trna2_abundance` | A-site tRNA abundance | `0.9` |
| `pausing_potential` | Ribosomal pausing score | `0.3` |
| `wobble_pairs` | Wobble base pairs detected | `P-site_G-U` |
| `trna_score` | Overall tRNA interaction score | `0.85` |

#### 2. BED Output (`*.prf_candidates.bed`)

Genomic coordinates for visualization and analysis:

```
NC_045512	13462	13469	-1_TTTAAAC_s7	0	+
NC_045512	76	83	-1_UUUAAAA_s5	0	+
```

### Advanced Features

#### RNA Secondary Structure Prediction

The PRF scanner uses multiple tools for comprehensive structure analysis:

- **RNAfold**: Minimum free energy (MFE) structure prediction
- **ProbKnot**: Pseudoknot detection and prediction
- **PKNOTS**: Alternative pseudoknot prediction
- **Simple Detection**: Fallback methods when external tools unavailable

#### tRNA Interaction Validation

Validates biological feasibility of frameshifting:

- **tRNA Abundance**: Uses organism-specific tRNA abundance data
- **Pausing Potential**: Calculates ribosomal pausing likelihood
- **Wobble Pairs**: Detects non-canonical base pairs
- **Interaction Scores**: Overall tRNA interaction assessment

#### Frame Context Analysis

When GFF annotations are provided:

- **CDS Overlap**: Identifies coding sequence overlaps
- **Frame Classification**: Determines reading frame context
- **Multi-frame Detection**: Handles overlapping ORFs

### PRF Scanner Scripts and Examples

VARIANT includes several helper scripts and examples for PRF scanning:

#### Quick PRF Scanning Script

```bash
# Scan all supported viruses for PRF sites
./scripts/run_prf_scan.sh

# This script automatically:
# - Activates the conda environment
# - Scans all virus reference genomes
# - Generates comprehensive PRF candidate files
# - Provides summary statistics
```

#### Python Examples Script

```bash
# Run comprehensive PRF scanning examples
python examples/prf_scanner_examples.py

# This script demonstrates:
# - Basic PRF scanning for each virus
# - Advanced parameter configurations
# - Different organism settings
# - Performance comparisons
```

#### Example Output Files

After running PRF scanning, you'll get:

```
result/prf_scan_results/
├── sars_cov2_prf.prf_candidates.csv    # SARS-CoV-2 PRF candidates
├── sars_cov2_prf.prf_candidates.bed    # Genomic coordinates
├── hiv1_prf.prf_candidates.csv         # HIV-1 PRF candidates
├── hiv1_prf.prf_candidates.bed         # Genomic coordinates
├── chikungunya_prf.prf_candidates.csv  # Chikungunya PRF candidates
├── chikungunya_prf.prf_candidates.bed  # Genomic coordinates
└── ... (additional virus results)
```

### Technical Implementation

The frameshift detection uses:
- **Pattern Matching**: Identifies slippery sequences and shifty stops
- **tRNA Validation**: Ensures biological feasibility of frameshifting
- **Structural Analysis**: Detects stem-loop structures that promote frameshifting
- **Frame Analysis**: Calculates reading frame shifts
- **CSV Output**: Provides structured data for further analysis

### CSV Output (genome_id_row_hot_mutations_date.csv)
Structured data for statistical analysis:
- Genome ID
- Mutation type
- Position(s)
- Nucleotide changes
- Affected protein
- Amino acid changes
- Biological classification

### Mutation Summary Output (genome_id_mutation_summary_date.csv)
Automatically generated for all processed genomes:
- Comprehensive mutation summary in CSV format
- Includes all mutation types with protein impact
- Structured for statistical analysis and visualization
- Generated automatically with `--process-all` or single genome analysis

## Visualization & Plotting

VARIANT includes comprehensive visualization capabilities for generating publication-ready figures from your mutation analysis data.

### 🎨 **Genome Organization Visualization**

Generate professional genome organization plots showing protein positions and mutation hotspots:

```bash
# Generate genome organization plots for all viruses
python plot/figure1_genome_analysis.py

# This creates interactive HTML and static PDF files:
# - plot/SARS-CoV-2_EPI_ISL_XXXXX_genome_organization.html
# - plot/HIV-1_MW881698.1_genome_organization.html
# - plot/Chikungunya_XXXXX_genome_organization.html
# - plot/ZaireEbola_XXXXX_genome_organization.html
```

**Features:**
- **Auto-detection**: Automatically finds proteome and mutation CSV files
- **Multi-virus support**: Works with SARS-CoV-2, HIV-1, Chikungunya, ZaireEbola
- **Mutation overlays**: Shows mutation hotspots on genome organization
- **Interactive plots**: HTML files with hover details and zoom capabilities
- **Publication-ready**: High-resolution PDF output for scientific papers

### 📊 **Available Plot Types**

| **Figure Type** | **Description** | **Data Source** | **Output Format** |
|-----------------|-----------------|-----------------|-------------------|
| **Genome Organization** | Protein positions + mutation hotspots | `mutation_summary.csv` | HTML + PDF |
| **Mutation Distribution** | Mutation type frequencies | `mutation_summary.csv` | HTML + PDF |
| **Protein Impact** | Amino acid change analysis | `mutation_summary.csv` | HTML + PDF |
| **Hot Mutations** | High-frequency mutation patterns | `row_hot_mutations.csv` | HTML + PDF |

### 🛠️ **Plotting Dependencies**

The visualization features require additional dependencies that are included in the updated environment:

```bash
# Core plotting dependencies (auto-installed)
plotly>=6.1.1      # Interactive plotting
kaleido>=1.0.0     # Static image export (PDF)
matplotlib>=3.5.0  # Additional plotting support
seaborn>=0.11.0    # Statistical visualization
```

### 📈 **Scientific Paper Integration**

These visualizations are designed for scientific publications:

- **Nature Reviews style**: Professional color schemes and typography
- **High resolution**: 1200x800px output suitable for print
- **Interactive elements**: HTML versions for online supplements
- **Real data integration**: Uses actual mutation analysis results
- **Multi-format output**: Both interactive (HTML) and static (PDF) versions

### 🔧 **Customization**

All plotting scripts support customization:

```python
# Example: Custom genome organization plot
from plot.figure1_genome_analysis import GenomeOrganizationVisualizer

# Create custom visualization
visualizer = GenomeOrganizationVisualizer('SARS-CoV-2', 'EPI_ISL_16127650')
visualizer.create_genome_organization_chart('my_custom_plot.html')
```

📖 **[Complete Visualization Guide](docs/plot_visualization_guide.md)** | 📊 **[Visualization Examples](docs/COMPREHENSIVE_VISUALIZATION_SUMMARY.md)**

## Recent Updates

### Version 2025.01.15 (Major Refactoring)
- **🔧 Complete main.py refactoring**: Eliminated code duplication, improved logic flow, and enhanced maintainability
- **🤖 Automatic mutation summary generation**: Removed `--generate-summaries` flag - summaries are now generated automatically
- **🚀 Enhanced `--process-all` functionality**: Now processes all virus types seamlessly with comprehensive output
- **🧬 Fixed complex protein coordinate parsing**: Handles `join()` coordinates correctly (e.g., RNA-dependent-polymerase)
- **📊 Improved mutation detection**: Better validation and error handling for all mutation types
- **🏗️ Clean architecture**: Introduced MutationProcessor class with proper separation of concerns
- **⚡ Better performance**: Optimized processing pipeline for both single and multi-segment viruses

### Version 2025.08.18
- **Streamlined `--process-all` command**: Removed generation of unnecessary files (root variants, majority variants, SNP records)
- **Clean output**: `--process-all` now only prints summary without creating additional files
- **Enhanced mutation summary generation**: Added `--generate-summaries` flag for comprehensive CSV output
- **Improved code organization**: Removed unused `_process_root_variants` methods for cleaner codebase

### Version 2025.08.07
- **Removed RNA editing functionality**: Focused analysis on PRF detection only
- **Simplified PRF output**: Removed confidence scores and impact analysis for cleaner results
- **CSV output format**: Changed from text to CSV for easier data analysis
- **Optimized processing**: SNP analysis is skipped when only PRF detection is requested
- **Enhanced tRNA validation**: Improved wobble base pairing rules for better detection accuracy

## Adding New Viruses

### Quick Setup (Recommended)

Use the setup script for easy virus configuration:

```bash
# Add virus to configuration
python scripts/setup_virus_dataset.py --add-virus "MyVirus" --reference "my_ref.fasta"

# Set up complete dataset
python scripts/setup_virus_dataset.py --setup "MyVirus" \
    --reference-path "/path/to/your/reference.fasta" \
    --proteome-path "/path/to/your/proteome.fasta" \
    --msa-files "/path/to/msa1.txt,/path/to/msa2.txt"

# Validate setup
python scripts/setup_virus_dataset.py --validate "MyVirus"
```

### Manual Setup

1. Create directory structure:
   ```bash
   # For single-segment virus
   mkdir -p data/NewVirus/{clustalW,refs,fasta}
   mkdir -p result/NewVirus
   
   # For multi-segment virus
   for i in {1..N}; do
     mkdir -p data/NewVirus/segment_$i/{clustalW,refs,fasta}
     mkdir -p result/NewVirus/segment_$i
   done
   ```

2. Add configuration to `virus_config.yaml`:
   ```yaml
   viruses:
     NewVirus:
       reference_genome: "my_ref.fasta"
       proteome_file: "my_proteome.fasta"
       codon_table_id: 1
       description: "My custom virus"
       default_msa_file: "my_msa.txt"
   ```

3. Place your files:
   ```bash
   # Copy reference genome
   cp /path/to/your/reference.fasta data/NewVirus/refs/my_ref.fasta
   
   # Copy proteome file
   cp /path/to/your/proteome.fasta data/NewVirus/refs/my_proteome.fasta
   
   # Copy MSA files
   cp /path/to/your/msa*.txt data/NewVirus/clustalW/
   ```

### File Requirements

#### Reference Genome File (Required)
- **Format**: FASTA format
- **Content**: Complete reference genome sequence
- **Example** (SARS-CoV-2 reference): [SARS-CoV-2 genome reference sequence](data/SARS-CoV-2/refs/NC_045512.fasta)

#### Proteome File (Required)
- **Format**: FASTA format
- **Content**: Protein sequences with headers
- **Example** (SARS-CoV-2 proteome reference): [SARS-CoV-2 proteome reference sequence](data/SARS-CoV-2/refs/SARS-CoV-2_proteome.fasta)

#### Variants Sequences File (Optional)
- **Format**: FASTA format
- **Content**: Reference genome sequence together with other variant genome sequences that user wanna decoded
- **Example** (SARS-CoV-2 proteome reference): [SARS-CoV-2 all variants sequences](data/SARS-CoV-2/fasta/sequences_w_reference.fasta)

#### Multiple Sequence Alignment (MSA) File (Required)
- **Format**: Text file with aligned sequences
- **Content**: Multiple aligned genome sequences obtained from [Clustal Omega](https://www.ebi.ac.uk/jdispatcher/msa/clustalo?stype=protein). The input file is [SARS-CoV-2 all variants sequences](data/SARS-CoV-2/fasta/sequences_w_reference.fasta)
- **Example** (SARS-CoV-2 MSA): [SARS-CoV-2 MSA](data/SARS-CoV-2/data/clustalW/gisaid_hcov-19_20221201_20221230_China_0_msa.txt)

### Codon Table IDs

Different organisms use different genetic codes:

- **1**: Standard genetic code (most viruses, bacteria, eukaryotes)
- **2**: Vertebrate mitochondrial code
- **3**: Yeast mitochondrial code
- **4**: Mold, protozoan, and coelenterate mitochondrial code
- **5**: Invertebrate mitochondrial code
- **6**: Ciliate, dasycladacean and hexamita nuclear code
- **9**: Echinoderm and flatworm mitochondrial code
- **10**: Euplotid nuclear code
- **11**: Bacterial, archaeal and plant plastid code
- **12**: Alternative yeast nuclear code
- **13**: Ascidian mitochondrial code
- **14**: Alternative flatworm mitochondrial code
- **16**: Chlorophycean mitochondrial code
- **21**: Trematode mitochondrial code
- **22**: Scenedesmus obliquus mitochondrial code
- **23**: Thraustochytrium mitochondrial code
- **24**: Pterobranchia mitochondrial code
- **25**: Candidate division SR1 and gracilibacteria code
- **26**: Pachysolen tannophilus nuclear code
- **27**: Karyorelict nuclear code
- **28**: Condylostoma nuclear code
- **29**: Mesodinium nuclear code
- **30**: Peritrich nuclear code
- **31**: Blastocrithidia nuclear code
- **33**: Cephalodiscidae mitochondrial code

### Example Setup

#### Example 1: Adding HIV Dataset
```bash
# 1. Add HIV to configuration
python scripts/setup_virus_dataset.py --add-virus "HIV" \
    --reference "HIV_reference.fasta" \
    --proteome "HIV_proteome.fasta" \
    --description "Human Immunodeficiency Virus"

# 2. Set up dataset with files
python scripts/setup_virus_dataset.py --setup "HIV" \
    --reference-path "/path/to/HIV_reference.fasta" \
    --proteome-path "/path/to/HIV_proteome.fasta" \
    --msa-files "/path/to/hiv_msa1.txt,/path/to/hiv_msa2.txt"

# 3. Run analysis
python main.py --virus HIV --process-all
```

#### Example 2: Manual Setup
```bash
# 1. Create directories
mkdir -p data/MyCustomVirus/{clustalW,refs,fasta}
mkdir -p result/MyCustomVirus

# 2. Add to virus_config.yaml
# Add the configuration shown above

# 3. Copy files
cp /path/to/reference.fasta data/MyCustomVirus/refs/my_ref.fasta
cp /path/to/proteome.fasta data/MyCustomVirus/refs/my_proteome.fasta
cp /path/to/msa.txt data/MyCustomVirus/clustalW/

# 4. Run analysis
python main.py --virus MyCustomVirus --genome-id GENOME_ID
```

## Troubleshooting

### Common Issues

1. **"Virus not found in configuration"**
   ```bash
   # Add the virus to configuration
   python scripts/setup_virus_dataset.py --add-virus "MyVirus" --reference "my_ref.fasta"
   ```

2. **"MSA file not found"**
   ```bash
   # Check available files
   ls data/MyVirus/clustalW/
   
   # Use a specific MSA file
   python main.py --virus MyVirus --msa-file my_msa.txt
   ```

3. **"Reference genome file not found"**
   ```bash
   # Check if file exists
   ls data/MyVirus/refs/
   
   # Copy your reference file
   cp /path/to/your/reference.fasta data/MyVirus/refs/my_ref.fasta
   ```

4. **"Proteome file not found"**
   - This is a warning, not an error
   - Analysis will continue with genome-level mutations only
   - Protein mutation analysis will be skipped

5. **Multi-segment Processing Issues**
   - Verify segment directory structure
   - Check segment-specific configurations
   - Ensure all required files exist for each segment

6. **PRF Detection Issues**
   - Verify reference genome is in DNA format (ATGC, not AUGC)
   - Check that MSA file contains the reference sequence
   - Ensure proper directory structure for virus data

### Validation

```bash
# Validate your dataset
python scripts/setup_virus_dataset.py --validate "MyVirus"

# List all configured viruses
python main.py --list-viruses

# Check file structure
tree data/MyVirus/
```

### Getting Help

```bash
# Setup script help
python scripts/setup_virus_dataset.py --help

# Main script help
python main.py --help

# List available viruses
python main.py --list-viruses
```

<!-- ## Citation

If you use this software in your research, please cite:

```bibtex
@software{wang2025variant,
  title={VARIANT: A Comprehensive Python Toolkit for Decoding and Analyzing Viral Mutations at Genome and Protein Levels},
  author={Wang, Rui},
  year={2025},
  url={https://github.com/wangru25/VARIANT}
}
``` -->

## Contact

- Author: Rui Wang
- Email: rw3594@nyu.edu
- Institution: New York University