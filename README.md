# Viralytics-mut: Advanced Viral Mutation Analysis Framework

A comprehensive Python framework for analyzing viral mutations, supporting both single-segment and multi-segment viruses. The framework provides detailed analysis of nucleotide changes, protein impacts, and mutation classifications including missense, nonsense, and frameshift mutations.

## Features

- Multi-segment virus support with automatic structure detection
- Comprehensive mutation analysis (point mutations, insertions, deletions)
- Advanced mutation classification (missense, nonsense, frameshift)
- Protein-level mutation impact analysis
- Antisense strand mutation support (complement notation)
- Structured output in both text and CSV formats
- Configurable per-virus settings

## Installation

### Prerequisites
- Python 3.8 or higher
- Required packages:
  - BioPython >= 1.79
  - NumPy >= 1.21.0
  - Pandas >= 1.3.0
  - PyYAML >= 6.0
  - fuzzysearch == 0.7.3
  - more-itertools >= 8.12.0

### Setup

```bash
# Clone the repository
git clone https://github.com/yourusername/viralytics-mut.git
cd viralytics-mut

# Install dependencies
pip install -r requirements.txt
```

## Usage

### Basic Commands

```bash
# Analyze a single genome
python main.py --virus SARS-CoV-2 --genome-id EPI_ISL_16327572

# Process all genomes in the MSA file
python main.py --virus SARS-CoV-2 --process-all

# Process specific segment of a multi-segment virus
python main.py --virus H3N2 --segment segment_1

# List available viruses
python main.py --list-viruses
```

### Command Line Arguments

- `--virus`: Name of the virus to analyze
- `--genome-id`: Specific genome ID to process
- `--process-all`: Process all genomes in the MSA file
- `--segment`: For multi-segment viruses, specify segment to analyze
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

### CSV Output (genome_id_row_hot_mutations_date.csv)
Structured data for statistical analysis:
- Genome ID
- Mutation type
- Position(s)
- Nucleotide changes
- Affected protein
- Amino acid changes
- Biological classification

## Adding New Viruses

1. Create directory structure:
   ```bash
   # For single-segment virus
   mkdir -p data/NewVirus/{clustalW,refs,fasta}
   
   # For multi-segment virus
   for i in {1..N}; do
     mkdir -p data/NewVirus/segment_$i/{clustalW,refs,fasta}
   done
   ```

2. Add configuration to virus_config.yaml
3. Place reference files in appropriate directories
4. Run analysis with new virus configuration

## Troubleshooting

Common issues and solutions:

1. **MSA File Not Found**
   - Verify file exists in correct directory
   - Check file name matches configuration

2. **Reference Genome/Proteome Issues**
   - Ensure files are in correct refs directory
   - Verify FASTA format is correct
   - Check file names match configuration

3. **Multi-segment Processing**
   - Verify segment directory structure
   - Check segment-specific configurations
   - Ensure all required files exist for each segment

## Citation

If you use this software in your research, please cite:

```bibtex
@software{wang2025viralyticsmut,
  title={Viralytics-mut: A Comprehensive Framework for Multi-Scale Viral Mutation Analysis},
  author={Wang, Rui},
  year={2025},
  url={https://github.com/yourusername/viralytics-mut}
}
```

## Contact

- Author: Rui Wang
- Email: rw3594@nyu.edu
- Institution: New York University