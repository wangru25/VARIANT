# Viralytics-mut: Comprehensive Viral Mutation Analysis Framework

A professional-grade Python package for analyzing viral mutations at multiple scales, from nucleotide changes to protein impacts, with integrated programmed ribosomal frameshifting detection for SARS-CoV-2 and other viruses.

## Quick Start

### Installation

```bash
# Install dependencies
pip install -r requirements.txt

# Or install dependencies individually
pip install biopython>=1.79 numpy>=1.21.0 pandas>=1.3.0 scikit-learn>=1.0.0 fuzzysearch==0.7.3 more-itertools>=8.12.0 pyyaml>=6.0
```

### How to Run

The project comes with sample SARS-CoV-2 data already configured. You can run it immediately:

#### **Option A: Analyze a Single Genome**
```bash
python main.py --virus SARS-CoV-2 --genome-id EPI_ISL_16327572
```

#### **Option B: Process All Genomes in the MSA File**
```bash
python main.py --virus SARS-CoV-2 --process-all
```

#### **Option C: Use a Different MSA File**
```bash
python main.py --virus SARS-CoV-2 --msa-file test_msa_1.txt
```

#### **Option D: List Available Viruses**
```bash
python main.py --list-viruses
```

### Command Line Arguments

- `--virus`: Virus name (default: "SARS-CoV-2")
- `--genome-id`: Specific genome ID to process (default: "EPI_ISL_16327572")
- `--msa-file`: MSA file name (uses default from config if not specified)
- `--process-all`: Process all genomes in the MSA file instead of just one
- `--list-viruses`: List all available viruses in configuration
- `--config`: Path to virus configuration file (default: "virus_config.yaml")

### Example Commands

```bash
# Quick test with a single genome
python main.py --virus SARS-CoV-2 --genome-id EPI_ISL_16327572

# Process all genomes (this will take longer)
python main.py --virus SARS-CoV-2 --process-all

# Use test data
python main.py --virus SARS-CoV-2 --msa-file test_msa_1.txt --genome-id EPI_ISL_16327572
```

### Output Files

Results are saved in:
- `result/SARS-CoV-2/` directory
- Individual genome files: `{genome_id}_{date}.txt`
- CSV files for row/hot mutations: `{genome_id}_row_hot_mutations_{date}.csv`

### Troubleshooting

If you get import errors, make sure you're running from the project root directory:
```bash
cd /path/to/viralytics-mut
python main.py --virus SARS-CoV-2 --genome-id EPI_ISL_16327572
```

## Features

- ** Multi-virus support**: Analyze mutations for SARS-CoV-2, HIV, Pox, and custom viruses
- ** Virus-specific organization**: Each virus has its own data and result directories
- **⚙ Flexible configuration**: Easy setup for new virus types via YAML configuration
- ** Comprehensive analysis**: Genome-level and protein-level mutation detection
- ** Hot mutation detection**: Advanced detection of non-consecutive nucleotide changes
- ** CSV output**: Structured results for further analysis
- ** User-friendly setup**: Helper scripts for dataset management

## Prerequisites

- **Python**: 3.8 or higher
- **Dependencies**: Automatically installed with the package
  - BioPython ≥ 1.79
  - NumPy ≥ 1.21.0
  - Pandas ≥ 1.3.0
  - PyYAML ≥ 6.0
  - Additional scientific computing libraries

## Installation Options

### Option 1: Direct Installation (Recommended)

```bash
# Clone the repository
git clone https://github.com/yourusername/viralytics-mut.git
cd viralytics-mut

# Install dependencies
pip install -r requirements.txt
```

### Option 2: Manual Dependency Installation

```bash
pip install biopython>=1.79 numpy>=1.21.0 pandas>=1.3.0 scikit-learn>=1.0.0 fuzzysearch==0.7.3 more-itertools>=8.12.0 pyyaml>=6.0
```

### Option 3: Conda Installation

```bash
conda install biopython numpy pandas scikit-learn pyyaml
pip install fuzzysearch more-itertools
```

## Quick Examples

### Analyze a Single Genome

```bash
python main.py --virus SARS-CoV-2 --genome-id EPI_ISL_16327572
```

### Process All Genomes in an MSA File

```bash
python main.py --virus SARS-CoV-2 --process-all
```

### Use a Custom MSA File

```bash
python main.py --virus SARS-CoV-2 --msa-file test_msa_1.txt
```

### List Available Viruses

```bash
python main.py --list-viruses
```

## 📁 Adding Your Own Virus Dataset

### Method 1: Manual Setup (Recommended)

```bash
# 1. Add virus to configuration in virus_config.yaml
# 2. Create directory structure
mkdir -p data/MyVirus/{clustalW,refs,fasta}
mkdir -p result/MyVirus

# 3. Place your files in the appropriate directories
# 4. Run analysis
python main.py --virus MyVirus
```

### Method 2: Manual Setup

1. **Create directory structure**:
   ```
   data/
   └── YourVirus/
       ├── clustalW/     # MSA files
       ├── refs/         # Reference genome and proteome
       └── fasta/        # FASTA files
   result/
   └── YourVirus/        # Output files
   ```

2. **Add virus configuration** to `virus_config.yaml`:
   ```yaml
   viruses:
     YourVirus:
       reference_genome: "your_reference.fasta"
       proteome_file: "your_proteome.fasta"
       codon_table_id: 1
       description: "Description of your virus"
       default_msa_file: "your_msa_file.txt"
   ```

3. **Place your files**:
   - Reference genome: `data/YourVirus/refs/your_reference.fasta`
   - Proteome file: `data/YourVirus/refs/your_proteome.fasta`
   - MSA files: `data/YourVirus/clustalW/`

4. **Run analysis**:
   ```bash
   python main.py --virus YourVirus
   ```

## 📄 File Requirements

### Reference Genome File
- **Format**: FASTA format
- **Content**: Complete reference genome sequence
- **Example**:
  ```
  >NC_045512.2
  ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATCTTCTGCAGGCTGCTTACGGTTTCGTCCGTGTTGCAGCCGATCATCAGCACATCTAGGTTTCGTCCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAACACACGTCCAACTCAGTTTGCCTGTTTTACAGGTTCGCGACGTGCTCGTACGTGGCTTTGGAGACTCCGTGGAGGAGGTCTTATCAGAGGCACGTCAACATCTTAAAGATGGCACTTGTGGCTTAGTAGAAGTTGAAAAAGGCGTTTTGCCTCAACTTGAACAGCCCTATGTGTTCATCAAACGTTCGGATGCTCGAACTGCACCTCATGGTCATGTTATGGTTGAGCTGGTAGCAGAACTCGAAGGCATTCAGTACGGTCGTAGTGGTGAGACACTTGGTGTCCTTGTCCCTCATGTGGGCGAAATACCAGTGGCTTACCGCAAGGTTCTTCTTCGTAAGAACGGTAATAAAGGAGCTGGTGGCCATAGTTACGGCGCCGATCTAAAGTCATTTGACTTAGGCGACGAGCTTGGCACTGATCCTTATGAAGATTTTCAAGAAAACTGGAACACTAAACATAGCAGTGGTGTTACCCGTGAACTCATGCGTGAGCTTAACGGAGGGGCATACACTCGCTATGTCGATAACAACTTCTGTGGCCCTGATGGCTACCCTCTTGAGTGCATTAAAGACCTTCTAGCACGTGCTGGTAAAGCTTCATGCACTTTGTCCGAACAACTGGACTTTATTGACACTAAGAGGGGTGTATACTGCTGCCGTGAACATGAGCATGAAATTGCTTGGTACACGGAACGTTCTGAAAAGAGCTATGAATTGCAGACACCTTTTGAAATTAAATTGGCAAAGAAATTTGACACCTTCAATGGGGAATGTCCAAATTTTGTATTTCCCTTAAATTCCATAATCAAGACTATTCAACCAAGGGTTGAAAAGAAAAAGCTTGATGGCTTTATGGGTAGAATTCGATCTGTCTATCCAGTTGCGTCACCAAATGAATGCAACCAAATGTGCCTTTCAACTCTCATGAAGTGTGATCATTGTGGTGAAACTTCATGGCAGACGGGCGATTTTGTTAAAGCCACTTGCGAATTTTGTGGCACTGAGAATTTGACTAAAGAAGGTGCCACTACTTGTGGTTACTTACCCCAAAATGCTGTTGTTAAAATTTATTGTCCAGCATGTCACAATTCAGAAGTAGGACCTGAGCATAGTCTTGCCGAATACCATAATGAATCTGGCTTGAAAACCATTCTTCGTAAGGGTGGTCGCACTATTGCCTTTGGAGGCTGTGTGTTCTCTTATGTTGGTTGCCATAACAAGTGTGCCTATTGGGTTCCACGTGCTAGCGCTAACATAGGTTGTAACCATACAGGTGTTGTTGGAGAAGGTTCCGAAGGTCTTAATGACAACCTTCTTGAAATACTCCAAAAAGAGAAAGTCAACATCAATATTGTTGGTGACTTTAAACTTAATGAAGAGATCGCCATTATTTTGGCATCTTTTTCTGCTTCCACAAGTGCTTTTGTGGAAACTGTGAAAGGTTTGGATTATAAAGCATTCAAACAAATTGTTGAATCCTGTGGTAATTTTAAAGTTACAAAAGGAAAAGCTAAAAAAGGTGCCTGGAATATTGGTGAACAGAAATCAATACTGAGTCCTCTTTATGCATTTGCATCAGAGGCTGCTCGTGTTGTACGATCAATACTTATCCTACTGATGGAAGCAAGTCATTCGAAGGACAGCTTAGTGATATTAAT
  ```

### Proteome File
- **Format**: FASTA format
- **Content**: Protein sequences with headers
- **Example**:
  ```
  >sp|P0DTD1|R1AB_SARS2 Replicase polyprotein 1ab
  MESLVPGFNEKTHVLLLDRSELANPVALNAGLMTHLPQDRTPYGSFKTCLRKGFSGQVYKDHVKLRGARDLQPGLNPVLPPVLTKLTVGARQLVDFLKRFLYALKRKGLSPVFLENKCNGVLGTFVLTKLSDVGVFNVESLKSFVLGKGFSKKIFKQVVDFTSEVLDFVGTDGSSEVAVKMFGPVTTVSNYVLTFPHGAGLSSYVKSIQYLTNSTFEYPITNDLVVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSHRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT
  ```

### MSA (Multiple Sequence Alignment) File
- **Format**: Text file with aligned sequences
- **Content**: Multiple aligned genome sequences
- **Example**:
  ```
  >hCoV-19/Wuhan/WH01/2019|EPI_ISL_406801|2019-12-26
  ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATCTTCTGCAGGCTGCTTACGGTTTCGTCCGTGTTGCAGCCGATCATCAGCACATCTAGGTTTCGTCCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAACACACGTCCAACTCAGTTTGCCTGTTTTACAGGTTCGCGACGTGCTCGTACGTGGCTTTGGAGACTCCGTGGAGGAGGTCTTATCAGAGGCACGTCAACATCTTAAAGATGGCACTTGTGGCTTAGTAGAAGTTGAAAAAGGCGTTTTGCCTCAACTTGAACAGCCCTATGTGTTCATCAAACGTTCGGATGCTCGAACTGCACCTCATGGTCATGTTATGGTTGAGCTGGTAGCAGAACTCGAAGGCATTCAGTACGGTCGTAGTGGTGAGACACTTGGTGTCCTTGTCCCTCATGTGGGCGAAATACCAGTGGCTTACCGCAAGGTTCTTCTTCGTAAGAACGGTAATAAAGGAGCTGGTGGCCATAGTTACGGCGCCGATCTAAAGTCATTTGACTTAGGCGACGAGCTTGGCACTGATCCTTATGAAGATTTTCAAGAAAACTGGAACACTAAACATAGCAGTGGTGTTACCCGTGAACTCATGCGTGAGCTTAACGGAGGGGCATACACTCGCTATGTCGATAACAACTTCTGTGGCCCTGATGGCTACCCTCTTGAGTGCATTAAAGACCTTCTAGCACGTGCTGGTAAAGCTTCATGCACTTTGTCCGAACAACTGGACTTTATTGACACTAAGAGGGGTGTATACTGCTGCCGTGAACATGAGCATGAAATTGCTTGGTACACGGAACGTTCTGAAAAGAGCTATGAATTGCAGACACCTTTTGAAATTAAATTGGCAAAGAAATTTGACACCTTCAATGGGGAATGTCCAAATTTTGTATTTCCCTTAAATTCCATAATCAAGACTATTCAACCAAGGGTTGAAAAGAAAAAGCTTGATGGCTTTATGGGTAGAATTCGATCTGTCTATCCAGTTGCGTCACCAAATGAATGCAACCAAATGTGCCTTTCAACTCTCATGAAGTGTGATCATTGTGGTGAAACTTCATGGCAGACGGGCGATTTTGTTAAAGCCACTTGCGAATTTTGTGGCACTGAGAATTTGACTAAAGAAGGTGCCACTACTTGTGGTTACTTACCCCAAAATGCTGTTGTTAAAATTTATTGTCCAGCATGTCACAATTCAGAAGTAGGACCTGAGCATAGTCTTGCCGAATACCATAATGAATCTGGCTTGAAAACCATTCTTCGTAAGGGTGGTCGCACTATTGCCTTTGGAGGCTGTGTGTTCTCTTATGTTGGTTGCCATAACAAGTGTGCCTATTGGGTTCCACGTGCTAGCGCTAACATAGGTTGTAACCATACAGGTGTTGTTGGAGAAGGTTCCGAAGGTCTTAATGACAACCTTCTTGAAATACTCCAAAAAGAGAAAGTCAACATCAATATTGTTGGTGACTTTAAACTTAATGAAGAGATCGCCATTATTTTGGCATCTTTTTCTGCTTCCACAAGTGCTTTTGTGGAAACTGTGAAAGGTTTGGATTATAAAGCATTCAAACAAATTGTTGAATCCTGTGGTAATTTTAAAGTTACAAAAGGAAAAGCTAAAAAAGGTGCCTGGAATATTGGTGAACAGAAATCAATACTGAGTCCTCTTTATGCATTTGCATCAGAGGCTGCTCGTGTTGTACGATCAATACTTATCCTACTGATGGAAGCAAGTCATTCGAAGGACAGCTTAGTGATATTAAT
  ```

## ⚙️ Configuration

### Virus Configuration File (`virus_config.yaml`)

The configuration file allows you to define virus-specific settings:

```yaml
viruses:
  SARS-CoV-2:
    reference_genome: "NC_045512.fasta"
    proteome_file: "SARS-CoV-2_proteome.fasta"
    codon_table_id: 1
    description: "SARS-CoV-2 coronavirus"
    default_msa_file: "gisaid_hcov-19_20221201_20221230_China_0_msa.txt"
    
  YourVirus:
    reference_genome: "your_reference.fasta"
    proteome_file: "your_proteome.fasta"
    codon_table_id: 1
    description: "Description of your virus"
    default_msa_file: "your_msa_file.txt"
```

### Codon Table IDs

Different organisms may use different genetic codes. Common codon table IDs:

- **1**: Standard genetic code (most viruses, bacteria, eukaryotes)
- **2**: Vertebrate mitochondrial code
- **3**: Yeast mitochondrial code
- **4**: Mold, protozoan, and coelenterate mitochondrial code
- **5**: Invertebrate mitochondrial code

## Output

Results are saved in virus-specific directories:

```
result/
└── YourVirus/
    ├── genome_id_20250623.txt
    ├── genome_id_row_hot_mutations_20250623.csv
    └── ...
```

### Output Files

1. **Text files** (`.txt`): Detailed mutation analysis
   - Mutation type (point, deletion, insertion, etc.)
   - Position in the genome
   - Nucleotide change
   - Protein mutation (if proteome file is available)

2. **CSV files** (`.csv`): Structured data for analysis
   - Genome ID
   - Mutation type
   - Position
   - Nucleotide change
   - Protein affected
   - Amino acid change (formatted as `"['S375F', 'T376A']"`)
   - Biological classification

## 🔧 Troubleshooting

### Common Issues

1. **"MSA file not found"**
   - Check that your MSA file is in the correct directory: `data/YourVirus/clustalW/`
   - Verify the filename matches your configuration

2. **"Reference genome file not found"**
   - Ensure your reference genome is in: `data/YourVirus/refs/`
   - Check the filename in your virus configuration

3. **"Proteome file not found"**
   - This is a warning, not an error
   - Analysis will continue with genome-level mutations only
   - Protein mutation analysis will be skipped

4. **"Virus not found in configuration"**
   - Add the virus to `virus_config.yaml` or use the setup script
   - Run: `viralytics-mut setup --add-virus "YourVirus"`

### Getting Help

1. **Validate your dataset**:
   ```bash
   viralytics-mut setup --validate "YourVirus"
   ```

2. **Check available files**:
   ```bash
   ls data/YourVirus/clustalW/
   ls data/YourVirus/refs/
   ```

3. **List configured viruses**:
   ```bash
   viralytics-mut list-viruses
   ```

## Contributing

We welcome contributions! To add support for new virus types or improve the tool:

1. Fork the repository
2. Create a feature branch
3. Add your virus configuration to `virus_config.yaml`
4. Test with your dataset
5. Submit a pull request

## 📄 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 📚 Citation

If you use Viralytics-mut in your research, please cite:

```bibtex
@software{wang2025viralyticsmut,
  title={Viralytics-mut: A Comprehensive Framework for Multi-Scale Viral Mutation Analysis},
  author={Wang, Rui},
  year={2025},
  url={https://github.com/yourusername/viralytics-mut}
}
```

## Contact

- **Author**: Rui Wang
- **Email**: rw3594@nyu.edu
- **Institution**: New York University
- **GitHub**: [@yourusername](https://github.com/yourusername)

## Acknowledgments

- BioPython community for sequence analysis tools
- SARS-CoV-2 research community for reference genomes
- Contributors and users of this package

## Command-Line Usage

After installing or cloning this repository, you can use the `viralytics-mut` command-line tool for all major functions.

### Basic Usage

```sh
# Make sure the script is executable
chmod +x viralytics-mut

# Show help
./viralytics-mut --help

# Analyze mutations for a virus
./viralytics-mut analyze --virus SARS-CoV-2 --msa data/SARS-CoV-2/clustalW/msa.txt

# Setup a new virus dataset
./viralytics-mut setup --virus HIV --config virus_config.yaml

# Analyze PRF sites
./viralytics-mut prf --genome data/refs/NC_045512.fasta
```

You can also run it with Python:
```sh
python viralytics-mut analyze --virus SARS-CoV-2 --msa data/SARS-CoV-2/clustalW/msa.txt
```

If you distribute your package via pip, the `viralytics-mut` command will be available globally after install. If running from source, use `./viralytics-mut` from the project root.
