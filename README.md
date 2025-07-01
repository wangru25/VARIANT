# Viralytics-mut: Comprehensive Viral Mutation Analysis Framework

A professional-grade Python package for analyzing viral mutations at multiple scales, from nucleotide changes to protein impacts, with integrated programmed ribosomal frameshifting detection for SARS-CoV-2 and other viruses.

## 🚀 Quick Start

### Installation

```bash
# Install from PyPI (recommended)
pip install viralytics-mut

# Or install from source
git clone https://github.com/yourusername/viralytics-mut.git
cd viralytics-mut
pip install -e .
```

### Basic Usage

```bash
# List available viruses
viralytics-mut list-viruses

# Analyze mutations for SARS-CoV-2
viralytics-mut analyze --virus SARS-CoV-2 --genome-id EPI_ISL_16327572

# Process all genomes in an MSA file
viralytics-mut analyze --virus SARS-CoV-2 --process-all
```

## ✨ Features

- **🔬 Multi-virus support**: Analyze mutations for SARS-CoV-2, HIV, Pox, and custom viruses
- **📁 Virus-specific organization**: Each virus has its own data and result directories
- **⚙️ Flexible configuration**: Easy setup for new virus types via YAML configuration
- **🧬 Comprehensive analysis**: Genome-level and protein-level mutation detection
- **🔄 Hot mutation detection**: Advanced detection of non-consecutive nucleotide changes
- **📊 CSV output**: Structured results for further analysis
- **🛠️ User-friendly setup**: Helper scripts for dataset management

## 📋 Prerequisites

- **Python**: 3.8 or higher
- **Dependencies**: Automatically installed with the package
  - BioPython ≥ 1.79
  - NumPy ≥ 1.21.0
  - Pandas ≥ 1.3.0
  - PyYAML ≥ 6.0
  - Additional scientific computing libraries

## 🎯 Installation Options

### Option 1: PyPI Installation (Recommended)

```bash
pip install viralytics-mut
```

### Option 2: Development Installation

```bash
git clone https://github.com/yourusername/viralytics-mut.git
cd viralytics-mut
pip install -e .
```

### Option 3: Conda Installation

```bash
conda install -c conda-forge viralytics-mut
```

## 🏃‍♂️ Quick Examples

### Analyze a Single Genome

```bash
mutparser analyze --virus SARS-CoV-2 --genome-id EPI_ISL_16327572
```

### Process All Genomes in an MSA File

```bash
mutparser analyze --virus SARS-CoV-2 --process-all
```

### Use a Custom MSA File

```bash
mutparser analyze --virus SARS-CoV-2 --msa-file custom_msa.txt
```

### Setup a New Virus Dataset

```bash
mutparser setup --virus HIV --config virus_config.yaml
```

## 📁 Adding Your Own Virus Dataset

### Method 1: Using the Setup Script (Recommended)

```bash
# 1. Add virus to configuration
viralytics-mut setup --add-virus "MyVirus" --reference "my_ref.fasta"

# 2. Set up complete dataset with files
viralytics-mut setup --setup "MyVirus" \
    --reference-path "path/to/your/reference.fasta" \
    --proteome-path "path/to/your/proteome.fasta" \
    --msa-files "path/to/msa1.txt,path/to/msa2.txt"

# 3. Validate the setup
viralytics-mut setup --validate "MyVirus"

# 4. Run analysis
viralytics-mut analyze --virus MyVirus
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
   viralytics-mut analyze --virus YourVirus
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

## 📊 Output

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

## 🤝 Contributing

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

## 📞 Contact

- **Author**: Rui Wang
- **Email**: rw3594@nyu.edu
- **Institution**: New York University
- **GitHub**: [@yourusername](https://github.com/yourusername)

## 🙏 Acknowledgments

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
