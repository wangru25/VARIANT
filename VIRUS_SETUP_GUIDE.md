# Virus Dataset Setup Guide

This guide explains how to add your own virus datasets to MutParser and run analysis smoothly.

## Overview

The new MutParser system supports multiple virus types with virus-specific organization:

```
data/
├── SARS-CoV-2/          # SARS-CoV-2 data
│   ├── clustalW/        # MSA files
│   ├── refs/            # Reference genome & proteome
│   └── fasta/           # FASTA files
├── HIV/                 # HIV data
│   ├── clustalW/
│   ├── refs/
│   └── fasta/
└── YourVirus/           # Your custom virus
    ├── clustalW/
    ├── refs/
    └── fasta/

result/
├── SARS-CoV-2/          # SARS-CoV-2 results
├── HIV/                 # HIV results
└── YourVirus/           # Your virus results
```

## Quick Setup (Recommended)

### Step 1: Add Your Virus to Configuration

```bash
python setup_virus_dataset.py --add-virus "MyVirus" --reference "my_ref.fasta"
```

This will:
- Add your virus to `virus_config.yaml`
- Create the directory structure
- Set default file names

### Step 2: Set Up Your Dataset

```bash
python setup_virus_dataset.py --setup "MyVirus" \
    --reference-path "path/to/your/reference.fasta" \
    --proteome-path "path/to/your/proteome.fasta" \
    --msa-files "path/to/msa1.txt,path/to/msa2.txt"
```

This will:
- Copy your files to the correct locations
- Validate the setup
- Make your dataset ready for analysis

### Step 3: Run Analysis

```bash
python main.py --virus MyVirus
```

## Detailed Setup

### Method 1: Using the Setup Script

#### 1. Add Virus Configuration

```bash
# Basic setup
python setup_virus_dataset.py --add-virus "MyVirus" --reference "my_ref.fasta"

# Advanced setup with all options
python setup_virus_dataset.py --add-virus "MyVirus" \
    --reference "my_ref.fasta" \
    --proteome "my_proteome.fasta" \
    --codon-table 1 \
    --description "My custom virus for analysis" \
    --default-msa "my_msa.txt"
```

#### 2. Set Up Complete Dataset

```bash
python setup_virus_dataset.py --setup "MyVirus" \
    --reference-path "/path/to/your/reference.fasta" \
    --proteome-path "/path/to/your/proteome.fasta" \
    --msa-files "/path/to/msa1.txt,/path/to/msa2.txt,/path/to/msa3.txt"
```

#### 3. Validate Your Setup

```bash
python setup_virus_dataset.py --validate "MyVirus"
```

### Method 2: Manual Setup

#### 1. Create Directory Structure

```bash
mkdir -p data/MyVirus/clustalW
mkdir -p data/MyVirus/refs
mkdir -p data/MyVirus/fasta
mkdir -p result/MyVirus
```

#### 2. Add Configuration

Edit `virus_config.yaml`:

```yaml
viruses:
  MyVirus:
    reference_genome: "my_ref.fasta"
    proteome_file: "my_proteome.fasta"
    codon_table_id: 1
    description: "My custom virus"
    default_msa_file: "my_msa.txt"
```

#### 3. Place Your Files

```bash
# Copy reference genome
cp /path/to/your/reference.fasta data/MyVirus/refs/my_ref.fasta

# Copy proteome file
cp /path/to/your/proteome.fasta data/MyVirus/refs/my_proteome.fasta

# Copy MSA files
cp /path/to/your/msa*.txt data/MyVirus/clustalW/
```

## File Requirements

### Reference Genome File

**Format**: FASTA format
**Content**: Complete reference genome sequence

**Example**:
```
>NC_045512.2
ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATCTTCTGCAGGCTGCTTACGGTTTCGTCCGTGTTGCAGCCGATCATCAGCACATCTAGGTTTCGTCCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAACACACGTCCAACTCAGTTTGCCTGTTTTACAGGTTCGCGACGTGCTCGTACGTGGCTTTGGAGACTCCGTGGAGGAGGTCTTATCAGAGGCACGTCAACATCTTAAAGATGGCACTTGTGGCTTAGTAGAAGTTGAAAAAGGCGTTTTGCCTCAACTTGAACAGCCCTATGTGTTCATCAAACGTTCGGATGCTCGAACTGCACCTCATGGTCATGTTATGGTTGAGCTGGTAGCAGAACTCGAAGGCATTCAGTACGGTCGTAGTGGTGAGACACTTGGTGTCCTTGTCCCTCATGTGGGCGAAATACCAGTGGCTTACCGCAAGGTTCTTCTTCGTAAGAACGGTAATAAAGGAGCTGGTGGCCATAGTTACGGCGCCGATCTAAAGTCATTTGACTTAGGCGACGAGCTTGGCACTGATCCTTATGAAGATTTTCAAGAAAACTGGAACACTAAACATAGCAGTGGTGTTACCCGTGAACTCATGCGTGAGCTTAACGGAGGGGCATACACTCGCTATGTCGATAACAACTTCTGTGGCCCTGATGGCTACCCTCTTGAGTGCATTAAAGACCTTCTAGCACGTGCTGGTAAAGCTTCATGCACTTTGTCCGAACAACTGGACTTTATTGACACTAAGAGGGGTGTATACTGCTGCCGTGAACATGAGCATGAAATTGCTTGGTACACGGAACGTTCTGAAAAGAGCTATGAATTGCAGACACCTTTTGAAATTAAATTGGCAAAGAAATTTGACACCTTCAATGGGGAATGTCCAAATTTTGTATTTCCCTTAAATTCCATAATCAAGACTATTCAACCAAGGGTTGAAAAGAAAAAGCTTGATGGCTTTATGGGTAGAATTCGATCTGTCTATCCAGTTGCGTCACCAAATGAATGCAACCAAATGTGCCTTTCAACTCTCATGAAGTGTGATCATTGTGGTGAAACTTCATGGCAGACGGGCGATTTTGTTAAAGCCACTTGCGAATTTTGTGGCACTGAGAATTTGACTAAAGAAGGTGCCACTACTTGTGGTTACTTACCCCAAAATGCTGTTGTTAAAATTTATTGTCCAGCATGTCACAATTCAGAAGTAGGACCTGAGCATAGTCTTGCCGAATACCATAATGAATCTGGCTTGAAAACCATTCTTCGTAAGGGTGGTCGCACTATTGCCTTTGGAGGCTGTGTGTTCTCTTATGTTGGTTGCCATAACAAGTGTGCCTATTGGGTTCCACGTGCTAGCGCTAACATAGGTTGTAACCATACAGGTGTTGTTGGAGAAGGTTCCGAAGGTCTTAATGACAACCTTCTTGAAATACTCCAAAAAGAGAAAGTCAACATCAATATTGTTGGTGACTTTAAACTTAATGAAGAGATCGCCATTATTTTGGCATCTTTTTCTGCTTCCACAAGTGCTTTTGTGGAAACTGTGAAAGGTTTGGATTATAAAGCATTCAAACAAATTGTTGAATCCTGTGGTAATTTTAAAGTTACAAAAGGAAAAGCTAAAAAAGGTGCCTGGAATATTGGTGAACAGAAATCAATACTGAGTCCTCTTTATGCATTTGCATCAGAGGCTGCTCGTGTTGTACGATCAATACTTATCCTACTGATGGAAGCAAGTCATTCGAAGGACAGCTTAGTGATATTAAT
```

### Proteome File (Optional)

**Format**: FASTA format
**Content**: Protein sequences with headers

**Example**:
```
>sp|P0DTD1|R1AB_SARS2 Replicase polyprotein 1ab
MESLVPGFNEKTHVLLLDRSELANPVALNAGLMTHLPQDRTPYGSFKTCLRKGFSGQVYKDHVKLRGARDLQPGLNPVLPPVLTKLTVGARQLVDFLKRFLYALKRKGLSPVFLENKCNGVLGTFVLTKLSDVGVFNVESLKSFVLGKGFSKKIFKQVVDFTSEVLDFVGTDGSSEVAVKMFGPVTTVSNYVLTFPHGAGLSSYVKSIQYLTNSTFEYPITNDLVVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSHRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT
```

### MSA (Multiple Sequence Alignment) File

**Format**: Text file with aligned sequences
**Content**: Multiple aligned genome sequences

**Example**:
```
>hCoV-19/Wuhan/WH01/2019|EPI_ISL_406801|2019-12-26
ATTAAAGGTTTATACCTTCCCAGGTAACAAACCAACCAACTTTCGATCTCTTGTAGATCTGTTCTCTAAACGAACTTTAAAATCTGTGTGGCTGTCACTCGGCTGCATGCTTAGTGCACTCACGCAGTATAATTAATAACTAATTACTGTCGTTGACAGGACACGAGTAACTCGTCTATCTTCTGCAGGCTGCTTACGGTTTCGTCCGTGTTGCAGCCGATCATCAGCACATCTAGGTTTCGTCCGGGTGTGACCGAAAGGTAAGATGGAGAGCCTTGTCCCTGGTTTCAACGAGAAAACACACGTCCAACTCAGTTTGCCTGTTTTACAGGTTCGCGACGTGCTCGTACGTGGCTTTGGAGACTCCGTGGAGGAGGTCTTATCAGAGGCACGTCAACATCTTAAAGATGGCACTTGTGGCTTAGTAGAAGTTGAAAAAGGCGTTTTGCCTCAACTTGAACAGCCCTATGTGTTCATCAAACGTTCGGATGCTCGAACTGCACCTCATGGTCATGTTATGGTTGAGCTGGTAGCAGAACTCGAAGGCATTCAGTACGGTCGTAGTGGTGAGACACTTGGTGTCCTTGTCCCTCATGTGGGCGAAATACCAGTGGCTTACCGCAAGGTTCTTCTTCGTAAGAACGGTAATAAAGGAGCTGGTGGCCATAGTTACGGCGCCGATCTAAAGTCATTTGACTTAGGCGACGAGCTTGGCACTGATCCTTATGAAGATTTTCAAGAAAACTGGAACACTAAACATAGCAGTGGTGTTACCCGTGAACTCATGCGTGAGCTTAACGGAGGGGCATACACTCGCTATGTCGATAACAACTTCTGTGGCCCTGATGGCTACCCTCTTGAGTGCATTAAAGACCTTCTAGCACGTGCTGGTAAAGCTTCATGCACTTTGTCCGAACAACTGGACTTTATTGACACTAAGAGGGGTGTATACTGCTGCCGTGAACATGAGCATGAAATTGCTTGGTACACGGAACGTTCTGAAAAGAGCTATGAATTGCAGACACCTTTTGAAATTAAATTGGCAAAGAAATTTGACACCTTCAATGGGGAATGTCCAAATTTTGTATTTCCCTTAAATTCCATAATCAAGACTATTCAACCAAGGGTTGAAAAGAAAAAGCTTGATGGCTTTATGGGTAGAATTCGATCTGTCTATCCAGTTGCGTCACCAAATGAATGCAACCAAATGTGCCTTTCAACTCTCATGAAGTGTGATCATTGTGGTGAAACTTCATGGCAGACGGGCGATTTTGTTAAAGCCACTTGCGAATTTTGTGGCACTGAGAATTTGACTAAAGAAGGTGCCACTACTTGTGGTTACTTACCCCAAAATGCTGTTGTTAAAATTTATTGTCCAGCATGTCACAATTCAGAAGTAGGACCTGAGCATAGTCTTGCCGAATACCATAATGAATCTGGCTTGAAAACCATTCTTCGTAAGGGTGGTCGCACTATTGCCTTTGGAGGCTGTGTGTTCTCTTATGTTGGTTGCCATAACAAGTGTGCCTATTGGGTTCCACGTGCTAGCGCTAACATAGGTTGTAACCATACAGGTGTTGTTGGAGAAGGTTCCGAAGGTCTTAATGACAACCTTCTTGAAATACTCCAAAAAGAGAAAGTCAACATCAATATTGTTGGTGACTTTAAACTTAATGAAGAGATCGCCATTATTTTGGCATCTTTTTCTGCTTCCACAAGTGCTTTTGTGGAAACTGTGAAAGGTTTGGATTATAAAGCATTCAAACAAATTGTTGAATCCTGTGGTAATTTTAAAGTTACAAAAGGAAAAGCTAAAAAAGGTGCCTGGAATATTGGTGAACAGAAATCAATACTGAGTCCTCTTTATGCATTTGCATCAGAGGCTGCTCGTGTTGTACGATCAATACTTATCCTACTGATGGAAGCAAGTCATTCGAAGGACAGCTTAGTGATATTAAT
```

## Running Analysis

### Basic Analysis

```bash
# Analyze a single genome
python main.py --virus MyVirus --genome-id GENOME_ID

# Process all genomes in an MSA file
python main.py --virus MyVirus --process-all

# Use a specific MSA file
python main.py --virus MyVirus --msa-file my_custom_msa.txt
```

### Advanced Options

```bash
# Use a custom configuration file
python main.py --virus MyVirus --config my_config.yaml

# List all available viruses
python main.py --list-viruses

# Get help
python main.py --help
```

## Configuration Options

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

### Example Configuration

```yaml
viruses:
  MyVirus:
    reference_genome: "my_ref.fasta"
    proteome_file: "my_proteome.fasta"
    codon_table_id: 1
    description: "My custom virus for analysis"
    default_msa_file: "my_msa.txt"
```

## Troubleshooting

### Common Issues

1. **"Virus not found in configuration"**
   ```bash
   # Add the virus to configuration
   python setup_virus_dataset.py --add-virus "MyVirus" --reference "my_ref.fasta"
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

### Validation

```bash
# Validate your dataset
python setup_virus_dataset.py --validate "MyVirus"

# List all configured viruses
python main.py --list-viruses

# Check file structure
tree data/MyVirus/
```

### Getting Help

```bash
# Setup script help
python setup_virus_dataset.py --help

# Main script help
python main.py --help

# List available viruses
python main.py --list-viruses
```

## Output

Results are saved in virus-specific directories:

```
result/
└── MyVirus/
    ├── genome_id_20250623.txt
    ├── another_genome_20250623.txt
    └── ...
```

Each output file contains:
- Mutation type (point, deletion, insertion, etc.)
- Position in the genome
- Nucleotide change
- Protein mutation (if proteome file is available)

## Examples

### Example 1: Adding HIV Dataset

```bash
# 1. Add HIV to configuration
python setup_virus_dataset.py --add-virus "HIV" \
    --reference "HIV_reference.fasta" \
    --proteome "HIV_proteome.fasta" \
    --description "Human Immunodeficiency Virus"

# 2. Set up dataset with files
python setup_virus_dataset.py --setup "HIV" \
    --reference-path "/path/to/HIV_reference.fasta" \
    --proteome-path "/path/to/HIV_proteome.fasta" \
    --msa-files "/path/to/hiv_msa1.txt,/path/to/hiv_msa2.txt"

# 3. Run analysis
python main.py --virus HIV --process-all
```

### Example 2: Adding Custom Virus

```bash
# 1. Add custom virus
python setup_virus_dataset.py --add-virus "MyCustomVirus" \
    --reference "custom_ref.fasta" \
    --codon-table 2 \
    --description "My custom virus with mitochondrial genetic code"

# 2. Set up dataset
python setup_virus_dataset.py --setup "MyCustomVirus" \
    --reference-path "/path/to/custom_ref.fasta" \
    --msa-files "/path/to/custom_msa.txt"

# 3. Validate
python setup_virus_dataset.py --validate "MyCustomVirus"

# 4. Run analysis
python main.py --virus MyCustomVirus --genome-id CUSTOM_GENOME_ID
```

## Tips

1. **Use descriptive virus names**: Avoid spaces and special characters
2. **Keep file names consistent**: Use the same names in configuration and actual files
3. **Validate before running**: Always validate your dataset before analysis
4. **Check file formats**: Ensure your files are in the correct format
5. **Backup your data**: Keep backups of your original files
6. **Use version control**: Track changes to your configuration files

## Support

If you encounter issues:

1. Check the troubleshooting section above
2. Validate your dataset using the setup script
3. Check file permissions and paths
4. Ensure all required dependencies are installed
5. Review the main README.md for additional information 