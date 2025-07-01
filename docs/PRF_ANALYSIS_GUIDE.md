<!--
 * @Author: Rui Wang
 * @Date: 2025-06-30 08:51:21
 * @LastModifiedBy: Rui Wang
 * @LastEditTime: 2025-06-30 08:54:47
 * @Email: wang.rui@nyu.edu
 * @FilePath: /7_MutParser/docs/PRF_ANALYSIS_GUIDE.md
 * @Description: 
-->

# Programmed Ribosomal Frameshifting (PRF) Analysis Guide

## Overview

The PRF Analyzer module provides comprehensive tools for analyzing Programmed Ribosomal Frameshifting in SARS-CoV-2. This mechanism is critical for viral replication as it allows the virus to produce two different proteins (ORF1a and ORF1ab) from the same mRNA.

## Biological Background

### What is Programmed Ribosomal Frameshifting?

Programmed Ribosomal Frameshifting (PRF) is a translational recoding mechanism where the ribosome shifts reading frames during translation, allowing the synthesis of multiple proteins from a single mRNA.

### SARS-CoV-2 PRF Site

- **Location**: Nucleotides 13,468-13,474
- **Slippery Sequence**: `UUUAAAC`
- **Frameshift Type**: -1 frameshift
- **Efficiency**: ~20-30%

### Protein Products

- **ORF1a**: 4,402 amino acids (nucleotides 266-13,468)
- **ORF1ab**: 7,096 amino acids (nucleotides 266-21,555)

## Installation and Setup

```bash
# Ensure you have the required dependencies
pip install biopython numpy

# The PRF analyzer is part of the mutparser package
# No additional installation required
```

## Basic Usage

### 1. Import the Module

```python
from src.utils.prf_analyzer import PRFAnalyzer, analyze_prf_in_sequence
```

### 2. Analyze a Genome

```python
# Create analyzer with genome sequence
analyzer = PRFAnalyzer(genome_sequence)

# Detect PRF site
prf_site = analyzer.detect_prf_site()

# Get results
print(f"Position: {prf_site.position}")
print(f"Sequence: {prf_site.sequence}")
print(f"Type: {prf_site.site_type.value}")
print(f"Efficiency: {prf_site.efficiency_score}")
```

### 3. Analyze Mutations

```python
# Define mutations
mutations = [
    {"pos": "13468", "type": "pointMutation", "SNP": "U->A"},
    {"pos": "13470:13472", "type": "deletion", "SNP": "UAA"}
]

# Analyze PRF mutations
prf_mutations = analyzer.analyze_prf_mutations(mutations)

# Check results
for mutation in prf_mutations:
    print(f"Affects PRF: {mutation['affects_prf']}")
    print(f"Impact: {mutation['prf_impact']['severity']}")
```

### 4. Predict Protein Ratios

```python
# Predict protein ratios based on PRF efficiency
ratios = analyzer.predict_protein_ratios(0.25)  # 25% efficiency

print(f"ORF1a ratio: {ratios['orf1a_ratio']}")
print(f"ORF1ab ratio: {ratios['orf1ab_ratio']}")
print(f"RdRp availability: {ratios['rdrp_availability']}")
```

### 5. One-Line Analysis

```python
report = analyze_prf_in_sequence(genome_sequence, mutations)
```

## Advanced Features

### 1. Comprehensive Analysis

```python
# Generate complete PRF report
report = analyzer.generate_prf_report(prf_site, mutations)

# Report includes:
# - PRF site details
# - Protein ratios
# - Mutation analysis
# - Summary statistics
```

### 2. Convenience Function

```python
# One-line analysis
report = analyze_prf_in_sequence(genome_sequence, mutations)
```

### 3. Custom Analysis

```python
# Analyze specific positions
prf_site = analyzer.detect_prf_site(position=13468)

# Check slippery site variants
is_variant = analyzer._is_slippery_site_variant("UUUAAAA")

# Score downstream structure
structure_score = analyzer._score_downstream_structure(downstream_sequence)
```

## API Reference

### PRFAnalyzer Class

#### Constructor
```python
PRFAnalyzer(ref_genome: str)
```
- `ref_genome`: Reference genome sequence

#### Methods

##### `detect_prf_site(position: Optional[int] = None) -> PRFSite`
Detect PRF site at specified position (default: 13468).

##### `analyze_prf_mutations(mutations: List[Dict]) -> List[Dict]`
Analyze mutations for PRF impact.

##### `predict_protein_ratios(prf_efficiency: float) -> Dict[str, float]`
Predict protein ratios based on PRF efficiency.

##### `generate_prf_report(prf_site: PRFSite, mutations: List[Dict] = None) -> Dict`
Generate comprehensive PRF analysis report.

### PRFSite Dataclass

```python
@dataclass
class PRFSite:
    position: int
    sequence: str
    site_type: PRFSiteType
    efficiency_score: float
    downstream_structure: Optional[str] = None
    mutations: List[Dict] = None
```

### PRFSiteType Enum

```python
class PRFSiteType(Enum):
    CANONICAL = "canonical"
    VARIANT = "variant"
    ABSENT = "absent"
```

## Example Workflows

### 1. Reference Genome Analysis

```python
from src.utils.sequence_utils import read_fasta
from src.utils.prf_analyzer import PRFAnalyzer

# Load reference genome
records = read_fasta("./data/refs/NC_045512.fasta")
ref_genome = records[0]["seq"]

# Analyze PRF
analyzer = PRFAnalyzer(ref_genome)
prf_site = analyzer.detect_prf_site()

print(f"PRF Site: {prf_site.sequence}")
print(f"Efficiency: {prf_site.efficiency_score:.3f}")
```

### 2. Mutation Impact Analysis

```python
# Load mutations from your analysis
mutations = load_mutations_from_file("mutations.json")

# Analyze PRF impact
prf_mutations = analyzer.analyze_prf_mutations(mutations)

# Filter high-impact mutations
high_impact = [m for m in prf_mutations if m['prf_impact']['severity'] == 'high']

print(f"High-impact PRF mutations: {len(high_impact)}")
```

### 3. Variant Comparison

```python
# Compare PRF efficiency across variants
variants = {
    "Reference": ref_genome,
    "Alpha": alpha_genome,
    "Beta": beta_genome
}

for name, genome in variants.items():
    analyzer = PRFAnalyzer(genome)
    prf_site = analyzer.detect_prf_site()
    ratios = analyzer.predict_protein_ratios(prf_site.efficiency_score)
    
    print(f"{name}: ORF1ab ratio = {ratios['orf1ab_ratio']:.3f}")
```

## Biological Interpretation

### Efficiency Categories

- **High (≥0.8)**: Optimal frameshifting, normal viral replication
- **Moderate (0.5-0.8)**: Reduced frameshifting, may affect replication
- **Low (0.2-0.5)**: Significantly reduced frameshifting
- **Very Low (<0.2)**: Minimal frameshifting, likely defective

### Mutation Impact

- **High Severity**: Deletions/insertions in slippery site
- **Moderate Severity**: Point mutations in slippery site
- **Low Severity**: Mutations in downstream region

### Protein Ratio Implications

- **ORF1ab ratio < 0.2**: Insufficient RdRp, defective replication
- **ORF1ab ratio 0.2-0.4**: Reduced replication efficiency
- **ORF1ab ratio > 0.4**: Normal replication capacity

## Troubleshooting

### Common Issues

1. **Genome too short**: Ensure genome is at least 15,000 nucleotides
2. **Invalid nucleotides**: Check for non-ATCGU characters
3. **PRF site not found**: Verify position 13468 contains expected sequence

### Error Messages

- `"Genome sequence too short"`: Increase genome length
- `"Invalid nucleotides found"`: Clean genome sequence
- `"Position out of range"`: Check mutation positions

## Integration with Existing Code

The PRF analyzer integrates seamlessly with your existing mutation analysis pipeline:

```python
# In your main analysis script
from src.utils.prf_analyzer import PRFAnalyzer

# After processing mutations
prf_analyzer = PRFAnalyzer(ref_genome)
prf_report = prf_analyzer.generate_prf_report(prf_site, processed_mutations)

# Add PRF analysis to your results
results["prf_analysis"] = prf_report
```

## Future Enhancements

Planned features for future versions:

1. **RNA Structure Prediction**: Integration with RNAfold
2. **Machine Learning Models**: Improved efficiency prediction
3. **Comparative Analysis**: Cross-species PRF comparison
4. **Visualization Tools**: PRF site visualization
5. **Database Integration**: PRF mutation database

## References

1. Kelly, J. A., et al. (2021). Structural and functional conservation of the programmed -1 ribosomal frameshift signal of SARS coronavirus 2 (SARS-CoV-2). *Journal of Biological Chemistry*, 296, 100410.

2. Plant, E. P., et al. (2005). A three-stemmed mRNA pseudoknot in the SARS coronavirus frameshift signal. *PLoS Biology*, 3(6), e172.

3. Atkins, J. F., et al. (2016). Recoding: Expansion of decoding rules enriches gene expression. *Nucleic Acids Research*, 44(15), 7227-7235.
