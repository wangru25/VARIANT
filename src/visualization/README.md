# VARIANT Visualization Module

This module generates publication-ready figures from VARIANT analysis results. It is
driven from the repository root by `plot.py`, which routes to the implementations here.

## Files

- `figure1_mutation_analysis.py` — combined genome + protein mutation analysis (`SimpleCombinedAnalyzer`).
- `figure2_row_hot_mutations.py` — row and hot-mutation views (`RowHotMutationVisualizer`).
- `figure3_PRF.py` — potential PRF-region visualization from `*.prf_candidates.*` result files (`PRFVisualizer`).
- `__init__.py` — module initialization and exports.

## Usage

### Command line (from the repository root)

```bash
python plot.py --list-viruses
python plot.py --type mutation --virus SARS-CoV-2 --genome-id EPI_ISL_16127650
python plot.py --type row-hot  --virus SARS-CoV-2
python plot.py --type prf      --virus SARS-CoV-2

# A figure script can also be run directly:
python src/visualization/figure3_PRF.py --help
```

`--type` selects the figure: `mutation`, `row-hot`, or `prf`.

### Python API

```python
from src.visualization import SimpleCombinedAnalyzer
```

## Output

Figures are written as interactive **HTML** and static **PDF** under
`imgs/visualizations/<Virus>/` (unless `--output` overrides the path), suitable for
scientific publications, web applications, and presentations.
