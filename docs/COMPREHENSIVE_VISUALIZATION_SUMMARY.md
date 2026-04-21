# Visualization overview (VARIANT)

This document describes the **current** visualization stack in the VARIANT repository. Implementation code lives under `src/visualization/`; the root script `plot.py` is the supported command-line entry point.

## Architecture

| Role | Location |
|------|----------|
| Unified CLI | `plot.py` |
| Combined genome + protein mutation analysis | `src/visualization/figure1_mutation_analysis.py` |
| Row and hot-mutation plots | `src/visualization/figure2_row_hot_mutations.py` |
| PRF candidate region plots | `src/visualization/figure3_PRF.py` |

The web application invokes the same modules via `plot.py`-style routing; keep CLI and web paths aligned when changing arguments or defaults.

## Data expectations

Modules read analysis outputs from the `result/<Virus>/` tree (mutation summaries, row/hot mutation tables, PRF candidate files, etc.). Run the main VARIANT analysis pipeline first so those inputs exist; see the root `README.md` for the full workflow.

## Outputs

Unless `--output` is given, figures are saved as paired **HTML** (interactive) and **PDF** (static) under:

`imgs/visualizations/<Virus>/`

Naming follows each module’s conventions (often including genome or analysis identifiers in the basename).

## Quick examples

```bash
python plot.py --list-viruses
python plot.py --type mutation --virus SARS-CoV-2 --genome-id EPI_ISL_123456
python plot.py --type row-hot --virus SARS-CoV-2
python plot.py --type prf --virus SARS-CoV-2
```

For module-specific options, use `--help` on `plot.py` or on the individual script under `src/visualization/`.

## Related documentation

- [plot_visualization_guide.md](plot_visualization_guide.md) — concise plotting guide and links.
- [PROTEIN_VISUALIZATION_SUMMARY.md](PROTEIN_VISUALIZATION_SUMMARY.md) — protein-level visualization context.
