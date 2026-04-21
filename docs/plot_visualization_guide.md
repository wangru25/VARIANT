# VARIANT visualization guide

## Where the code lives

- **CLI entry point:** `plot.py` at the repository root routes to implementations in `src/visualization/`.
- **Modules:** `figure1_mutation_analysis.py`, `figure2_row_hot_mutations.py`, and `figure3_PRF.py` contain the current Plotly workflows used by the web app and command line.

## Requirements

Install Plotly (and other project dependencies) as described in the root `README.md` / `environment.yaml`.

## Command-line usage

From the repository root:

```bash
python plot.py --list-viruses
python plot.py --type mutation --virus SARS-CoV-2 --genome-id EPI_ISL_123456
python plot.py --type row-hot --virus SARS-CoV-2
python plot.py --type prf --virus SARS-CoV-2
```

Optional flags: `--genome-id` to pin a sample, `--output` to set an explicit output path.

You can also run the underlying modules directly (see each file’s `argparse` / `main()` for flags), for example:

```bash
python src/visualization/figure3_PRF.py --help
```

## Outputs

By default, HTML and PDF figures are written under `imgs/visualizations/<Virus>/` unless you pass `--output`.

## Further reading

- [COMPREHENSIVE_VISUALIZATION_SUMMARY.md](COMPREHENSIVE_VISUALIZATION_SUMMARY.md) — high-level overview of visualization types and data expectations.
- [PROTEIN_VISUALIZATION_SUMMARY.md](PROTEIN_VISUALIZATION_SUMMARY.md) — protein-focused notes where applicable.
