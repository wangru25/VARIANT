# User Guide: Existing-Dual-Search Pipeline

This guide walks through installing dependencies, preparing inputs, and running the three-step pipeline that converts RNA dot-bracket notation into dual graph IDs.

---

## Overview

The pipeline takes an RNA secondary structure in dot-bracket format and assigns a dual graph ID to each structural domain. It runs in three steps:

1. **`convert_dbn_to_dssr.py`** — converts a `.dbn` file to X3DNA-DSSR format.
2. **`PDBto2D.py`** — parses DSSR output into 2D structure summaries and CT files.
3. **`Dual_Library.py`** — assigns dual graph IDs from CT files.

---

## Prerequisites

### 1. Conda

Install [Miniconda](https://docs.conda.io/en/latest/miniconda.html) or [Anaconda](https://www.anaconda.com/) if you do not already have it.

### 2. RNAstructure (provides `dot2ct` and `data_tables`)

RNAstructure is available through Bioconda:

```bash
conda install -c bioconda rnastructure
```

This installs the `dot2ct` command-line tool and the required `data_tables` directory into your conda environment.

### 3. Python packages

The scripts in this folder require:

- `python-igraph`
- `matplotlib`
- `numpy`

Install them with:

```bash
conda install -c conda-forge python-igraph matplotlib numpy
```

Or, if you are using the project-level environment file from the repository root:

```bash
conda env create -f ../environment.yaml
conda activate variant
conda install -c bioconda rnastructure        # add dot2ct
conda install -c conda-forge python-igraph    # add igraph
```

---

## Setting `DATAPATH`

`dot2ct` requires the `DATAPATH` environment variable to point to RNAstructure's `data_tables` directory. Run this once per terminal session before executing the pipeline:

```bash
export DATAPATH=$(find $(conda info --base) -name "data_tables" -type d | head -1)
```

To make this permanent, add the line to your shell profile (`~/.bashrc`, `~/.zshrc`, etc.):

```bash
echo 'export DATAPATH=$(find $(conda info --base) -name "data_tables" -type d | head -1)' >> ~/.zshrc
source ~/.zshrc
```

Verify the setup:

```bash
echo $DATAPATH   # should print a path ending in /data_tables
dot2ct --help    # should print usage information
```

---

## Preparing Input

### 1. Create a `.dbn` file

Place your RNA structure file in `PDB_DBN/` using the naming convention `<ID>-2D.dbn`.

A `.dbn` file must contain two lines: the RNA sequence and the dot-bracket notation:

```text
>my_structure
AUGCAUGCAUGC
(((...)))...
```

Dot-bracket symbols:
- `(` — opening base pair
- `)` — closing base pair
- `.` — unpaired nucleotide
- `&` — chain break (for multi-chain structures)

**If you do not already have a dot-bracket structure**, generate one from a PDB file using either of these tools:
- [X3DNA-DSSR web interface](https://x3dna.org/)
- [RNApdbee](http://rnapdbee.cs.put.poznan.pl/)

Extract the sequence and dot-bracket lines from the output and save them into a `.dbn` file.

### 2. Register the ID in `list_file.txt`

Open `list_file.txt` and add your structure ID as a comma-separated entry:

```text
7LYJ, MY_ID,
```

Each ID must match the base name of your `.dbn` file (e.g., `MY_ID` for `PDB_DBN/MY_ID-2D.dbn`).

---

## Running the Pipeline

All commands must be run from inside the `Existing-Dual-Search/` directory.

```bash
cd Existing-Dual-Search

# Step 0: set DATAPATH (required each session unless added to shell profile)
export DATAPATH=$(find $(conda info --base) -name "data_tables" -type d | head -1)

# Step 1: convert .dbn to DSSR format  (repeat for each new ID)
python3 convert_dbn_to_dssr.py PDB_DBN/<ID>-2D.dbn PDB_DSSR/<ID>.out <ID>

# Step 2: parse DSSR output into 2D summaries and CT files (processes all IDs in list_file.txt)
python3 PDBto2D.py

# Step 3: assign dual graph IDs (processes all IDs in list_file.txt)
python3 Dual_Library.py
```

Replace `<ID>` with your structure identifier (e.g., `7LYJ`).

---

## Outputs

After a successful run the following files are produced:

| Path | Description |
|---|---|
| `PDB_DSSR/<ID>.out` | DSSR-format structure output |
| `PDB_DSSR_2D/<ID>.txt` | 2D structure summary |
| `PDB_DSSR_CT/<ID>.ct` | CT file for the full structure |
| `PDB_DSSR_CT/<ID>_<chains>.ct` | CT files for each substructure |
| `PDB_DSSR_Dual/<ID>_<chains>.txt` | Dual graph assignment(s) |

---

## Troubleshooting

| Symptom | Likely cause | Fix |
|---|---|---|
| `dot2ct: command not found` | RNAstructure not installed or not activated | `conda install -c bioconda rnastructure` and activate the environment |
| `dot2ct` runs but produces empty output | `DATAPATH` not set | Run `export DATAPATH=$(find $(conda info --base) -name "data_tables" -type d | head -1)` |
| `ModuleNotFoundError: igraph` | `python-igraph` not installed | `conda install -c conda-forge python-igraph` |
| `ERROR: DSSR file not found` | Step 1 was skipped or the ID in `list_file.txt` does not match the filename | Run Step 1 first; check spelling in `list_file.txt` |
| `No base pair!` | The dot-bracket string contains no paired nucleotides | Verify the `.dbn` file content |
