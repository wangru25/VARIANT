# Existing-Dual-Search

This folder contains a runnable workflow for processing RNA secondary structure data and assigning dual graph IDs.

## Workflow Overview

The pipeline runs in three steps:

1. Convert dot-bracket input into DSSR-style output.
2. Convert DSSR output into 2D structure and CT files.
3. Assign dual graph IDs from CT files.

## Required Input

- A DBN file with sequence and dot-bracket notation.
- One or more IDs in `list_file.txt` for downstream steps.

Dot-bracket symbols:
- `(` opening base pair
- `)` closing base pair
- `.` unpaired nucleotide
- `&` chain break

## How To Get Dot-Bracket Notation

If you do not already have sequence + dot-bracket data:

1. Use X3DNA-DSSR web interface: [https://x3dna.org/](https://x3dna.org/)
2. Or use RNApdbee: [http://rnapdbee.cs.put.poznan.pl/](http://rnapdbee.cs.put.poznan.pl/)
3. Extract sequence and dot-bracket and place them into a `.dbn` file.

Example `.dbn`:

```text
>my_structure
AUGCAUGCAUGC
(((...)))...
```

## Run The Pipeline

From `Existing-Dual-Search`:

```bash
# Set DATAPATH for RNAstructure
export DATAPATH=$(find $(conda info --base) -name "data_tables" -type d | head -1)

python3 convert_dbn_to_dssr.py PDB_DBN/<ID>-2D.dbn PDB_DSSR/<ID>.out <ID>
python3 PDBto2D.py
python3 Dual_Library.py
```

`PDBto2D.py` and `Dual_Library.py` process all IDs from `list_file.txt`.

## Outputs

- `PDB_DSSR/{id}.out` DSSR-style structure output
- `PDB_DSSR_2D/{id}.txt` 2D structure summary
- `PDB_DSSR_CT/{id}.ct` CT file for whole structure
- `PDB_DSSR_CT/{id}_{chains}.ct` CT files for substructures
- `PDB_DSSR_Dual/{id}_{chains}.txt` dual graph assignments

## PDBto2D Complete Version Notes

The included complete version restores key logic from the original implementation:

- NT array parsing and relabel mapping
- Subchain handling
- Chain interaction detection
- Dimer detection and weak interaction filtering
- Base pair cleanup in dot-bracket notation
- More robust DSSR parsing

These improve chain mapping accuracy and data integrity compared with simplified variants.

## Dependencies

- `DATAPATH` should point to RNAstructure `data_tables`
- `dot2ct` available in `PATH`
- Python dependencies used by scripts in this folder (`python-igraph`, `matplotlib`, `numpy`)
