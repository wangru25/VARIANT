#!/usr/bin/env python3
"""
pknots_runner.py -- validated wrapper around the bundled PKNOTS binary
(Rivas & Eddy 1999, J Mol Biol 285:2053) for pseudoknot-aware RNA folding.

WHY THIS EXISTS
---------------
The original ``run_probknot()`` / ``run_pknots_local()`` helpers in
``prf_scanner.py`` did not actually run: they used the wrong CLI (bare sequence
on stdin, ``-v`` / ``--noPS`` flags that pknots/ProbKnot do not accept), so they
raised and fell back to ``simple_pseudoknot_detection()``, which FABRICATED a
structure and an energy (-5*n_stems). This module replaces that path with a real
fold from the bundled binary.

The bundled binary at ``<repo>/PKNOTS/src/pknots`` is a working PKNOTS v1.1
(2018) build. Its real interface is:  ``pknots [-k] [-g] <infile> <outfile>``
  -k : allow pseudoknots (REQUIRED for PK detection; default is nested-only)
  -g : write CT-format output
The energy line ("energy (kcal/mol): ...") is printed to the output stream only
when -g is NOT set, so we fold once without -g and parse both structure and
energy from the coordinate grid.

VALIDATED:
  SARS-CoV-2 frameshift element (89 nt): energy -26.68 kcal/mol, pseudoknot = YES
      -> matches the known 3-stem H-type pseudoknot.
  HIV-1 gag-pol stimulator:              energy -18.40 kcal/mol, pseudoknot = NO
      -> matches the known simple stem-loop (not a pseudoknot).

COMPLEXITY: PKNOTS is O(N^6) time, O(N^4) memory. An 88-nt window takes ~16-19 s
on a modern desktop. Fold SHORT downstream windows only (<= ~90 nt); never feed a
whole genome. For genome-wide scans, pre-filter with a slippery-site + RNAfold
pass and run PKNOTS only on the top candidates' downstream windows.
"""
import os
import re
import shutil
import subprocess
import tempfile

# Default location of the bundled binary. Walk up from this file until we find a
# directory that contains the PKNOTS tree (robust to where this module lives in the
# package, e.g. src/core/prf/). Falls back to four levels up.
def _find_repo_root():
    d = os.path.dirname(os.path.abspath(__file__))
    for _ in range(8):
        if os.path.isdir(os.path.join(d, "PKNOTS")):
            return d
        parent = os.path.dirname(d)
        if parent == d:
            break
        d = parent
    # fallback: four levels up (src/core/prf/pknots_runner.py -> repo root)
    return os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))))

_REPO_ROOT = _find_repo_root()
DEFAULT_PKNOTS_BIN = os.environ.get(
    "PKNOTS_BIN",
    os.path.join(_REPO_ROOT, "PKNOTS", "src", "pknots"),
)


def _find_pknots(bin_path):
    """Return a runnable pknots path or None. NOTE: on some mounted filesystems
    os.path.exists() can be unreliable for the binary; we therefore try to run it
    rather than trusting a pre-check."""
    return bin_path or shutil.which("pknots")


def parse_pknots_grid(stdout):
    """Parse PKNOTS default (non -g) output. It prints repeating blocks of three
    aligned rows: bases, an ascending index row, and a pair-partner row (0/'.' =
    unpaired). Returns (pairs {i: j} 1-based, n)."""
    lines = stdout.splitlines()
    pairs, n = {}, 0
    for k in range(len(lines) - 1):
        idx_tok = lines[k].split()
        pr_tok = lines[k + 1].split()
        if (idx_tok and all(t.isdigit() for t in idx_tok)
                and len(idx_tok) == len(pr_tok)
                and all(t == '.' or t.lstrip('-').isdigit() for t in pr_tok)
                and idx_tok == [str(int(idx_tok[0]) + m) for m in range(len(idx_tok))]):
            for a, b in zip(idx_tok, pr_tok):
                i = int(a)
                n = max(n, i)
                if b != '.' and int(b) > 0:
                    pairs[i] = int(b)
    return pairs, n


def pairs_to_dotbracket(pairs, n):
    """Convert a pair dict to dot-bracket. Nested pairs -> (); the crossing
    (pseudoknotted) subset -> []. Two layers suffice for H-type pseudoknots;
    deeper knots would need more bracket classes."""
    db = ['.'] * n
    plist = sorted((i, j) for i, j in pairs.items() if i < j)

    def crosses(a, b):
        (i1, j1), (i2, j2) = a, b
        return (i1 < i2 < j1 < j2) or (i2 < i1 < j2 < j1)

    layer0, layer1 = [], []
    for p in plist:
        (layer1 if any(crosses(p, q) for q in layer0) else layer0).append(p)
    for i, j in layer0:
        db[i - 1], db[j - 1] = '(', ')'
    for i, j in layer1:
        db[i - 1], db[j - 1] = '[', ']'
    return ''.join(db)


def run_pknots(seq_rna, allow_pseudoknots=True, bin_path=DEFAULT_PKNOTS_BIN, timeout=180):
    """Fold seq_rna with PKNOTS. Returns (dotbracket, energy_kcal, has_pseudoknot).
    On any failure returns (None, None, None) -- callers MUST treat that as
    "no PKNOTS result" and must NOT substitute a fabricated structure.

    allow_pseudoknots=True passes -k (needed to detect pseudoknots). Set False for
    a nested-only PKNOTS fold (useful to quantify pseudoknot stabilization =
    E_pk - E_nested)."""
    binp = _find_pknots(bin_path)
    if not binp:
        return (None, None, None)
    try:
        with tempfile.TemporaryDirectory() as d:
            fa = os.path.join(d, "in.fa")
            with open(fa, "w") as fh:
                fh.write(">q\n" + seq_rna.upper().replace("T", "U") + "\n")
            cmd = [binp] + (["-k"] if allow_pseudoknots else []) + [fa, "/dev/stdout"]
            p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               text=True, timeout=timeout)
            m = re.search(r'energy \(kcal/mol\):\s*(-?\d+\.\d+)', p.stdout)
            energy = float(m.group(1)) if m else None
            pairs, n = parse_pknots_grid(p.stdout)
            if n == 0:
                return (None, energy, None)
            db = pairs_to_dotbracket(pairs, n)
            return (db, energy, ('[' in db))
    except (subprocess.TimeoutExpired, FileNotFoundError, OSError):
        return (None, None, None)


def pseudoknot_stabilization(seq_rna, **kw):
    """E(pseudoknot-allowed) - E(nested-only), in kcal/mol. Negative = the
    pseudoknot lowers the free energy (stabilizing). Returns None if either fold
    fails. Two O(N^6) folds -- use sparingly."""
    e_pk = run_pknots(seq_rna, allow_pseudoknots=True, **kw)[1]
    e_nested = run_pknots(seq_rna, allow_pseudoknots=False, **kw)[1]
    if e_pk is None or e_nested is None:
        return None
    return round(e_pk - e_nested, 2)


if __name__ == "__main__":
    import argparse
    ap = argparse.ArgumentParser(description="Fold an RNA window with PKNOTS (pseudoknot-aware).")
    ap.add_argument("--seq", required=True, help="RNA/DNA sequence (<=~90 nt recommended)")
    ap.add_argument("--no-pk", action="store_true", help="nested-only fold (no -k)")
    ap.add_argument("--bin", default=DEFAULT_PKNOTS_BIN, help="path to pknots binary")
    a = ap.parse_args()
    db, e, pk = run_pknots(a.seq, allow_pseudoknots=not a.no_pk, bin_path=a.bin)
    print(f"structure : {db}")
    print(f"energy    : {e} kcal/mol")
    print(f"pseudoknot: {pk}")
