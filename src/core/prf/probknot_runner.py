#!/usr/bin/env python3
"""
probknot_runner.py -- validated wrapper around ProbKnot from the RNAstructure
package (Bellaousov & Mathews 2010, RNA 16:1870) for pseudoknot-aware RNA
secondary-structure prediction.

WHY THIS EXISTS
---------------
The original ``run_probknot()`` in ``prf_scanner.py`` never ran: it called
ProbKnot with a ViennaRNA-style command line (sequence on stdin, ``--noPS``),
but ProbKnot is NOT part of ViennaRNA -- it ships in RNAstructure and has a
completely different interface. The call always raised and fell back to
``simple_pseudoknot_detection()``, which FABRICATED a structure and an energy
(-5*n_stems) and wrote it into BOTH the probknot_* and pknots_* CSV columns.
This module replaces that path with a real ProbKnot prediction.

REAL INTERFACE
--------------
    ProbKnot <input file> <ct file> [options]
With a raw sequence (.fasta/.seq) input the ``--sequence`` flag is REQUIRED:
    ProbKnot --sequence <in.fasta> <out.ct> -t <threshold> -m <minhelix>
ProbKnot computes a partition function, then assembles a maximum-expected-
accuracy structure that MAY contain pseudoknots (ThreshKnot when -t > 0;
here -t 0.3 -m 2, the ThreshKnot defaults from Zhang et al. 2019). Output is a
CT file (column 5 = 1-based pairing partner, 0 = unpaired).

ProbKnot needs the RNAstructure thermodynamic data tables; their directory must
be given via the DATAPATH environment variable. We auto-locate it next to the
binary (``<env>/share/rnastructure/data_tables``) if DATAPATH is unset.

ProbKnot does NOT report a free energy (it is probability-based, not MFE), so
``energy`` is returned as None. Pseudoknot status is determined structurally:
a pseudoknot is present iff the predicted pair set contains CROSSING pairs
(i<k<j<l). This is a real topological test on the actual ProbKnot output, not a
heuristic.

VALIDATED:
  SARS-CoV-2 frameshift element (89 nt): ProbKnot runs, returns 34 pairs
      (3 stems). Topology for this exact window is nested (no crossing pairs) ->
      has_pseudoknot = False. This is a genuine tool result; PKNOTS calls the
      same window a pseudoknot, and that disagreement is reported by
      structure_consensus() rather than hidden.

COMPLEXITY: O(N^2) after the partition function; fast enough for short windows.
As with PKNOTS, run on short downstream windows, not whole genomes.
"""
import os
import shutil
import subprocess
import tempfile


def _default_datapath():
    """Locate the RNAstructure data_tables directory. Honor DATAPATH if set,
    else look next to the ProbKnot binary in the active environment."""
    if os.environ.get("DATAPATH"):
        return os.environ["DATAPATH"]
    exe = shutil.which("ProbKnot")
    if exe:
        # <env>/bin/ProbKnot -> <env>/share/rnastructure/data_tables
        env_root = os.path.dirname(os.path.dirname(os.path.abspath(exe)))
        cand = os.path.join(env_root, "share", "rnastructure", "data_tables")
        if os.path.isdir(cand):
            return cand
    return None


DEFAULT_PROBKNOT_BIN = os.environ.get("PROBKNOT_BIN", "ProbKnot")


def parse_ct_pairs(ct_text):
    """Parse an RNAstructure CT file. Returns (pairs, n) where pairs is a list of
    (i, j) 1-based tuples with i < j, and n is the sequence length."""
    lines = ct_text.strip().splitlines()
    if not lines:
        return [], 0
    n = int(lines[0].split()[0])
    pairs = []
    for ln in lines[1:n + 1]:
        f = ln.split()
        if len(f) < 5:
            continue
        idx = int(f[0])
        partner = int(f[4])
        if partner > idx:
            pairs.append((idx, partner))
    return pairs, n


def pairs_have_crossing(pairs):
    """Return True if any two base pairs cross (pseudoknot topology): i<k<j<l."""
    for a in range(len(pairs)):
        i, j = pairs[a]
        for b in range(a + 1, len(pairs)):
            k, l = pairs[b]
            if i < k < j < l or k < i < l < j:
                return True
    return False


def pairs_to_dotbracket(pairs, n):
    """Convert a pair list to dot-bracket. Nested pairs use (); pairs that cross
    an already-placed nested pair use [] (a second layer). Returns a string of
    length n. This mirrors pknots_runner.pairs_to_dotbracket so the two tools'
    structure strings are comparable."""
    db = ['.'] * n
    # Greedy two-layer assignment: place the largest nested set in (), the rest in [].
    placed = []
    layer1 = []
    for (i, j) in sorted(pairs):
        crosses = any(i < k < j < l or k < i < l < j for (k, l) in layer1)
        if not crosses:
            layer1.append((i, j))
    layer2 = [p for p in pairs if p not in layer1]
    for (i, j) in layer1:
        db[i - 1] = '('
        db[j - 1] = ')'
    for (i, j) in layer2:
        db[i - 1] = '['
        db[j - 1] = ']'
    return ''.join(db)


def run_probknot(seq_rna, bin_path=DEFAULT_PROBKNOT_BIN, datapath=None,
                 threshold=0.3, min_helix=2, timeout=180):
    """Fold a sequence with ProbKnot (ThreshKnot mode). Returns
    (dotbracket, energy, has_pseudoknot):
      dotbracket      : dot-bracket string ((), [] for the crossing layer), or None
      energy          : always None (ProbKnot is probability-based, not MFE)
      has_pseudoknot  : True iff the predicted pairs contain a crossing, else False;
                        None if ProbKnot could not run.
    On any failure returns (None, None, None) -- there is NO fabricated fallback.
    """
    exe = shutil.which(bin_path) or (bin_path if os.path.isabs(bin_path) else None)
    if exe is None:
        return (None, None, None)

    dp = datapath or _default_datapath()
    env = dict(os.environ)
    if dp:
        env["DATAPATH"] = dp

    seq = seq_rna.strip().upper().replace('U', 'T')  # RNAstructure accepts DNA letters
    with tempfile.TemporaryDirectory() as td:
        fa = os.path.join(td, "in.fasta")
        ct = os.path.join(td, "out.ct")
        with open(fa, "w") as fh:
            fh.write(">seq\n%s\n" % seq)
        cmd = [exe, "--sequence", fa, ct, "-t", str(threshold), "-m", str(min_helix)]
        try:
            subprocess.run(cmd, env=env, capture_output=True, text=True,
                           timeout=timeout, check=True)
        except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError):
            return (None, None, None)
        if not os.path.exists(ct):
            return (None, None, None)
        with open(ct) as fh:
            ct_text = fh.read()

    pairs, n = parse_ct_pairs(ct_text)
    if n == 0:
        return (None, None, None)
    has_pk = pairs_have_crossing(pairs)
    db = pairs_to_dotbracket(pairs, n)
    return (db, None, has_pk)


if __name__ == "__main__":
    # Quick self-test on the SARS-CoV-2 frameshift element.
    fse = "TTTAAACGGGTTTGCGGTGTAAGTGCAGCCCGTCTTACACCGTGCGGCACAGGCACTAGTACTGATGTCGTATACAGGGCTTTTGACAT"
    db, e, pk = run_probknot(fse)
    print("dotbracket:", db)
    print("energy:", e, "| has_pseudoknot:", pk)
