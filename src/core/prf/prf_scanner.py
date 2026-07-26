#!/usr/bin/env python3
"""
PRF Scanner: detect candidate programmed ribosomal frameshifting (-1) sites.

Pipeline:
1) Scan for -1 PRF slippery heptamers of the form X XXY YYZ (spacing marks the
   0-frame codons). Strict (default): X = any base (three identical), Y in {A, U}
   (three identical), Z in {A, C, U} (canonically not G). Relaxed (--relaxed):
   any Y, any Z. Canonical examples: U U U A A A C (SARS-CoV-2), U U U U U U A
   (HIV-1 gag-pol). The heptamer is matched directly (regex with backreferences),
   so novel sites are found, not only a hardcoded list.
2) For each site, fold the downstream window across spacer lengths (default 5..9 nt).
3) Structure evidence (--use-rnafold and/or --use-pknots):
     RNAfold (ViennaRNA)  -> nested minimum-free-energy fold. Does NOT model
                             pseudoknots; MFE is a proxy for local stability.
     PKNOTS (Rivas & Eddy) -> pseudoknot-aware fold via src/core/prf/pknots_runner.py
                             using the bundled binary. O(N^6): run once per site on
                             the most stable spacer window (<= --pknots-window nt).
   A transparent consensus is reported (structure_call / structure_confidence):
   a pseudoknot is called only when PKNOTS finds crossing base pairs; tool
   agreement raises confidence and disagreement is reported, not hidden. There is
   NO fabricated-structure fallback: if a tool does not run, its fields are blank.
4) Report candidates as CSV and BED.

Notes:
- tRNA columns are a PROXY for A-/P-site decoding availability derived from an
  optional user-supplied codon-abundance CSV. Two non-redundant quantities are
  reported: trna_score = sqrt(a1*a2) (geometric mean, tAI-consistent overall tRNA
  supply) and pausing_potential = 1 - min(a1,a2) (limiting/"hungry"-codon pausing).
  They are NOT measured tRNA concentrations, kinetics, or ribosome-profiling data.
- If you have Ribo-seq or conservation tracks, join those externally to prioritize candidates.

Usage:
    python src/core/prf/prf_scanner.py --fasta genome.fasta --out out_prefix \
        [--relaxed --spacer-min 5 --spacer-max 9 --window 120 \
         --use-rnafold --use-pknots --pknots-window 85 \
         --trna trna_abundance.csv]

"""

import argparse
import re
import sys
import os
from pathlib import Path
from typing import List, Tuple, Optional, Dict
import csv
import shutil
import subprocess
import json

# Real pseudoknot-aware folding via the bundled PKNOTS binary. Import must work in
# BOTH modes this file is used in: (a) imported as a module, and (b) run as a script
# via `python src/core/prf/prf_scanner.py` (how mutation_processor invokes it as a
# subprocess). In script mode sys.path[0] is the caller's cwd (the repo root), NOT
# this file's directory, so a bare `import pknots_runner` fails. Add this file's own
# directory to sys.path first so the sibling module always resolves.
_THIS_DIR = os.path.dirname(os.path.abspath(__file__))
if _THIS_DIR not in sys.path:
    sys.path.insert(0, _THIS_DIR)
try:
    from pknots_runner import run_pknots as _run_pknots
except ImportError:
    try:
        from .pknots_runner import run_pknots as _run_pknots
    except ImportError:
        _run_pknots = None

# ProbKnot (RNAstructure, Bellaousov & Mathews 2010) -- a second, independent
# pseudoknot-capable predictor. Same tolerant import as PKNOTS. Available only if
# the RNAstructure package (ProbKnot binary + data tables) is installed.
try:
    from probknot_runner import run_probknot as _run_probknot
except ImportError:
    try:
        from .probknot_runner import run_probknot as _run_probknot
    except ImportError:
        _run_probknot = None

# NUPACK (Zadeh et al. 2011) -- a third pseudoknot-capable predictor. Reports both a
# dot-parens-plus structure AND an MFE. Same tolerant import. Available only inside
# the user's docker image (rw3594/dual_rag_if, NUPACK 3.2.2 at /workspace/nupack3.2.2).
try:
    from nupack_runner import run_nupack as _run_nupack, nupack_available as _nupack_available
except ImportError:
    try:
        from .nupack_runner import run_nupack as _run_nupack, nupack_available as _nupack_available
    except ImportError:
        _run_nupack = None
        _nupack_available = None

DNA2RNA = str.maketrans({'T':'U','t':'u'})
RNA2DNA = str.maketrans({'U':'T','u':'t'})

def read_fasta(fp: Path) -> Dict[str,str]:
    seqs = {}
    name = None
    buf = []
    with fp.open() as f:
        for raw in f:
            line = raw.strip()
            if not line:
                continue
            if line.startswith('>'):
                # Flush previous
                if name is not None:
                    seqs[name] = ''.join(buf)
                hdr = line[1:].strip()
                if hdr:
                    # Take first token as ID, sanitize spaces
                    name = hdr.split()[0]
                else:
                    # Robust fallback for malformed headers like '>'
                    name = f"seq{len(seqs)+1}"
                buf = []
            else:
                # Sequence line
                buf.append(line.upper())
        if name is not None:
            seqs[name] = ''.join(buf)
    return seqs

def load_trna_data(trna_path: Optional[Path]) -> Optional[Dict]:
    """Load tRNA abundance data from file.
    
    Expected format: CSV with columns: codon,organism,abundance
    Example:
        codon,organism,abundance
        AAA,human,0.8
        UUU,human,0.9
        CCC,human,0.6
        GGG,human,0.7
        AAA,sars_cov2,0.3
        UUU,sars_cov2,0.4
    """
    if not trna_path or not trna_path.exists():
        return None
    
    trna_data = {}
    try:
        with trna_path.open() as f:
            # Skip blank lines and '#' comment lines so the CSV can carry a
            # human-readable header/provenance note without breaking parsing.
            data_lines = (ln for ln in f if ln.strip() and not ln.lstrip().startswith('#'))
            reader = csv.DictReader(data_lines)
            for row in reader:
                codon = row['codon']
                organism = row['organism']
                abundance = float(row['abundance'])
                
                if organism not in trna_data:
                    trna_data[organism] = {}
                trna_data[organism][codon] = abundance
                
        print(f"[INFO] Loaded tRNA data for {len(trna_data)} organisms")
        return trna_data
    except Exception as e:
        print(f"[WARN] Failed to load tRNA data: {e}", file=sys.stderr)
        return None

# -1 PRF slippery heptamer grammar  X XXY YYZ  (nucleotides X X X Y Y Y Z):
#   X : any base, three identical  -> capture group 1 with backreferences \1\1
#   Y : A or U,   three identical  -> strict; group 2 with \2\2
#   Z : A, C or U (canonically not G)
# Backreferences enforce the "three identical" requirement for ANY base, so novel
# sites are found without enumerating a hardcoded triplet list. This correctly
# matches X == Y sites (e.g. HIV-1 U UUU UUA) too, since \1 and \2 are independent.
SLIP_STRICT_RX  = re.compile(r"(?=(([ACGU])\2\2([AU])\3\3([ACU])))")   # Y in {A,U}, Z != G
SLIP_RELAXED_RX = re.compile(r"(?=(([ACGU])\2\2([ACGU])\3\3([ACGU])))") # any Y, any Z

def build_minus1_regex(relaxed: bool = False) -> re.Pattern:
    """Return the slippery-heptamer regex. Strict (default): Y in {A,U}, Z in
    {A,C,U}. Relaxed: any Y, any Z. Both require X, Y each be three identical
    bases via backreferences (the biologically defined X XXY YYZ pattern)."""
    return SLIP_RELAXED_RX if relaxed else SLIP_STRICT_RX

def find_minus1_slippery_sites(seq_rna: str, relaxed: bool = False) -> List[Tuple[int,str]]:
    """Return list of (0-based start, heptamer) for -1 PRF slippery sites."""
    pat = build_minus1_regex(relaxed)
    out = []
    for m in pat.finditer(seq_rna):
        out.append((m.start(), m.group(1)))
    return out

def run_rnafold(seq_rna: str) -> Tuple[Optional[str], Optional[float]]:
    """Call RNAfold if available on PATH. Return (structure, mfe_kcal)."""
    if shutil.which('RNAfold') is None:
        return None, None
    try:
        p = subprocess.run(
            ['RNAfold','--noPS'],
            input=(seq_rna+'\n').encode(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True
        )
        out = p.stdout.decode().strip().splitlines()
        # Expected format:
        # sequence
        # structure ( ... ) ( -xx.xx )
        # e.g.  ..((....)).. ( -12.30)
        if len(out) >= 2:
            struct_line = out[1]
            # The dot-bracket string is the leading run of . ( ) characters; the
            # energy is the last parenthesised signed float on the line. Parse each
            # explicitly rather than rsplit('(') (which breaks on '(' inside the
            # structure or on '( -0.00)' spacing variants).
            mstruct = re.match(r'^([.()]+)', struct_line)
            dotbr = mstruct.group(1) if mstruct else None
            menergy = re.search(r'\(\s*(-?\d+\.\d+)\s*\)\s*$', struct_line)
            mfe = float(menergy.group(1)) if menergy else None
            if dotbr is not None:
                return dotbr, mfe
        return None, None
    except Exception:
        return None, None

def run_pknots_local(seq_rna: str) -> Tuple[Optional[str], Optional[float], Optional[bool]]:
    """Fold seq_rna with the bundled PKNOTS binary (pseudoknot-aware, -k) via
    pknots_runner.run_pknots. Returns (dotbracket, energy_kcal, has_pseudoknot).

    Returns (None, None, None) if the wrapper/binary is unavailable or the fold
    fails. There is NO fabricated fallback: callers must treat None as "no result".

    NOTE: PKNOTS is O(N^6); fold short windows only (<= ~90 nt) and once per site.
    """
    if _run_pknots is None:
        return (None, None, None)
    return _run_pknots(seq_rna, allow_pseudoknots=True)


def run_probknot_local(seq_rna: str) -> Tuple[Optional[str], Optional[float], Optional[bool]]:
    """Fold seq_rna with ProbKnot (RNAstructure) via probknot_runner.run_probknot.
    Returns (dotbracket, energy, has_pseudoknot). ProbKnot is probability-based
    (ThreshKnot), so energy is always None; has_pseudoknot is a real topological
    test (crossing pairs) on the ProbKnot output.

    Returns (None, None, None) if the wrapper/binary is unavailable or the fold
    fails. There is NO fabricated fallback: callers must treat None as "no result".
    """
    if _run_probknot is None:
        return (None, None, None)
    return _run_probknot(seq_rna)


def run_nupack_local(seq_rna: str) -> Tuple[Optional[str], Optional[float], Optional[bool]]:
    """Fold seq_rna with NUPACK's pseudoknot-aware mfe via nupack_runner.run_nupack.
    Returns (dotbracket, energy_kcal, has_pseudoknot). NUPACK reports both a
    structure and an MFE; has_pseudoknot is a real crossing-pair test.

    Returns (None, None, None) if the wrapper/binary is unavailable or the fold
    fails. There is NO fabricated fallback: callers must treat None as "no result".
    NUPACK is only present in the user's docker image (rw3594/dual_rag_if).
    """
    if _run_nupack is None:
        return (None, None, None)
    return _run_nupack(seq_rna)


def structure_consensus(rnafold_struct, rnafold_mfe, pknots_struct, pknots_energy,
                        pknots_has_pk, mfe_flag=-15.0,
                        probknot_struct=None, probknot_has_pk=None,
                        nupack_struct=None, nupack_energy=None, nupack_has_pk=None):
    """Combine RNAfold (nested MFE), PKNOTS (pseudoknot-aware, MFE), ProbKnot
    (pseudoknot-aware, probability-based) and NUPACK (pseudoknot-aware, MFE)
    evidence into a transparent, non-fabricated call. Returns
    (structure_call, confidence, note).

    Rules (structure is weighted EVIDENCE, not a gate):
      * Pseudoknot voters = every pseudoknot-capable tool that actually ran
        (PKNOTS, ProbKnot, NUPACK), each voting via a real crossing-pair test.
        RNAfold cannot represent pseudoknots and never votes on PK presence.
      * A pseudoknot is CALLED when at least one PK-capable tool reports crossing
        pairs. When PK tools split, the MAJORITY sets the call and the split is
        reported in the note (ties -> call PK, confidence low).
      * Stability is judged from the best (lowest) available MFE among the tools
        that report an energy (RNAfold, PKNOTS, NUPACK) vs mfe_flag. ProbKnot is
        probability-based and reports no energy, so it contributes topology only.
      * Confidence:
          - 'high'   : >=2 PK-capable tools ran and ALL agree on PK presence; or
                       (no PK-tool ran) RNAfold+PKNOTS agree on stability.
          - 'medium' : exactly one PK-capable tool ran; or only RNAfold ran.
          - 'low'    : PK-capable tools disagree, or RNAfold/PKNOTS disagree on
                       stability. The disagreement is REPORTED, never hidden.
      * No tool result -> ('none', 'none', ...).
    """
    have_rnafold = rnafold_struct is not None and rnafold_mfe is not None
    have_pknots = pknots_struct is not None and pknots_energy is not None
    have_probknot = probknot_struct is not None
    have_nupack = nupack_struct is not None
    if not any((have_rnafold, have_pknots, have_probknot, have_nupack)):
        return ('none', 'none', 'no structure tool produced a result')

    # Stability from every tool that reports an energy.
    energies = [e for e in (rnafold_mfe, pknots_energy, nupack_energy) if e is not None]
    best_energy = min(energies) if energies else None
    stable = best_energy is not None and best_energy <= mfe_flag

    # --- pseudoknot voting across every PK-capable tool that ran ----------------
    pk_votes = []  # (tool_name, bool)
    if have_pknots:
        pk_votes.append(('PKNOTS', bool(pknots_has_pk)))
    if have_probknot:
        pk_votes.append(('ProbKnot', bool(probknot_has_pk)))
    if have_nupack:
        pk_votes.append(('NUPACK', bool(nupack_has_pk)))
    n_pk_tools = len(pk_votes)
    yes = [t for t, v in pk_votes if v]
    no = [t for t, v in pk_votes if not v]
    n_pk_yes = len(yes)

    # Majority of the PK-capable tools that ran sets the PK call; a tie (equal
    # yes/no) is resolved in favour of calling PK and flagged 'low' below. A lone
    # dissenting "yes" in a >=3-tool minority does NOT trigger a PK call (it is a
    # minority) but is still reported in the note.
    n_pk_no = n_pk_tools - n_pk_yes
    pk_called = n_pk_yes > 0 and n_pk_yes >= n_pk_no
    if pk_called:
        call = 'pseudoknot' if stable else 'weak_pseudoknot'
    elif stable:
        call = 'stem_loop'
    else:
        call = 'unstructured'

    # --- confidence + transparent note -----------------------------------------
    notes = []
    if n_pk_tools >= 2:
        if n_pk_yes == n_pk_tools:
            confidence = 'high'
            notes.append('%s agree: pseudoknot' % ', '.join(yes))
        elif n_pk_yes == 0:
            # all PK tools agree: no crossing. Confidence from stability agreement.
            if have_rnafold and have_pknots:
                confidence = 'high' if (rnafold_mfe <= mfe_flag) == (pknots_energy <= mfe_flag) else 'low'
            else:
                confidence = 'medium'
            notes.append('%s agree: no pseudoknot' % ', '.join([t for t, _ in pk_votes]))
        else:
            confidence = 'low'
            notes.append('pseudoknot disagreement: %s call PK, %s do not'
                         % (', '.join(yes), ', '.join(no)))
    elif n_pk_tools == 1:
        # Only one PK-capable tool ran (commonly PKNOTS). Its topology is trusted
        # for the call, but confidence still derives from RNAfold+PKNOTS stability
        # agreement so the historical RNAfold+PKNOTS-only case keeps 'high'.
        t, v = pk_votes[0]
        if have_rnafold and have_pknots and (rnafold_mfe <= mfe_flag) == (pknots_energy <= mfe_flag):
            confidence = 'high'
        else:
            confidence = 'medium'
        notes.append('single pseudoknot tool (%s): %s' % (t, 'PK' if v else 'no PK'))
    else:
        confidence = 'medium'
        notes.append('no pseudoknot-capable tool ran')

    # stability disagreement between RNAfold and PKNOTS is always worth reporting
    if have_rnafold and have_pknots and (rnafold_mfe <= mfe_flag) != (pknots_energy <= mfe_flag):
        confidence = 'low'
        notes.append('stability disagreement: RNAfold %.1f vs PKNOTS %.1f vs threshold %.1f'
                     % (rnafold_mfe, pknots_energy, mfe_flag))

    ran = [n for n, h in (('RNAfold', have_rnafold), ('PKNOTS', have_pknots),
                          ('ProbKnot', have_probknot), ('NUPACK', have_nupack)) if h]
    notes.append('tools run: ' + ', '.join(ran))
    return (call, confidence, '; '.join(notes))

def classify_structure_type(structure: str) -> str:
    """Classify RNA structure type based on dot-bracket notation."""
    if not structure:
        return 'unknown'
    
    # Check for pseudoknot indicators
    if '[' in structure and ']' in structure:
        return 'pseudoknot'
    elif '<' in structure and '>' in structure:
        return 'pseudoknot'
    elif '{' in structure and '}' in structure:
        return 'pseudoknot'
    
    # Check for stem-loop patterns
    open_parens = structure.count('(')
    close_parens = structure.count(')')
    dots = structure.count('.')
    
    if open_parens > 0 and close_parens > 0:
        if open_parens == close_parens:
            return 'stem_loop'
        else:
            return 'complex_nested'
    
    # Check for simple structures
    if dots > len(structure) * 0.8:
        return 'unstructured'
    
    return 'other'

# NOTE: The former simple_pseudoknot_detection() / simple_stem_loop_detection()
# helpers were removed. They FABRICATED a dot-bracket string and an energy
# (energy = -5 * n_stems) whenever the real tools did not run, and filled the
# probknot_* and pknots_* columns with that fake data. Structure evidence now
# comes only from RNAfold (nested MFE) and PKNOTS (pseudoknot-aware); if a tool
# does not run, its columns are left blank rather than fabricated.

def validate_trna_interaction(slippery_motif: str, trna_data: Optional[Dict] = None, organism: str = 'human') -> Dict:
    """Validate tRNA interaction potential for a slippery site.
    
    Args:
        slippery_motif: 7-nt motif (XXXYYYZ)
        trna_data: tRNA abundance data dictionary, or None if not available
        organism: organism for tRNA abundance data
    
    Returns:
        Dict with tRNA validation results
    """
    # Extract codons from the slippery site
    # Format: XXXYYYZ where X and Y are identical triplets
    codon1 = slippery_motif[:3]  # First triplet (P-site)
    codon2 = slippery_motif[3:6]  # Second triplet (A-site)
    
    # Check if tRNA data is available
    if trna_data is None:
        return {
            'codon1': codon1,
            'codon2': codon2,
            'trna1_abundance': 'NA',
            'trna2_abundance': 'NA',
            'pausing_potential': 'NA',
            'wobble_pairs': [],
            'trna_score': 'NA'
        }
    
    # Get tRNA abundance for the codons (user-supplied relative abundance proxy in
    # [0,1]; 0.5 = neutral default when a codon is absent from the table).
    trna1_abundance = trna_data.get(organism, {}).get(codon1, 0.5)
    trna2_abundance = trna_data.get(organism, {}).get(codon2, 0.5)

    # Ribosomal pausing potential -- LIMITING-CODON model.
    # Mechanistically, a ribosomal pause at the slippery site is driven by the
    # SLOWEST-decoded codon (the "hungry codon" effect: a single scarce cognate
    # tRNA sets the dwell time), not by the average of the two. We therefore define
    # pausing from the least-abundant of the P-/A-site codons:
    #     pausing_potential = 1 - min(a1, a2)      range [0,1]
    # (Replaces the earlier additive (1-a1)+(1-a2), which was a pure linear
    # transform of trna_score below and thus carried no independent information.)
    pausing_potential = 1.0 - min(trna1_abundance, trna2_abundance)

    # Wobble base pairing validation
    # Check for G-U wobble pairs in the P-site and A-site
    wobble_pairs = []
    if 'G' in codon1 and 'U' in codon1:
        wobble_pairs.append('P-site_G-U')
    if 'G' in codon2 and 'U' in codon2:
        wobble_pairs.append('A-site_G-U')

    # Overall tRNA interaction score -- GEOMETRIC MEAN (tAI-consistent).
    # The tRNA Adaptation Index aggregates per-codon weights as a geometric mean
    # (dos Reis et al. 2004), which is dominated by the smaller term and so is
    # more faithful to decoding kinetics than an arithmetic mean. Higher = both
    # codons are well supplied with tRNA (faster decoding).
    #     trna_score = sqrt(a1 * a2)               range [0,1]
    # This is NOT a deterministic function of pausing_potential above: a fast+slow
    # pair and two medium codons can share a trna_score yet differ in pausing.
    trna_score = (trna1_abundance * trna2_abundance) ** 0.5
    
    return {
        'codon1': codon1,
        'codon2': codon2,
        'trna1_abundance': trna1_abundance,
        'trna2_abundance': trna2_abundance,
        'pausing_potential': round(pausing_potential, 3),
        'wobble_pairs': wobble_pairs,
        'trna_score': round(trna_score, 3)
    }

def scan_seq(name: str,
             seq: str,
             spacer_min: int = 5,
             spacer_max: int = 9,
             window: int = 120,
             use_rnafold: bool = False,
             use_pknots: bool = False,
             use_nupack: bool = False,
             pknots_window: int = 85,
             pknots_top_n: int = 5,
             relaxed: bool = False,
             mfe_flag: float = -15.0,
             organism: str = 'human',
             trna_data: Optional[Dict] = None) -> List[dict]:
    """Scan one sequence for -1 PRF slippery heptamers and score downstream
    structure. One row is emitted per slippery site (the most stable spacer
    window is reported).

    Two-stage design (important for cost and correctness):
      1. RNAfold (nested MFE, fast) is run over every spacer in [spacer_min,
         spacer_max] for EVERY slippery site when use_rnafold=True; the most
         stable window is kept and used to RANK sites.
      2. PKNOTS (pseudoknot-aware, O(N^6), slow) is a CONFIRMATION step run on
         the short (<= pknots_window nt) best window of ONLY the top
         `pknots_top_n` sites by RNAfold MFE (set pknots_top_n<=0 or None to run
         it on all sites -- not recommended for genomes). PKNOTS is never run on
         the whole genome; it only ever folds these short PRF-like windows.

    There is no fabricated fallback: a tool that does not run leaves its columns
    blank (empty string).
    """
    # Normalize to RNA alphabet
    seq_rna = seq.upper().translate(DNA2RNA)
    hits = find_minus1_slippery_sites(seq_rna, relaxed=relaxed)
    results = []
    for pos0, motif in hits:
        spacer_start = pos0 + 7  # end of heptamer (0-based, exclusive)

        # --- RNAfold over the spacer range; keep the most stable window --------
        best = None  # (mfe, spacer, fold_start, fold_end, subseq, struct)
        for spacer in range(spacer_min, spacer_max + 1):
            fold_start = spacer_start + spacer
            fold_end = min(fold_start + window, len(seq_rna))
            if fold_start >= len(seq_rna) or fold_end - fold_start < 20:
                continue
            subseq = seq_rna[fold_start:fold_end]
            r_struct, r_mfe = run_rnafold(subseq) if use_rnafold else (None, None)
            # Rank by RNAfold MFE when available; otherwise keep the first valid
            # window so the site is still reported (structure columns blank).
            key = r_mfe if r_mfe is not None else float('inf')
            if best is None or key < best[0]:
                best = (key, spacer, fold_start, fold_end, subseq, r_struct, r_mfe)

        if best is None:
            continue
        _, spacer, fold_start, fold_end, subseq, rnafold_struct, rnafold_mfe = best

        trna_validation = validate_trna_interaction(motif, trna_data, organism)
        rnafold_type = classify_structure_type(rnafold_struct) if rnafold_struct else 'unknown'

        # PKNOTS is deferred: we build every row with RNAfold first, then run the
        # expensive PKNOTS confirmation only on the top-N by MFE (below). Columns
        # are left blank here and filled in the second pass for selected rows.
        results.append({
            'seqid': name,
            'site_start_1based': pos0 + 1,
            'slippery_motif': motif,
            'type': '-1',
            'spacer_nt': spacer,
            'fold_window_start': fold_start + 1,
            'fold_window_end': fold_end,
            'fold_window_seq': subseq,
            'rnafold_structure': rnafold_struct if rnafold_struct is not None else '',
            'rnafold_mfe': rnafold_mfe if rnafold_mfe is not None else '',
            'rnafold_type': rnafold_type,
            'pknots_structure': '',
            'pknots_energy': '',
            'pknots_has_pseudoknot': '',
            'pknots_type': 'unknown',
            'probknot_structure': '',
            'probknot_has_pseudoknot': '',
            'probknot_type': 'unknown',
            'nupack_structure': '',
            'nupack_energy': '',
            'nupack_has_pseudoknot': '',
            'nupack_type': 'unknown',
            'structure_call': '',
            'structure_confidence': '',
            'structure_note': '',
            'codon1': trna_validation['codon1'],
            'codon2': trna_validation['codon2'],
            'trna1_abundance': trna_validation['trna1_abundance'],
            'trna2_abundance': trna_validation['trna2_abundance'],
            'pausing_potential': trna_validation['pausing_potential'],
            'wobble_pairs': ','.join(trna_validation['wobble_pairs']) if trna_validation['wobble_pairs'] else 'none',
            'trna_score': trna_validation['trna_score']
        })

    # --- second pass: PKNOTS confirmation on the top-N sites by RNAfold MFE -----
    # Select which rows get the expensive pseudoknot-aware fold. Rows with an
    # RNAfold MFE are ranked (most stable first); if no RNAfold was run, fall back
    # to the first pknots_top_n rows so PKNOTS still targets a bounded set.
    if (use_pknots or use_nupack) and results:
        def _mfe(r):
            v = r['rnafold_mfe']
            return v if isinstance(v, (int, float)) else float('inf')
        ranked = sorted(range(len(results)), key=lambda i: _mfe(results[i]))
        selected = ranked if (pknots_top_n is None or pknots_top_n <= 0) else ranked[:pknots_top_n]
        for i in selected:
            r = results[i]
            window = r['fold_window_seq'][:pknots_window]
            # PKNOTS (bundled) -- only when requested (O(N^6), slow).
            pk_struct, pk_energy, pk_has = run_pknots_local(window) if use_pknots else (None, None, None)
            # ProbKnot (RNAstructure) -- independent pseudoknot-capable predictor.
            # Auto-runs on the same short window when the tools pass is active;
            # blank if RNAstructure is not installed.
            probk_struct, _probk_energy, probk_has = run_probknot_local(window)
            # NUPACK -- third pseudoknot-capable predictor (structure + MFE). Only
            # when requested; present only in the rw3594/dual_rag_if docker image.
            nup_struct, nup_energy, nup_has = run_nupack_local(window) if use_nupack else (None, None, None)
            call, conf, note = structure_consensus(
                r['rnafold_structure'] or None,
                r['rnafold_mfe'] if isinstance(r['rnafold_mfe'], (int, float)) else None,
                pk_struct, pk_energy, pk_has, mfe_flag=mfe_flag,
                probknot_struct=probk_struct, probknot_has_pk=probk_has,
                nupack_struct=nup_struct, nupack_energy=nup_energy, nupack_has_pk=nup_has)
            r['pknots_structure'] = pk_struct if pk_struct is not None else ''
            r['pknots_energy'] = pk_energy if pk_energy is not None else ''
            r['pknots_has_pseudoknot'] = ('' if pk_has is None else bool(pk_has))
            r['pknots_type'] = classify_structure_type(pk_struct) if pk_struct else 'unknown'
            r['probknot_structure'] = probk_struct if probk_struct is not None else ''
            r['probknot_has_pseudoknot'] = ('' if probk_has is None else bool(probk_has))
            r['probknot_type'] = classify_structure_type(probk_struct) if probk_struct else 'unknown'
            r['nupack_structure'] = nup_struct if nup_struct is not None else ''
            r['nupack_energy'] = nup_energy if nup_energy is not None else ''
            r['nupack_has_pseudoknot'] = ('' if nup_has is None else bool(nup_has))
            r['nupack_type'] = classify_structure_type(nup_struct) if nup_struct else 'unknown'
            r['structure_call'] = call
            r['structure_confidence'] = conf
            r['structure_note'] = note

    # For rows that did NOT get PKNOTS, still provide an RNAfold-only consensus so
    # structure_call is never blank when RNAfold ran.
    for r in results:
        if r['structure_call'] == '':
            call, conf, note = structure_consensus(
                r['rnafold_structure'] or None,
                r['rnafold_mfe'] if isinstance(r['rnafold_mfe'], (int, float)) else None,
                None, None, None, mfe_flag=mfe_flag)
            r['structure_call'], r['structure_confidence'], r['structure_note'] = call, conf, note

    return results

def write_csv(rows: List[dict], out_csv: Path):
    if not rows:
        with out_csv.open('w', newline='') as f:
            f.write('')  # empty file
        return
    fields = list(rows[0].keys())
    with out_csv.open('w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        for r in rows:
            w.writerow(r)

def write_bed(rows: List[dict], out_bed: Path):
    # BED6: chrom, chromStart, chromEnd, name, score, strand.
    # BED is 0-based, half-open. site_start_1based is 1-based inclusive, so the
    # heptamer spans 1-based [s, s+6]; in BED that is [s-1, s+6) -> chromStart = s-1,
    # chromEnd = s-1 + 7 = s + 6. (The previous code wrote the 1-based start directly
    # and used end = start + 6, shifting the interval by one nt and making it 6 nt.)
    with out_bed.open('w') as f:
        for r in rows:
            chrom = r['seqid']
            s = r['site_start_1based']
            chrom_start = s - 1        # 0-based
            chrom_end = s - 1 + 7      # half-open end (heptamer length 7)
            name = f"{r['type']}_{r['slippery_motif']}_s{r['spacer_nt']}"
            score = '0'
            strand = '+'
            f.write(f"{chrom}\t{chrom_start}\t{chrom_end}\t{name}\t{score}\t{strand}\n")

def main():
    ap = argparse.ArgumentParser(description="Scan genome for candidate PRF sites (-1)")
    ap.add_argument("--fasta", required=True, help="Genome FASTA")
    ap.add_argument("--out", required=True, help="Output prefix")
    ap.add_argument("--spacer-min", type=int, default=5)
    ap.add_argument("--spacer-max", type=int, default=9)
    ap.add_argument("--window", type=int, default=120, help="Downstream fold window size (nt) for the RNAfold ranking stage (default 120). This is deliberately generous: the exact 3' boundary of a -1 PRF stimulatory element is not known a priori, and typical stimulatory structures (pseudoknot/stem-loop) run ~50-90 nt (e.g. SARS-CoV-2 FSE ~88 nt, HIV-1 FSE ~40-80 nt), so 120 nt contains them with margin and avoids truncating a structure. RNAfold is O(N^3) (cheap), so a large window costs little here. NOTE: this is separate from --pknots-window; the expensive pseudoknot fold uses that smaller window. Lower it (e.g. 55) to target a specific short element.")
    ap.add_argument("--use-rnafold", action="store_true", help="Fold downstream window with RNAfold (ViennaRNA) if available on PATH")
    ap.add_argument("--use-pknots", action="store_true", help="Also fold the best window with PKNOTS (pseudoknot-aware, bundled binary). O(N^6): slow")
    ap.add_argument("--use-nupack", dest="use_nupack", action="store_true", default=None, help="Force NUPACK on (pseudoknot-aware, structure + MFE). NUPACK 3.2.2 must be on PATH (present in the rw3594/dual_rag_if docker image). Runs on the same top-N windows as PKNOTS. DEFAULT is auto: NUPACK runs automatically when its binary + parameter tables are detected, and is skipped silently otherwise.")
    ap.add_argument("--no-nupack", dest="use_nupack", action="store_false", help="Force NUPACK off even if it is installed (overrides the auto default).")
    ap.add_argument("--pknots-window", type=int, default=85, help="Max nt of the best window passed to the pseudoknot-aware folds PKNOTS/ProbKnot/NUPACK (default 85). Kept SMALLER than --window because PKNOTS is O(N^6): 120 nt would be ~(120/85)^6 ~ 8x slower than 85 nt. 85 nt still contains the short PRF-like stimulatory windows; larger is very slow.")
    ap.add_argument("--pknots-top-n", type=int, default=5, help="Run PKNOTS only on the N most stable sites by RNAfold MFE (default 5). PKNOTS is a confirmation step on short PRF-like windows, never a genome-wide fold. Use <=0 to run on all sites (not recommended for genomes).")
    ap.add_argument("--relaxed", action="store_true", help="Relaxed slippery grammar: allow any Y and any Z (default: Y in {A,U}, Z in {A,C,U})")
    ap.add_argument("--mfe-flag", type=float, default=-15.0, help="MFE (kcal/mol) at or below which a downstream structure counts as stable (default -15.0)")
    ap.add_argument("--organism", type=str, default="human", help="Organism label for tRNA proxy data (default: human)")
    ap.add_argument("--trna", type=str, default=None, help="Optional CSV file with tRNA abundance proxy data")
    args = ap.parse_args()

    # NUPACK is AUTO by default: run it when the binary + parameter tables are
    # present (e.g. inside the rw3594/dual_rag_if image), skip it silently
    # otherwise. --use-nupack / --no-nupack force the decision explicitly.
    if args.use_nupack is None:
        args.use_nupack = bool(_nupack_available and _nupack_available())
        if args.use_nupack:
            print("[INFO] NUPACK detected on PATH -> pseudoknot-aware NUPACK folding enabled (auto).")
    elif args.use_nupack and not (_nupack_available and _nupack_available()):
        print("[WARN] --use-nupack set but NUPACK 'mfe' + parameters not found; NUPACK columns will be blank.", file=sys.stderr)

    fasta = Path(args.fasta)
    if not fasta.exists():
        sys.stderr.write(f"[ERR] FASTA not found: {fasta}\n")
        sys.exit(2)
    seqs = read_fasta(fasta)

    trna_data = load_trna_data(Path(args.trna)) if args.trna else None

    all_rows = []
    for name, seq in seqs.items():
        rows = scan_seq(
            name, seq,
            spacer_min=args.spacer_min,
            spacer_max=args.spacer_max,
            window=args.window,
            use_rnafold=args.use_rnafold,
            use_pknots=args.use_pknots,
            use_nupack=args.use_nupack,
            pknots_window=args.pknots_window,
            pknots_top_n=args.pknots_top_n,
            relaxed=args.relaxed,
            mfe_flag=args.mfe_flag,
            organism=args.organism,
            trna_data=trna_data
        )
        all_rows.extend(rows)

    out_csv = Path(args.out + ".prf_candidates.csv")
    out_bed = Path(args.out + ".prf_candidates.bed")
    write_csv(all_rows, out_csv)
    write_bed(all_rows, out_bed)

    print(f"[OK] Candidates: {out_csv}")
    print(f"[OK] BED: {out_bed}")
    if args.use_rnafold and shutil.which('RNAfold') is None:
        print("[WARN] --use-rnafold set but RNAfold not found on PATH; structures left blank.", file=sys.stderr)

if __name__ == "__main__":
    main()
