# `src/core/prf` — Programmed Ribosomal Frameshifting (PRF) detection

This subpackage detects candidate **−1** and **+1** programmed ribosomal
frameshifting (PRF) sites in viral genomes. It has two independent detectors and a
set of RNA-folding wrappers.

| File | Role |
|------|------|
| `prf_scanner.py` | **−1 PRF** two-stage scanner (regex slippery search → RNAfold ranking → PKNOTS/ProbKnot/NUPACK pseudoknot confirmation → transparent consensus). Runnable as a CLI or importable (`scan_seq`). |
| `frameshift_detector.py` | **+1 PRF** codon-context detector (`FrameshiftDetector`), plus its own −1 slippery calls. Wired into the `variant prf` CLI. |
| `pknots_runner.py` | Wrapper around the bundled **PKNOTS** binary (`PKNOTS/src/pknots`) — pseudoknot-aware, MFE, O(N⁶). |
| `probknot_runner.py` | Wrapper around **ProbKnot** (RNAstructure) — pseudoknot-aware, probability-based (ThreshKnot, no MFE). Optional. |
| `nupack_runner.py` | Wrapper around **NUPACK 3.2.2** — pseudoknot-aware, MFE. Reached via a Docker image. Optional. |

> **No fabricated structures.** Whenever a folding tool is unavailable or fails, its
> columns are left blank (or `NA`). The pipeline never invents a structure or energy.

---

## 1. −1 PRF: `prf_scanner.py`

### Algorithm (two stages)

1. **Slippery-site search.** A regular expression with backreferences matches the
   canonical heptamer grammar `X XXY YYZ` (nucleotides X X X Y Y Y Z), so *every*
   conforming site is found — not a hardcoded list.
   - Strict (default): `X` = any base ×3, `Y ∈ {A,U}` ×3, `Z ∈ {A,C,U}` (not G).
     Implemented pattern: `(?=(([ACGU])\2\2([AU])\3\3([ACU])))`.
   - Relaxed (`--relaxed`): any `Y`, any `Z`.
   - Canonical elements (SARS-CoV-2 `UUUAAAC`, HIV-1 `UUUUUUA`) are recovered as
     special cases.
2. **RNAfold ranking (fast).** For each site, RNAfold (ViennaRNA) folds the downstream
   window over spacer lengths `--spacer-min..--spacer-max` and keeps the most stable
   window; the nested MFE ranks all sites.
3. **Pseudoknot confirmation (slow).** On the top `--pknots-top-n` short
   (≤ `--pknots-window` nt) windows only, three pseudoknot-aware predictors run:
   **PKNOTS** (MFE), **ProbKnot** (ThreshKnot), **NUPACK** (MFE). A pseudoknot is
   determined by a real topological test for crossing base pairs on each predicted
   structure.
4. **Transparent consensus.** `structure_consensus()` combines the evidence:
   a pseudoknot is *called* by majority vote among the pseudoknot-capable tools that
   ran; ties favour a pseudoknot at low confidence; disagreements are reported in
   `structure_note`, never hidden. RNAfold cannot represent pseudoknots and never
   votes on pseudoknot presence.

### CLI

```bash
python src/core/prf/prf_scanner.py --fasta <genome.fasta> --out <prefix> \
    --use-rnafold --use-pknots \
    [--pknots-top-n N] [--pknots-window 85] \
    [--use-nupack | --no-nupack] [--relaxed] \
    [--spacer-min 5 --spacer-max 9 --window 120] [--mfe-flag -15.0] \
    [--trna <table.csv> --organism human]
```

Key flags:

| Flag | Default | Meaning |
|------|---------|---------|
| `--use-rnafold` | off | Fold + rank with RNAfold (required for ranking). |
| `--use-pknots` | off | Run the PKNOTS/ProbKnot/NUPACK confirmation pass. |
| `--pknots-top-n` | 5 | Confirm only the N most stable sites. **`<= 0` folds every candidate window.** |
| `--window` | 120 | Downstream fold window (nt) for the RNAfold **ranking** stage (see below). |
| `--pknots-window` | 85 | Max nt passed to the pseudoknot folds — smaller because PKNOTS is O(N⁶). |
| `--spacer-min` / `--spacer-max` | 5 / 9 | Spacer range (nt) between the heptamer and the fold window. |
| `--use-nupack` / `--no-nupack` | auto | Force NUPACK on/off. Default: auto-enable when reachable. |
| `--relaxed` | off | Relaxed slippery grammar (any Y, any Z). |
| `--trna` / `--organism` | none / human | Optional tRNA proxy table (see §3). |

**Two fold windows, and why the defaults are 120 / 85.** For each slippery site the
downstream stimulatory region is folded twice, with *different* window sizes:

- **`--window` (120 nt, RNAfold ranking).** The 3′ boundary of a −1 PRF stimulatory
  element is not known a priori, and real stimulatory structures run ~50–90 nt
  (SARS-CoV-2 FSE ~88 nt; HIV-1 FSE ~40–80 nt), so a generous 120-nt window contains
  them with margin and avoids truncating a structure. RNAfold is O(N³) (cheap), so a
  large window costs almost nothing at this stage. It is a conservative "don't cut the
  structure" default, **not** a biological constant — lower it (e.g. `--window 55`) to
  target a specific short element. (This is why `fold_window_seq` in the output is
  120 nt even for a ~53-nt canonical FSE construct.)
- **`--pknots-window` (85 nt, pseudoknot folds).** PKNOTS is O(N⁶), so 120 nt would be
  ~(120/85)⁶ ≈ 8× slower; 85 nt still covers the short PRF-like windows while keeping
  the expensive fold affordable.

### Output columns (`*.prf_candidates.csv`)

`seqid, site_start_1based, slippery_motif, type, spacer_nt, fold_window_start,
fold_window_end, fold_window_seq,` then per tool
`rnafold_structure/mfe/type`, `pknots_structure/energy/has_pseudoknot/type`,
`probknot_structure/has_pseudoknot/type`, `nupack_structure/energy/has_pseudoknot/type`,
the consensus `structure_call / structure_confidence / structure_note`, and the tRNA
proxy columns `codon1, codon2, trna1_abundance, trna2_abundance, pausing_potential,
wobble_pairs, trna_score`. A `*.prf_candidates.bed` (BED6) gives coordinates.

`probknot_has_pseudoknot` has no energy (ProbKnot is probability-based).
NUPACK columns are blank when NUPACK is not reachable.

---

## 2. +1 PRF: `frameshift_detector.py`

+1 PRF is codon-context driven, so `FrameshiftDetector.detect_frameshift_sites()`
scans **in frame**: for each (P-site, A-site) codon pair it flags, in priority order,

1. **Curated, literature-validated sites** (exact P-/A-site pair): antizyme OAZ1
   `UCC-UGA`, antizyme `GCG-UGA`, Ty1 `CUU-AGG`, Ty3 `GCG-AGU`, *E. coli* prfB
   `CUU-UGA`. (`mechanism = known_plus1_site` — high-confidence subset.)
2. **Proline P-site + stop A-site** (`CCN` then a stop): forward P-site slippage,
   e.g. the Ebola `CCU-UAA`-type site. (`mechanism = proline_psite_slip`.)
3. **Leaky "shifty" stops** in the A-site (UGA > UAG > UAA leakiness). Frame-agnostic
   and therefore **low-specificity** — treat as candidates. (`mechanism = shifty_stop`.)

The same class also emits its own −1 slippery-heptamer calls (`type = -1 PRF`).

---

## 3. tRNA proxy columns (optional, not measured data)

With `--trna <table.csv>`, the P-site (`a1`) and A-site (`a2`) codons are looked up in
a user-supplied relative-abundance table (values in `[0,1]`) and two **non-redundant**
proxies are reported:

- `trna_score = sqrt(a1 · a2)` — geometric mean, tRNA-adaptation-index-consistent
  (dos Reis et al., 2004, *Nucleic Acids Res.* 32:5036).
- `pausing_potential = 1 − min(a1, a2)` — limiting-codon ("hungry-codon") pausing.

These are **interpretable, sequence-derived proxies — not measured tRNA
concentrations, kinetics, or ribosome-profiling data.** Without `--trna` they are
reported as `NA`; a codon missing from a supplied table defaults to a neutral `0.5`.
The bundled `data/tRNA_abundance_sample.csv` is an **illustrative placeholder** (4
codons) — do not use it for real analysis. See `data/tRNA_abundance_README.md`.

---

## 4. Tool availability

| Tool | How it is provided | Notes |
|------|--------------------|-------|
| RNAfold | ViennaRNA (conda `variant` env) | Nested MFE; ranking only. |
| PKNOTS | Bundled at `PKNOTS/src/pknots` | No install needed. |
| ProbKnot | RNAstructure — `conda install -n variant rnastructure` | Optional; needs `DATAPATH` to its `data_tables` (auto-located next to the binary). |
| NUPACK | Docker image `rw3594/dual_rag_if:1.0` (~6.5 GB) | Optional; VARIANT auto-shells-out to the container (you do NOT run inside it). Auto-enabled when reachable, else blank. Full setup + verification: **root README → "Enabling NUPACK via Docker"**. |

Validation on the SARS-CoV-2 88-nt frameshift element: PKNOTS −26.68 kcal/mol
(pseudoknot = yes), NUPACK ≈ −27.6 to −29.5 kcal/mol (pseudoknot = yes), ProbKnot
calls the same window nested — a genuine tool disagreement reported (not hidden) by
the consensus logic.

---

## 5. How each entry point uses this subpackage

`src/core/mutation_processor.py::_run_frameshift_analysis` decides what runs, keyed
on `enable_structure_prediction`:

- **Always** runs the fast `FrameshiftDetector` (pure codon/regex scan, **no RNA
  folding**) → `<Virus>_frameshift_detector.csv`. This is the file the PRF figure
  (`src/visualization/figure3_PRF.py`) draws from.
- **Only when `enable_structure_prediction` is True** (local CLI default) it *also*
  runs the −1 two-stage folding scanner (`prf_scanner.py`, RNAfold + PKNOTS/ProbKnot/
  NUPACK) → `<Virus>.prf_candidates.csv/.bed`.

Per entry point:
- `python main.py --virus <V> --detect-frameshifts` → detector CSV **+** the −1
  folding scanner. Add `--skip-structure-prediction` to get only the detector CSV.
- `variant prf --genome <fasta> --output <csv>` → the **+1/−1 `FrameshiftDetector`**
  directly, via `src/cli/commands.py`.
- `web_app.py` → sets `enable_structure_prediction=False`, so it runs **only the
  `FrameshiftDetector`** (no secondary-structure prediction at all) and the web PRF
  figure is drawn from `<Virus>_frameshift_detector.csv`.
