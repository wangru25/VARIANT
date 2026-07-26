# tRNA abundance proxy — input format and how to supply real data

The PRF scanner can annotate each candidate frameshift site with two optional,
codon-level quantities derived from a **user-supplied** tRNA-abundance proxy
table. These are interpretable annotations, **not** measured tRNA concentrations
or kinetics.

## Is this required?

**No.** It is entirely optional.

- If you do **not** pass `--trna <file>`, the scanner runs normally and reports
  the tRNA/pausing columns as `NA`. Sequence-grammar detection and the
  RNAfold/PKNOTS/ProbKnot/NUPACK structure consensus are unaffected.
- If you pass a table but a particular codon is missing from it, that codon
  defaults to a neutral abundance of `0.5` (no bias in either direction).

The scanner never fabricates tRNA values.

## File format

A CSV with exactly three columns:

```
codon,organism,abundance
UUU,human,0.42
AAA,human,0.31
...
```

| column     | meaning |
|------------|---------|
| `codon`    | codon in RNA letters (e.g. `UUU`) |
| `organism` | organism tag; must match the scanner's `--organism` value (e.g. `human`, `sars_cov2`, `hiv_1`) |
| `abundance`| relative decoding-tRNA supply in `[0,1]`. **Higher = plentiful cognate tRNA, fast decoding; lower = scarce, slow decoding (more pausing).** |

Lines beginning with `#` and blank lines are ignored, so you may keep a
provenance note at the top of the file.

## How the values are used

For a slippery site, the P-site codon (`a1`) and A-site codon (`a2`) are looked
up and combined into two **non-redundant** quantities:

- `trna_score = sqrt(a1 * a2)` — geometric mean; overall aminoacyl-tRNA supply
  at the site. Geometric (not arithmetic) mean for consistency with the tRNA
  Adaptation Index (dos Reis et al., 2004, *Nucleic Acids Res* 32:5036-5044).
- `pausing_potential = 1 - min(a1, a2)` — limiting-codon pausing: the pause is
  set by the slowest-decoded of the two codons ("hungry-codon" effect).

## The bundled sample is a PLACEHOLDER

`tRNA_abundance_sample.csv` in this folder contains only 4 codons
(AAA/UUU/CCC/GGG) with hand-picked round numbers. It exists to demonstrate the
**format only** and must not be used for real analysis or in a publication.

## Supplying real data

Replace the sample with a full table covering the 61 sense codons for **your
host organism**. Recommended sources:

- **GtRNAdb** — tRNA gene-copy numbers per anticodon (gtrnadb.ucsc.edu); map
  to codons and normalise to `[0,1]`.
- **tRNA Adaptation Index (tAI) weights** — per-codon `w_i` values (dos Reis
  et al., 2004); already in `[0,1]` and directly usable as the `abundance`
  column.
- **Ribosome-profiling dwell times** — invert and normalise so that long dwell
  (slow codon) maps to low abundance.

Whatever the source, document it in a `#` header at the top of your CSV so the
provenance travels with the file.
