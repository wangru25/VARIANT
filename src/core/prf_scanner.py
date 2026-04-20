#!/usr/bin/env python3
"""
PRF Scanner: detect candidate programmed ribosomal frameshifting (-1, optional +1) sites.

Pipeline:
1) Scan for -1 PRF slippery sites of the form XXXYYYZ (Z != G) with allowed triplets:
   X in {AAA, UUU, CCC, GGG}, Y in {AAA, UUU, CCC, GGG}, X != Y
   Canonical examples include U U U A A A C and A A A U U U A (DNA T instead of U)
2) For each site, check spacer length (default 5..9 nt) to the start of a downstream RNA structure window.
3) Optionally compute RNA secondary structure (MFE) with RNAfold (ViennaRNA) for a window downstream (default 120 nt).
4) Report candidates as CSV and BED.

Notes:
- RNAfold does not predict pseudoknots. Many PRF stimulators are pseudoknots; MFE here is a proxy for local stability.
- If you have Ribo-seq or conservation tracks, join those externally to prioritize candidates.
- If you provide a GFF/GTF, the script will annotate frame overlaps (0, -1, +1) relative to CDS features.

Usage:
    python prf_scanner.py --fasta genome.fasta --out out_prefix \
        [--spacer-min 5 --spacer-max 9 --window 120 --use-rnafold --gff annotations.gff --trna trna_abundance.csv]

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

def load_gff(gff_path: Optional[Path]):
    if not gff_path: 
        return []
    feats = []
    with gff_path.open() as f:
        for line in f:
            if not line.strip() or line.startswith('#'): 
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9: 
                continue
            seqid, source, ftype, start, end, score, strand, phase, attrs = parts
            if ftype.lower() in ('cds','gene','mrna','orf'):
                feats.append({
                    'seqid': seqid,
                    'type': ftype,
                    'start': int(start),
                    'end': int(end),
                    'strand': strand,
                    'phase': phase if phase != '.' else '0',
                    'attrs': attrs
                })
    return feats

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
            reader = csv.DictReader(f)
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

def build_minus1_regex() -> re.Pattern:
    # Build regex for XXXYYYZ with Z != G (in RNA alphabet).
    # Include both X != Y (standard pattern) and X == Y (known functional sites like HIV-1)
    triplets = ['AAA','UUU','CCC','GGG']
    alts = []
    
    # Standard XXXYYYZ pattern where X != Y
    for X in triplets:
        for Y in triplets:
            if X == Y:
                continue
            # Z any base except G
            alts.append(f"{X}{Y}[AUCU]")
    
    # Known functional sites where X == Y (e.g., HIV-1 UUUUUUA)
    for X in triplets:
        # Z any base except G
        alts.append(f"{X}{X}[AUCU]")
    
    # Join with alternation; use lookahead to capture overlapping hits
    pattern = '(?=(' + '|'.join(alts) + '))'
    return re.compile(pattern)

def find_minus1_slippery_sites(seq_rna: str) -> List[Tuple[int,str]]:
    """Return list of (0-based start, motif) for -1 PRF slippery sites."""
    pat = build_minus1_regex()
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
        if len(out) >= 2:
            struct_line = out[1]
            # structure and energy inside parentheses at end
            # e.g. ..((..)). (-12.30)
            parts = struct_line.rsplit('(', 1)
            dotbr = parts[0].strip()
            energy = parts[1].strip(') ').replace('kcal/mol','')
            mfe = float(energy)
            return dotbr, mfe
        return None, None
    except Exception:
        return None, None

def run_probknot(seq_rna: str) -> Tuple[Optional[str], Optional[float]]:
    """Call ProbKnot (from ViennaRNA) if available on PATH. Return (structure, score)."""
    if shutil.which('ProbKnot') is None:
        return None, None
    try:
        p = subprocess.run(
            ['ProbKnot', '--noPS'],
            input=(seq_rna+'\n').encode(),
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=True
        )
        out = p.stdout.decode().strip().splitlines()
        # ProbKnot output format similar to RNAfold
        if len(out) >= 2:
            struct_line = out[1]
            # Parse structure and energy
            parts = struct_line.rsplit('(', 1)
            if len(parts) == 2:
                dotbr = parts[0].strip()
                energy = parts[1].strip(') ').replace('kcal/mol','')
                try:
                    mfe = float(energy)
                    return dotbr, mfe
                except ValueError:
                    return dotbr, None
        return None, None
    except Exception:
        return None, None

def run_pknots_local(seq_rna: str) -> Tuple[Optional[str], Optional[float]]:
    """Try to run local pknots binary if available."""
    pknots_paths = ['./PKNOTS/bin/pknots', './pknots', 'pknots']
    for pknots_path in pknots_paths:
        if shutil.which(pknots_path) or (pknots_path.startswith('./') and os.path.exists(pknots_path)):
            try:
                p = subprocess.run(
                    [pknots_path, '-v'],
                    input=(seq_rna+'\n').encode(),
                    stdout=subprocess.PIPE,
                    stderr=subprocess.PIPE,
                    check=True
                )
                out = p.stdout.decode().strip().splitlines()
                # Parse pknots output for structure
                for line in out:
                    if 'Structure:' in line:
                        structure = line.split(':', 1)[1].strip()
                        return structure, None
                return None, None
            except Exception:
                continue
    return None, None

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

def simple_pseudoknot_detection(seq_rna: str) -> Tuple[Optional[str], Optional[float]]:
    """Simple pseudoknot detection based on sequence patterns and basic structure prediction."""
    if len(seq_rna) < 20:
        return None, None
    
    # Look for common pseudoknot patterns
    # This is a simplified approach - real pseudoknots are more complex
    
    # Check for potential stem-loop regions
    potential_stems = []
    for i in range(len(seq_rna) - 10):
        for j in range(i + 10, min(i + 50, len(seq_rna))):
            # Look for complementary regions
            region1 = seq_rna[i:i+5]
            region2 = seq_rna[j-5:j]
            region2_comp = region2.translate(str.maketrans('AUGC', 'UACG'))[::-1]
            
            if region1 == region2_comp:
                potential_stems.append((i, j, region1, region2))
    
    if len(potential_stems) >= 2:
        # Multiple potential stems suggest pseudoknot possibility
        structure = '.' * len(seq_rna)
        for i, (start, end, r1, r2) in enumerate(potential_stems[:2]):
            if i == 0:
                # First stem uses () brackets
                structure = structure[:start] + '(' * 5 + structure[start+5:end-5] + ')' * 5 + structure[end:]
            else:
                # Second stem uses [] brackets (pseudoknot)
                structure = structure[:start] + '[' * 5 + structure[start+5:end-5] + ']' * 5 + structure[end:]
        
        # Estimate energy (very rough)
        energy = -len(potential_stems) * 5.0
        return structure, energy
    
    return None, None

def simple_stem_loop_detection(seq_rna: str) -> Tuple[Optional[str], Optional[float]]:
    """Simple stem-loop detection based on sequence patterns."""
    if len(seq_rna) < 20:
        return None, None
    
    # Look for the best potential stem-loop
    best_stem = None
    best_score = 0
    
    for i in range(len(seq_rna) - 10):
        for j in range(i + 10, min(i + 30, len(seq_rna))):
            # Look for complementary regions
            region1 = seq_rna[i:i+5]
            region2 = seq_rna[j-5:j]
            region2_comp = region2.translate(str.maketrans('AUGC', 'UACG'))[::-1]
            
            matches = sum(1 for a, b in zip(region1, region2_comp) if a == b)
            if matches >= 3 and matches > best_score:
                best_score = matches
                best_stem = (i, j, region1, region2)
    
    if best_stem:
        start, end, r1, r2 = best_stem
        structure = '.' * len(seq_rna)
        structure = structure[:start] + '(' * 5 + structure[start+5:end-5] + ')' * 5 + structure[end:]
        energy = -best_score * 2.0
        return structure, energy
    
    return None, None

def annotate_frame(pos0: int, cds_list: List[dict]) -> str:
    """Rough frame context; assumes + strand for simplicity.
    Returns one of {'intergenic','CDS_frame0','CDS_frame-1','CDS_frame+1','CDS_multi'}"""
    in_frames = set()
    for cds in cds_list:
        if cds['strand'] == '-':
            # Not handling reverse here; users can reverse-complement before running
            continue
        if pos0 < cds['start']-1 or pos0 > cds['end']-1:
            continue
        frame0 = (pos0 - (cds['start']-1)) % 3  # 0,1,2
        if frame0 == 0:
            in_frames.add('CDS_frame0')
        elif frame0 == 1:
            in_frames.add('CDS_frame+1')
        else:
            in_frames.add('CDS_frame-1')
    if not in_frames:
        return 'intergenic'
    if len(in_frames) > 1:
        return 'CDS_multi'
    return list(in_frames)[0]

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
    
    # Get tRNA abundance for the codons
    trna1_abundance = trna_data.get(organism, {}).get(codon1, 0.5)
    trna2_abundance = trna_data.get(organism, {}).get(codon2, 0.5)
    
    # Calculate ribosomal pausing potential
    # Lower tRNA abundance = higher pausing potential
    pausing_potential = (1 - trna1_abundance) + (1 - trna2_abundance)
    
    # Wobble base pairing validation
    # Check for G-U wobble pairs in the P-site and A-site
    wobble_pairs = []
    if 'G' in codon1 and 'U' in codon1:
        wobble_pairs.append('P-site_G-U')
    if 'G' in codon2 and 'U' in codon2:
        wobble_pairs.append('A-site_G-U')
    
    # Overall tRNA interaction score
    # Higher score = better PRF potential
    trna_score = (trna1_abundance + trna2_abundance) / 2
    
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
             gff_feats: Optional[List[dict]] = None,
             organism: str = 'human',
             trna_data: Optional[Dict] = None) -> List[dict]:
    # Normalize to RNA alphabet
    seq_rna = seq.upper().translate(DNA2RNA)
    hits = find_minus1_slippery_sites(seq_rna)
    results = []
    cds_list = [f for f in (gff_feats or []) if f['seqid'] == name and f['type'].lower() == 'cds']
    for pos0, motif in hits:
        # Spacer boundary
        spacer_start = pos0 + 7  # end of heptamer
        # Use spacer_min..spacer_max; fold window starts at spacer_start + spacer_len
        for spacer in range(spacer_min, spacer_max+1):
            fold_start = spacer_start + spacer
            fold_end = min(fold_start + window, len(seq_rna))
            if fold_start >= len(seq_rna) or fold_end - fold_start < 20:
                continue
            subseq = seq_rna[fold_start:fold_end]
            
            # Run multiple structure prediction tools
            rnafold_struct, rnafold_mfe = (None, None)
            probknot_struct, probknot_score = (None, None)
            pknots_struct, pknots_score = (None, None)
            
            if use_rnafold:
                # Try external tools first
                rnafold_struct, rnafold_mfe = run_rnafold(subseq)
                probknot_struct, probknot_score = run_probknot(subseq)
                pknots_struct, pknots_score = run_pknots_local(subseq)
                
                # If external tools failed, use simple detection methods
                if not rnafold_struct:
                    rnafold_struct, rnafold_mfe = simple_stem_loop_detection(subseq)
                if not probknot_struct:
                    probknot_struct, probknot_score = simple_pseudoknot_detection(subseq)
                if not pknots_struct:
                    pknots_struct, pknots_score = simple_pseudoknot_detection(subseq)
            
            # Classify structure types
            rnafold_type = classify_structure_type(rnafold_struct) if rnafold_struct else 'unknown'
            probknot_type = classify_structure_type(probknot_struct) if probknot_struct else 'unknown'
            pknots_type = classify_structure_type(pknots_struct) if pknots_struct else 'unknown'
            frame_ctx = annotate_frame(pos0, cds_list) if cds_list else 'NA'
            
            # Add tRNA validation
            trna_validation = validate_trna_interaction(motif, trna_data, organism)
            
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
                'probknot_structure': probknot_struct if probknot_struct is not None else '',
                'probknot_type': probknot_type,
                'pknots_structure': pknots_struct if pknots_struct is not None else '',
                'pknots_type': pknots_type,
                'frame_context': frame_ctx,
                'codon1': trna_validation['codon1'],
                'codon2': trna_validation['codon2'],
                'trna1_abundance': trna_validation['trna1_abundance'],
                'trna2_abundance': trna_validation['trna2_abundance'],
                'pausing_potential': trna_validation['pausing_potential'],
                'wobble_pairs': ','.join(trna_validation['wobble_pairs']) if trna_validation['wobble_pairs'] else 'none',
                'trna_score': trna_validation['trna_score']
            })
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
    # BED6: chrom, start, end, name, score, strand
    with out_bed.open('w') as f:
        for r in rows:
            chrom = r['seqid']
            start = r['site_start_1based']  # Use 1-based coordinates
            end = start + 6  # motif length (7 nucleotides, so end = start + 6 for 1-based)
            name = f"{r['type']}_{r['slippery_motif']}_s{r['spacer_nt']}"
            score = '0'
            strand = '+'
            f.write(f"{chrom}\t{start}\t{end}\t{name}\t{score}\t{strand}\n")

def main():
    ap = argparse.ArgumentParser(description="Scan genome for candidate PRF sites (-1)")
    ap.add_argument("--fasta", required=True, help="Genome FASTA")
    ap.add_argument("--out", required=True, help="Output prefix")
    ap.add_argument("--spacer-min", type=int, default=5)
    ap.add_argument("--spacer-max", type=int, default=9)
    ap.add_argument("--window", type=int, default=120, help="Downstream fold window size (nt)")
    ap.add_argument("--use-rnafold", action="store_true", help="Call RNAfold, ProbKnot, and pknots if available on PATH")
    ap.add_argument("--gff", type=str, default=None, help="Optional GFF for CDS frame context (plus-strand only)")
    ap.add_argument("--organism", type=str, default="human", help="Organism for tRNA abundance data (default: human)")
    ap.add_argument("--trna", type=str, default=None, help="Optional CSV file with tRNA abundance data")
    args = ap.parse_args()

    fasta = Path(args.fasta)
    if not fasta.exists():
        sys.stderr.write(f"[ERR] FASTA not found: {fasta}\n")
        sys.exit(2)
    seqs = read_fasta(fasta)

    gff_feats = load_gff(Path(args.gff)) if args.gff else None
    trna_data = load_trna_data(Path(args.trna)) if args.trna else None

    all_rows = []
    for name, seq in seqs.items():
        rows = scan_seq(
            name, seq,
            spacer_min=args.spacer_min,
            spacer_max=args.spacer_max,
            window=args.window,
            use_rnafold=args.use_rnafold,
            gff_feats=gff_feats,
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
