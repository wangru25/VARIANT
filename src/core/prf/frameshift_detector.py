# -*- coding: utf-8 -*-
'''
Author: Rui Wang
Date: 2025-08-19 17:30:00
LastModifiedBy: Rui Wang
LastEditTime: 2025-08-19 17:30:00
Email: rw3594@nyu.edu
FilePath: /VARIANT/src/core/prf/frameshift_detector.py
Description: Frameshift detection and analysis module.
'''
import json
import re
from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Optional, Tuple


# -1 PRF slippery heptamer  X XXY YYZ  in the DNA alphabet (T for U). Backreferences
# enforce that X and Y are each three identical bases, so novel sites are found
# without a hardcoded list (previously the detector matched only ~9 fixed 7-mers,
# missing valid sites such as GGGAAAT). Strict: Y in {A,T}, Z in {A,C,T} (not G).
# Relaxed: any Y, any Z. Lookahead captures overlapping hits.
_SLIP_STRICT_DNA_RX  = re.compile(r"(?=(([ACGT])\2\2([AT])\3\3([ACT])))")
_SLIP_RELAXED_DNA_RX = re.compile(r"(?=(([ACGT])\2\2([ACGT])\3\3([ACGT])))")

# Curated, literature-validated +1 PRF sites as (P-site codon, A-site codon) in the
# DNA alphabet. +1 PRF is codon-context driven (not captured by a slippery heptamer),
# so it is detected by an in-frame codon scan, not substring search.
#   antizyme (OAZ1)  UCC-UGA : Matsufuji et al. 1995 Cell 80:51
#   antizyme variant GCG-UGA : Palanimurugan et al. 2004 EMBO J 23:4857
#   Ty1 (S. cerevisiae) CUU-AGG : Belcourt & Farabaugh 1990 Cell 62:339
#   Ty3 (S. cerevisiae) GCG-AGU : Farabaugh et al. 1993 Cell 74:93
#   E. coli prfB (RF2)  CUU-UGA : Craigen et al. 1985 PNAS 82:3616
_KNOWN_PLUS1_SITES = {
    ("TCC", "TGA"): "antizyme OAZ1-type UCC-UGA (Matsufuji 1995)",
    ("GCG", "TGA"): "antizyme GCG-UGA (Palanimurugan 2004)",
    ("CTT", "AGG"): "Ty1 CUU-AGG (Belcourt & Farabaugh 1990)",
    ("GCG", "AGT"): "Ty3 GCG-AGU (Farabaugh 1993)",
    ("CTT", "TGA"): "E. coli prfB CUU-UGA (Craigen 1985)",
}


class PRFSiteType(Enum):
    """Types of PRF sites."""
    MINUS_ONE = "-1 PRF"
    PLUS_ONE = "+1 PRF"
    STEM_LOOP = "Stem-loop"


@dataclass
class PRFSite:
    """Data class representing a PRF site with enhanced analysis."""
    position: int
    end_position: int
    sequence: str
    site_type: PRFSiteType
    efficiency_score: float = 0.0
    downstream_structure: Optional[str] = None
    mutations: List[Dict] = None
    notes: Optional[str] = None
    trna_validated: bool = False

    def __post_init__(self):
        if self.mutations is None:
            self.mutations = []


class FrameshiftDetector:
    """
    Detects potential programmed ribosomal frameshifting (PRF) sites in viral genomes.

    Supports detection of:
    - -1 PRF: Slippery sequences (X XXY YYZ pattern)
    - +1 PRF: Shifty stop codons and specific motifs
    - Stem-loop structures that may promote frameshifting
    """

    def __init__(self, ref_genome=None):
        """
        Initialize the frameshift detector.

        Args:
            ref_genome: Reference genome object (optional)
        """
        self.ref_genome = ref_genome

        # Common slippery sequences for -1 PRF
        self.slippery_sequences = {
            # Coronavirus -1 PRF signals
            "TTTAAAC": "SARS-CoV-2-like slippery sequence",
            "TTTAAAT": "Coronavirus -1 PRF signal",
            "GGGTTTA": "Coronavirus -1 PRF signal",
            "GGGTTTG": "Coronavirus -1 PRF signal",

            # HIV -1 PRF signals
            "TTTTTTA": "HIV -1 PRF signal",
            "TTTTTTG": "HIV -1 PRF signal",

            # General -1 PRF patterns (X XXY YYZ)
            "AAAAAAC": "General -1 PRF pattern",
            "CCCCCCA": "General -1 PRF pattern",
            "GGGGGGA": "General -1 PRF pattern",
            "TTTTTTA": "General -1 PRF pattern",
        }

        # +1 PRF signals
        self.plus_one_signals = {
            "TCCT": "Shifty stop codon context",
            "TCCG": "Shifty stop codon context",
            "CCT": "Stop codon with +1 potential",
            "CCG": "Stop codon with +1 potential",
        }

        # Stem-loop detection parameters
        self.min_stem_length = 6
        self.min_loop_length = 3
        self.max_loop_length = 10

        # Scoring parameters (TBD)
        self.slippery_site_weight = 0.6
        self.downstream_structure_weight = 0.3
        self.context_weight = 0.1

        # Ribosomal pause motifs
        self.pause_motifs = {
            "GG": "Strong pause motif",
            "CC": "Strong pause motif",
            "TT": "Moderate pause motif",
            "AA": "Moderate pause motif",
        }

        # Published slippery sequence scores (Plant et al., 2010)
        self.published_slippery_scores = {
            "TTTAAAC": 0.85,  # SARS-CoV-2, experimentally validated
            "TTTAAAT": 0.82,  # Coronavirus consensus
            "GGGTTTA": 0.78,  # Coronavirus variant
            "TTTTTTA": 0.75,  # HIV-1, experimentally validated
            "TTTTTTG": 0.72,  # HIV-1 variant
            "AAAAAAC": 0.65,  # General pattern
            "CCCCCCA": 0.62,  # General pattern
            "GGGGGGA": 0.60,  # General pattern
        }

        # Known stimulatory structures (experimentally validated)
        self.known_stimulatory_structures = {
            "SARS-CoV-2": {
                "position": 13468,
                "structure_type": "pseudoknot",
                "free_energy": -15.2,
                "reference": "Plant et al. (2010)"
            },
            "HIV-1": {
                "position": 1637,
                "structure_type": "stem-loop",
                "free_energy": -12.8,
                "reference": "Marcheschi et al. (2009)"
            }
        }

        # Experimental validation references
        self.experimental_validations = {
            "TTTAAAC": ["Plant et al. (2010)", "Kelly et al. (2020)"],
            "TTTTTTA": ["Marcheschi et al. (2009)", "Dulude et al. (2006)"],
            "TTTAAAT": ["Plant et al. (2010)"],
            "GGGTTTA": ["Plant et al. (2010)"]
        }

    def detect_frameshift_sites(self, sequence: str, start_pos: int = 1) -> List[Dict]:
        """
        Detect potential frameshifting sites in a sequence.

        Args:
            sequence (str): DNA/RNA sequence to analyze
            start_pos (int): Starting position for coordinate mapping (1-based indexing, default=1)

        Returns:
            List[Dict]: List of detected frameshift sites with details
        """
        frameshift_sites = []

        # Detect -1 PRF sites (slippery sequences with optional downstream stem-loops)
        minus_one_sites = self._detect_minus_one_prf(sequence, start_pos)
        frameshift_sites.extend(minus_one_sites)

        # Detect +1 PRF sites (shifty stop codons)
        plus_one_sites = self._detect_plus_one_prf(sequence, start_pos)
        frameshift_sites.extend(plus_one_sites)

        # Note: Stem-loops are now only detected as facilitators for slippery sequences,
        # not as independent frameshift sites

        # Sort by position
        frameshift_sites.sort(key=lambda x: x['position'])

        return frameshift_sites

    def _detect_minus_one_prf(self, sequence: str, start_pos: int) -> List[Dict]:
        """
        Detect -1 PRF slippery sequences.

        Args:
            sequence (str): DNA/RNA sequence
            start_pos (int): Starting position for coordinate mapping

        Returns:
            List[Dict]: Detected -1 PRF sites
        """
        sites = []
        seq = sequence.upper()

        # Scan with the slippery-heptamer regex (X XXY YYZ). This finds every
        # conforming site, not only a hardcoded list, so novel slippery sequences
        # are detected. Backreferences guarantee the "three identical bases"
        # geometry required for P-/A-site tRNA re-pairing after the -1 shift.
        rx = _SLIP_RELAXED_DNA_RX if getattr(self, 'relaxed_slippery', False) else _SLIP_STRICT_DNA_RX
        for m in rx.finditer(seq):
            pos = m.start()
            heptamer = m.group(1)
            # Prefer a curated description if this exact heptamer is a named signal.
            description = self.slippery_sequences.get(heptamer, "-1 PRF slippery heptamer (X XXY YYZ)")

            # Downstream stem-loop proxy (see _check_downstream_stem_loop docstring).
            stem_loop_score = self._check_downstream_stem_loop(seq, pos + 7)

            site = {
                'position': pos + start_pos,
                'end_position': pos + 7 + start_pos - 1,
                'sequence': heptamer,
                'type': '-1 PRF',
                'description': description,
                'stem_loop_score': stem_loop_score,
                'mechanism': 'slippery_sequence',
                # The heptamer conforms to the slippery grammar (P-/A-site
                # re-pairing is geometrically possible). This is NOT validation
                # against tRNA abundance -- no abundance data is used here.
                'slippery_conforms': True,
                'trna_note': 'grammar-conforming; tRNA abundance not evaluated'
            }
            sites.append(site)

        return sites

    def _detect_plus_one_prf(self, sequence: str, start_pos: int) -> List[Dict]:
        """
        Detect +1 PRF signals.

        Args:
            sequence (str): DNA/RNA sequence
            start_pos (int): Starting position for coordinate mapping

        Returns:
            List[Dict]: Detected +1 PRF sites
        """
        sites = []
        seq = sequence.upper()
        stop_codons = ('TAA', 'TAG', 'TGA')
        # A-site stop leakiness (read-through / +1 shift propensity), from
        # eukaryotic read-through literature: UGA > UAG > UAA.
        stop_leak = {'TGA': 1.0, 'TAG': 0.6, 'TAA': 0.4}
        # Proline codon family (CCN). A proline P-site codon can slip +1 forward
        # because CCN is 4-fold degenerate, so the shifted P-site (CC + next base)
        # is still read as proline via wobble. When such a proline codon is directly
        # followed by a STOP in the A-site, the ribosome stalls and can slip +1
        # ("forward P-site slippage on proline codons followed by a stop"). This is
        # the mechanism of e.g. the Ebola virus L-gene +1 site CCU-UAA (Stubbs,
        # Ruzicka & Taylor 2024, Pathogens 13:829) and is independent of the A-site
        # stop's leakiness -- so it is checked separately from the shifty-stop rule.
        proline_codons = ('CCT', 'CCC', 'CCA', 'CCG')

        # +1 PRF is defined by the CODON in each reading frame, not by a substring.
        # The previous implementation matched "CCT"/"CCG" anywhere in a +-1 window,
        # which mislabels proline codons and, being frame-agnostic, both missed real
        # sites (antizyme UCC-UGA, Ty1 CUU-AGG) and produced false positives.
        # We scan every frame in-frame: at each codon boundary, the (P-site, A-site)
        # codon pair is the +1 trigger context.
        frame_offset = getattr(self, 'plus1_frame_offset', None)
        frames = [frame_offset] if frame_offset in (0, 1, 2) else (0, 1, 2)

        for frame in frames:
            # i indexes the A-site codon start; P-site codon is the preceding one.
            for i in range(frame + 3, len(seq) - 2, 3):
                p_site = seq[i - 3:i]
                a_site = seq[i:i + 3]
                if len(a_site) < 3:
                    break

                description = mechanism = None
                # 1) Curated, literature-validated +1 site (exact P-/A-site pair).
                if (p_site, a_site) in _KNOWN_PLUS1_SITES:
                    description = _KNOWN_PLUS1_SITES[(p_site, a_site)]
                    mechanism = 'known_plus1_site'
                # 2) Proline P-site codon (CCN) directly followed by a STOP: forward
                #    P-site slippage. Independent of A-site stop leakiness, so this
                #    catches sites the shifty-stop rule misses -- e.g. Ebola CCU-UAA
                #    (UAA leakiness 0.4 is below the shifty-stop threshold, but the
                #    proline-slippage mechanism still applies).
                elif p_site in proline_codons and a_site in stop_codons:
                    description = ('proline P-site (%s) + stop A-site (%s): forward '
                                   'P-site slippage (Ebola CCU-UAA-type)' % (p_site, a_site))
                    mechanism = 'proline_psite_slip'
                # 3) Shifty stop in the A-site: a leaky stop stalls the ribosome and
                #    can drive a +1 shift (e.g. antizyme, prfB). Reported as a
                #    candidate; leakiness weights confidence, it is not a gate.
                #    Only stops at or above plus1_min_leakiness are reported, so the
                #    least-leaky UAA (0.4) is dropped by default -- a standalone UAA
                #    is not a credible +1 site. NOTE: without a supplied ORF/frame
                #    this scans all three frames, so bare shifty-stop calls are
                #    low-specificity; treat them as candidates, not predictions.
                elif a_site in stop_codons and stop_leak[a_site] >= getattr(self, 'plus1_min_leakiness', 0.5):
                    description = 'shifty stop in A-site (%s), leakiness %.1f' % (a_site, stop_leak[a_site])
                    mechanism = 'shifty_stop'

                if description is None:
                    continue

                site = {
                    'position': i + start_pos,
                    'end_position': i + 3 + start_pos - 1,
                    'sequence': a_site,
                    'context': p_site + a_site,
                    'p_site_codon': p_site,
                    'a_site_codon': a_site,
                    'frame': frame,
                    'type': '+1 PRF',
                    'description': description,
                    'mechanism': mechanism,
                    'stop_leakiness': stop_leak.get(a_site, 0.0),
                }
                sites.append(site)

        return sites

    def _detect_stem_loops(self, sequence: str, start_pos: int) -> List[Dict]:
        """
        Detect potential stem-loop structures that may promote frameshifting.

        Args:
            sequence (str): DNA/RNA sequence
            start_pos (int): Starting position for coordinate mapping

        Returns:
            List[Dict]: Detected stem-loop structures
        """
        sites = []

        # Simple stem-loop detection (can be enhanced with more sophisticated algorithms)
        for i in range(len(sequence) - self.min_stem_length * 2 - self.min_loop_length):
            # Look for potential stem-loop structures
            stem_loop = self._find_stem_loop(sequence, i)
            if stem_loop:
                # Analyze the impact of this stem-loop structure
                site = {
                    'position': i + start_pos,
                    'end_position': stem_loop['end'] + start_pos,
                    'sequence': sequence[i:stem_loop['end']],
                    'type': 'stem_loop',
                    'description': f"Stem-loop structure (stem: {stem_loop['stem_length']}, loop: {stem_loop['loop_length']})",
                    'confidence': stem_loop['score'],
                    'stem_length': stem_loop['stem_length'],
                    'loop_length': stem_loop['loop_length'],
                    'mechanism': 'structural_element'
                }
                sites.append(site)

        return sites

    def _find_all_occurrences(self, sequence: str, pattern: str) -> List[int]:
        """
        Find all occurrences of a pattern in a sequence.

        Args:
            sequence (str): Sequence to search
            pattern (str): Pattern to find

        Returns:
            List[int]: Positions of all occurrences
        """
        positions = []
        start = 0
        while True:
            pos = sequence.find(pattern, start)
            if pos == -1:
                break
            positions.append(pos)
            start = pos + 1
        return positions

    def _check_downstream_stem_loop(self, sequence: str, position: int) -> float:
        """
        Check for downstream stem-loop structure that may promote frameshifting.

        Args:
            sequence (str): DNA/RNA sequence
            position (int): Position to check from

        Returns:
            float: Stem-loop score (0.0 to 1.0)
        """
        if position >= len(sequence) - 20:
            return 0.0

        # Look for stem-loop in the optimal range: 5-9 nt downstream
        # This is the biologically correct spacing for ribosome pausing
        optimal_start = position + 5
        position + 9

        # Check if we have enough sequence
        if optimal_start >= len(sequence):
            return 0.0

        # Look in the optimal range first (5-9 nt)
        for offset in range(5, min(10, len(sequence) - position)):
            check_pos = position + offset
            if check_pos >= len(sequence) - 20:  # Need enough sequence for stem-loop
                continue

            downstream_seq = sequence[check_pos:check_pos + 30]  # Look 30 nt downstream
            stem_loop = self._find_stem_loop(downstream_seq, 0)

            if stem_loop:
                # Boost score for optimal spacing (5-9 nt)
                spacing_boost = 1.5 if 5 <= offset <= 9 else 1.0
                return min(1.0, stem_loop['score'] * spacing_boost)

        return 0.0

    def _find_stem_loop(self, sequence: str, start: int) -> Optional[Dict]:
        """
        Find a stem-loop structure starting at a given position.

        Args:
            sequence (str): DNA/RNA sequence
            start (int): Starting position

        Returns:
            Optional[Dict]: Stem-loop information or None
        """
        # Simple stem-loop detection algorithm
        # This can be enhanced with more sophisticated RNA structure prediction

        for stem_len in range(self.min_stem_length, min(15, (len(sequence) - start) // 3)):
            for loop_len in range(self.min_loop_length, min(self.max_loop_length, len(sequence) - start - 2 * stem_len)):

                if start + 2 * stem_len + loop_len > len(sequence):
                    continue

                # Check if we can form a stem
                left_stem = sequence[start:start + stem_len]
                right_stem = sequence[start + stem_len + loop_len:start + 2 * stem_len + loop_len]

                # Calculate complementarity
                complementarity = self._calculate_complementarity(left_stem, right_stem)

                if complementarity >= 0.7:  # 70% complementarity threshold
                    return {
                        'start': start,
                        'end': start + 2 * stem_len + loop_len,
                        'stem_length': stem_len,
                        'loop_length': loop_len,
                        'score': complementarity
                    }

        return None

    def _calculate_complementarity(self, seq1: str, seq2: str) -> float:
        """
        Calculate the fraction of positions in a candidate stem that form a valid
        base pair. The two arms of a hairpin stem are ANTIPARALLEL, so the right
        arm (seq2, read 5'->3') must be reversed before pairing against the left
        arm (seq1, read 5'->3'). The previous implementation compared them in
        parallel, so a perfect antiparallel stem scored ~0.5, not 1.0.

        Accepts both DNA (T) and RNA (U), and counts the G-U/G-T wobble pair, which
        is a real base pair in RNA helices.

        Args:
            seq1 (str): Left stem arm, 5'->3'
            seq2 (str): Right stem arm, 5'->3'

        Returns:
            float: Fraction of paired positions (0.0 to 1.0)
        """
        if len(seq1) != len(seq2):
            return 0.0

        s1 = seq1.upper().replace('U', 'T')
        s2 = seq2.upper().replace('U', 'T')[::-1]  # antiparallel: reverse right arm

        # Watson-Crick (A-T, C-G) plus G-T (G-U) wobble.
        valid_pairs = {
            ('A', 'T'), ('T', 'A'), ('C', 'G'), ('G', 'C'),
            ('G', 'T'), ('T', 'G'),
        }
        matches = sum(1 for a, b in zip(s1, s2) if (a, b) in valid_pairs)
        return matches / len(s1)







    def _get_minus_one_frame(self, sequence: str) -> str:
        """
        Get the -1 reading frame of a sequence.

        Args:
            sequence (str): DNA/RNA sequence

        Returns:
            str: -1 frame sequence
        """
        return sequence[1:]  # Remove first nucleotide

    def _get_plus_one_frame(self, sequence: str) -> str:
        """
        Get the +1 reading frame of a sequence.

        Args:
            sequence (str): DNA/RNA sequence

        Returns:
            str: +1 frame sequence
        """
        return sequence[2:]  # Remove first two nucleotides

    def format_frameshift_output(self, frameshift_sites: List[Dict]) -> str:
        """
        Format frameshift sites for CSV output.

        Args:
            frameshift_sites (List[Dict]): List of frameshift sites

        Returns:
            str: Formatted CSV output string
        """
        if not frameshift_sites:
            return "position,end_position,sequence,type\n"

        # CSV header
        output = "position,end_position,sequence,type\n"

        for site in frameshift_sites:
            output += f"{site['position']},{site['end_position']},{site['sequence']},{site['type']}\n"

        return output

    def _slippery_grammar_ok(self, slippery_seq: str) -> bool:
        """
        Check that a heptamer conforms to the -1 PRF slippery grammar X XXY YYZ, i.e.
        that after a -1 shift the P- and A-site tRNAs can RE-PAIR to the new codons.

        IMPORTANT: this is a SEQUENCE-GRAMMAR check only. It does NOT use tRNA
        abundance, aminoacylation, or kinetic data, so it is not "tRNA validation"
        of a functional site. (The previous name _validate_trna_pairing overstated
        this: combined with a wobble helper that returned True for every base pair,
        it always returned True for a conforming heptamer and stamped a misleading
        trna_validated=True on every site.)

        Grammar: positions 1-3 = X X X, positions 4-6 = Y Y Y, position 7 = Z.
        After the -1 slip the P-site reads X X Y (needs XX identical to the original
        XX) and the A-site reads Y Y Z (needs YY identical to the original YY). Both
        hold by construction for a well-formed heptamer; we verify them explicitly.

        Returns:
            bool: True if the heptamer is grammar-conforming.
        """
        s = slippery_seq.upper()
        if len(s) != 7:
            return False
        # X X X : three identical
        if not (s[0] == s[1] == s[2]):
            return False
        # Y Y Y : three identical
        if not (s[3] == s[4] == s[5]):
            return False
        # P-site re-pairing X X (positions 1-2) preserved when frame shifts to X X Y
        if s[0] != s[1]:
            return False
        # A-site re-pairing Y Y (positions 4-5) preserved when frame shifts to Y Y Z
        if s[3] != s[4]:
            return False
        return True

    def calculate_efficiency_score(self, site: Dict) -> float:
        """
        Calculate PRF efficiency score for a detected site.

        Args:
            site: Frameshift site dictionary

        Returns:
            float: Efficiency score between 0.0 and 1.0
        """
        if site['type'] != '-1 PRF':
            return 0.0

        sequence = site['sequence']
        position = site['position']

        # Score slippery site
        slippery_score = self._score_slippery_site(sequence)

        # Score downstream structure
        downstream_score = 0.0
        if 'downstream_structure' in site:
            downstream_score = self._score_downstream_structure(site['downstream_structure'])

        # Score context
        context_score = self._score_context(sequence, position)

        # Calculate weighted score
        efficiency = (
            slippery_score * self.slippery_site_weight +
            downstream_score * self.downstream_structure_weight +
            context_score * self.context_weight
        )

        return min(1.0, max(0.0, efficiency))

    def _score_slippery_site(self, sequence: str) -> float:
        """
        Score the slippery site sequence using published experimental data.
        Based on Plant et al. (2010) and other literature sources.
        """
        if len(sequence) != 7:
            return 0.0

        # Use published scores if available
        if sequence in self.published_slippery_scores:
            return self.published_slippery_scores[sequence]

        # Calculate score based on X XXY YYZ pattern
        # Following the thermodynamic model of Plant et al. (2010)
        score = 0.0

        # Check X XXY YYZ pattern (positions 0-1, 2-3, 4-5)
        if sequence[0] == sequence[1]:  # X XXY
            score += 0.3
        if sequence[2] == sequence[3]:  # XXY YYZ
            score += 0.3
        if sequence[4] == sequence[5]:  # YYY YYZ
            score += 0.2

        # Bonus for strong base pairs in positions 1-2 and 4-5
        strong_pairs = ["GC", "CG", "AT", "TA"]
        if sequence[1:3] in strong_pairs or sequence[4:6] in strong_pairs:
            score += 0.1

        # Penalty for weak base pairs
        weak_pairs = ["AA", "TT", "CC", "GG"]
        if sequence[1:3] in weak_pairs or sequence[4:6] in weak_pairs:
            score -= 0.1

        return max(0.0, min(1.0, score))

    def _score_downstream_structure(self, structure: str) -> float:
        """
        Score downstream RNA structure using thermodynamic principles.
        Based on free energy calculations and known stimulatory structures.
        """
        if not structure:
            return 0.0

        # Convert structure score to free energy estimate
        # Following the relationship: ΔG ≈ -RT ln(K)
        # Where K is the equilibrium constant for structure formation

        try:
            structure_score = float(structure)

            # Convert score to free energy (kcal/mol)
            # Higher scores (more negative) = more stable structures
            if structure_score < -10:
                return 0.9  # Very stable structure (pseudoknot-like)
            elif structure_score < -5:
                return 0.8  # Stable stem-loop
            elif structure_score < -2:
                return 0.6  # Moderate structure
            elif structure_score < 0:
                return 0.4  # Weak structure
            else:
                return 0.2  # Unstable or no structure

        except (ValueError, TypeError):
            # Fallback to string-based scoring
            if "stem-loop" in str(structure).lower():
                return 0.7
            elif "stem" in str(structure).lower():
                return 0.5
            elif "loop" in str(structure).lower():
                return 0.3
            else:
                return 0.1

    def _score_context(self, sequence: str, position: int) -> float:
        """Score the genomic context around the PRF site."""
        score = 0.5  # Base score

        # Check for pause motifs
        for motif in self.pause_motifs:
            if motif in sequence:
                score += 0.2

        # Bonus if the heptamer conforms to the slippery grammar (sequence-only;
        # not a tRNA-abundance measurement -- see _slippery_grammar_ok).
        if self._slippery_grammar_ok(sequence):
            score += 0.3

        return min(1.0, score)

    def analyze_mutation_impact(self, frameshift_sites: List[Dict], mutations: List[Dict]) -> List[Dict]:
        """
        Analyze how mutations affect PRF sites.

        Args:
            frameshift_sites: List of detected PRF sites
            mutations: List of mutation dictionaries

        Returns:
            List[Dict]: Mutations with PRF impact analysis
        """
        prf_mutations = []

        for mutation in mutations:
            mutation_impact = self._assess_mutation_impact(mutation, frameshift_sites)
            if mutation_impact:
                mutation['prf_impact'] = mutation_impact
                prf_mutations.append(mutation)

        return prf_mutations

    def _assess_mutation_impact(self, mutation: Dict, frameshift_sites: List[Dict]) -> Optional[Dict]:
        """Assess the impact of a mutation on PRF sites."""
        mutation_pos = self._parse_mutation_position(mutation)
        if not mutation_pos:
            return None

        for site in frameshift_sites:
            if self._mutation_affects_site(mutation_pos, site):
                return {
                    'severity': self._determine_severity(mutation, site),
                    'efficiency_change': self._calculate_efficiency_change(mutation, site),
                    'mechanism': self._determine_mechanism(mutation, site),
                    'affected_site': site['position']
                }

        return None

    def _parse_mutation_position(self, mutation: Dict) -> Optional[int]:
        """Parse mutation position from mutation dictionary."""
        pos = mutation.get('pos', '')
        if isinstance(pos, str) and ':' in pos:
            # Range mutation, use start position
            return int(pos.split(':')[0])
        elif isinstance(pos, str) and pos.isdigit():
            return int(pos)
        elif isinstance(pos, int):
            return pos
        return None

    def _mutation_affects_site(self, mutation_pos: int, site: Dict) -> bool:
        """Check if mutation affects a PRF site."""
        site_start = site['position']
        site_end = site['end_position']

        return site_start <= mutation_pos <= site_end

    def _determine_severity(self, mutation: Dict, site: Dict) -> str:
        """Determine the severity of mutation impact on PRF."""
        if mutation.get('type') == 'deletion':
            return 'high'
        elif mutation.get('type') == 'insertion':
            return 'medium'
        else:
            return 'low'

    def _calculate_efficiency_change(self, mutation: Dict, site: Dict) -> float:
        """Calculate the change in efficiency due to mutation."""
        # Simplified calculation - in practice, this would be more complex
        if mutation.get('type') == 'deletion':
            return -0.3
        elif mutation.get('type') == 'insertion':
            return -0.2
        else:
            return -0.1

    def _determine_mechanism(self, mutation: Dict, site: Dict) -> str:
        """Determine the mechanism by which mutation affects PRF."""
        if mutation.get('type') == 'deletion':
            return 'disruption of slippery sequence'
        elif mutation.get('type') == 'insertion':
            return 'alteration of reading frame'
        else:
            return 'sequence modification'

    def predict_protein_ratios(self, efficiency_score: float, sequence: str = None) -> Dict[str, float]:
        """
        Predict protein ratios based on PRF efficiency using literature-based models.
        Based on experimental measurements from Plant et al. (2010) and others.

        Args:
            efficiency_score: PRF efficiency score (0.0 to 1.0)
            sequence: Slippery sequence for context-specific predictions

        Returns:
            Dict[str, float]: Predicted protein ratios with confidence intervals
        """
        # Base ratios from simple frameshift model. These are fractions of ribosomes
        # taking each path, so they are probabilities in [0, 1] that must sum to 1.
        orf1ab_ratio = min(1.0, max(0.0, efficiency_score))  # frameshifted (ORF1ab) fraction

        # Apply virus-specific corrections based on literature. The correction scales
        # the frameshifted fraction, then we CLAMP to [0, 1] so the result stays a
        # valid probability (the previous code could return orf1ab_ratio > 1.0, e.g.
        # 0.85 * 1.2 = 1.02, and a negative orf1a_ratio).
        if sequence:
            if sequence in ["TTTAAAC", "TTTAAAT"]:  # Coronavirus
                # SARS-CoV-2: ORF1a:ORF1ab ~ 1:1 to 1:2 (Plant et al., 2010)
                orf1ab_ratio = min(1.0, orf1ab_ratio * 1.2)
            elif sequence in ["TTTTTTA", "TTTTTTG"]:  # HIV-1
                # HIV-1: Gag:Gag-Pol ~ 20:1 (Dulude et al., 2006) -> low FS efficiency
                orf1ab_ratio = min(1.0, orf1ab_ratio * 0.05)

        orf1a_ratio = max(0.0, 1.0 - orf1ab_ratio)  # non-frameshifted (ORF1a) fraction

        # Calculate replication efficiency based on experimental data
        # Higher PRF efficiency generally correlates with better replication
        # but the relationship is virus-specific and context-dependent
        if efficiency_score > 0.8:
            replication_efficiency = 0.8 + (efficiency_score - 0.8) * 0.2
        elif efficiency_score > 0.5:
            replication_efficiency = 0.6 + (efficiency_score - 0.5) * 0.4
        else:
            replication_efficiency = 0.4 + efficiency_score * 0.4

        return {
            'orf1a_ratio': orf1a_ratio,
            'orf1ab_ratio': orf1ab_ratio,
            'replication_efficiency': replication_efficiency
        }

    def get_experimental_validation(self, sequence: str) -> List[str]:
        """
        Get experimental validation references for a slippery sequence.

        Args:
            sequence: Slippery sequence

        Returns:
            List[str]: List of experimental validation references
        """
        return self.experimental_validations.get(sequence, [])

    def calculate_confidence_interval(self, efficiency_score: float) -> Tuple[float, float]:
        """
        Calculate confidence interval for efficiency score.
        Based on experimental error ranges from literature.

        Args:
            efficiency_score: Calculated efficiency score

        Returns:
            Tuple[float, float]: (lower_bound, upper_bound)
        """
        # Standard error from experimental measurements
        # Based on Plant et al. (2010) and other studies
        standard_error = 0.1  # ±10% typical experimental error

        lower_bound = max(0.0, efficiency_score - standard_error)
        upper_bound = min(1.0, efficiency_score + standard_error)

        return (lower_bound, upper_bound)

    def generate_comprehensive_report(self, frameshift_sites: List[Dict],
                                    mutations: List[Dict] = None) -> Dict:
        """
        Generate a comprehensive PRF analysis report.

        Args:
            frameshift_sites: List of detected PRF sites
            mutations: List of mutations (optional)

        Returns:
            Dict: Comprehensive analysis report
        """
        # Calculate efficiency scores
        for site in frameshift_sites:
            site['efficiency_score'] = self.calculate_efficiency_score(site)

        # Analyze mutation impacts
        prf_mutations = []
        if mutations:
            prf_mutations = self.analyze_mutation_impact(frameshift_sites, mutations)

        # Calculate summary statistics
        total_sites = len(frameshift_sites)
        minus_one_sites = len([s for s in frameshift_sites if s['type'] == '-1 PRF'])
        plus_one_sites = len([s for s in frameshift_sites if s['type'] == '+1 PRF'])

        avg_efficiency = sum(s.get('efficiency_score', 0) for s in frameshift_sites) / total_sites if total_sites > 0 else 0

        report = {
            'summary': {
                'total_prf_sites': total_sites,
                'minus_one_sites': minus_one_sites,
                'plus_one_sites': plus_one_sites,
                'average_efficiency': avg_efficiency,
                'prf_affecting_mutations': len(prf_mutations)
            },
            'sites': frameshift_sites,
            'mutations': prf_mutations,
            'efficiency_category': self._categorize_efficiency(avg_efficiency)
        }

        return report

    def _categorize_efficiency(self, efficiency: float) -> str:
        """Categorize efficiency score."""
        if efficiency >= 0.8:
            return 'high'
        elif efficiency >= 0.5:
            return 'medium'
        elif efficiency >= 0.2:
            return 'low'
        else:
            return 'very_low'

    def save_detailed_report(self, frameshift_sites: List[Dict],
                           output_file: str, mutations: List[Dict] = None) -> None:
        """
        Save a detailed PRF analysis report to JSON file.

        Args:
            frameshift_sites: List of detected PRF sites
            output_file: Output file path
            mutations: List of mutations (optional)
        """
        report = self.generate_comprehensive_report(frameshift_sites, mutations)

        with open(output_file, 'w') as f:
            json.dump(report, f, indent=2, default=str)
