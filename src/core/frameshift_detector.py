"""
Author: Rui Wang
Date: 2025-01-07
LastModifiedBy: Rui Wang
LastEditTime: 2025-01-07
Email: rw3594@nyu.edu
FilePath: /muvitrack/src/core/frameshift_detector.py
Description: 
Frameshift detection module for MuViTrack. Detects potential programmed ribosomal 
frameshifting (PRF) sites including +1 and -1 frameshifting mechanisms.
"""

from typing import Dict, List, Tuple, Optional
from Bio.Seq import Seq
from Bio.SeqUtils import nt_search
import re


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
        
    def detect_frameshift_sites(self, sequence: str, start_pos: int = 0) -> List[Dict]:
        """
        Detect potential frameshifting sites in a sequence.
        
        Args:
            sequence (str): DNA/RNA sequence to analyze
            start_pos (int): Starting position for coordinate mapping
            
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
        
        for slippery_seq, description in self.slippery_sequences.items():
            # Find all occurrences of the slippery sequence
            positions = self._find_all_occurrences(sequence, slippery_seq)
            
            for pos in positions:
                # Validate tRNA pairing compatibility (CRITICAL for biological accuracy)
                if not self._validate_trna_pairing(slippery_seq):
                    continue  # Skip if tRNA pairing is incompatible
                
                # Check if this position has a downstream stem-loop
                stem_loop_score = self._check_downstream_stem_loop(sequence, pos + len(slippery_seq))
                
                site = {
                    'position': pos + start_pos,
                    'end_position': pos + len(slippery_seq) + start_pos,
                    'sequence': slippery_seq,
                    'type': '-1 PRF',
                    'description': description,
                    'stem_loop_score': stem_loop_score,
                    'mechanism': 'slippery_sequence',
                    'trna_validated': True  # Mark as tRNA-validated
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
        
        # Find stop codons that might be shifty
        stop_codons = ['TAA', 'TAG', 'TGA']
        
        for i in range(len(sequence) - 2):
            codon = sequence[i:i+3]
            if codon in stop_codons:
                # Check for shifty context
                context = sequence[max(0, i-1):i+4]
                
                for signal, description in self.plus_one_signals.items():
                    if signal in context:
                        # Analyze the impact
                        site = {
                            'position': i + start_pos,
                            'end_position': i + 3 + start_pos,
                            'sequence': codon,
                            'context': context,
                            'type': '+1 PRF',
                            'description': description,
                            'mechanism': 'shifty_stop'
                        }
                        sites.append(site)
                        break
        
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
        optimal_end = position + 9
        
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
        Calculate complementarity between two sequences.
        
        Args:
            seq1 (str): First sequence (DNA format ATGC)
            seq2 (str): Second sequence (DNA format ATGC)
            
        Returns:
            float: Complementarity score (0.0 to 1.0)
        """
        if len(seq1) != len(seq2):
            return 0.0
        
        # DNA base pairing rules (A-T, C-G)
        complement_pairs = {
            'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'
        }
        
        matches = 0
        for i in range(len(seq1)):
            if seq2[i] == complement_pairs.get(seq1[i], ''):
                matches += 1
        
        return matches / len(seq1)
    

    

    

    
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

    def _validate_trna_pairing(self, slippery_seq: str) -> bool:
        """
        Validate tRNA pairing compatibility for -1 PRF.
        
        Args:
            slippery_seq (str): Slippery sequence (7 nt)
            
        Returns:
            bool: True if tRNA pairing is compatible
        """
        if len(slippery_seq) != 7:
            return False
            
        # Split into P-site (positions 1-3) and A-site (positions 4-6)
        # Position 7 is the spacer
        p_site = slippery_seq[0:3]  # XXX
        a_site = slippery_seq[3:6]  # YYY
        spacer = slippery_seq[6]    # Z
        
        # After -1 shift:
        # P-site becomes positions 2-4 (XXY)
        # A-site becomes positions 5-7 (YYZ)
        p_site_shifted = slippery_seq[1:4]  # XXY
        a_site_shifted = slippery_seq[4:7]  # YYZ
        
        # Check P-site compatibility
        # Original: XXX, Shifted: XXY
        # First two positions must be identical for tRNA re-pairing
        if p_site[0] != p_site_shifted[0] or p_site[1] != p_site_shifted[1]:
            return False
            
        # Check A-site compatibility  
        # Original: YYY, Shifted: YYZ
        # First two positions must be identical for tRNA re-pairing
        if a_site[0] != a_site_shifted[0] or a_site[1] != a_site_shifted[1]:
            return False
            
        # Check third position wobble compatibility for A-site
        # Original third: Y, Shifted third: Z
        # These should be wobble-compatible
        wobble_compatible = self._check_wobble_compatibility(a_site[2], a_site_shifted[2])
        
        return wobble_compatible
    
    def _check_wobble_compatibility(self, base1: str, base2: str) -> bool:
        """
        Check if two bases are wobble-compatible.
        
        Args:
            base1 (str): First base
            base2 (str): Second base
            
        Returns:
            bool: True if wobble-compatible
        """
        # Wobble pairing rules (including non-canonical pairs observed in viral PRF)
        wobble_pairs = {
            'G': ['U', 'C', 'T', 'A'],  # G can pair with U, C, T, A (observed in viral PRF)
            'U': ['A', 'G', 'T'],       # U can pair with A, G, T
            'T': ['A', 'G', 'U'],       # T can pair with A, G, U (DNA format)
            'I': ['A', 'U', 'C'],       # Inosine (I) can pair with A, U, or C
            'A': ['C', 'U', 'T', 'G'],  # A can pair with C, U, T, G (observed in viral PRF)
            'C': ['A', 'G', 'U']        # C can pair with A, G, U
        }
        
        # Check if base1 can wobble-pair with base2
        if base1 in wobble_pairs and base2 in wobble_pairs[base1]:
            return True
            
        # Check reverse direction
        if base2 in wobble_pairs and base1 in wobble_pairs[base2]:
            return True
            
        # Special case: identical bases are always compatible
        if base1 == base2:
            return True
            
        return False
