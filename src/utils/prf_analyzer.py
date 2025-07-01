"""
Programmed Ribosomal Frameshifting (PRF) Analyzer for SARS-CoV-2.

This module provides tools for analyzing the programmed ribosomal frameshifting
mechanism in SARS-CoV-2, including site detection, efficiency calculation,
mutation impact assessment, and RNA structure prediction.
"""

from dataclasses import dataclass
from enum import Enum
from typing import Dict, List, Optional


class PRFSiteType(Enum):
    """Types of PRF sites."""

    CANONICAL = "canonical"
    VARIANT = "variant"
    ABSENT = "absent"
    REFERENCE_VARIANT = "reference_variant"  # NC_045512.2 variant


@dataclass
class PRFSite:
    """Data class representing a PRF site."""

    position: int
    sequence: str
    site_type: PRFSiteType
    efficiency_score: float
    downstream_structure: Optional[str] = None
    mutations: List[Dict] = None
    notes: Optional[str] = None

    def __post_init__(self):
        if self.mutations is None:
            self.mutations = []


class PRFAnalyzer:
    """
    Analyzer for Programmed Ribosomal Frameshifting in SARS-CoV-2.

    This class provides comprehensive tools for analyzing PRF sites,
    calculating efficiency, and assessing the impact of mutations.
    """

    # SARS-CoV-2 PRF site constants
    CANONICAL_PRF_POSITION = 13468
    CANONICAL_SLIPPERY_SITE = "UUUAAAC"
    REFERENCE_SLIPPERY_SITE = "CGGGTTT"  # Actual sequence in NC_045512.2
    DOWNSTREAM_REGION_LENGTH = 50

    # Efficiency scoring weights
    SLIPPERY_SITE_WEIGHT = 0.6
    DOWNSTREAM_STRUCTURE_WEIGHT = 0.3
    CONTEXT_WEIGHT = 0.1

    def __init__(self, ref_genome: str):
        """
        Initialize the PRF analyzer.

        Args:
            ref_genome: Reference genome sequence
        """
        self.ref_genome = ref_genome.upper()
        self._validate_genome()

    def _validate_genome(self) -> None:
        """Validate the reference genome sequence."""
        if (
            len(self.ref_genome)
            < self.CANONICAL_PRF_POSITION + self.DOWNSTREAM_REGION_LENGTH
        ):
            raise ValueError(
                f"Genome sequence too short. Expected at least "
                f"{self.CANONICAL_PRF_POSITION + self.DOWNSTREAM_REGION_LENGTH} nucleotides"
            )

        # Check for valid nucleotides
        invalid_chars = set(self.ref_genome) - set("ATCGU")
        if invalid_chars:
            raise ValueError(f"Invalid nucleotides found: {invalid_chars}")

    def detect_prf_site(self, position: Optional[int] = None) -> PRFSite:
        """
        Detect the programmed ribosomal frameshifting site.

        Args:
            position: Position to analyze (default: canonical PRF position)

        Returns:
            PRFSite object with analysis results
        """
        if position is None:
            position = self.CANONICAL_PRF_POSITION

        # Extract sequences
        slippery_site = self.ref_genome[position - 1 : position + 6]
        downstream_seq = self.ref_genome[
            position + 6 : position + 6 + self.DOWNSTREAM_REGION_LENGTH
        ]

        # Determine site type
        if slippery_site == self.CANONICAL_SLIPPERY_SITE:
            site_type = PRFSiteType.CANONICAL
            notes = "Optimal slippery site for frameshifting"
        elif slippery_site == self.REFERENCE_SLIPPERY_SITE:
            site_type = PRFSiteType.REFERENCE_VARIANT
            notes = "Reference genome variant (NC_045512.2) - reduced efficiency"
        elif self._is_slippery_site_variant(slippery_site):
            site_type = PRFSiteType.VARIANT
            notes = "Known slippery site variant"
        else:
            site_type = PRFSiteType.ABSENT
            notes = "No recognizable slippery site"

        # Calculate efficiency
        efficiency_score = self._calculate_prf_efficiency(slippery_site, downstream_seq)

        # Predict downstream structure
        downstream_structure = self._predict_downstream_structure(downstream_seq)

        return PRFSite(
            position=position,
            sequence=slippery_site,
            site_type=site_type,
            efficiency_score=efficiency_score,
            downstream_structure=downstream_structure,
            notes=notes,
        )

    def _is_slippery_site_variant(self, sequence: str) -> bool:
        """
        Check if a sequence is a variant of the canonical slippery site.

        Args:
            sequence: 7-nucleotide sequence to check

        Returns:
            True if it's a known slippery site variant
        """
        # Known slippery site variants that can still support frameshifting
        slippery_variants = {
            "UUUAAAC",  # Canonical
            "CGGGTTT",  # Reference genome variant (NC_045512.2)
            "UUUAAAA",  # Common variant
            "UUUAAAU",  # Common variant
            "UUUAAAG",  # Less common
            "UUUAAUU",  # Less common
        }

        return sequence in slippery_variants

    def _calculate_prf_efficiency(
        self, slippery_site: str, downstream_seq: str
    ) -> float:
        """
        Calculate predicted PRF efficiency based on sequence features.

        Args:
            slippery_site: 7-nucleotide slippery site sequence
            downstream_seq: Downstream sequence for structure analysis

        Returns:
            Efficiency score between 0.0 and 1.0
        """
        # Slippery site score
        slippery_score = self._score_slippery_site(slippery_site)

        # Downstream structure score
        structure_score = self._score_downstream_structure(downstream_seq)

        # Context score
        context_score = self._score_context(slippery_site, downstream_seq)

        # Weighted combination
        total_score = (
            slippery_score * self.SLIPPERY_SITE_WEIGHT
            + structure_score * self.DOWNSTREAM_STRUCTURE_WEIGHT
            + context_score * self.CONTEXT_WEIGHT
        )

        return min(1.0, max(0.0, total_score))

    def _score_slippery_site(self, sequence: str) -> float:
        """
        Score the slippery site based on conservation and known variants.

        Args:
            sequence: 7-nucleotide slippery site sequence

        Returns:
            Score between 0.0 and 1.0
        """
        if sequence == self.CANONICAL_SLIPPERY_SITE:
            return 1.0

        # Score known variants
        variant_scores = {
            "CGGGTTT": 0.3,  # Reference genome variant - low efficiency
            "UUUAAAA": 0.9,
            "UUUAAAU": 0.85,
            "UUUAAAG": 0.7,
            "UUUAAUU": 0.6,
        }

        return variant_scores.get(sequence, 0.0)

    def _score_downstream_structure(self, sequence: str) -> float:
        """
        Score downstream sequence for potential structure formation.

        Args:
            sequence: Downstream sequence

        Returns:
            Score between 0.0 and 1.0
        """
        # Simple heuristics for structure prediction
        # In practice, you'd use RNAfold or similar tools

        # Check for GC content (higher GC = more stable structure)
        gc_count = sequence.count("G") + sequence.count("C")
        gc_content = gc_count / len(sequence)

        # Check for potential stem-loop patterns
        stem_loop_score = self._detect_stem_loop_patterns(sequence)

        # Combine scores
        structure_score = gc_content * 0.6 + stem_loop_score * 0.4

        return min(1.0, structure_score)

    def _detect_stem_loop_patterns(self, sequence: str) -> float:
        """
        Detect potential stem-loop patterns in the sequence.

        Args:
            sequence: Sequence to analyze

        Returns:
            Score between 0.0 and 1.0
        """
        # Simple pattern detection
        # Look for complementary regions that could form stems

        score = 0.0

        # Check for palindromic sequences
        for length in range(6, 12):
            for i in range(len(sequence) - length):
                segment = sequence[i : i + length]
                if self._is_palindromic(segment):
                    score += 0.2

        # Check for GC-rich regions
        gc_rich_count = 0
        for i in range(0, len(sequence) - 5, 5):
            window = sequence[i : i + 5]
            if window.count("G") + window.count("C") >= 3:
                gc_rich_count += 1

        score += (gc_rich_count / (len(sequence) // 5)) * 0.3

        return min(1.0, score)

    def _is_palindromic(self, sequence: str) -> bool:
        """
        Check if a sequence is palindromic (can form hairpin).

        Args:
            sequence: Sequence to check

        Returns:
            True if palindromic
        """
        complement_map = {"A": "U", "U": "A", "G": "C", "C": "G"}
        complement = "".join(complement_map.get(base, base) for base in sequence)
        return sequence == complement[::-1]

    def _score_context(self, slippery_site: str, downstream_seq: str) -> float:
        """
        Score the context around the PRF site.

        Args:
            slippery_site: Slippery site sequence
            downstream_seq: Downstream sequence

        Returns:
            Score between 0.0 and 1.0
        """
        # Check for optimal spacing and context
        score = 0.5  # Base score

        # Optimal spacing after slippery site
        if len(downstream_seq) >= 10:
            # Check for optimal spacing patterns
            if downstream_seq[0:3] in ["GGC", "GGA", "GGU"]:
                score += 0.2

        # Check for stop codons too close (would terminate translation)
        for i in range(0, min(30, len(downstream_seq)), 3):
            codon = downstream_seq[i : i + 3]
            if codon in ["UAA", "UAG", "UGA"]:
                score -= 0.3
                break

        return max(0.0, min(1.0, score))

    def _predict_downstream_structure(self, sequence: str) -> str:
        """
        Predict RNA secondary structure for downstream sequence.

        Args:
            sequence: Downstream sequence

        Returns:
            Structure prediction (simplified)
        """
        # This is a simplified prediction
        # In practice, use RNAfold or ViennaRNA

        structure = ""
        for i, _base in enumerate(sequence):
            if i < len(sequence) // 2:
                structure += "("  # Potential stem formation
            else:
                structure += ")"  # Potential stem formation

        return structure

    def analyze_prf_mutations(self, mutations: List[Dict]) -> List[Dict]:
        """
        Analyze mutations that might affect PRF efficiency.

        Args:
            mutations: List of mutation dictionaries

        Returns:
            List of mutations with PRF impact analysis
        """
        prf_region_start = self.CANONICAL_PRF_POSITION - 10
        prf_region_end = self.CANONICAL_PRF_POSITION + 20

        prf_mutations = []

        for mutation in mutations:
            # Check if mutation is in PRF region
            if self._mutation_in_prf_region(mutation, prf_region_start, prf_region_end):
                mutation_copy = mutation.copy()
                mutation_copy["affects_prf"] = True
                mutation_copy["prf_impact"] = self._assess_prf_impact(mutation)
                mutation_copy["prf_region"] = (
                    "slippery_site"
                    if self._in_slippery_site(mutation)
                    else "downstream"
                )
                prf_mutations.append(mutation_copy)

        return prf_mutations

    def _mutation_in_prf_region(self, mutation: Dict, start: int, end: int) -> bool:
        """
        Check if a mutation is in the PRF region.

        Args:
            mutation: Mutation dictionary
            start: Start of PRF region
            end: End of PRF region

        Returns:
            True if mutation is in PRF region
        """
        if "pos" not in mutation:
            return False

        pos_str = str(mutation["pos"])
        if ":" in pos_str:
            # Range mutation
            pos_start, pos_end = map(int, pos_str.split(":"))
            return pos_start <= end and pos_end >= start
        else:
            # Point mutation
            pos = int(pos_str)
            return start <= pos <= end

    def _in_slippery_site(self, mutation: Dict) -> bool:
        """
        Check if mutation is in the slippery site.

        Args:
            mutation: Mutation dictionary

        Returns:
            True if in slippery site
        """
        slippery_start = self.CANONICAL_PRF_POSITION
        slippery_end = self.CANONICAL_PRF_POSITION + 6

        if "pos" not in mutation:
            return False

        pos_str = str(mutation["pos"])
        if ":" in pos_str:
            pos_start, pos_end = map(int, pos_str.split(":"))
            return pos_start <= slippery_end and pos_end >= slippery_start
        else:
            pos = int(pos_str)
            return slippery_start <= pos <= slippery_end

    def _assess_prf_impact(self, mutation: Dict) -> Dict:
        """
        Assess the impact of a mutation on PRF efficiency.

        Args:
            mutation: Mutation dictionary

        Returns:
            Impact assessment dictionary
        """
        impact = {
            "severity": "unknown",
            "efficiency_change": 0.0,
            "mechanism": "unknown",
        }

        # Analyze based on mutation type and position
        if mutation.get("type") == "deletion":
            impact["severity"] = "high"
            impact["efficiency_change"] = -0.5
            impact["mechanism"] = "sequence disruption"
        elif mutation.get("type") == "insertion":
            impact["severity"] = "high"
            impact["efficiency_change"] = -0.4
            impact["mechanism"] = "spacing disruption"
        elif mutation.get("type") == "pointMutation":
            impact["severity"] = "moderate"
            impact["efficiency_change"] = -0.2
            impact["mechanism"] = "sequence alteration"

        return impact

    def predict_protein_ratios(self, prf_efficiency: float) -> Dict[str, float]:
        """
        Predict ORF1a:ORF1ab protein ratios based on PRF efficiency.

        Args:
            prf_efficiency: PRF efficiency (0.0 to 1.0)

        Returns:
            Dictionary with protein ratios
        """
        orf1a_ratio = 1.0 - prf_efficiency
        orf1ab_ratio = prf_efficiency

        return {
            "orf1a_ratio": orf1a_ratio,
            "orf1ab_ratio": orf1ab_ratio,
            "rdrp_availability": orf1ab_ratio,  # RdRp is in ORF1ab
            "replication_efficiency": min(1.0, orf1ab_ratio * 2.0),  # Simplified model
        }

    def generate_prf_report(
        self, prf_site: PRFSite, mutations: List[Dict] = None
    ) -> Dict:
        """
        Generate a comprehensive PRF analysis report.

        Args:
            prf_site: PRF site analysis results
            mutations: Optional list of mutations

        Returns:
            Comprehensive PRF report
        """
        if mutations is None:
            mutations = []

        prf_mutations = self.analyze_prf_mutations(mutations)
        protein_ratios = self.predict_protein_ratios(prf_site.efficiency_score)

        report = {
            "prf_site": {
                "position": prf_site.position,
                "sequence": prf_site.sequence,
                "site_type": prf_site.site_type.value,
                "efficiency_score": prf_site.efficiency_score,
                "downstream_structure": prf_site.downstream_structure,
                "notes": prf_site.notes,
            },
            "protein_ratios": protein_ratios,
            "prf_mutations": prf_mutations,
            "summary": {
                "total_mutations": len(mutations),
                "prf_affecting_mutations": len(prf_mutations),
                "efficiency_category": self._categorize_efficiency(
                    prf_site.efficiency_score
                ),
            },
        }

        return report

    def _categorize_efficiency(self, efficiency: float) -> str:
        """
        Categorize PRF efficiency.

        Args:
            efficiency: Efficiency score

        Returns:
            Efficiency category
        """
        if efficiency >= 0.8:
            return "high"
        elif efficiency >= 0.5:
            return "moderate"
        elif efficiency >= 0.2:
            return "low"
        else:
            return "very_low"


def analyze_prf_in_sequence(sequence: str, mutations: List[Dict] = None) -> Dict:
    """
    Convenience function to analyze PRF in a sequence.

    Args:
        sequence: Genome sequence
        mutations: Optional list of mutations

    Returns:
        PRF analysis report
    """
    analyzer = PRFAnalyzer(sequence)
    prf_site = analyzer.detect_prf_site()
    return analyzer.generate_prf_report(prf_site, mutations or [])
 