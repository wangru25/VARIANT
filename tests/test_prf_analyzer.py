"""
Tests for the PRF Analyzer module.
"""

import pytest

from src.utils.prf_analyzer import (
    PRFAnalyzer,
    PRFSite,
    PRFSiteType,
    analyze_prf_in_sequence,
)


class TestPRFAnalyzer:
    """Test cases for PRFAnalyzer class."""

    @pytest.fixture
    def sample_genome(self):
        """Create a sample genome with canonical PRF site."""
        # Create a genome with the canonical PRF site at position 13468
        genome = "A" * 13467 + "UUUAAAC" + "GGC" * 20 + "A" * 1000
        return genome

    @pytest.fixture
    def analyzer(self, sample_genome):
        """Create a PRFAnalyzer instance."""
        return PRFAnalyzer(sample_genome)

    def test_init_valid_genome(self, sample_genome):
        """Test initialization with valid genome."""
        analyzer = PRFAnalyzer(sample_genome)
        assert analyzer.ref_genome == sample_genome

    def test_init_invalid_genome(self):
        """Test initialization with invalid genome."""
        short_genome = "A" * 1000  # Too short
        with pytest.raises(ValueError, match="Genome sequence too short"):
            PRFAnalyzer(short_genome)

        invalid_genome = "A" * 15000 + "X" + "A" * 1000  # Invalid nucleotide
        with pytest.raises(ValueError, match="Invalid nucleotides found"):
            PRFAnalyzer(invalid_genome)

    def test_detect_prf_site_canonical(self, analyzer):
        """Test detection of canonical PRF site."""
        prf_site = analyzer.detect_prf_site()

        assert prf_site.position == 13468
        assert prf_site.sequence == "UUUAAAC"
        assert prf_site.site_type == PRFSiteType.CANONICAL
        assert prf_site.efficiency_score > 0.8

    def test_detect_prf_site_variant(self):
        """Test detection of variant PRF site."""
        variant_genome = "A" * 13467 + "UUUAAAA" + "GGC" * 20 + "A" * 1000
        analyzer = PRFAnalyzer(variant_genome)
        prf_site = analyzer.detect_prf_site()

        assert prf_site.sequence == "UUUAAAA"
        assert prf_site.site_type == PRFSiteType.VARIANT
        assert prf_site.efficiency_score > 0.5

    def test_detect_prf_site_absent(self):
        """Test detection when PRF site is absent."""
        no_prf_genome = "A" * 13467 + "AAAAAAA" + "GGC" * 20 + "A" * 1000
        analyzer = PRFAnalyzer(no_prf_genome)
        prf_site = analyzer.detect_prf_site()

        assert prf_site.site_type == PRFSiteType.ABSENT
        assert prf_site.efficiency_score < 0.3  # Relaxed threshold

    def test_slippery_site_variants(self, analyzer):
        """Test slippery site variant detection."""
        assert analyzer._is_slippery_site_variant("UUUAAAC")  # Canonical
        assert analyzer._is_slippery_site_variant("UUUAAAA")  # Known variant
        assert analyzer._is_slippery_site_variant("UUUAAAU")  # Known variant
        assert not analyzer._is_slippery_site_variant("AAAAAAA")  # Not slippery

    def test_palindromic_detection(self, analyzer):
        """Test palindromic sequence detection."""
        assert analyzer._is_palindromic("GCGCGC")
        assert analyzer._is_palindromic("AUAUAU")
        assert not analyzer._is_palindromic("AAAAAA")
        assert not analyzer._is_palindromic("GCATGC")

    def test_protein_ratios_prediction(self, analyzer):
        """Test protein ratio prediction."""
        ratios = analyzer.predict_protein_ratios(0.25)  # 25% efficiency

        assert ratios["orf1a_ratio"] == 0.75
        assert ratios["orf1ab_ratio"] == 0.25
        assert ratios["rdrp_availability"] == 0.25
        assert ratios["replication_efficiency"] == 0.5

    def test_mutation_analysis(self, analyzer):
        """Test PRF mutation analysis."""
        mutations = [
            {"pos": "13468", "type": "pointMutation", "SNP": "U->A"},
            {"pos": "13470:13472", "type": "deletion", "SNP": "UAA"},
            {
                "pos": "20000",
                "type": "pointMutation",
                "SNP": "A->G",
            },  # Outside PRF region
        ]

        prf_mutations = analyzer.analyze_prf_mutations(mutations)

        assert len(prf_mutations) == 2  # Only mutations in PRF region
        assert all(m["affects_prf"] for m in prf_mutations)
        assert prf_mutations[0]["prf_region"] == "slippery_site"
        assert prf_mutations[1]["prf_region"] == "slippery_site"


class TestPRFSite:
    """Test cases for PRFSite dataclass."""

    def test_prf_site_creation(self):
        """Test PRFSite object creation."""
        prf_site = PRFSite(
            position=13468,
            sequence="UUUAAAC",
            site_type=PRFSiteType.CANONICAL,
            efficiency_score=0.9,
        )

        assert prf_site.position == 13468
        assert prf_site.sequence == "UUUAAAC"
        assert prf_site.site_type == PRFSiteType.CANONICAL
        assert prf_site.efficiency_score == 0.9
        assert prf_site.mutations == []  # Default empty list


class TestConvenienceFunctions:
    """Test cases for convenience functions."""

    def test_analyze_prf_in_sequence(self):
        """Test the convenience function."""
        genome = "A" * 13467 + "UUUAAAC" + "GGC" * 20 + "A" * 1000
        mutations = [{"pos": "13468", "type": "pointMutation", "SNP": "U->A"}]

        report = analyze_prf_in_sequence(genome, mutations)

        assert "prf_site" in report
        assert "protein_ratios" in report
        assert "prf_mutations" in report
        assert "summary" in report
        assert report["prf_site"]["sequence"] == "UUUAAAC"
        assert len(report["prf_mutations"]) == 1


if __name__ == "__main__":
    pytest.main([__file__])
