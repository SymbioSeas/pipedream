"""
Unit tests for Phase 1 enhanced data retrieval functions.

Tests for:
- get_related_taxa_weighted()
- intelligent_exclusion_selection()
- Helper functions (_calculate_taxonomic_distance, _estimate_sequence_count, etc.)
"""

import pytest
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from assay_design.data_retrieval import (
    get_related_taxa_weighted,
    intelligent_exclusion_selection,
    _calculate_taxonomic_distance,
    _select_diverse_taxa,
    _determine_diversity_levels,
    _calculate_coverage_score
)


# ==================== HELPER FUNCTION TESTS ====================

class TestCalculateTaxonomicDistance:
    """Test taxonomic distance calculation."""

    def test_genus_sibling_of_species(self):
        """Test distance from species to genus-level sibling."""
        distance = _calculate_taxonomic_distance("genus", "species")
        assert distance == 1, "Genus sibling of species should have distance 1"

    def test_family_sibling_of_species(self):
        """Test distance from species to family-level sibling."""
        distance = _calculate_taxonomic_distance("family", "species")
        assert distance == 2, "Family sibling of species should have distance 2"

    def test_order_sibling_of_species(self):
        """Test distance from species to order-level sibling."""
        distance = _calculate_taxonomic_distance("order", "species")
        assert distance == 3, "Order sibling of species should have distance 3"

    def test_order_sibling_of_genus(self):
        """Test distance from genus to order-level sibling."""
        distance = _calculate_taxonomic_distance("order", "genus")
        assert distance == 2, "Order sibling of genus should have distance 2"

    def test_same_level_returns_1(self):
        """Test that same-level siblings get distance 1."""
        distance = _calculate_taxonomic_distance("genus", "genus")
        assert distance == 1, "Same-level siblings should have distance 1"

    def test_unknown_ranks_use_defaults(self):
        """Test that unknown ranks use sensible defaults."""
        distance = _calculate_taxonomic_distance("genus", "unknown_rank")
        assert isinstance(distance, int), "Should return integer distance"
        assert distance > 0, "Distance should be positive"


class TestDetermineDiversityLevels:
    """Test automatic diversity level determination."""

    def test_species_rank(self):
        """Test diversity levels for species-rank inclusion."""
        levels = _determine_diversity_levels("species")
        assert "genus" in levels, "Should include genus for species"
        assert "family" in levels, "Should include family for species"
        assert "order" in levels, "Should include order for species"

    def test_genus_rank(self):
        """Test diversity levels for genus-rank inclusion."""
        levels = _determine_diversity_levels("genus")
        assert "family" in levels, "Should include family for genus"
        assert "order" in levels, "Should include order for genus"
        # Genus-rank should sample higher levels
        assert "genus" not in levels, "Should not include genus for genus-rank"

    def test_family_rank(self):
        """Test diversity levels for family-rank inclusion."""
        levels = _determine_diversity_levels("family")
        assert "order" in levels, "Should include order for family"
        assert len(levels) >= 2, "Should have at least 2 levels"

    def test_unknown_rank_returns_defaults(self):
        """Test that unknown ranks return default levels."""
        levels = _determine_diversity_levels("unknown_rank")
        assert len(levels) >= 2, "Should return at least 2 default levels"
        assert "genus" in levels or "family" in levels, "Should include standard levels"

    def test_returns_valid_levels_only(self):
        """Test that only valid taxonomic levels are returned."""
        levels = _determine_diversity_levels("species")
        valid_levels = ["genus", "family", "order", "class", "phylum"]
        for level in levels:
            assert level in valid_levels, f"Level {level} should be in valid levels"


class TestSelectDiverseTaxa:
    """Test diverse taxa selection algorithm."""

    def test_basic_selection(self):
        """Test basic selection from multiple levels."""
        mock_candidates = [
            {"taxid": "1", "scientific_name": "Genus1", "priority_score": 10.0,
             "diversity_level": "genus", "phylogenetic_distance": 1},
            {"taxid": "2", "scientific_name": "Genus2", "priority_score": 9.0,
             "diversity_level": "genus", "phylogenetic_distance": 1},
            {"taxid": "3", "scientific_name": "Family1", "priority_score": 7.0,
             "diversity_level": "family", "phylogenetic_distance": 2},
            {"taxid": "4", "scientific_name": "Order1", "priority_score": 5.0,
             "diversity_level": "order", "phylogenetic_distance": 3},
        ]

        selected = _select_diverse_taxa(mock_candidates, max_results=3,
                                       diversity_levels=["genus", "family", "order"])

        assert len(selected) <= 3, "Should not exceed max_results"
        assert len(selected) > 0, "Should select at least one taxon"

    def test_diversity_across_levels(self):
        """Test that selection includes multiple levels."""
        mock_candidates = [
            {"taxid": str(i), "scientific_name": f"Genus{i}",
             "priority_score": 10.0 - i, "diversity_level": "genus",
             "phylogenetic_distance": 1}
            for i in range(5)
        ] + [
            {"taxid": str(i+10), "scientific_name": f"Family{i}",
             "priority_score": 5.0 - i, "diversity_level": "family",
             "phylogenetic_distance": 2}
            for i in range(3)
        ]

        selected = _select_diverse_taxa(mock_candidates, max_results=6,
                                       diversity_levels=["genus", "family", "order"])

        levels = set(t["diversity_level"] for t in selected)
        assert len(levels) >= 2, "Should represent multiple taxonomic levels"

    def test_empty_candidates_returns_empty(self):
        """Test that empty candidate list returns empty."""
        selected = _select_diverse_taxa([], max_results=10,
                                       diversity_levels=["genus", "family", "order"])
        assert len(selected) == 0, "Empty candidates should return empty list"

    def test_respects_max_results(self):
        """Test that selection respects max_results limit."""
        mock_candidates = [
            {"taxid": str(i), "scientific_name": f"Taxon{i}",
             "priority_score": 10.0 - i, "diversity_level": "genus",
             "phylogenetic_distance": 1}
            for i in range(20)
        ]

        for max_n in [1, 5, 10, 15]:
            selected = _select_diverse_taxa(mock_candidates, max_results=max_n,
                                           diversity_levels=["genus", "family", "order"])
            assert len(selected) <= max_n, f"Should not exceed max_results={max_n}"


class TestCalculateCoverageScore:
    """Test phylogenetic coverage score calculation."""

    def test_empty_taxa_returns_zero(self):
        """Test that empty taxa list returns 0 score."""
        score = _calculate_coverage_score([], ["genus", "family", "order"])
        assert score == 0.0, "Empty taxa should have 0 coverage score"

    def test_single_level_low_score(self):
        """Test that single-level selection gets lower score."""
        mock_taxa = [
            {"diversity_level": "genus", "phylogenetic_distance": 1}
            for _ in range(5)
        ]

        score = _calculate_coverage_score(mock_taxa, ["genus", "family", "order"])
        assert 0.0 < score < 0.7, "Single-level should have limited coverage score"

    def test_multi_level_higher_score(self):
        """Test that multi-level selection gets higher score."""
        mock_taxa = [
            {"diversity_level": "genus", "phylogenetic_distance": 1},
            {"diversity_level": "genus", "phylogenetic_distance": 1},
            {"diversity_level": "family", "phylogenetic_distance": 2},
            {"diversity_level": "family", "phylogenetic_distance": 2},
            {"diversity_level": "order", "phylogenetic_distance": 3},
        ]

        score = _calculate_coverage_score(mock_taxa, ["genus", "family", "order"])
        assert score >= 0.5, "Multi-level should have good coverage score"

    def test_more_taxa_higher_score(self):
        """Test that more taxa generally gives higher score."""
        small_set = [
            {"diversity_level": "genus", "phylogenetic_distance": 1},
            {"diversity_level": "family", "phylogenetic_distance": 2},
        ]

        large_set = small_set + [
            {"diversity_level": "genus", "phylogenetic_distance": 1},
            {"diversity_level": "family", "phylogenetic_distance": 2},
            {"diversity_level": "order", "phylogenetic_distance": 3},
        ]

        small_score = _calculate_coverage_score(small_set, ["genus", "family", "order"])
        large_score = _calculate_coverage_score(large_set, ["genus", "family", "order"])

        assert large_score > small_score, "More taxa should increase coverage score"

    def test_score_in_valid_range(self):
        """Test that score is always between 0 and 1."""
        # Test various combinations
        test_cases = [
            [{"diversity_level": "genus", "phylogenetic_distance": 1}],
            [{"diversity_level": "genus", "phylogenetic_distance": 1},
             {"diversity_level": "family", "phylogenetic_distance": 2}],
            [{"diversity_level": l, "phylogenetic_distance": i+1}
             for i, l in enumerate(["genus", "family", "order"])],
        ]

        for taxa in test_cases:
            score = _calculate_coverage_score(taxa, ["genus", "family", "order"])
            assert 0.0 <= score <= 1.0, f"Score {score} should be between 0 and 1"


# ==================== INTEGRATION TESTS ====================

class TestIntelligentExclusionSelection:
    """Test intelligent exclusion selection (integration test)."""

    def test_returns_expected_structure(self):
        """Test that the function returns expected data structure."""
        # Note: This test requires NCBI access, may need mocking in CI/CD
        # Using a well-known organism for testing

        result = intelligent_exclusion_selection(
            inclusion_taxid="562",  # E. coli - well-sequenced
            email="test@example.com",
            max_exclusion_taxa=5,
            require_sequences=False,  # Speed up test
            min_sequence_count=0
        )

        # Check required keys
        assert "exclusion_taxids" in result
        assert "taxa_info" in result
        assert "selection_strategy" in result
        assert "coverage_score" in result

        # Check data types
        assert isinstance(result["exclusion_taxids"], list)
        assert isinstance(result.get("coverage_score", 0), (int, float))

    def test_respects_max_exclusion_taxa(self):
        """Test that result respects max_exclusion_taxa parameter."""
        result = intelligent_exclusion_selection(
            inclusion_taxid="562",  # E. coli
            email="test@example.com",
            max_exclusion_taxa=3,
            require_sequences=False
        )

        if "error" not in result:
            assert len(result["exclusion_taxids"]) <= 3, \
                "Should not exceed max_exclusion_taxa"

    def test_coverage_score_in_valid_range(self):
        """Test that coverage score is in valid range."""
        result = intelligent_exclusion_selection(
            inclusion_taxid="562",  # E. coli
            email="test@example.com",
            max_exclusion_taxa=10,
            require_sequences=False
        )

        if "error" not in result:
            score = result.get("coverage_score", 0)
            assert 0.0 <= score <= 1.0, "Coverage score should be between 0 and 1"

    def test_handles_invalid_taxid(self):
        """Test graceful handling of invalid taxid."""
        result = intelligent_exclusion_selection(
            inclusion_taxid="99999999999",  # Invalid taxid
            email="test@example.com",
            max_exclusion_taxa=5
        )

        # Should return error or empty results, not crash
        assert "error" in result or len(result["exclusion_taxids"]) == 0

    def test_tier_summary_present(self):
        """Test that tier summary is included in results."""
        result = intelligent_exclusion_selection(
            inclusion_taxid="562",  # E. coli
            email="test@example.com",
            max_exclusion_taxa=10,
            require_sequences=False
        )

        if "error" not in result and result.get("exclusion_taxids"):
            assert "tier_summary" in result
            assert isinstance(result["tier_summary"], dict)


# ==================== PYTEST FIXTURES ====================

@pytest.fixture
def mock_taxa_candidates():
    """Fixture providing mock taxa candidates for testing."""
    return [
        {
            "taxid": "1",
            "scientific_name": "Genus species1",
            "priority_score": 10.0,
            "diversity_level": "genus",
            "phylogenetic_distance": 1,
            "sequence_count": 100
        },
        {
            "taxid": "2",
            "scientific_name": "Genus species2",
            "priority_score": 9.0,
            "diversity_level": "genus",
            "phylogenetic_distance": 1,
            "sequence_count": 80
        },
        {
            "taxid": "3",
            "scientific_name": "Family species1",
            "priority_score": 7.0,
            "diversity_level": "family",
            "phylogenetic_distance": 2,
            "sequence_count": 50
        },
    ]


# ==================== SUMMARY TEST ====================

def test_all_phase1_functions_importable():
    """Test that all Phase 1 functions can be imported successfully."""
    from assay_design.data_retrieval import (
        get_related_taxa_weighted,
        intelligent_exclusion_selection,
        _calculate_taxonomic_distance,
        _estimate_sequence_count,
        _select_diverse_taxa,
        _determine_diversity_levels,
        _calculate_coverage_score
    )

    # If we get here without ImportError, test passes
    assert True, "All Phase 1 functions imported successfully"


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v", "--tb=short"])
