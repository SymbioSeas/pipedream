"""
Unit tests for gene selection module (Phase 1.3 & 1.4).

Tests for:
- GeneMetadata dataclass
- GeneSelectionCriteria class
- Gene database completeness
- Ranking algorithms
- Recommendation system
"""

import pytest
import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from assay_design.gene_selection import (
    GeneMetadata,
    GeneSelectionCriteria
)


# ==================== GENE DATABASE TESTS ====================

class TestGeneDatabases:
    """Test gene database completeness and consistency."""

    def test_bacteria_single_copy_genes_exist(self):
        """Test that bacterial single-copy genes are defined."""
        genes = GeneSelectionCriteria.BACTERIA_SINGLE_COPY
        assert len(genes) > 0, "Should have bacterial single-copy genes"
        assert "rpoB" in genes, "Should include rpoB"
        assert "gyrB" in genes, "Should include gyrB"
        assert "recA" in genes, "Should include recA"

    def test_bacteria_marker_genes_exist(self):
        """Test that bacterial marker genes are defined."""
        genes = GeneSelectionCriteria.BACTERIA_MARKER
        assert len(genes) > 0, "Should have bacterial marker genes"
        assert "16S ribosomal RNA" in genes, "Should include 16S rRNA"

    def test_archaea_genes_exist(self):
        """Test that archaeal genes are defined."""
        single_copy = GeneSelectionCriteria.ARCHAEA_SINGLE_COPY
        marker = GeneSelectionCriteria.ARCHAEA_MARKER

        assert len(single_copy) > 0, "Should have archaeal single-copy genes"
        assert len(marker) > 0, "Should have archaeal marker genes"

    def test_fungi_genes_exist(self):
        """Test that fungal genes are defined."""
        single_copy = GeneSelectionCriteria.FUNGI_SINGLE_COPY
        marker = GeneSelectionCriteria.FUNGI_MARKER

        assert len(single_copy) > 0, "Should have fungal single-copy genes"
        assert len(marker) > 0, "Should have fungal marker genes"

    def test_gene_metadata_complete(self):
        """Test that all genes have complete metadata."""
        all_genes = []
        all_genes.extend(GeneSelectionCriteria.BACTERIA_SINGLE_COPY.values())
        all_genes.extend(GeneSelectionCriteria.BACTERIA_MARKER.values())
        all_genes.extend(GeneSelectionCriteria.ARCHAEA_SINGLE_COPY.values())
        all_genes.extend(GeneSelectionCriteria.ARCHAEA_MARKER.values())
        all_genes.extend(GeneSelectionCriteria.FUNGI_SINGLE_COPY.values())
        all_genes.extend(GeneSelectionCriteria.FUNGI_MARKER.values())

        for gene_meta in all_genes:
            assert isinstance(gene_meta, GeneMetadata), "Should be GeneMetadata instance"
            assert gene_meta.gene_name, "Should have gene name"
            assert gene_meta.copy_number in ["single", "low", "multi"], \
                f"Invalid copy number: {gene_meta.copy_number}"
            assert gene_meta.hgt_propensity in ["low", "medium", "high"], \
                f"Invalid HGT propensity: {gene_meta.hgt_propensity}"
            assert gene_meta.evolutionary_rate in ["slow", "moderate", "fast"], \
                f"Invalid evolutionary rate: {gene_meta.evolutionary_rate}"
            assert len(gene_meta.typical_length_bp) == 2, "Should have length range"
            assert gene_meta.typical_length_bp[0] < gene_meta.typical_length_bp[1], \
                "Min length should be less than max length"
            assert gene_meta.functional_category, "Should have functional category"
            assert gene_meta.description, "Should have description"

    def test_single_copy_genes_marked_correctly(self):
        """Test that single-copy databases only contain single/low copy genes."""
        single_copy_dbs = [
            GeneSelectionCriteria.BACTERIA_SINGLE_COPY,
            GeneSelectionCriteria.ARCHAEA_SINGLE_COPY,
            GeneSelectionCriteria.FUNGI_SINGLE_COPY
        ]

        for db in single_copy_dbs:
            for gene_name, metadata in db.items():
                assert metadata.copy_number in ["single", "low"], \
                    f"{gene_name} in single-copy DB should be single/low copy, got {metadata.copy_number}"

    def test_marker_genes_can_be_multi_copy(self):
        """Test that marker gene databases can include multi-copy genes."""
        marker_dbs = [
            GeneSelectionCriteria.BACTERIA_MARKER,
            GeneSelectionCriteria.ARCHAEA_MARKER,
            GeneSelectionCriteria.FUNGI_MARKER
        ]

        # At least one marker should be multi-copy
        has_multi_copy = False
        for db in marker_dbs:
            for metadata in db.values():
                if metadata.copy_number == "multi":
                    has_multi_copy = True
                    break

        assert has_multi_copy, "Marker genes should include multi-copy genes like rRNA"


# ==================== CLASS METHOD TESTS ====================

class TestGetGenesForDomain:
    """Test get_genes_for_domain() class method."""

    def test_bacteria_single_copy_only(self):
        """Test bacteria with single-copy preference."""
        genes = GeneSelectionCriteria.get_genes_for_domain('bacteria', prefer_single_copy=True)

        assert len(genes) > 0, "Should return genes"
        # All should be single/low copy
        for metadata in genes.values():
            assert metadata.copy_number in ["single", "low"], \
                "With prefer_single_copy, should only get single/low copy genes"

    def test_bacteria_include_markers(self):
        """Test bacteria with markers included."""
        genes = GeneSelectionCriteria.get_genes_for_domain('bacteria', prefer_single_copy=False)

        # Should include 16S rRNA
        assert "16S ribosomal RNA" in genes, "Should include 16S when markers allowed"

    def test_archaea_domain(self):
        """Test archaea domain selection."""
        genes = GeneSelectionCriteria.get_genes_for_domain('archaea', prefer_single_copy=True)

        assert len(genes) > 0, "Should return archaeal genes"
        # Check for known archaeal genes
        gene_names = list(genes.keys())
        assert any('rps' in g or 'rpl' in g or 'ef' in g.lower() for g in gene_names), \
            "Should include ribosomal protein or elongation factor genes"

    def test_fungi_domain(self):
        """Test fungi domain selection."""
        genes = GeneSelectionCriteria.get_genes_for_domain('fungi', prefer_single_copy=True)

        assert len(genes) > 0, "Should return fungal genes"
        # Check for known fungal genes
        gene_names = list(genes.keys())
        assert any('RPB' in g or 'TEF' in g for g in gene_names), \
            "Should include RNA polymerase or translation elongation factor genes"

    def test_eukaryote_alias_for_fungi(self):
        """Test that 'eukaryote' works as alias for fungi."""
        fungi_genes = GeneSelectionCriteria.get_genes_for_domain('fungi')
        euk_genes = GeneSelectionCriteria.get_genes_for_domain('eukaryote')

        assert fungi_genes.keys() == euk_genes.keys(), \
            "Eukaryote should return same genes as fungi"

    def test_unknown_domain_defaults_to_bacteria(self):
        """Test that unknown domain defaults to bacteria."""
        genes = GeneSelectionCriteria.get_genes_for_domain('unknown_domain')

        assert len(genes) > 0, "Should return genes even with unknown domain"
        # Should be similar to bacteria
        assert "rpoB" in genes or "gyrB" in genes, \
            "Default should include bacterial genes"

    def test_case_insensitive(self):
        """Test that domain names are case-insensitive."""
        lower = GeneSelectionCriteria.get_genes_for_domain('bacteria')
        upper = GeneSelectionCriteria.get_genes_for_domain('BACTERIA')
        mixed = GeneSelectionCriteria.get_genes_for_domain('BaCteRia')

        assert lower.keys() == upper.keys() == mixed.keys(), \
            "Domain names should be case-insensitive"


class TestRankGenesByCriteria:
    """Test rank_genes_by_criteria() class method."""

    def test_returns_sorted_list(self):
        """Test that ranking returns sorted list of tuples."""
        genes = GeneSelectionCriteria.get_genes_for_domain('bacteria')
        ranked = GeneSelectionCriteria.rank_genes_by_criteria(genes)

        assert isinstance(ranked, list), "Should return list"
        assert len(ranked) == len(genes), "Should rank all genes"

        # Check tuple structure
        for gene_name, score in ranked:
            assert isinstance(gene_name, str), "First element should be gene name"
            assert isinstance(score, (int, float)), "Second element should be score"

    def test_scores_in_valid_range(self):
        """Test that all scores are between 0 and 1."""
        genes = GeneSelectionCriteria.get_genes_for_domain('bacteria')
        ranked = GeneSelectionCriteria.rank_genes_by_criteria(genes)

        for gene_name, score in ranked:
            assert 0.0 <= score <= 1.0, f"Score {score} for {gene_name} should be in [0, 1]"

    def test_sorted_descending(self):
        """Test that results are sorted in descending order."""
        genes = GeneSelectionCriteria.get_genes_for_domain('bacteria')
        ranked = GeneSelectionCriteria.rank_genes_by_criteria(genes)

        scores = [score for _, score in ranked]
        assert scores == sorted(scores, reverse=True), "Should be sorted descending"

    def test_single_copy_ranked_higher(self):
        """Test that single-copy genes generally rank higher."""
        # Get mix of single and multi-copy genes
        genes = GeneSelectionCriteria.get_genes_for_domain('bacteria', prefer_single_copy=False)
        ranked = GeneSelectionCriteria.rank_genes_by_criteria(genes)

        # Top genes should tend to be single-copy
        top_5_genes = [name for name, score in ranked[:5]]
        single_copy_count = sum(
            1 for name in top_5_genes
            if genes[name].copy_number == 'single'
        )

        assert single_copy_count >= 3, "At least 3 of top 5 should be single-copy"

    def test_custom_weights_change_ranking(self):
        """Test that custom weights can affect scores (order might not always change)."""
        genes = GeneSelectionCriteria.get_genes_for_domain('bacteria')

        # Default weights
        default_ranked = GeneSelectionCriteria.rank_genes_by_criteria(genes)

        # Custom weights heavily favoring length
        custom_weights = {
            'single_copy': 0.1,
            'hgt_resistance': 0.1,
            'evolutionary_rate': 0.1,
            'length': 0.7
        }
        custom_ranked = GeneSelectionCriteria.rank_genes_by_criteria(genes, custom_weights)

        # Scores should differ (even if order might be similar for very similar genes)
        default_scores = [score for _, score in default_ranked]
        custom_scores = [score for _, score in custom_ranked]

        # At least some scores should be different
        assert default_scores != custom_scores, "Custom weights should affect scores"


class TestGetGeneMetadata:
    """Test get_gene_metadata() class method."""

    def test_retrieves_known_gene(self):
        """Test retrieving metadata for known gene."""
        metadata = GeneSelectionCriteria.get_gene_metadata('rpoB')

        assert metadata is not None, "Should find rpoB"
        assert isinstance(metadata, GeneMetadata), "Should return GeneMetadata"
        assert metadata.gene_name == 'rpoB', "Gene name should match"

    def test_unknown_gene_returns_none(self):
        """Test that unknown gene returns None."""
        metadata = GeneSelectionCriteria.get_gene_metadata('unknown_fake_gene')
        assert metadata is None, "Unknown gene should return None"

    def test_domain_specific_search(self):
        """Test searching within specific domain."""
        # rpoB is bacterial
        metadata = GeneSelectionCriteria.get_gene_metadata('rpoB', domain='bacteria')
        assert metadata is not None, "Should find rpoB in bacteria"

        # RPB1 is fungal
        metadata = GeneSelectionCriteria.get_gene_metadata('RPB1', domain='fungi')
        assert metadata is not None, "Should find RPB1 in fungi"

    def test_searches_all_domains_if_none(self):
        """Test that search covers all domains when domain=None."""
        # Known genes from different domains
        test_genes = ['rpoB', 'rps3', 'RPB1', 'ITS']

        for gene in test_genes:
            metadata = GeneSelectionCriteria.get_gene_metadata(gene, domain=None)
            assert metadata is not None, f"Should find {gene} across all domains"


class TestListAllGenes:
    """Test list_all_genes() class method."""

    def test_returns_all_genes_without_filter(self):
        """Test listing all genes across all domains."""
        all_genes = GeneSelectionCriteria.list_all_genes()

        assert isinstance(all_genes, list), "Should return list"
        assert len(all_genes) > 0, "Should have genes"

        # Should include genes from all domains
        gene_set = set(all_genes)
        assert 'rpoB' in gene_set, "Should include bacterial genes"
        assert 'ITS' in gene_set or 'RPB1' in gene_set, "Should include fungal genes"

    def test_filters_by_domain(self):
        """Test filtering by domain."""
        bacteria_genes = GeneSelectionCriteria.list_all_genes('bacteria')
        fungi_genes = GeneSelectionCriteria.list_all_genes('fungi')

        # Should be different
        assert set(bacteria_genes) != set(fungi_genes), \
            "Different domains should have different genes"

        # Bacterial genes should include rpoB
        assert 'rpoB' in bacteria_genes, "Bacteria should include rpoB"

        # Fungal genes should include ITS or RPB1
        assert 'ITS' in fungi_genes or 'RPB1' in fungi_genes, \
            "Fungi should include typical fungal genes"

    def test_returns_sorted_list(self):
        """Test that gene list is sorted."""
        genes = GeneSelectionCriteria.list_all_genes()
        assert genes == sorted(genes), "Gene list should be sorted alphabetically"


class TestGetRecommendedGenes:
    """Test get_recommended_genes() class method."""

    def test_quantification_use_case(self):
        """Test recommendations for quantification use case."""
        recs = GeneSelectionCriteria.get_recommended_genes(
            'bacteria', use_case='quantification', max_genes=5
        )

        assert len(recs) <= 5, "Should not exceed max_genes"
        assert len(recs) > 0, "Should return recommendations"

        # Check structure
        for gene_name, metadata, score in recs:
            assert isinstance(gene_name, str)
            assert isinstance(metadata, GeneMetadata)
            assert isinstance(score, (int, float))

        # For quantification, should prefer single-copy
        single_copy_count = sum(
            1 for _, meta, _ in recs
            if meta.copy_number == 'single'
        )
        assert single_copy_count >= len(recs) * 0.6, \
            "Quantification should favor single-copy genes"

    def test_phylogeny_use_case(self):
        """Test recommendations for phylogeny use case."""
        recs = GeneSelectionCriteria.get_recommended_genes(
            'bacteria', use_case='phylogeny', max_genes=5
        )

        assert len(recs) > 0, "Should return recommendations"

        # For phylogeny, should have low HGT propensity
        low_hgt_count = sum(
            1 for _, meta, _ in recs
            if meta.hgt_propensity == 'low'
        )
        assert low_hgt_count >= len(recs) * 0.6, \
            "Phylogeny should favor HGT-resistant genes"

    def test_detection_use_case(self):
        """Test recommendations for detection use case."""
        recs = GeneSelectionCriteria.get_recommended_genes(
            'bacteria', use_case='detection', max_genes=5
        )

        assert len(recs) > 0, "Should return recommendations"

        # Detection may include marker genes (multi-copy for sensitivity)
        # Just verify it returns valid results
        for _, metadata, _ in recs:
            assert metadata.copy_number in ['single', 'low', 'multi']

    def test_different_domains(self):
        """Test recommendations work for different domains."""
        for domain in ['bacteria', 'archaea', 'fungi']:
            recs = GeneSelectionCriteria.get_recommended_genes(
                domain, use_case='quantification', max_genes=3
            )
            assert len(recs) > 0, f"Should return recommendations for {domain}"

    def test_unknown_use_case_uses_defaults(self):
        """Test that unknown use case doesn't crash."""
        recs = GeneSelectionCriteria.get_recommended_genes(
            'bacteria', use_case='unknown_case', max_genes=3
        )

        assert len(recs) > 0, "Should return recommendations even with unknown use case"

    def test_recommendations_sorted_by_score(self):
        """Test that recommendations are sorted by score."""
        recs = GeneSelectionCriteria.get_recommended_genes(
            'bacteria', use_case='quantification', max_genes=10
        )

        scores = [score for _, _, score in recs]
        assert scores == sorted(scores, reverse=True), \
            "Recommendations should be sorted by score (descending)"


# ==================== INTEGRATION TESTS ====================

class TestGeneSelectionWorkflow:
    """Test complete gene selection workflow."""

    def test_complete_workflow_bacteria(self):
        """Test complete workflow for bacterial gene selection."""
        # Step 1: Get genes for domain
        genes = GeneSelectionCriteria.get_genes_for_domain('bacteria', prefer_single_copy=True)
        assert len(genes) > 0

        # Step 2: Rank genes
        ranked = GeneSelectionCriteria.rank_genes_by_criteria(genes)
        assert len(ranked) > 0

        # Step 3: Get top gene metadata
        top_gene_name = ranked[0][0]
        metadata = GeneSelectionCriteria.get_gene_metadata(top_gene_name)
        assert metadata is not None

        # Step 4: Verify top gene is single-copy
        assert metadata.copy_number == 'single', "Top gene should be single-copy"

    def test_recommended_genes_are_in_database(self):
        """Test that recommended genes are actually in the database."""
        recs = GeneSelectionCriteria.get_recommended_genes('bacteria', 'quantification', 5)

        all_genes = GeneSelectionCriteria.list_all_genes('bacteria')

        for gene_name, _, _ in recs:
            assert gene_name in all_genes, f"{gene_name} should be in database"


# ==================== SUMMARY TEST ====================

def test_gene_selection_module_complete():
    """Test that gene selection module is complete and functional."""
    # Test that all major components work
    domains = ['bacteria', 'archaea', 'fungi']

    for domain in domains:
        # Get genes
        genes = GeneSelectionCriteria.get_genes_for_domain(domain)
        assert len(genes) > 0, f"Should have genes for {domain}"

        # Rank genes
        ranked = GeneSelectionCriteria.rank_genes_by_criteria(genes)
        assert len(ranked) > 0, f"Should rank genes for {domain}"

        # Get recommendations
        recs = GeneSelectionCriteria.get_recommended_genes(domain, 'quantification', 3)
        assert len(recs) > 0, f"Should get recommendations for {domain}"

    print("✓ Gene selection module is complete and functional")


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v", "--tb=short"])
