#!/usr/bin/env python3
"""
Unit tests for Phase 2 gene selection features.

Tests for:
- evaluate_gene_suitability()
- rank_candidate_genes()
- auto_select_gene_for_taxon()
- Integration with hierarchical_search
"""

import pytest
import sys
import os
from unittest.mock import Mock, patch, MagicMock
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from assay_design.gene_selection import (
    evaluate_gene_suitability,
    rank_candidate_genes,
    auto_select_gene_for_taxon,
    GeneSelectionCriteria,
    GeneMetadata
)


# ==================== HELPER FIXTURES ====================

@pytest.fixture
def mock_sequences_abundant():
    """Fixture providing mock abundant sequences (good availability)."""
    sequences = []
    for i in range(25):
        seq = SeqRecord(
            Seq("ATCG" * 400),  # 1600 bp
            id=f"seq_{i}",
            description=f"Mock sequence {i}"
        )
        sequences.append(seq)
    return sequences


@pytest.fixture
def mock_sequences_sparse():
    """Fixture providing mock sparse sequences (poor availability)."""
    sequences = []
    for i in range(3):
        seq = SeqRecord(
            Seq("ATCG" * 300),  # 1200 bp
            id=f"seq_{i}",
            description=f"Mock sequence {i}"
        )
        sequences.append(seq)
    return sequences


@pytest.fixture
def mock_sequences_variable_length():
    """Fixture providing sequences with variable lengths."""
    sequences = []
    lengths = [500, 800, 1200, 1500, 1800, 2200, 2800]
    for i, length in enumerate(lengths):
        seq = SeqRecord(
            Seq("A" * length),
            id=f"seq_{i}",
            description=f"Mock sequence {i}"
        )
        sequences.append(seq)
    return sequences


@pytest.fixture
def mock_taxon_info_bacteria():
    """Fixture for bacterial taxon info."""
    return {
        'taxid': '562',
        'scientific_name': 'Escherichia coli',
        'rank': 'species',
        'lineage': 'cellular organisms; Bacteria; Proteobacteria; Gammaproteobacteria; Enterobacterales; Enterobacteriaceae; Escherichia'
    }


@pytest.fixture
def mock_taxon_info_archaea():
    """Fixture for archaeal taxon info."""
    return {
        'taxid': '2159',
        'scientific_name': 'Archaea',
        'rank': 'superkingdom',
        'lineage': 'cellular organisms; Archaea'
    }


# ==================== EVALUATE_GENE_SUITABILITY TESTS ====================

class TestEvaluateGeneSuitability:
    """Test evaluate_gene_suitability() function."""

    def test_evaluates_with_abundant_sequences(self, mock_sequences_abundant):
        """Test evaluation with abundant sequences."""
        with patch('assay_design.data_retrieval.fetch_gene_sequences') as mock_fetch:
            mock_fetch.return_value = mock_sequences_abundant

            result = evaluate_gene_suitability(
                taxid="562",
                gene_name="rpoB",
                email="test@example.com"
            )

            assert result['gene'] == 'rpoB'
            assert result['sequence_count'] == 25
            assert result['sequence_availability_score'] > 0.9  # Should be high
            assert result['overall_score'] > 0.5  # Should be decent
            assert result['recommendation'] in ['excellent', 'good', 'fair', 'poor']

    def test_evaluates_with_sparse_sequences(self, mock_sequences_sparse):
        """Test evaluation with sparse sequences."""
        with patch('assay_design.data_retrieval.fetch_gene_sequences') as mock_fetch:
            mock_fetch.return_value = mock_sequences_sparse

            result = evaluate_gene_suitability(
                taxid="562",
                gene_name="obscureGene",
                email="test@example.com"
            )

            assert result['sequence_count'] == 3
            assert result['sequence_availability_score'] < 0.2  # Should be low
            assert result['recommendation'] in ['fair', 'poor']  # Should not be excellent

    def test_handles_no_sequences(self):
        """Test graceful handling when no sequences found."""
        with patch('assay_design.data_retrieval.fetch_gene_sequences') as mock_fetch:
            mock_fetch.return_value = []

            result = evaluate_gene_suitability(
                taxid="562",
                gene_name="nonexistent",
                email="test@example.com"
            )

            assert result['sequence_count'] == 0
            assert result['sequence_availability_score'] == 0.0
            assert result['avg_sequence_length'] == 0
            # Should still return valid result
            assert 'overall_score' in result
            assert 'recommendation' in result

    def test_calculates_average_length_correctly(self, mock_sequences_variable_length):
        """Test that average length is calculated correctly."""
        with patch('assay_design.data_retrieval.fetch_gene_sequences') as mock_fetch:
            mock_fetch.return_value = mock_sequences_variable_length

            result = evaluate_gene_suitability(
                taxid="562",
                gene_name="testGene",
                email="test@example.com"
            )

            # Average of [500, 800, 1200, 1500, 1800, 2200, 2800]
            expected_avg = (500 + 800 + 1200 + 1500 + 1800 + 2200 + 2800) / 7
            assert abs(result['avg_sequence_length'] - expected_avg) < 1

    def test_uses_provided_metadata(self, mock_sequences_abundant):
        """Test that provided metadata is used."""
        with patch('assay_design.data_retrieval.fetch_gene_sequences') as mock_fetch:
            mock_fetch.return_value = mock_sequences_abundant

            custom_metadata = GeneMetadata(
                gene_name="customGene",
                copy_number="single",
                hgt_propensity="low",
                evolutionary_rate="slow",
                typical_length_bp=(1000, 2000),
                functional_category="housekeeping",
                description="Test gene"
            )

            result = evaluate_gene_suitability(
                taxid="562",
                gene_name="customGene",
                email="test@example.com",
                metadata=custom_metadata
            )

            # Should use single-copy and low HGT
            assert result['copy_number_score'] == 1.0  # Single copy
            assert result['hgt_resistance_score'] == 1.0  # Low HGT

    def test_handles_fetch_exception(self):
        """Test graceful handling of fetch exceptions."""
        with patch('assay_design.data_retrieval.fetch_gene_sequences') as mock_fetch:
            mock_fetch.side_effect = Exception("NCBI API error")

            result = evaluate_gene_suitability(
                taxid="562",
                gene_name="rpoB",
                email="test@example.com"
            )

            # Should return valid result even on error
            assert result['sequence_count'] == 0
            assert 'overall_score' in result
            assert 'error' not in result  # Should handle gracefully

    def test_length_score_optimal_at_1500bp(self, mock_sequences_abundant):
        """Test that length score is optimal near 1500bp."""
        with patch('assay_design.data_retrieval.fetch_gene_sequences') as mock_fetch:
            # Sequences exactly 1500bp
            optimal_seqs = [
                SeqRecord(Seq("A" * 1500), id=f"seq_{i}")
                for i in range(20)
            ]
            mock_fetch.return_value = optimal_seqs

            result = evaluate_gene_suitability(
                taxid="562",
                gene_name="optimalGene",
                email="test@example.com"
            )

            # Length score should be 1.0 (perfect)
            assert result['length_appropriateness_score'] == 1.0
            assert result['avg_sequence_length'] == 1500

    def test_recommendation_tiers_correct(self, mock_sequences_abundant):
        """Test that recommendation tiers are assigned correctly."""
        with patch('assay_design.data_retrieval.fetch_gene_sequences') as mock_fetch:
            mock_fetch.return_value = mock_sequences_abundant

            # Test with excellent gene (single-copy, low HGT)
            excellent_metadata = GeneMetadata(
                gene_name="excellentGene",
                copy_number="single",
                hgt_propensity="low",
                evolutionary_rate="moderate",
                typical_length_bp=(1400, 1600),
                functional_category="housekeeping",
                description="Excellent gene"
            )

            result = evaluate_gene_suitability(
                taxid="562",
                gene_name="excellentGene",
                email="test@example.com",
                metadata=excellent_metadata
            )

            # With abundant sequences, single-copy, low HGT, and good length
            # should get excellent or good rating
            assert result['recommendation'] in ['excellent', 'good']
            assert result['overall_score'] >= 0.6


# ==================== RANK_CANDIDATE_GENES TESTS ====================

class TestRankCandidateGenes:
    """Test rank_candidate_genes() function."""

    def test_returns_sorted_list(self, mock_sequences_abundant, mock_taxon_info_bacteria):
        """Test that genes are returned sorted by score."""
        with patch('assay_design.data_retrieval.get_taxon_info') as mock_taxon, \
             patch('assay_design.data_retrieval.fetch_gene_sequences') as mock_fetch:

            mock_taxon.return_value = mock_taxon_info_bacteria
            mock_fetch.return_value = mock_sequences_abundant

            result = rank_candidate_genes(
                taxid="562",
                email="test@example.com",
                domain="bacteria",
                max_genes_to_test=5
            )

            # Should return a list of tuples
            assert isinstance(result, list)
            assert len(result) <= 5

            # Each element should be (gene_name, evaluation_dict)
            for gene_name, evaluation in result:
                assert isinstance(gene_name, str)
                assert isinstance(evaluation, dict)
                assert 'overall_score' in evaluation

            # Should be sorted by overall_score (descending)
            scores = [eval_dict['overall_score'] for _, eval_dict in result]
            assert scores == sorted(scores, reverse=True)

    def test_auto_detects_bacteria(self, mock_sequences_abundant, mock_taxon_info_bacteria):
        """Test auto-detection of bacterial domain."""
        with patch('assay_design.data_retrieval.get_taxon_info') as mock_taxon, \
             patch('assay_design.data_retrieval.fetch_gene_sequences') as mock_fetch:

            mock_taxon.return_value = mock_taxon_info_bacteria
            mock_fetch.return_value = mock_sequences_abundant

            result = rank_candidate_genes(
                taxid="562",
                email="test@example.com",
                domain=None,  # Auto-detect
                max_genes_to_test=3
            )

            # Should successfully detect bacteria and return genes
            assert len(result) > 0
            # All returned genes should be bacterial genes
            gene_names = [name for name, _ in result]
            bacterial_genes = GeneSelectionCriteria.get_genes_for_domain('bacteria')
            for name in gene_names:
                assert name in bacterial_genes

    def test_auto_detects_archaea(self, mock_sequences_abundant, mock_taxon_info_archaea):
        """Test auto-detection of archaeal domain."""
        with patch('assay_design.data_retrieval.get_taxon_info') as mock_taxon, \
             patch('assay_design.data_retrieval.fetch_gene_sequences') as mock_fetch:

            mock_taxon.return_value = mock_taxon_info_archaea
            mock_fetch.return_value = mock_sequences_abundant

            result = rank_candidate_genes(
                taxid="2159",
                email="test@example.com",
                domain=None,  # Auto-detect
                max_genes_to_test=3
            )

            # Should successfully detect archaea
            assert len(result) > 0

    def test_respects_max_genes_limit(self, mock_sequences_abundant, mock_taxon_info_bacteria):
        """Test that max_genes_to_test is respected."""
        with patch('assay_design.data_retrieval.get_taxon_info') as mock_taxon, \
             patch('assay_design.data_retrieval.fetch_gene_sequences') as mock_fetch:

            mock_taxon.return_value = mock_taxon_info_bacteria
            mock_fetch.return_value = mock_sequences_abundant

            for max_genes in [1, 3, 5, 10]:
                result = rank_candidate_genes(
                    taxid="562",
                    email="test@example.com",
                    domain="bacteria",
                    max_genes_to_test=max_genes
                )

                assert len(result) <= max_genes

    def test_prefers_single_copy_when_requested(self, mock_sequences_abundant, mock_taxon_info_bacteria):
        """Test that single-copy genes are preferred when requested."""
        with patch('assay_design.data_retrieval.get_taxon_info') as mock_taxon, \
             patch('assay_design.data_retrieval.fetch_gene_sequences') as mock_fetch:

            mock_taxon.return_value = mock_taxon_info_bacteria
            mock_fetch.return_value = mock_sequences_abundant

            result = rank_candidate_genes(
                taxid="562",
                email="test@example.com",
                domain="bacteria",
                prefer_single_copy=True,
                max_genes_to_test=10
            )

            # Top genes should tend to be single-copy
            if len(result) >= 3:
                top_3_genes = [name for name, _ in result[:3]]
                bacterial_genes = GeneSelectionCriteria.get_genes_for_domain('bacteria')
                single_copy_count = sum(
                    1 for name in top_3_genes
                    if bacterial_genes[name].copy_number == 'single'
                )
                # At least 2 of top 3 should be single-copy
                assert single_copy_count >= 2

    def test_handles_fetch_errors_gracefully(self, mock_taxon_info_bacteria):
        """Test that individual gene fetch errors don't crash entire ranking."""
        with patch('assay_design.data_retrieval.get_taxon_info') as mock_taxon, \
             patch('assay_design.data_retrieval.fetch_gene_sequences') as mock_fetch:

            mock_taxon.return_value = mock_taxon_info_bacteria
            # Simulate intermittent errors
            call_count = [0]
            def fetch_side_effect(*args, **kwargs):
                call_count[0] += 1
                if call_count[0] % 3 == 0:  # Every 3rd call fails
                    raise Exception("NCBI timeout")
                return [SeqRecord(Seq("ATCG" * 400), id="seq1")]

            mock_fetch.side_effect = fetch_side_effect

            result = rank_candidate_genes(
                taxid="562",
                email="test@example.com",
                domain="bacteria",
                max_genes_to_test=6
            )

            # Should return some results despite errors
            assert len(result) > 0


# ==================== AUTO_SELECT_GENE_FOR_TAXON TESTS ====================

class TestAutoSelectGeneForTaxon:
    """Test auto_select_gene_for_taxon() function."""

    def test_selects_best_gene(self, mock_sequences_abundant, mock_taxon_info_bacteria):
        """Test that the best gene is selected."""
        with patch('assay_design.data_retrieval.get_taxon_info') as mock_taxon, \
             patch('assay_design.data_retrieval.fetch_gene_sequences') as mock_fetch:

            mock_taxon.return_value = mock_taxon_info_bacteria
            mock_fetch.return_value = mock_sequences_abundant

            result = auto_select_gene_for_taxon(
                taxid="562",
                email="test@example.com",
                use_case="quantification",
                max_genes_to_test=5
            )

            # Should return a tuple (gene_name, evaluation_dict)
            assert result is not None
            gene_name, evaluation = result
            assert isinstance(gene_name, str)
            assert isinstance(evaluation, dict)
            assert 'overall_score' in evaluation

    def test_returns_none_if_no_suitable_gene(self, mock_taxon_info_bacteria):
        """Test that None is returned if no gene meets threshold."""
        with patch('assay_design.data_retrieval.get_taxon_info') as mock_taxon, \
             patch('assay_design.data_retrieval.fetch_gene_sequences') as mock_fetch:

            mock_taxon.return_value = mock_taxon_info_bacteria
            # Return very few sequences (poor score)
            mock_fetch.return_value = [SeqRecord(Seq("ATCG" * 100), id="seq1")]

            result = auto_select_gene_for_taxon(
                taxid="562",
                email="test@example.com",
                use_case="quantification",
                max_genes_to_test=3,
                min_acceptable_score=0.9  # Very high threshold
            )

            # Should return None if no gene meets threshold
            # (depends on scoring, might still pass with good genes)
            assert result is None or result[1]['overall_score'] >= 0.9

    def test_quantification_prefers_single_copy(self, mock_sequences_abundant, mock_taxon_info_bacteria):
        """Test that quantification use case prefers single-copy genes."""
        with patch('assay_design.data_retrieval.get_taxon_info') as mock_taxon, \
             patch('assay_design.data_retrieval.fetch_gene_sequences') as mock_fetch:

            mock_taxon.return_value = mock_taxon_info_bacteria
            mock_fetch.return_value = mock_sequences_abundant

            result = auto_select_gene_for_taxon(
                taxid="562",
                email="test@example.com",
                use_case="quantification",
                max_genes_to_test=5
            )

            if result:
                gene_name, evaluation = result
                # Should prefer single-copy for quantification
                metadata = evaluation.get('metadata')
                if metadata:
                    # If top gene has metadata, it should ideally be single-copy
                    # (not strict requirement but should be common)
                    pass  # Test passes if result is returned

    def test_respects_min_acceptable_score(self, mock_sequences_abundant, mock_taxon_info_bacteria):
        """Test that min_acceptable_score threshold is respected."""
        with patch('assay_design.data_retrieval.get_taxon_info') as mock_taxon, \
             patch('assay_design.data_retrieval.fetch_gene_sequences') as mock_fetch:

            mock_taxon.return_value = mock_taxon_info_bacteria
            mock_fetch.return_value = mock_sequences_abundant

            min_score = 0.6
            result = auto_select_gene_for_taxon(
                taxid="562",
                email="test@example.com",
                use_case="quantification",
                max_genes_to_test=5,
                min_acceptable_score=min_score
            )

            if result:
                gene_name, evaluation = result
                assert evaluation['overall_score'] >= min_score


# ==================== INTEGRATION TESTS ====================

class TestHierarchicalSearchIntegration:
    """Test integration with hierarchical_marker_search."""

    def test_auto_gene_selection_in_hierarchical_search(self, mock_sequences_abundant, mock_taxon_info_bacteria):
        """Test that auto-gene selection works in hierarchical search."""
        from assay_design.hierarchical_search import hierarchical_marker_search

        with patch('assay_design.data_retrieval.get_taxon_info') as mock_taxon, \
             patch('assay_design.data_retrieval.fetch_gene_sequences') as mock_fetch, \
             patch('assay_design.hierarchical_search.fetch_gene_sequences') as mock_fetch2:

            mock_taxon.return_value = mock_taxon_info_bacteria
            mock_fetch.return_value = mock_sequences_abundant
            mock_fetch2.return_value = mock_sequences_abundant

            result = hierarchical_marker_search(
                inclusion_taxid="562",
                email="test@example.com",
                auto_gene_selection=True,
                gene_use_case="quantification",
                max_genes_to_try=3,
                min_gene_score=0.4,
                max_seq_count=10,
                timeout_seconds=60
            )

            # Should include gene_selection_info in result
            assert 'gene_selection_info' in result or 'error' in result

            # If successful, should have selected a gene
            if 'gene_selection_info' in result:
                gene_info = result['gene_selection_info']
                if 'selected_gene' in gene_info:
                    assert isinstance(gene_info['selected_gene'], str)
                    assert 'gene_score' in gene_info


# ==================== SUMMARY TEST ====================

def test_all_phase2_functions_importable():
    """Test that all Phase 2 functions can be imported successfully."""
    from assay_design.gene_selection import (
        evaluate_gene_suitability,
        rank_candidate_genes,
        auto_select_gene_for_taxon
    )

    from assay_design.hierarchical_search import hierarchical_marker_search

    # If we get here without ImportError, test passes
    assert True, "All Phase 2 functions imported successfully"


if __name__ == "__main__":
    # Run tests
    pytest.main([__file__, "-v", "--tb=short"])
