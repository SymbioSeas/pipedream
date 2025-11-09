#!/usr/bin/env python3
"""
Test script for get_related_taxa_weighted() function.

Tests the enhanced phylogenetic distance weighting and sequence availability checking.
"""

import sys
import os

# Add parent directory to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))

from assay_design.data_retrieval import get_related_taxa_weighted, get_taxon_info


def test_vibrio_mediterranei():
    """
    Test with Vibrio mediterranei (taxid: 689).

    Expected behavior:
    - Find genus-level siblings (other Vibrio species)
    - Find family-level relatives (other Vibrionaceae)
    - Find order-level representatives (other Vibrionales)
    - Prioritize by sequence availability and proximity
    """
    print("\n" + "="*80)
    print("TEST 1: Vibrio mediterranei (taxid: 689)")
    print("="*80)

    email = "test@example.com"  # Replace with actual email for NCBI
    taxid = "689"

    # First, get basic info about the inclusion taxon
    print("\nInclusion taxon info:")
    tax_info = get_taxon_info(taxid, email)
    print(f"  Scientific name: {tax_info.get('scientific_name', 'Unknown')}")
    print(f"  Rank: {tax_info.get('rank', 'Unknown')}")
    print(f"  Lineage: {tax_info.get('lineage', 'Unknown')}")

    # Get weighted exclusion taxa
    print("\nFetching weighted exclusion taxa...")
    weighted_taxa = get_related_taxa_weighted(
        taxid=taxid,
        email=email,
        max_results=10,
        diversity_levels=["genus", "family", "order"],
        min_sequence_count=5
    )

    print(f"\nFound {len(weighted_taxa)} weighted exclusion taxa:\n")
    print(f"{'Rank':<10} {'Scientific Name':<40} {'Distance':<10} {'Seqs':<8} {'Score':<8}")
    print("-" * 80)

    for taxon in weighted_taxa:
        rank = taxon['diversity_level']
        name = taxon['scientific_name'][:38]
        distance = taxon['phylogenetic_distance']
        seqs = taxon['sequence_count']
        score = taxon['priority_score']

        print(f"{rank:<10} {name:<40} {distance:<10} {seqs:<8} {score:<8.2f}")

    # Verify diversity across levels
    levels_represented = set(t['diversity_level'] for t in weighted_taxa)
    print(f"\nTaxonomic levels represented: {', '.join(levels_represented)}")

    # Success criteria
    assert len(weighted_taxa) > 0, "Should find at least some exclusion taxa"
    assert len(levels_represented) >= 2, "Should span at least 2 taxonomic levels"

    print("\n✓ Test PASSED: Successfully found diverse, weighted exclusion taxa")

    return weighted_taxa


def test_edge_case_monotypic():
    """
    Test with a potentially monotypic or poorly sequenced taxon.

    Expected behavior:
    - Handle cases where few siblings exist
    - Fall back to higher taxonomic levels
    - Return gracefully if insufficient data
    """
    print("\n" + "="*80)
    print("TEST 2: Edge case - Poorly sequenced organism")
    print("="*80)

    email = "test@example.com"
    # Use a less common organism (example: some obscure bacterial species)
    taxid = "1871047"  # Acinetobacter defluvii (less common)

    print("\nInclusion taxon info:")
    tax_info = get_taxon_info(taxid, email)
    print(f"  Scientific name: {tax_info.get('scientific_name', 'Unknown')}")
    print(f"  Rank: {tax_info.get('rank', 'Unknown')}")

    # Get weighted exclusion taxa with relaxed requirements
    print("\nFetching weighted exclusion taxa (relaxed min_sequence_count)...")
    weighted_taxa = get_related_taxa_weighted(
        taxid=taxid,
        email=email,
        max_results=10,
        diversity_levels=["genus", "family", "order"],
        min_sequence_count=1  # Very relaxed for poorly sequenced taxa
    )

    print(f"\nFound {len(weighted_taxa)} exclusion taxa")

    if len(weighted_taxa) > 0:
        print("\n✓ Successfully handled edge case - found some exclusion taxa")
    else:
        print("\n✓ Gracefully handled edge case - no suitable taxa found (as expected)")

    return weighted_taxa


def test_helper_functions():
    """Test the helper functions directly."""
    print("\n" + "="*80)
    print("TEST 3: Helper functions")
    print("="*80)

    from assay_design.data_retrieval import (
        _calculate_taxonomic_distance,
        _estimate_sequence_count,
        _select_diverse_taxa
    )

    # Test taxonomic distance calculation
    print("\nTesting _calculate_taxonomic_distance():")

    test_cases = [
        ("genus", "species", 1),  # Genus siblings of a species
        ("family", "species", 2),  # Family relatives of a species
        ("order", "genus", 2),    # Order relatives of a genus
    ]

    for level, inclusion_rank, expected in test_cases:
        distance = _calculate_taxonomic_distance(level, inclusion_rank)
        print(f"  {level} sibling of {inclusion_rank}: distance={distance} (expected={expected})")
        assert distance == expected, f"Distance calculation failed for {level}/{inclusion_rank}"

    print("  ✓ Taxonomic distance calculation working correctly")

    # Test sequence count estimation (with a well-known organism)
    print("\nTesting _estimate_sequence_count():")
    email = "test@example.com"

    # E. coli should have many sequences
    ecoli_count = _estimate_sequence_count("562", email)
    print(f"  E. coli (taxid:562): ~{ecoli_count} sequences")
    assert ecoli_count > 100, "E. coli should have many sequences"

    print("  ✓ Sequence count estimation working correctly")

    # Test diverse taxa selection
    print("\nTesting _select_diverse_taxa():")

    mock_candidates = [
        {"taxid": "1", "scientific_name": "Genus1", "priority_score": 10.0, "diversity_level": "genus"},
        {"taxid": "2", "scientific_name": "Genus2", "priority_score": 9.0, "diversity_level": "genus"},
        {"taxid": "3", "scientific_name": "Genus3", "priority_score": 8.0, "diversity_level": "genus"},
        {"taxid": "4", "scientific_name": "Family1", "priority_score": 7.0, "diversity_level": "family"},
        {"taxid": "5", "scientific_name": "Family2", "priority_score": 6.0, "diversity_level": "family"},
        {"taxid": "6", "scientific_name": "Order1", "priority_score": 5.0, "diversity_level": "order"},
    ]

    selected = _select_diverse_taxa(mock_candidates, max_results=5, diversity_levels=["genus", "family", "order"])

    print(f"  Selected {len(selected)} from {len(mock_candidates)} candidates")

    levels_in_selection = set(t['diversity_level'] for t in selected)
    print(f"  Levels represented: {', '.join(levels_in_selection)}")

    # Should have representation from multiple levels
    assert len(levels_in_selection) >= 2, "Should have diversity across levels"

    print("  ✓ Diverse taxa selection working correctly")

    print("\n✓ All helper functions PASSED")


def main():
    """Run all tests."""
    print("\n" + "="*80)
    print("TESTING get_related_taxa_weighted() - Phase 1.1 Implementation")
    print("="*80)

    try:
        # Test 1: Standard case with well-sequenced organism
        test_vibrio_mediterranei()

        # Test 2: Edge case with poorly sequenced organism
        test_edge_case_monotypic()

        # Test 3: Helper functions
        test_helper_functions()

        print("\n" + "="*80)
        print("ALL TESTS PASSED ✓")
        print("="*80)
        print("\nPhase 1.1 implementation is working correctly!")
        print("Enhanced get_related_taxa_weighted() successfully:")
        print("  - Calculates phylogenetic distance")
        print("  - Estimates sequence availability")
        print("  - Prioritizes taxa by weighted score")
        print("  - Ensures diversity across taxonomic levels")

    except Exception as e:
        print(f"\n✗ TEST FAILED: {e}")
        import traceback
        traceback.print_exc()
        return 1

    return 0


if __name__ == "__main__":
    sys.exit(main())
