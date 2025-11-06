#!/usr/bin/env python3

import logging
from assay_design.specificity_validation import (
    validate_primer_specificity,
    simulate_inclusion_hits,
    simulate_exclusion_hits
)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def test_validation():
    """Test primer specificity validation."""
    # Mock primers
    primers = [
        {
            "name": "Forward_Primer",
            "sequence": "ATGCATGCATGCAAAAAA",
            "properties": {
                "length": 18,
                "gc_content": 38.9,
                "tm": 55.4
            },
            "orientation": "forward"
        },
        {
            "name": "Reverse_Primer",
            "sequence": "TTTTTTGCATGCATGCAT",
            "properties": {
                "length": 18,
                "gc_content": 38.9,
                "tm": 55.4
            },
            "orientation": "reverse"
        }
    ]
    
    # Test taxids
    inclusion_taxid = "689"  # Vibrio mediterranei
    exclusion_taxids = ["666", "670"]  # Related Vibrio species
    
    try:
        # Test validation
        validation_results = validate_primer_specificity(
            primers, inclusion_taxid, exclusion_taxids
        )
        
        print("Specificity validation results:")
        print(f"  Specificity score: {validation_results.get('specificity_score', 0):.2f}")
        
        # Check inclusion hits
        inclusion_hits = validation_results.get('inclusion_hits', [])
        print(f"\nExpected amplification in target ({len(inclusion_hits)} hits):")
        for hit in inclusion_hits:
            print(f"  {hit.get('organism')} - Product size: {hit.get('product_size')} bp")
        
        # Check exclusion hits
        exclusion_hits = validation_results.get('exclusion_hits', [])
        print(f"\nPotential cross-reactivity ({len(exclusion_hits)} hits):")
        for hit in exclusion_hits:
            print(f"  {hit.get('organism')} - Product size: {hit.get('product_size')} bp")
        
        # Check potential issues
        issues = validation_results.get('potential_issues', [])
        print(f"\nPotential issues ({len(issues)}):")
        for issue in issues:
            print(f"  - {issue}")
            
        return True
    except Exception as e:
        logger.error(f"Error: {e}")
        return False

def test_simulation():
    """Test the simulation functions for in silico PCR."""
    forward_primer = "ATGCATGCATGCAAAAAA"
    reverse_primer = "TTTTTTGCATGCATGCAT"
    inclusion_taxid = "689"
    exclusion_taxids = ["666", "670"]
    
    try:
        # Test inclusion simulation
        inclusion_hits = simulate_inclusion_hits(forward_primer, reverse_primer, inclusion_taxid)
        print(f"Simulated {len(inclusion_hits)} inclusion hits")
        
        # Test exclusion simulation
        exclusion_hits = simulate_exclusion_hits(forward_primer, reverse_primer, exclusion_taxids)
        print(f"Simulated {len(exclusion_hits)} exclusion hits")
        
        return True
    except Exception as e:
        logger.error(f"Error: {e}")
        return False

if __name__ == "__main__":
    print("Testing specificity_validation module...")
    
    tests = [
        ("Primer specificity validation", test_validation),
        ("In silico PCR simulation", test_simulation),
    ]
    
    for name, test_func in tests:
        print(f"\n=== Testing: {name} ===")
        success = test_func()
        print(f"Result: {'SUCCESS' if success else 'FAILURE'}")