#!/usr/bin/env python3

import logging
from assay_design.primer_design import (
    calculate_primer_properties,
    design_primers,
    optimize_primer_pair
)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def test_primer_properties():
    """Test calculating primer properties."""
    test_primers = [
        "ATGCATGCATGCATGCATGC",
        "GGGGGGGGGGCCCCCCCCC",
        "ATATATATATATATATAT",
    ]
    
    try:
        for primer in test_primers:
            props = calculate_primer_properties(primer)
            print(f"Properties for primer {primer}:")
            print(f"  Length: {props.get('length')}")
            print(f"  GC Content: {props.get('gc_content'):.1f}%")
            print(f"  Tm: {props.get('tm'):.1f}°C")
            print(f"  Self-complementarity: {props.get('self_complementarity')}")
            print(f"  Has repeats: {props.get('has_repeats')}")
            print(f"  End stability: {props.get('end_stability'):.2f}")
        return True
    except Exception as e:
        logger.error(f"Error: {e}")
        return False

def test_primer_design():
    """Test designing primers."""
    # Create a mock marker info
    marker_info = {
        "marker_region_start": 10,
        "marker_region_end": 110,
        "marker_length": 100,
        "marker_sequence": "ATGCATGCATGCAAAAAAAAAAAATTGCATGCATGCATGCATGCAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAAATGCATGCATGCTT",
        "conservation_score": 0.95,
    }
    
    try:
        primers = design_primers(marker_info)
        print(f"Successfully designed {len(primers)} primers:")
        for i, primer in enumerate(primers, 1):
            print(f"Primer {i} ({primer.get('orientation', 'unknown')}):")
            print(f"  Sequence: {primer.get('sequence', '')}")
            properties = primer.get('properties', {})
            if properties:
                print(f"  Length: {properties.get('length', 0)}")
                print(f"  GC Content: {properties.get('gc_content', 0):.1f}%")
                print(f"  Tm: {properties.get('tm', 0):.1f}°C")
                
        # Test primer pair optimization
        primer_pair = optimize_primer_pair(primers)
        print("\nOptimized primer pair:")
        print(f"  Product size: {primer_pair.get('product_size')} bp")
        print(f"  Compatibility score: {primer_pair.get('compatibility_score'):.2f}")
        
        return True
    except Exception as e:
        logger.error(f"Error: {e}")
        return False

if __name__ == "__main__":
    print("Testing primer_design module...")
    
    tests = [
        ("Calculate primer properties", test_primer_properties),
        ("Design and optimize primers", test_primer_design),
    ]
    
    for name, test_func in tests:
        print(f"\n=== Testing: {name} ===")
        success = test_func()
        print(f"Result: {'SUCCESS' if success else 'FAILURE'}")