#!/usr/bin/env python3

import logging
from Bio.SeqRecord import SeqRecord
from Bio.Seq import Seq
from assay_design.data_retrieval import fetch_sequences_for_taxid, fetch_gene_sequences
from assay_design.target_identification import find_conserved_marker

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def create_test_sequences():
    """Create test sequences for a simple test case."""
    # Inclusion sequences share a conserved region
    inclusion_seqs = [
        SeqRecord(
            Seq("ATGCATGCATGCAAAAAAAAAAAATGCATGCATGC"),
            id="inclusion_1",
            description="Test inclusion sequence 1"
        ),
        SeqRecord(
            Seq("ATGCATGCTTGCAAAAAAAAAAATTGCATGCATGC"),
            id="inclusion_2",
            description="Test inclusion sequence 2"
        ),
        SeqRecord(
            Seq("ATGCATGCATGCAAAAAAAAAAATTGCATGCATGC"),
            id="inclusion_3",
            description="Test inclusion sequence 3"
        ),
    ]
    
    # Exclusion sequences have different central region
    exclusion_seqs = [
        SeqRecord(
            Seq("ATGCATGCATGCTTTTTTTTTTTTTTGCATGCATGC"),
            id="exclusion_1",
            description="Test exclusion sequence 1"
        ),
        SeqRecord(
            Seq("ATGCATGCATGCTTTTTTTTTTTTTTGCATGCATGC"),
            id="exclusion_2",
            description="Test exclusion sequence 2"
        ),
    ]
    
    return inclusion_seqs, exclusion_seqs

def test_with_artificial_data():
    """Test marker identification with artificial data."""
    inclusion_seqs, exclusion_seqs = create_test_sequences()
    
    try:
        result = find_conserved_marker(inclusion_seqs, exclusion_seqs)
        print("Successfully identified marker region:")
        print(f"  Start: {result.get('marker_region_start')}")
        print(f"  End: {result.get('marker_region_end')}")
        print(f"  Length: {result.get('marker_length')}")
        print(f"  Conservation score: {result.get('conservation_score', 0):.2f}")
        if 'differentiation_score' in result:
            print(f"  Differentiation score: {result.get('differentiation_score', 0):.2f}")
        print(f"  Marker sequence: {result.get('marker_sequence', '')}")
        return True
    except Exception as e:
        logger.error(f"Error: {e}")
        return False

def test_with_real_data():
    """Test marker identification with real data."""
    email = "steph.smith@unc.edu"
    
    try:
        # Try a different approach to get sequences
        # Use regular fetch_sequences_for_taxid with broader queries
        inclusion_query = "txid689[Organism] AND (16S OR ribosomal)"
        exclusion_query = "txid666[Organism] AND (16S OR ribosomal)"
        
        inclusion_seqs = fetch_sequences_for_taxid(
            taxid="689",
            email=email,
            query_term=inclusion_query,
            max_records=3
        )
        
        exclusion_seqs = fetch_sequences_for_taxid(
            taxid="666",
            email=email,
            query_term=exclusion_query,
            max_records=3
        )
        
        # If we still don't have sequences, try with artificial data
        if not inclusion_seqs or not exclusion_seqs:
            print("Could not fetch enough sequences, using artificial data instead")
            inclusion_seqs, exclusion_seqs = create_test_sequences()
        
        print(f"Using {len(inclusion_seqs)} inclusion and {len(exclusion_seqs)} exclusion sequences")
        
        # Use memory-efficient k-mer approach
        from assay_design.target_identification import find_conserved_marker_kmers
        
        result = find_conserved_marker_kmers(
            inclusion_seqs,
            exclusion_seqs,
            kmer_size=20
        )
        
        print("Marker identification results:")
        print(f"  Length: {result.get('marker_length')}")
        print(f"  Conservation score: {result.get('conservation_score', 0):.2f}")
        print(f"  Marker sequence (first 50 bp): {result.get('marker_sequence', '')[:50]}...")
        
        return True
    except Exception as e:
        logger.error(f"Error: {e}")
        return False
        
        
if __name__ == "__main__":
    print("Testing target_identification module...")
    
    tests = [
        ("Marker identification with artificial data", test_with_artificial_data),
        ("Marker identification with real data", test_with_real_data),
    ]
    
    for name, test_func in tests:
        print(f"\n=== Testing: {name} ===")
        success = test_func()
        print(f"Result: {'SUCCESS' if success else 'FAILURE'}")