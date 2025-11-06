#!/usr/bin/env python3

import logging
from assay_design.data_retrieval import (
    fetch_sequence_by_accession,
    fetch_sequences_for_taxid,
    get_taxon_info,
    suggest_marker_genes
)

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def test_fetch_sequence():
    """Test fetching a single sequence by accession."""
    email = "your.email@example.com"  # Replace with your email
    accession = "NC_045512.2"  # SARS-CoV-2 reference genome
    
    try:
        seq_record = fetch_sequence_by_accession(accession, email)
        print(f"Successfully fetched {seq_record.id} of length {len(seq_record.seq)}")
        return True
    except Exception as e:
        logger.error(f"Error: {e}")
        return False

def test_fetch_sequences_for_taxid():
    """Test fetching sequences for a specific taxid."""
    email = "your.email@example.com"  # Replace with your email
    taxid = "562"  # E. coli
    
    try:
        sequences = fetch_sequences_for_taxid(taxid, email, max_records=3)
        print(f"Successfully fetched {len(sequences)} sequences for taxid {taxid}")
        for seq in sequences:
            print(f"  {seq.id}: {len(seq.seq)} bp")
        return True
    except Exception as e:
        logger.error(f"Error: {e}")
        return False

def test_get_taxon_info():
    """Test getting taxonomy information."""
    email = "your.email@example.com"  # Replace with your email
    taxid = "689"  # Vibrio mediterranei
    
    try:
        info = get_taxon_info(taxid, email)
        print(f"Successfully retrieved taxonomy info for {info.get('scientific_name')}")
        print(f"  Rank: {info.get('rank')}")
        print(f"  Lineage: {info.get('lineage')}")
        return True
    except Exception as e:
        logger.error(f"Error: {e}")
        return False

def test_suggest_marker_genes():
    """Test marker gene suggestions."""
    email = "your.email@example.com"  # Replace with your email
    taxid = "689"  # Vibrio mediterranei
    
    try:
        genes = suggest_marker_genes(taxid, email)
        print(f"Successfully suggested {len(genes)} marker genes:")
        for gene in genes:
            print(f"  {gene['gene']}: {gene['description']}")
        return True
    except Exception as e:
        logger.error(f"Error: {e}")
        return False

if __name__ == "__main__":
    print("Testing data_retrieval module...")
    
    tests = [
        ("Fetch sequence by accession", test_fetch_sequence),
        ("Fetch sequences for taxid", test_fetch_sequences_for_taxid),
        ("Get taxon info", test_get_taxon_info),
        ("Suggest marker genes", test_suggest_marker_genes),
    ]
    
    for name, test_func in tests:
        print(f"\n=== Testing: {name} ===")
        success = test_func()
        print(f"Result: {'SUCCESS' if success else 'FAILURE'}")