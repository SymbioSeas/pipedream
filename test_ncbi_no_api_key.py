#!/usr/bin/env python3
"""
Test NCBI API without using an API key to diagnose API key issues.
"""

import sys
from Bio import Entrez

def test_without_api_key(email):
    """Test NCBI API connection WITHOUT any API key."""

    print("="*60)
    print("NCBI API Test - NO API KEY")
    print("="*60)

    # ONLY set email, explicitly NO API key
    Entrez.email = email
    print(f"\n1. Email set to: {email}")
    print("2. NO API KEY - Testing at 3 requests/second rate limit")

    # Test 1: Taxonomy query
    print("\n" + "-"*60)
    print("Test 1: Fetching taxonomy info for taxid 1854 (Frankia)")
    print("-"*60)
    try:
        handle = Entrez.efetch(db="taxonomy", id="1854", rettype="full", retmode="xml")
        records = Entrez.read(handle)
        handle.close()

        if records:
            tax_info = records[0]
            print(f"✓ SUCCESS: Retrieved taxonomy for {tax_info.get('ScientificName', 'Unknown')}")
            print(f"  Rank: {tax_info.get('Rank', 'Unknown')}")
            print(f"  Division: {tax_info.get('Division', 'Unknown')}")
        else:
            print("✗ FAILED: No records returned")

    except Exception as e:
        print(f"✗ FAILED: {type(e).__name__}: {str(e)}")

    # Test 2: Search query
    print("\n" + "-"*60)
    print("Test 2: Searching for sequences with taxid 1854")
    print("-"*60)
    try:
        query = "txid1854[Organism] AND biomol_genomic[PROP] NOT wgs[Filter]"
        print(f"Query: {query}")

        search_handle = Entrez.esearch(db="nucleotide", term=query, retmax=5)
        search_results = Entrez.read(search_handle)
        search_handle.close()

        ids = search_results.get("IdList", [])
        count = search_results.get("Count", 0)
        print(f"✓ SUCCESS: Found {count} total records")
        print(f"  First 5 IDs: {ids}")

    except Exception as e:
        print(f"✗ FAILED: {type(e).__name__}: {str(e)}")

    # Test 3: Fetch sequence
    print("\n" + "-"*60)
    print("Test 3: Fetching a known sequence (NC_000913.3 - E. coli)")
    print("-"*60)
    try:
        handle = Entrez.efetch(db="nucleotide", id="NC_000913.3", rettype="fasta", retmode="text")
        fasta_data = handle.read()
        handle.close()

        lines = fasta_data.strip().split('\n')
        print(f"✓ SUCCESS: Retrieved sequence")
        print(f"  Header: {lines[0][:80]}...")
        print(f"  Sequence length: {len(''.join(lines[1:]))} bp")

    except Exception as e:
        print(f"✗ FAILED: {type(e).__name__}: {str(e)}")

    print("\n" + "="*60)
    print("Test complete")
    print("="*60)

if __name__ == "__main__":
    email = sys.argv[1] if len(sys.argv) > 1 else "test@example.com"
    test_without_api_key(email)
