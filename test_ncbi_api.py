#!/usr/bin/env python3
"""
Diagnostic script to test NCBI API connection and API key validity.
"""

import sys
import os
from Bio import Entrez

def test_ncbi_connection(email, api_key=None):
    """Test NCBI API connection with and without API key."""

    print("="*60)
    print("NCBI API Diagnostic Test")
    print("="*60)

    # Set up Entrez
    Entrez.email = email
    print(f"\n1. Email set to: {email}")

    if api_key:
        # Clean the API key (remove whitespace, newlines)
        api_key_cleaned = api_key.strip()
        print(f"2. API Key provided: {api_key_cleaned[:10]}... (first 10 chars)")
        print(f"   API Key length: {len(api_key_cleaned)}")
        print(f"   Has whitespace: {api_key != api_key_cleaned}")

        Entrez.api_key = api_key_cleaned
    else:
        print("2. No API key provided")

    # Test 1: Simple taxonomy query
    print("\n" + "-"*60)
    print("Test 1: Fetching taxonomy info for taxid 1854 (Frankia)")
    print("-"*60)
    try:
        # Note: taxonomy database doesn't require retmode parameter
        handle = Entrez.efetch(db="taxonomy", id="1854")
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
        print(f"\nFull error:")
        import traceback
        traceback.print_exc()

    # Test 2: Simple search query
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
        print(f"\nFull error:")
        import traceback
        traceback.print_exc()

    # Test 3: Fetch a known sequence
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
        print(f"\nFull error:")
        import traceback
        traceback.print_exc()

    print("\n" + "="*60)
    print("Diagnostic test complete")
    print("="*60)

if __name__ == "__main__":
    # Get email from command line or use default
    email = sys.argv[1] if len(sys.argv) > 1 else "test@example.com"

    # Get API key from environment or command line
    api_key = os.environ.get('NCBI_API_KEY')
    if len(sys.argv) > 2:
        api_key = sys.argv[2]

    test_ncbi_connection(email, api_key)
