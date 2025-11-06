#!/usr/bin/env python3

import logging
import os
from assay_design.data_retrieval import fetch_sequence_by_accession, load_local_fasta

# Configure basic logging
logging.basicConfig(level=logging.INFO)

def main():
    # Example: Fetch an E. coli K-12 MG1655 genome by accession:
    # You must provide your own valid email!
    my_email = "steph.smith@unc.edu"
    ecoli_accession = "NC_00913.3"
    
    # Remote fetch example
    seq_record = fetch_sequence_by_accession(ecoli_accession, my_email)
    print(f"Fetched {seq_record.id} of length {len(seq_record.seq)}")
    
    # Local FASTA example
    # Assumes a local FASTA in the same folder named ecoli_local.fasta
    local_path = os.path.join(os.path.dirname(__file__), "ecoli_local.fasta")
    # If the file exists, load it:
    if os.path.exists(local_path):
        local_sequences = load_local_fasta(local_path)
        print(f"Loaded {len(local_sequences)} local sequences from {local_path}")
    else:
        print(f"No local FASTA found at {local_path}. Skipping local load.")
        
if __name__ == "__main__":
    main()