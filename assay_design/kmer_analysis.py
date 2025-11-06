import multiprocessing as mp
import time
from typing import List, Dict, Set, Tuple
from Bio.SeqRecord import SeqRecord
try:
    from pybloom_live import BloomFilter
    BLOOM_FILTER_AVAILABLE = True
except ImportError:
    BLOOM_FILTER_AVAILABLE = False

def extract_kmers_from_chunk(chunk_data: Tuple) -> Set[str]:
    """Extract k-mers from a chunk of sequences using rolling hash."""
    sequences, kmer_size = chunk_data
    kmers = set()
    
    # Map bases to integers for faster hash computation
    base_map = {'A': 0, 'C': 1, 'G': 2, 'T': 3, 'N': 4, '-': 5}
    
    for seq_record in sequences:
        seq_str = str(seq_record.seq).upper()
        seq_len = len(seq_str)
        
        if seq_len < kmer_size:
            continue
            
        # Process k-mers with rolling hash
        rolling_hash = 0
        for i in range(kmer_size):
            base = seq_str[i]
            rolling_hash = (rolling_hash * 4 + base_map.get(base, 5)) % (2**31 - 1)
            
        # Get first k-mer
        first_kmer = seq_str[:kmer_size]
        if 'N' not in first_kmer and '-' not in first_kmer:
            kmers.add(first_kmer)
            
        # Use rolling hash for all other k-mers
        for i in range(1, seq_len - kmer_size + 1):
            # Update rolling hash
            out_base = seq_str[i-1]
            out_val = base_map.get(out_base, 5)
            rolling_hash = (rolling_hash - out_val * (4**(kmer_size-1))) % (2**31 - 1)
            
            in_base = seq_str[i+kmer_size-1]
            in_val = base_map.get(in_base, 5)
            rolling_hash = (rolling_hash * 4 + in_val) % (2**31 - 1)
            
            kmer = seq_str[i:i+kmer_size]
            if 'N' not in kmer and '-' not in kmer:
                kmers.add(kmer)
    
    return kmers

def process_inclusion_chunk(chunk_data: Tuple) -> Dict[str, int]:
    """
    Process a chunk of inclusion sequences to find k-mers not in exclusion filter.
    """
    sequences, kmer_size, exclusion_bloom = chunk_data
    local_kmer_counts = {}
    
    for seq in sequences:
        seq_str = str(seq.seq).upper()
        seq_kmers = set()  # Only count each k-mer once per sequence
        
        for i in range(0, len(seq_str) - kmer_size + 1, 2):  # Step by 2 for speed
            kmer = seq_str[i:i+kmer_size]
            if 'N' not in kmer and '-' not in kmer:
                # Check if kmer is not in exclusion set
                if kmer not in exclusion_bloom:
                    seq_kmers.add(kmer)
        
        # Update local counts
        for kmer in seq_kmers:
            local_kmer_counts[kmer] = local_kmer_counts.get(kmer, 0) + 1
    
    return local_kmer_counts

def kmer_analysis(
    inclusion_sequences: List[SeqRecord],
    exclusion_sequences: List[SeqRecord],
    kmer_size: int = 20,
    min_conservation: float = 0.8,
    n_processes: int = None
) -> Dict[str, float]:
    """
    Highly optimized k-mer analysis combining parallel processing,
    Bloom filters, and rolling hash.
    
    Args:
        inclusion_sequences: Target sequences
        exclusion_sequences: Non-target sequences
        kmer_size: Size of k-mers to analyze
        min_conservation: Minimum conservation level (0-1)
        n_processes: Number of processes to use (default: CPU count)
        
    Returns:
        Dict of specific conserved k-mers with conservation scores
    """
    
    # Check if dependencies are available
    if not BLOOM_FILTER_AVAILABLE:
        logger.warning("pybloom-live package not found. Falling back to standard k-mer analysis.")
        # Call the non-optimized version instead
        return find_conserved_kmers_parallel(
            inclusion_sequences,
            kmer_size,
            min_conservation,
            n_processes
        )
        
    start_time = time.time()
    
    # Determine number of processes
    if n_processes is None:
        n_processes = mp.cpu_count()
    
    # Step 1: Extract exclusion k-mers in parallel
    # Split exclusion sequences into chunks
    exclusion_chunks = []
    chunk_size = max(1, len(exclusion_sequences) // n_processes)
    
    for i in range(0, len(exclusion_sequences), chunk_size):
        chunk = exclusion_sequences[i:i+chunk_size]
        exclusion_chunks.append((chunk, kmer_size))
    
    # Process exclusion chunks in parallel
    with mp.Pool(processes=n_processes) as pool:
        exclusion_kmer_sets = pool.map(extract_kmers_from_chunk, exclusion_chunks)
    
    # Combine all exclusion k-mers
    all_exclusion_kmers = set()
    for kmer_set in exclusion_kmer_sets:
        all_exclusion_kmers.update(kmer_set)
    
    # Create bloom filter from exclusion k-mers
    exclusion_bloom = BloomFilter(capacity=len(all_exclusion_kmers), error_rate=0.001)
    for kmer in all_exclusion_kmers:
        exclusion_bloom.add(kmer)
    
    print(f"Exclusion Bloom filter created with {len(all_exclusion_kmers)} k-mers")
    
    # Step 2: Process inclusion sequences in parallel
    inclusion_chunks = []
    chunk_size = max(1, len(inclusion_sequences) // n_processes)
    
    for i in range(0, len(inclusion_sequences), chunk_size):
        chunk = inclusion_sequences[i:i+chunk_size]
        inclusion_chunks.append((chunk, kmer_size, exclusion_bloom))
    
    # Process inclusion chunks in parallel
    with mp.Pool(processes=n_processes) as pool:
        inclusion_results = pool.map(process_inclusion_chunk, inclusion_chunks)
    
    # Combine results from all processes
    combined_kmer_counts = {}
    for local_counts in inclusion_results:
        for kmer, count in local_counts.items():
            combined_kmer_counts[kmer] = combined_kmer_counts.get(kmer, 0) + count
    
    # Calculate conservation scores and filter
    conserved_specific_kmers = {}
    threshold = min_conservation * len(inclusion_sequences)
    
    for kmer, count in combined_kmer_counts.items():
        if count >= threshold:
            conserved_specific_kmers[kmer] = count / len(inclusion_sequences)
    
    elapsed_time = time.time() - start_time
    print(f"K-mer analysis completed in {elapsed_time:.2f} seconds")
    print(f"Found {len(conserved_specific_kmers)} conserved specific k-mers")
    
    return conserved_specific_kmers