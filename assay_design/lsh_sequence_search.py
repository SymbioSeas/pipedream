# assay_design/lsh_sequence_search.py

import logging
import time
from typing import List, Dict, Any, Set, Tuple
import numpy as np
from datasketch import MinHash, MinHashLSH
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)

from .kmer_analysis import extract_kmers_from_chunk

class LSHSequenceSearch:
    """
    Locality-Sensitive Hashing for efficient sequence comparison.
    Uses MinHash to rapidly identify similar regions between sequences.
    """
    
    def __init__(
        self,
        kmer_size: int = 10,
        window_size: int = 200,
        step_size: int = 50,
        num_perm: int = 128,
        threshold: float = 0.7
    ):
        """
        Initialize the LSH sequence search.
        
        Args:
            kmer_size: Size of k-mers for hashing
            window_size: Size of sliding windows
            step_size: Step size between windows
            num_perm: Number of permutations for MinHash
            threshold: Similarity threshold (0-1)
        """
        self.kmer_size = kmer_size
        self.window_size = window_size
        self.step_size = step_size
        self.num_perm = num_perm
        self.threshold = threshold
        self.lsh = MinHashLSH(threshold=threshold, num_perm=num_perm)
        self.exclusion_windows = []
        
    def _sequence_to_kmers(self, sequence: str) -> Set[str]:
        """Convert a sequence to a set of k-mers."""
        return extract_kmers_from_chunk([(sequence, self.kmer_size)])
    
    def _create_minhash(self, kmers: Set[str]) -> MinHash:
        """Create a MinHash object from a set of k-mers."""
        minhash = MinHash(num_perm=self.num_perm)
        for kmer in kmers:
            minhash.update(kmer.encode('utf-8'))
        return minhash
    
    def index_exclusion_sequences(self, exclusion_sequences: List[SeqRecord]) -> None:
        """
        Index exclusion sequences for rapid similarity search.
        
        Args:
            exclusion_sequences: List of sequences to exclude
        """
        start_time = time.time()
        logger.info(f"Indexing {len(exclusion_sequences)} exclusion sequences with LSH")
        
        self.exclusion_windows = []
        window_count = 0
        
        # Process each exclusion sequence
        for seq_idx, seq_record in enumerate(exclusion_sequences):
            seq_str = str(seq_record.seq).upper()
            
            # Skip if sequence is too short
            if len(seq_str) < self.window_size:
                continue
                
            # Create sliding windows
            for start in range(0, len(seq_str) - self.window_size + 1, self.step_size):
                window = seq_str[start:start + self.window_size]
                
                # Convert window to k-mers
                kmers = self._sequence_to_kmers(window)
                if not kmers:
                    continue
                    
                # Create MinHash and add to LSH index
                minhash = self._create_minhash(kmers)
                window_id = f"excl_{seq_idx}_{start}"
                self.lsh.insert(window_id, minhash)
                
                # Store window information
                self.exclusion_windows.append({
                    "id": window_id,
                    "sequence_id": seq_record.id,
                    "start": start,
                    "end": start + self.window_size,
                    "minhash": minhash
                })
                
                window_count += 1
                
                # Log progress periodically
                if window_count % 1000 == 0:
                    elapsed = time.time() - start_time
                    logger.info(f"Processed {window_count} windows in {elapsed:.2f} seconds")
        
        elapsed = time.time() - start_time
        logger.info(f"Completed LSH indexing of {window_count} windows in {elapsed:.2f} seconds")
    
    def find_specific_regions(
        self,
        inclusion_sequences: List[SeqRecord],
        min_region_size: int = 50,
        max_region_size: int = 200,
        min_conservation: float = 0.8,
        max_results: int = 10,
        timeout_seconds: int = 60
    ) -> List[Dict[str, Any]]:
        """
        Find regions in inclusion sequences that don't match exclusion sequences.
        
        Args:
            inclusion_sequences: Target sequences
            min_region_size: Minimum region size to return
            max_region_size: Maximum region size to return
            min_conservation: Minimum conservation within inclusion sequences
            max_results: Maximum number of results to return
            timeout_seconds: Maximum runtime in seconds
            
        Returns:
            List of specific regions with their information
        """
        start_time = time.time()
        logger.info(f"Searching for specific regions in {len(inclusion_sequences)} inclusion sequences")
        
        # Ensure exclusion sequences have been indexed
        if not self.exclusion_windows:
            logger.warning("No exclusion windows indexed. All regions will be considered specific.")
        
        specific_regions = []
        candidate_regions = []
        
        # Process each inclusion sequence
        for seq_idx, seq_record in enumerate(inclusion_sequences):
            seq_str = str(seq_record.seq).upper()
            
            # Skip if sequence is too short
            if len(seq_str) < self.window_size:
                continue
            
            # Create sliding windows with variable sizes
            for window_size in range(min_region_size, min(max_region_size, len(seq_str)), 10):
                # Check timeout
                if time.time() - start_time > timeout_seconds:
                    logger.info(f"Timeout reached after {timeout_seconds} seconds")
                    break
                    
                # Process windows of this size
                for start in range(0, len(seq_str) - window_size + 1, self.step_size):
                    window = seq_str[start:start + window_size]
                    
                    # Convert window to k-mers
                    kmers = self._sequence_to_kmers(window)
                    if not kmers:
                        continue
                    
                    # Create MinHash
                    minhash = self._create_minhash(kmers)
                    
                    # Check for matches in exclusion sequences
                    matches = self.lsh.query(minhash)
                    
                    # If no matches, this region is potentially specific
                    if not matches:
                        candidate_regions.append({
                            "sequence_id": seq_record.id,
                            "start": start,
                            "end": start + window_size,
                            "length": window_size,
                            "sequence": window,
                            "minhash": minhash
                        })
                        
                        # Log progress periodically
                        if len(candidate_regions) % 100 == 0:
                            elapsed = time.time() - start_time
                            logger.info(f"Found {len(candidate_regions)} candidate regions in {elapsed:.2f} seconds")
                            
                        # Early termination if we have enough candidates
                        if len(candidate_regions) >= max_results * 10:
                            logger.info(f"Found {len(candidate_regions)} candidate regions, stopping search")
                            break
            
            # Check timeout after processing each sequence
            if time.time() - start_time > timeout_seconds:
                break
        
        # Verify conservation across inclusion sequences
        logger.info(f"Verifying conservation of {len(candidate_regions)} candidate regions")
        
        for region in candidate_regions:
            # Check if we've found enough regions
            if len(specific_regions) >= max_results:
                break
                
            # Check timeout
            if time.time() - start_time > timeout_seconds:
                logger.info(f"Timeout reached during conservation verification")
                break
                
            # Check conservation across inclusion sequences
            conservation_score = self._calculate_conservation(
                region["sequence"],
                inclusion_sequences
            )
            
            if conservation_score >= min_conservation:
                region["conservation_score"] = conservation_score
                region["specificity_score"] = 1.0  # Fully specific by LSH definition
                specific_regions.append(region)
                
                logger.info(f"Found specific region: {region['sequence_id']}:{region['start']}-{region['end']}, " +
                           f"conservation: {conservation_score:.2f}")
        
        # Sort by conservation score
        specific_regions.sort(key=lambda x: x["conservation_score"], reverse=True)
        
        elapsed = time.time() - start_time
        logger.info(f"Found {len(specific_regions)} specific regions in {elapsed:.2f} seconds")
        
        return specific_regions[:max_results]
    
    def _calculate_conservation(
        self,
        region_sequence: str,
        inclusion_sequences: List[SeqRecord],
        max_check: int = 10,
        similarity_threshold: float = 0.8
    ) -> float:
        """
        Calculate conservation of a region across inclusion sequences.
        
        Args:
            region_sequence: Sequence of the region
            inclusion_sequences: List of inclusion sequences
            max_check: Maximum number of sequences to check
            similarity_threshold: Threshold for considering a match
            
        Returns:
            Conservation score (0-1)
        """
        region_kmers = self._sequence_to_kmers(region_sequence)
        region_minhash = self._create_minhash(region_kmers)
        
        match_count = 0
        checked_count = 0
        
        # Check each inclusion sequence
        for seq_record in inclusion_sequences:
            seq_str = str(seq_record.seq).upper()
            
            # Skip if sequence is too short
            if len(seq_str) < len(region_sequence):
                continue
                
            checked_count += 1
            if checked_count > max_check:
                break
                
            # Search for similar regions
            has_match = False
            
            # Create sliding windows
            for start in range(0, len(seq_str) - len(region_sequence) + 1, self.step_size):
                window = seq_str[start:start + len(region_sequence)]
                
                # Calculate direct similarity (faster than MinHash for verification)
                similarity = self._calculate_direct_similarity(region_sequence, window)
                
                if similarity >= similarity_threshold:
                    has_match = True
                    break
            
            if has_match:
                match_count += 1
        
        # Return conservation score
        return match_count / checked_count if checked_count > 0 else 0
    
    def _calculate_direct_similarity(self, seq1: str, seq2: str) -> float:
        """
        Calculate direct sequence similarity using shared k-mers.
        
        Args:
            seq1: First sequence
            seq2: Second sequence
            
        Returns:
            Similarity score (0-1)
        """
        kmers1 = self._sequence_to_kmers(seq1)
        kmers2 = self._sequence_to_kmers(seq2)
        
        if not kmers1 or not kmers2:
            return 0.0
            
        # Jaccard similarity: |intersection| / |union|
        intersection = kmers1.intersection(kmers2)
        union = kmers1.union(kmers2)
        
        return len(intersection) / len(union)
    
    def format_results(self, specific_regions: List[Dict[str, Any]]) -> Dict[str, Any]:
        """
        Format specific regions into a standard marker info format.
        
        Args:
            specific_regions: List of specific regions
            
        Returns:
            Dict with marker information
        """
        if not specific_regions:
            return {
                "error": "No specific regions found",
                "marker_sequence": "",
                "marker_length": 0
            }
        
        # Take the best region
        best_region = specific_regions[0]
        
        # Format as marker info
        return {
            "marker_sequence": best_region["sequence"],
            "marker_length": best_region["length"],
            "conservation_score": best_region["conservation_score"],
            "specificity_score": best_region["specificity_score"],
            "description": f"LSH-identified specific region from {best_region['sequence_id']}:{best_region['start']}-{best_region['end']}",
            "forward_primer_region": best_region["sequence"][:25] if len(best_region["sequence"]) >= 25 else "",
            "probe_region": best_region["sequence"][len(best_region["sequence"])//3:len(best_region["sequence"])//3+25] if len(best_region["sequence"]) >= 75 else "",
            "reverse_primer_region": best_region["sequence"][-25:] if len(best_region["sequence"]) >= 25 else ""
        }


def lsh_find_markers(
    inclusion_sequences: List[SeqRecord],
    exclusion_sequences: List[SeqRecord],
    min_region_size: int = 100,
    max_region_size: int = 200,
    min_conservation: float = 0.8,
    timeout_seconds: int = 60
) -> Dict[str, Any]:
    """
    Find marker regions using Locality-Sensitive Hashing.
    
    Args:
        inclusion_sequences: Sequences from target organisms
        exclusion_sequences: Sequences from non-target organisms
        min_region_size: Minimum region size to return
        max_region_size: Maximum region size to return
        min_conservation: Minimum conservation within inclusion sequences
        timeout_seconds: Maximum runtime in seconds
        
    Returns:
        Dict containing marker information
    """
    start_time = time.time()
    logger.info("Finding marker regions with LSH...")
    
    # Create LSH searcher with appropriate parameters
    lsh_searcher = LSHSequenceSearch(
        kmer_size=10,
        window_size=max(200, min_region_size),
        step_size=25,
        num_perm=128,
        threshold=0.7
    )
    
    # Index exclusion sequences
    lsh_searcher.index_exclusion_sequences(exclusion_sequences)
    
    # Find specific regions
    specific_regions = lsh_searcher.find_specific_regions(
        inclusion_sequences=inclusion_sequences,
        min_region_size=min_region_size,
        max_region_size=max_region_size,
        min_conservation=min_conservation,
        max_results=5,
        timeout_seconds=timeout_seconds
    )
    
    # Format results
    result = lsh_searcher.format_results(specific_regions)
    
    elapsed = time.time() - start_time
    logger.info(f"LSH marker identification completed in {elapsed:.2f} seconds")
    
    return result