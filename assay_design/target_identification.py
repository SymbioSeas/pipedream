# assay_design/target_identification.py

import logging
import time
from typing import List, Dict, Any, Optional, Set
from collections import Counter
from Bio.SeqRecord import SeqRecord

logger = logging.getLogger(__name__)

from .kmer_analysis import kmer_analysis, extract_kmers_from_chunk
from .primer_probe_design import find_conserved_kmers

def process_long_sequences(sequences, window_size=100000, step_size=50000, max_windows=20):
    """
    Process long sequences using a sliding window approach.
    
    Args:
        sequences: Original long sequences
        window_size: Size of each window to analyze
        step_size: Distance between consecutive windows
        max_windows: Maximum number of windows to process per sequence
        
    Returns:
        List of sequence windows for analysis
    """
    processed_sequences = []
    
    for seq_record in sequences:
        seq_length = len(seq_record.seq)
        
        # Only process if sequence is longer than window_size
        if seq_length > window_size:
            window_count = 0
            
            # Create overlapping windows
            for start in range(0, seq_length - window_size + 1, step_size):
                if window_count >= max_windows:
                    break
                    
                window_seq = seq_record[start:start + window_size]
                window_seq.id = f"{seq_record.id}_window_{start}-{start+window_size}"
                processed_sequences.append(window_seq)
                window_count += 1
        else:
            # Keep short sequences as is
            processed_sequences.append(seq_record)
    
    return processed_sequences
    
def extract_gene_regions(sequences, gene_name="16S", max_length=10000):
    """
    Extract specific gene regions from sequences based on annotation.
    
    Args:
        sequences: Original sequences
        gene_name: Target gene to extract
        max_length: Maximum length to consider if gene is too long
        
    Returns:
        List of extracted gene regions
    """
    extracted_regions = []
    
    for seq_record in sequences:
        # Look for gene features in annotations
        if hasattr(seq_record, 'features'):
            for feature in seq_record.features:
                if feature.type == "gene" and gene_name.lower() in str(feature.qualifiers.get('gene', '')).lower():
                    # Extract the gene region
                    start = feature.location.start
                    end = feature.location.end
                    
                    # Apply length limit if needed
                    if end - start > max_length:
                        mid_point = (start + end) // 2
                        start = max(0, mid_point - max_length // 2)
                        end = min(len(seq_record), mid_point + max_length // 2)
                    
                    gene_region = seq_record[start:end]
                    gene_region.id = f"{seq_record.id}_{gene_name}_{start}-{end}"
                    extracted_regions.append(gene_region)
    
    return extracted_regions
    
def adaptive_preprocessing(sequences, target_base_count=10000000):
    """
    Adaptively determine sequence length limit based on number of sequences.
    
    Args:
        sequences: Original sequences
        target_base_count: Target total bases to analyze
        
    Returns:
        Processed sequences with adaptive length limits
    """
    sequence_count = len(sequences)
    if sequence_count == 0:
        return []
    
    # Calculate target length per sequence
    target_length = target_base_count // sequence_count
    
    # Set minimum and maximum bounds
    min_length = 1000  # Ensure at least 1kb per sequence
    max_length = 1000000  # Cap at 1Mb per sequence
    
    sequence_length = max(min_length, min(max_length, target_length))
    
    # Process sequences with calculated length
    processed_sequences = []
    for seq_record in sequences:
        if len(seq_record.seq) > sequence_length:
            # Take sequence from the middle to avoid ends
            mid_point = len(seq_record.seq) // 2
            start = max(0, mid_point - sequence_length // 2)
            end = min(len(seq_record.seq), mid_point + sequence_length // 2)
            
            trimmed_seq = seq_record[start:end]
            trimmed_seq.id = f"{seq_record.id}_trimmed_{start}-{end}"
            processed_sequences.append(trimmed_seq)
        else:
            processed_sequences.append(seq_record)
    
    return processed_sequences
    
def find_markers_multi_pass(inclusion_sequences, exclusion_sequences):
    """
    Multi-pass approach to find markers, starting with small regions.
    
    Args:
        inclusion_sequences: Target sequences
        exclusion_sequences: Sequences to exclude
        
    Returns:
        Best markers found
    """
    # Try with small windows first (faster)
    small_inc_seqs = process_sequences(inclusion_sequences, max_length=5000)
    small_exc_seqs = process_sequences(exclusion_sequences, max_length=5000)
    
    result = find_discriminative_kmers(small_inc_seqs, small_exc_seqs)
    
    # If no good markers found, try with medium windows
    if result.get("conservation_score", 0) < 0.9:
        medium_inc_seqs = process_sequences(inclusion_sequences, max_length=20000)
        medium_exc_seqs = process_sequences(exclusion_sequences, max_length=20000)
        
        result = find_discriminative_kmers(medium_inc_seqs, medium_exc_seqs)
    
    # If still no good markers, try with large windows
    if result.get("conservation_score", 0) < 0.8:
        large_inc_seqs = process_sequences(inclusion_sequences, max_length=100000)
        large_exc_seqs = process_sequences(exclusion_sequences, max_length=100000)
        
        result = find_discriminative_kmers(large_inc_seqs, large_exc_seqs)
    
    return result
    
def find_discriminative_kmers(
    inclusion_sequences: List[SeqRecord],
    exclusion_sequences: List[SeqRecord],
    kmer_size: int = 40,
    min_conservation: float = 0.95,
    max_sequences: int = 100,
    max_sequence_length: int = 100000,
    use_optimized:  bool = True
) -> Dict[str, Any]:
    """
    Find conserved k-mers present in inclusion sequences but not in exclusion sequences.
    
    Args:
        inclusion_sequences: Sequences from target organisms
        exclusion_sequences: Sequences from non-target organisms
        kmer_size: Size of k-mers to analyze
        min_conservation: Minimum conservation level (0-1)
        max_sequences: Maximum number of sequences to process
        max_sequence_length: Maximum length of sequence to consider
        use_optimized: Whether to use the optimized k-mer analysis (do it)
        
    Returns:
        Dict containing marker information
    """
    start_time = time.time()
    logger.info(f"Finding discriminative k-mers of size {kmer_size}...")
    
    # Preprocess sequences to limit dataset size
    inc_seqs = preprocess_sequences(inclusion_sequences, max_sequences, max_sequence_length)
    exc_seqs = preprocess_sequences(exclusion_sequences, max_sequences, max_sequence_length)
    
    if not inc_seqs:
        return {"error": "No valid inclusion sequences after preprocessing"}
    
    logger.info(f"Analyzing {len(inc_seqs)} inclusion and {len(exc_seqs)} exclusion sequences")
    
    if use_optimized:
        conserved_specific_kmers = kmer_analysis(
            inc_seqs,
            exc_seqs,
            kmer_size=kmer_size,
            min_conservation=min_conservation
        )
    else:
    
        # Extract k-mers from inclusion sequences with counts
        inclusion_kmers = Counter()
        for seq in inc_seqs:
            seq_str = str(seq.seq).upper()
            for i in range(len(seq_str) - kmer_size + 1):
                kmer = seq_str[i:i+kmer_size]
                if 'N' not in kmer and '-' not in kmer:  # Skip k-mers with gaps or ambiguous bases
                    inclusion_kmers[kmer] += 1
                
    # Added check for empty exclusion_sequences
    if not exclusion_sequences:
        logger.info("No exclusion sequences provided - finding only conserved regions")
        exclusion_kmers = set() # Empty set of exclusion k-mers
    else:
        # Extract k-mers from exclusion sequences
        exclusion_kmers = set()
        for seq in exc_seqs:
            seq_str = str(seq.seq).upper()
            for i in range(len(seq_str) - kmer_size + 1):
                kmer = seq_str[i:i+kmer_size]
                if 'N' not in kmer and '-' not in kmer:
                    exclusion_kmers.add(kmer)
        
    # Find conserved k-mers specific to inclusion sequences
    conserved_specific_kmers = {}
    conservation_threshold = min_conservation * len(inc_seqs)
    
    for kmer, count in inclusion_kmers.items():
        if count >= conservation_threshold and kmer not in exclusion_kmers:
            conservation = count / len(inc_seqs)
            conserved_specific_kmers[kmer] = conservation
    
    logger.info(f"Found {len(conserved_specific_kmers)} conserved specific k-mers")
    
    # No discriminative k-mers found, try the best conserved k-mer
    if not conserved_specific_kmers:
        logger.info("No specific k-mers found, using most conserved k-mer")
        best_conserved = None
        best_conservation = 0
        
        for kmer, count in inclusion_kmers.items():
            conservation = count / len(inc_seqs)
            if conservation > best_conservation:
                best_conserved = kmer
                best_conservation = conservation
        
        if best_conserved:
            return {
                "marker_sequence": best_conserved,
                "marker_length": kmer_size,
                "conservation_score": best_conservation,
                "specificity_score": 0,
                "description": f"Conserved but non-specific marker (present in {best_conservation:.2f} of inclusion sequences)"
            }
        else:
            return {"error": "No suitable markers found"}
    
    # Find the best k-mer (highest conservation)
    best_kmer = max(conserved_specific_kmers.items(), key=lambda x: x[1])
    kmer, conservation = best_kmer
    
    # Check if we should extend the k-mer
    if kmer_size < 100:  # Only try extension for smaller k-mers
        logger.info(f"Extending best k-mer (conservation: {conservation:.2f})")
        extended_kmer = extend_kmer(kmer, inc_seqs, exc_seqs)
        if len(extended_kmer) > len(kmer):
            kmer = extended_kmer
    
    # Calculate specificity as 1 - (presence in exclusion / total exclusion)
    specificity = 1.0  # Fully specific
    
    elapsed_time = time.time() - start_time
    logger.info(f"Found marker in {elapsed_time:.2f} seconds")
    
    return {
        "marker_sequence": kmer,
        "marker_length": len(kmer),
        "conservation_score": conservation,
        "specificity_score": specificity,
        "description": f"Specific marker region (present in {conservation:.2f} of inclusion sequences, absent in exclusion)"
    }

def preprocess_sequences(
    sequences: List[SeqRecord], 
    max_sequences: int = 50, 
    max_length: int = 2000000
) -> List[SeqRecord]:
    """
    Preprocess sequences to limit dataset size.
    
    Args:
        sequences: Input sequences
        max_sequences: Maximum number of sequences to keep
        max_length: Maximum length of each sequence
        
    Returns:
        Processed sequences
    """
    if not sequences:
        return []
    
    # Limit number of sequences
    if len(sequences) > max_sequences:
        logger.info(f"Limiting from {len(sequences)} to {max_sequences} sequences")
        sequences = sequences[:max_sequences]
    
    # Process each sequence
    result = []
    for seq in sequences:
        # Trim sequence if too long
        if len(seq.seq) > max_length:
            logger.info(f"Trimming sequence {seq.id} from {len(seq.seq)} to {max_length} bp")
            from Bio.SeqRecord import SeqRecord
            from Bio.Seq import Seq
            
            trimmed_seq = SeqRecord(
                Seq(str(seq.seq[:max_length])),
                id=seq.id,
                description=f"{seq.description} (trimmed)"
            )
            result.append(trimmed_seq)
        else:
            result.append(seq)
    
    return result

def find_conserved_without_exclusion(
    inclusion_sequences: List[SeqRecord],
    min_conservation: float = 0.8,
    max_amplicon_length: int = 200,
    min_region_size: int = 50,
    timeout_seconds: int = 120
) -> Dict[str, Any]:
    """
    Find conserved regions across diverse taxa without exclusion criteria.
    Optimized for conserved functional gene assay design.
    
    Args:
        inclusion_sequences: Sequences from diverse taxa for the target gene
        min_conservation: Minimum conservation level (0-1)
        max_amplicon_length: Maximum amplicon length
        min_region_size: Minimum region size
        timeout_seconds: Maximum runtime in seconds
        
    Returns:
        Dict containing marker information
    """
    start_time = time.time()
    logger.info(f"Finding conserved regions across {len(inclusion_sequences)} sequences")
    
    # If very few sequences, adjust conservation threshold
    if len(inclusion_sequences) < 5:
        min_conservation = max(0.7, min_conservation * 0.9)
        logger.info(f"Small sample size, adjusting conservation threshold to {min_conservation}")
    
    # Prepare sequences for analysis (alignments would be better here, but for simplicity...)
    processed_sequences = preprocess_sequences(
        inclusion_sequences, 
        max_sequences=100,  # Allow more sequences for better conservation analysis
        max_length=10000   # Keep sequences relatively short
    )
    
    # Different strategies for finding conserved regions
    conserved_regions = []
    kmer_sizes = [30, 25, 20, 18, 15, 12]  # Try different k-mer sizes
    
    for kmer_size in kmer_sizes:
        if time.time() - start_time > timeout_seconds:
            logger.warning(f"Timeout reached after {timeout_seconds} seconds")
            break
            
        # Find conserved k-mers
        conserved_kmers = find_conserved_kmers(
            processed_sequences,
            min_conservation=min_conservation,
            kmer_size=kmer_size
        )
        
        logger.info(f"Found {len(conserved_kmers)} conserved {kmer_size}-mers")
        
        # Sort by conservation level
        sorted_kmers = sorted(
            [(kmer, conservation) for kmer, conservation in conserved_kmers.items()],
            key=lambda x: x[1],
            reverse=True
        )
        
        # Take the top conserved k-mers
        top_kmers = sorted_kmers[:50]  # Limit to top 50 for processing
        
        for kmer, conservation in top_kmers:
            # Try to extend this k-mer
            extended_kmer = extend_conserved_kmer(
                kmer,
                processed_sequences,
                max_extension=max_amplicon_length - kmer_size,
                min_conservation=min_conservation * 0.95  # Slightly relaxed for extension
            )
            
            # Add to candidate regions
            conserved_regions.append({
                "sequence": extended_kmer,
                "length": len(extended_kmer),
                "conservation": conservation,
                "original_kmer": kmer
            })
    
    # If no conserved regions found, fallback to using the most conserved k-mer directly
    if not conserved_regions:
        # Try with a smaller k-mer size and lower conservation threshold
        conserved_kmers = find_conserved_kmers(
            processed_sequences,
            min_conservation=min_conservation * 0.8,  # Lower threshold for fallback
            kmer_size=20
        )
        
        if conserved_kmers:
            best_kmer, conservation = max(conserved_kmers.items(), key=lambda x: x[1])
            logger.info(f"Using fallback approach with conservation {conservation:.2f}")
            
            return {
                "marker_sequence": best_kmer,
                "marker_length": len(best_kmer),
                "conservation_score": conservation,
                "specificity_score": 1.0,  # No exclusion criteria to evaluate specificity
                "description": f"Conserved region across diverse taxa (conservation: {conservation:.2f})"
            }
    
    # Select the best region based on length and conservation
    if conserved_regions:
        # Sort by a combined score favoring longer regions with high conservation
        conserved_regions.sort(
            key=lambda x: (
                x["conservation"] * 0.7 +  # Weight for conservation
                min(1.0, x["length"] / max_amplicon_length) * 0.3  # Weight for length
            ),
            reverse=True
        )
        
        best_region = conserved_regions[0]
        
        return {
            "marker_sequence": best_region["sequence"],
            "marker_length": best_region["length"],
            "conservation_score": best_region["conservation"],
            "specificity_score": 1.0,  # No exclusion criteria
            "description": f"Conserved region across diverse taxa (conservation: {best_region['conservation']:.2f})"
        }
    
    # If all else fails
    return {
        "error": "Could not identify sufficiently conserved regions",
        "marker_sequence": "",
        "marker_length": 0
    }

def extend_conserved_kmer(
    kmer: str,
    sequences: List[SeqRecord],
    max_extension: int = 150,
    min_conservation: float = 0.8
) -> str:
    """
    Extend a conserved k-mer to create a longer conserved marker.
    
    Args:
        kmer: Starting k-mer
        sequences: Sequences to analyze
        max_extension: Maximum length to extend
        min_conservation: Minimum conservation for extension
        
    Returns:
        Extended conserved sequence
    """
    # Find positions of k-mer in sequences
    positions = []
    for seq in sequences:
        seq_str = str(seq.seq).upper()
        pos = seq_str.find(kmer)
        if pos >= 0:
            positions.append((seq, pos))
    
    if not positions:
        return kmer
        
    # Start with original k-mer
    extended = kmer
    extensions_made = 0
    
    # Try extending to the right
    while extensions_made < max_extension:
        # Look at next position
        next_pos = len(extended)
        
        # Count bases at this position
        bases_count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        total_valid = 0
        
        for seq, start_pos in positions:
            if start_pos + next_pos < len(seq.seq):
                base = seq.seq[start_pos + next_pos].upper()
                if base in bases_count:
                    bases_count[base] += 1
                    total_valid += 1
        
        if total_valid == 0:
            break
            
        # Find most common base
        most_common = max(bases_count.items(), key=lambda x: x[1])
        base, count = most_common
        
        # Check conservation
        conservation = count / total_valid
        if conservation >= min_conservation:
            extended += base
            extensions_made += 1
        else:
            break
    
    # Try extending to the left
    left_extensions = 0
    left_extended = extended
    
    while left_extensions < max_extension:
        # Look at previous position
        prev_pos = -1 - left_extensions
        
        # Count bases at this position
        bases_count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        total_valid = 0
        
        for seq, start_pos in positions:
            if start_pos + prev_pos >= 0:
                base = seq.seq[start_pos + prev_pos].upper()
                if base in bases_count:
                    bases_count[base] += 1
                    total_valid += 1
        
        if total_valid == 0:
            break
            
        # Find most common base
        most_common = max(bases_count.items(), key=lambda x: x[1])
        base, count = most_common
        
        # Check conservation
        conservation = count / total_valid
        if conservation >= min_conservation:
            left_extended = base + left_extended
            left_extensions += 1
        else:
            break
    
    return left_extended
    
def extend_conserved_kmer(
    kmer: str,
    sequences: List[SeqRecord],
    max_extension: int = 150,
    min_conservation: float = 0.8
) -> str:
    """
    Extend a conserved k-mer to create a longer conserved marker.
    
    Args:
        kmer: Starting k-mer
        sequences: Sequences to analyze
        max_extension: Maximum length to extend
        min_conservation: Minimum conservation for extension
        
    Returns:
        Extended conserved sequence
    """
    # Find positions of k-mer in sequences
    positions = []
    for seq in sequences:
        seq_str = str(seq.seq).upper()
        pos = seq_str.find(kmer)
        if pos >= 0:
            positions.append((seq, pos))
    
    if not positions:
        return kmer
        
    # Start with original k-mer
    extended = kmer
    extensions_made = 0
    
    # Try extending to the right
    while extensions_made < max_extension:
        # Look at next position
        next_pos = len(extended)
        
        # Count bases at this position
        bases_count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        total_valid = 0
        
        for seq, start_pos in positions:
            if start_pos + next_pos < len(seq.seq):
                base = seq.seq[start_pos + next_pos].upper()
                if base in bases_count:
                    bases_count[base] += 1
                    total_valid += 1
        
        if total_valid == 0:
            break
            
        # Find most common base
        most_common = max(bases_count.items(), key=lambda x: x[1])
        base, count = most_common
        
        # Check conservation
        conservation = count / total_valid
        if conservation >= min_conservation:
            extended += base
            extensions_made += 1
        else:
            break
    
    # Try extending to the left
    left_extensions = 0
    left_extended = extended
    
    while left_extensions < max_extension:
        # Look at previous position
        prev_pos = -1 - left_extensions
        
        # Count bases at this position
        bases_count = {'A': 0, 'C': 0, 'G': 0, 'T': 0}
        total_valid = 0
        
        for seq, start_pos in positions:
            if start_pos + prev_pos >= 0:
                base = seq.seq[start_pos + prev_pos].upper()
                if base in bases_count:
                    bases_count[base] += 1
                    total_valid += 1
        
        if total_valid == 0:
            break
            
        # Find most common base
        most_common = max(bases_count.items(), key=lambda x: x[1])
        base, count = most_common
        
        # Check conservation
        conservation = count / total_valid
        if conservation >= min_conservation:
            left_extended = base + left_extended
            left_extensions += 1
        else:
            break
    
    return left_extended


def extend_kmer(
    kmer: str, 
    inclusion_seqs: List[SeqRecord], 
    exclusion_seqs: List[SeqRecord],
    max_extension: int = 100,
    min_conservation: float = 0.85
) -> str:
    """
    Extend k-mer to create a longer marker sequence.
    
    Args:
        kmer: Starting k-mer sequence
        inclusion_seqs: Inclusion sequences
        exclusion_seqs: Exclusion sequences
        max_extension: Maximum extension length
        min_conservation: Minimum conservation for extension
        
    Returns:
        Extended marker sequence
    """
    # Convert exclusion sequences to strings for faster lookup
    exclusion_strs = [str(seq.seq).upper() for seq in exclusion_seqs]
    
    # Find positions of the k-mer in inclusion sequences
    positions = []
    for seq in inclusion_seqs:
        seq_str = str(seq.seq).upper()
        # Find all occurrences of k-mer
        pos = seq_str.find(kmer)
        while pos >= 0:
            positions.append((seq, pos))
            pos = seq_str.find(kmer, pos + 1)
    
    if not positions:
        return kmer
    
    # Start with the original k-mer
    extended = kmer
    
    # Try to extend to the left
    extension_count = 0
    left_pos = 0
    
    while extension_count < max_extension:
        left_pos -= 1
        extension_count += 1
        
        # Count bases at this position
        bases = Counter()
        valid_positions = 0
        
        for seq, pos in positions:
            if pos + left_pos >= 0:
                base = seq.seq[pos + left_pos].upper()
                if base not in 'N-':  # Skip gaps and ambiguous bases
                    bases[base] += 1
                    valid_positions += 1
        
        if valid_positions == 0:
            break
        
        # Find most common base
        if not bases:
            break
            
        common_base, count = bases.most_common(1)[0]
        conservation = count / valid_positions if valid_positions > 0 else 0
        
        # Check if extension is conserved enough
        if conservation >= min_conservation:
            # Check if extended sequence appears in exclusion
            test_extension = common_base + extended
            is_specific = all(test_extension not in exc_str for exc_str in exclusion_strs)
            
            if is_specific:
                extended = test_extension
            else:
                break  # Stop if extension loses specificity
        else:
            break
    
    # Try to extend to the right
    extension_count = 0
    right_pos = len(extended)
    
    while extension_count < max_extension:
        extension_count += 1
        
        # Count bases at this position
        bases = Counter()
        valid_positions = 0
        
        for seq, pos in positions:
            if pos + right_pos < len(seq.seq):
                base = seq.seq[pos + right_pos].upper()
                if base not in 'N-':
                    bases[base] += 1
                    valid_positions += 1
        
        if valid_positions == 0:
            break
        
        # Find most common base
        if not bases:
            break
            
        common_base, count = bases.most_common(1)[0]
        conservation = count / valid_positions if valid_positions > 0 else 0
        
        # Check if extension is conserved enough
        if conservation >= min_conservation:
            # Check if extended sequence appears in exclusion
            test_extension = extended + common_base
            is_specific = all(test_extension not in exc_str for exc_str in exclusion_strs)
            
            if is_specific:
                extended = test_extension
                right_pos += 1
            else:
                break  # Stop if extension loses specificity
        else:
            break
    
    return extended

def find_optimal_marker(
    inclusion_sequences: List[SeqRecord],
    exclusion_sequences: List[SeqRecord],
    timeout_seconds: int = 60
) -> Dict[str, Any]:
    """
    Find the optimal marker with timeout.
    
    Args:
        inclusion_sequences: Sequences from target organisms
        exclusion_sequences: Sequences from non-target organisms
        timeout_seconds: Maximum execution time in seconds
        
    Returns:
        Dict containing marker information
    """
    start_time = time.time()
    logger.info("Finding optimal marker with timeout...")
    
    # Try different k-mer sizes to optimize results
    best_result = None
    best_score = 0
    
    # Try k-mer sizes in this order (prioritize medium k-mers)
    kmer_sizes = [25, 20, 30, 15, 35]
    
    for kmer_size in kmer_sizes:
        # Check timeout
        if time.time() - start_time > timeout_seconds:
            logger.info(f"Timeout reached after {timeout_seconds} seconds")
            break
        
        try:
            result = find_discriminative_kmers(
                inclusion_sequences,
                exclusion_sequences,
                kmer_size=kmer_size
            )
            
            if "error" in result:
                continue
                
            # Calculate combined score (conservation * specificity * length_factor)
            length_factor = min(1.0, result.get("marker_length", 0) / 100)
            score = (
                result.get("conservation_score", 0) * 
                result.get("specificity_score", 0) * 
                length_factor
            )
            
            if score > best_score:
                best_score = score
                best_result = result
                
            # If we found a good marker, stop early
            if score > 0.8:
                logger.info(f"Found high-quality marker with score {score:.2f}")
                break
                
        except Exception as e:
            logger.error(f"Error with k-mer size {kmer_size}: {e}")
    
    if best_result:
        elapsed_time = time.time() - start_time
        logger.info(f"Marker identification completed in {elapsed_time:.2f} seconds")
        return best_result
    else:
        # Fallback to a simple strategy
        logger.info("No optimal marker found, using fallback strategy")
        return fallback_marker_strategy(inclusion_sequences)

def fallback_marker_strategy(sequences: List[SeqRecord]) -> Dict[str, Any]:
    """
    Fallback strategy when optimal marker identification fails.
    
    Args:
        sequences: Inclusion sequences
        
    Returns:
        Dict containing marker information
    """
    if not sequences:
        return {"error": "No sequences provided"}
    
    # Use a region from the first sequence
    seq = sequences[0]
    length = min(150, len(seq.seq))
    
    marker = str(seq.seq[:length])
    
    return {
        "marker_sequence": marker,
        "marker_length": length,
        "conservation_score": 1.0 / len(sequences),
        "specificity_score": 0.5,  # Unknown specificity
        "description": "Fallback marker (limited evaluation)"
    }