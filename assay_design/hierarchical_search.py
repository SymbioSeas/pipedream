# assay_design/hierarchical_search.py

import logging
from typing import List, Dict, Any, Optional
from Bio.SeqRecord import SeqRecord

from .data_retrieval import (
    fetch_sequences_for_taxid,
    fetch_gene_sequences,
    get_related_taxa,
    get_taxon_info
)
from .target_identification import (
    find_discriminative_kmers,
    preprocess_sequences
)
from .lsh_sequence_search import lsh_find_markers
from .kmer_analysis import kmer_analysis
from .gene_selection import (
    auto_select_gene_for_taxon,
    rank_candidate_genes
)

logger = logging.getLogger(__name__)

def hierarchical_marker_search(
    inclusion_taxid: str,
    email: str,
    gene_name: Optional[str] = None,
    exclusion_taxids: Optional[List[str]] = None,
    max_amplicon_length: int = 200,
    max_seq_count: int = 50,
    max_seq_length: int = 20000,
    min_conservation: float = 0.8,
    min_specificity: float = 0.7,
    timeout_seconds: int = 180,
    use_lsh: bool = True,
    binding_sites_only: bool = True,
    std_kmer: bool = True,
    parallel_processes: Optional[int] = 8,
    inclusion_sequences: Optional[List[SeqRecord]] = None,
    exclusion_sequences: Optional[List[SeqRecord]] = None,
    auto_gene_selection: bool = False,
    gene_use_case: str = 'quantification',
    max_genes_to_try: int = 3,
    min_gene_score: float = 0.4,
    api_key: Optional[str] = None
) -> Dict[str, Any]:
    """
    Perform a hierarchical search for marker regions.

    Args:
        inclusion_taxid (str): NCBI taxonomy ID for target organism
        email (str): Email for NCBI queries
        gene_name (Optional[str]): Specific gene to target (if None and auto_gene_selection=True, will auto-select)
        exclusion_taxids (Optional[List[str]]): List of explicit exclusion taxids
        max_seq_count (int): Maximum sequences to process per taxid
        max_seq_length (int): Maximum sequence length to consider
        min_conservation (float): Minimum conservation within inclusion sequences
        min_specificity (float): Minimum specificity compared to exclusion sequences
        timeout_seconds (int): Maximum runtime in seconds
        use_lsh (bool): Use Locality-Sensitive Hashing for faster searching
        binding_sites_only (bool): Use primer and probe binding sites to determine specificity
        std_kmer (bool): Use standard k-mer analysis
        parallel_processes (Optional[int]): Number of parallel processes for k-mer analysis
        inclusion_sequences (Optional[List[SeqRecord]]): Pre-loaded inclusion sequences
        exclusion_sequences (Optional[List[SeqRecord]]): Pre-loaded exclusion sequences
        auto_gene_selection (bool): Automatically select optimal gene if gene_name not provided
        gene_use_case (str): Use case for gene selection ('quantification', 'phylogeny', 'detection')
        max_genes_to_try (int): Maximum number of genes to try in fallback strategy
        min_gene_score (float): Minimum acceptable gene suitability score (0-1)
        api_key (Optional[str]): NCBI API key for higher rate limits

    Returns:
        Dict[str, Any]: Marker region information including:
            - marker_sequence: The identified marker sequence
            - marker_length: Length of the marker
            - conservation_score: Conservation score (0-1)
            - specificity_score: Specificity score (0-1)
            - selected_gene (if auto_gene_selection used): The gene that was selected
            - gene_selection_info (if auto_gene_selection used): Details about gene selection
    """
    import time
    start_time = time.time()

    # Step 0: Auto-select gene if requested
    gene_selection_info = None
    candidate_genes = []

    if auto_gene_selection and gene_name is None:
        logger.info(f"Auto-selecting optimal gene for taxid {inclusion_taxid} (use_case: {gene_use_case})")

        try:
            # Get ranked list of candidate genes
            ranked_genes = rank_candidate_genes(
                taxid=inclusion_taxid,
                email=email,
                prefer_single_copy=(gene_use_case == 'quantification'),
                max_genes_to_test=max_genes_to_try,
                timeout_per_gene=30,
                api_key=api_key
            )

            if ranked_genes:
                # Filter genes that meet minimum score threshold
                suitable_genes = [
                    (gene, info) for gene, info in ranked_genes
                    if info.get('overall_score', 0) >= min_gene_score
                ]

                if suitable_genes:
                    # Use the top-scoring gene
                    gene_name, gene_info = suitable_genes[0]
                    candidate_genes = [g for g, _ in suitable_genes[:max_genes_to_try]]

                    gene_selection_info = {
                        'selected_gene': gene_name,
                        'gene_score': gene_info.get('overall_score', 0),
                        'recommendation_tier': gene_info.get('recommendation_tier', 'unknown'),
                        'fallback_genes': candidate_genes[1:] if len(candidate_genes) > 1 else [],
                        'selection_method': 'auto_rank',
                        'use_case': gene_use_case
                    }

                    logger.info(f"Auto-selected gene '{gene_name}' with score {gene_info.get('overall_score', 0):.2f} "
                               f"({gene_info.get('recommendation_tier', 'unknown')} tier)")
                    if candidate_genes[1:]:
                        logger.info(f"Fallback genes available: {', '.join(candidate_genes[1:])}")
                else:
                    logger.warning(f"No genes met minimum score threshold ({min_gene_score})")
                    gene_selection_info = {
                        'error': f'No genes met minimum score threshold ({min_gene_score})',
                        'selection_method': 'auto_rank',
                        'use_case': gene_use_case
                    }
            else:
                logger.warning("Auto-gene selection found no suitable genes")
                gene_selection_info = {
                    'error': 'No suitable genes found for this taxon',
                    'selection_method': 'auto_rank',
                    'use_case': gene_use_case
                }

        except Exception as e:
            logger.error(f"Error during auto-gene selection: {str(e)}")
            gene_selection_info = {
                'error': f'Gene selection failed: {str(e)}',
                'selection_method': 'auto_rank',
                'use_case': gene_use_case
            }

    # Step 1: Use pre-loaded inclusion sequences if provided, otherwise fetch them
    if inclusion_sequences is not None:
        logger.info(f"Using {len(inclusion_sequences)} pre-loaded inclusion sequences")
    else:
        # Step 1: Fetch inclusion sequences
        logger.info(f"Fetching sequences for inclusion taxid {inclusion_taxid}")
        
        try:
            if gene_name:
                inclusion_sequences = fetch_gene_sequences(
                    taxid=inclusion_taxid,
                    gene_name=gene_name,
                    email=email,
                    max_records=max_seq_count * 2,
                    api_key=api_key
                )
            else:
                inclusion_sequences = fetch_sequences_for_taxid(
                    taxid=inclusion_taxid,
                    email=email,
                    max_records=max_seq_count * 2,
                    api_key=api_key
                )
            
            if not inclusion_sequences:
                error_result = {
                    "error": f"No sequences found for inclusion taxid {inclusion_taxid}" +
                             (f" with gene {gene_name}" if gene_name else "")
                }
                if gene_selection_info:
                    error_result['gene_selection_info'] = gene_selection_info

                # Try fallback gene if available
                if candidate_genes and len(candidate_genes) > 1:
                    logger.info(f"Trying fallback gene: {candidate_genes[1]}")
                    return hierarchical_marker_search(
                        inclusion_taxid=inclusion_taxid,
                        email=email,
                        gene_name=candidate_genes[1],  # Try next gene
                        exclusion_taxids=exclusion_taxids,
                        max_amplicon_length=max_amplicon_length,
                        max_seq_count=max_seq_count,
                        max_seq_length=max_seq_length,
                        min_conservation=min_conservation,
                        min_specificity=min_specificity,
                        timeout_seconds=timeout_seconds,
                        use_lsh=use_lsh,
                        binding_sites_only=binding_sites_only,
                        std_kmer=std_kmer,
                        parallel_processes=parallel_processes,
                        auto_gene_selection=False,  # Don't auto-select again
                        max_genes_to_try=max_genes_to_try - 1,  # Decrement
                        min_gene_score=min_gene_score,
                        api_key=api_key
                    )

                return error_result
            
            logger.info(f"Retrieved {len(inclusion_sequences)} sequences for inclusion taxid")
        
        except Exception as e:
            logger.error(f"Error fetching inclusion sequences: {str(e)}")
            error_result = {"error": f"Error fetching inclusion sequences: {str(e)}"}
            if gene_selection_info:
                error_result['gene_selection_info'] = gene_selection_info
            return error_result
    
    # Step 2: Find specific oligo sites within inclusion sequences only
    logger.info("Finding specific regions for primers and probe")
    
    # Preprocess sequences to limit size
    processed_inclusion = preprocess_sequences(
        inclusion_sequences, 
        max_sequences=max_seq_count,
        max_length=max_seq_length
    )
    
    # If pre-loaded exclusion sequences are provided, use them
    all_exclusion_sequences = []
    if exclusion_sequences is not None:
        logger.info(f"Using {len(exclusion_sequences)} pre-loaded exclusion sequences")
        processed_exclusion = preprocess_sequences(
            exclusion_sequences,
            max_sequences=max_seq_count,
            max_length=max_seq_length
        )
        all_exclusion_sequences.extend(processed_exclusion)
        
    # Otherwise, if not exclusion sequences provided, we need to progressively build them
    elif not exclusion_taxids:            
        # Step 3a: Start with genus-level siblings
        logger.info("Finding genus-level siblings for progressive filtering")
        genus_siblings = get_related_taxa(
            taxid=inclusion_taxid,
            email=email,
            relationship="sibling",
            rank="genus",
            max_results=3,
            api_key=api_key
        )
        
        genus_taxids = [taxon.get("taxid") for taxon in genus_siblings 
                       if taxon.get("taxid")]
        
        if genus_taxids:
            logger.info(f"Using genus siblings as first exclusion set: {genus_taxids}")
            
            # Get sequences for these taxa
            genus_exclusion_seqs = []
            for taxid in genus_taxids:
                try:
                    seqs = fetch_sequences_for_taxid(
                        taxid=taxid,
                        email=email,
                        max_records=max_seq_count // 2,
                        api_key=api_key
                    )
                    genus_exclusion_seqs.extend(seqs)
                    logger.info(f"Retrieved {len(seqs)} sequences for exclusion taxid {taxid}")
                except Exception as e:
                    logger.warning(f"Error fetching sequences for {taxid}: {str(e)}")
            
            # Preprocess sequences
            processed_genus_exclusion = preprocess_sequences(
                genus_exclusion_seqs,
                max_sequences=max_seq_count,
                max_length=max_seq_length
            )
            
            all_exclusion_sequences.extend(processed_genus_exclusion)
    else:
        # Step 3b: Use user-provided exclusion taxids
        logger.info(f"Using user-provided exclusion taxids: {exclusion_taxids}")
        
        for taxid in exclusion_taxids:
            if gene_name:
                sequences = fetch_gene_sequences(
                    taxid=taxid,
                    gene_name=gene_name,
                    email=email,
                    max_records=max_seq_count // 2,
                    allow_wgs=True,
                    api_key=api_key
                )
            else:
                sequences = fetch_sequences_for_taxid(
                    taxid=taxid,
                    email=email,
                    max_records=max_seq_count // 2,
                    api_key=api_key
                )
                
                # Preprocess these sequences
                processed_seqs = preprocess_sequences(
                    seqs,
                    max_sequences=max_seq_count // len(exclusion_taxids),
                    max_length=max_seq_length
                )
                
                all_exclusion_sequences.extend(processed_seqs)
                
    # If binding_sites_only is True, use the discrete binding sites approach
    if binding_sites_only and all_exclusion_sequences:
        logger.info("Using discrete binding sites approach for assay design (default)")
        marker_result = find_discrete_binding_sites(
            processed_inclusion,
            all_exclusion_sequences,
            min_amplicon_length=50,
            max_amplicon_length=max_amplicon_length,
            min_conservation=min_conservation,
            min_specificity=min_specificity,
            std_kmer=std_kmer,
            parallel_processes=parallel_processes
        )
        
        if "error" not in marker_result:
            logger.info("Found suitable binding sites for assay design")
            if gene_selection_info:
                marker_result['gene_selection_info'] = gene_selection_info
            if gene_name and auto_gene_selection:
                marker_result['selected_gene'] = gene_name
            return marker_result
        else:
            logger.info("Could not find suitable binding sites, falling back to standard approaches")
    
    # Step 4: Filter conserved regions against exclusion set using LSH
    if all_exclusion_sequences:
        if use_lsh:
            logger.info(f"Using LSH to filter against {len(all_exclusion_sequences)} exclusion sequences")
        
            # Use LSH for faster marker identification
            marker_result = lsh_find_markers(
                processed_inclusion,
                all_exclusion_sequences,
                min_region_size=50,
                max_region_size=200,
                min_conservation=min_conservation,
                timeout_seconds=min(timeout_seconds, max(30, timeout_seconds - (time.time() - start_time)))
            )
        
            # If LSH found good markers, return them
            if "error" not in marker_result:
                logger.info(f"LSH found marker in {time.time() - start_time:.2f} seconds")
                if gene_selection_info:
                    marker_result['gene_selection_info'] = gene_selection_info
                if gene_name and auto_gene_selection:
                    marker_result['selected_gene'] = gene_name
                return marker_result
            else:
                logger.info("LSH search did not find suitable markers, falling back to traditional methods")
                
        # Either LSH is disabled or it failed to find good markers
        logger.info(f"Using traditional methods to filter against {len(all_exclusion_sequences)} exclusion sequences")
        
        # Second pass: Find multiple specific regions
        marker_result = find_specific_oligonucleotide_sites(
            processed_inclusion,
            all_exclusion_sequences,
            max_amplicon_length=max_amplicon_length,
            min_oligo_length=18,
            max_oligo_length=30,
            min_amplicon_length=50,
            min_conservation=min_conservation
        )

        if gene_selection_info:
            marker_result['gene_selection_info'] = gene_selection_info
        if gene_name and auto_gene_selection:
            marker_result['selected_gene'] = gene_name

        return marker_result
                                        
    # For case 2, try a different approach: family-level siblings
    if not exclusion_taxids and time.time() - start_time < timeout_seconds:
        logger.info("Trying family-level siblings for additional specificity")

        family_siblings = get_related_taxa(
            taxid=inclusion_taxid,
            email=email,
            relationship="sibling",
            rank="family",
            max_results=3,
            api_key=api_key
        )
        
        family_taxids = [taxon.get("taxid") for taxon in family_siblings 
                        if taxon.get("taxid")]
        
        if family_taxids:
            logger.info(f"Using family siblings as additional exclusion set: {family_taxids}")
            
            # Get sequences for these taxa
            family_exclusion_seqs = []
            for taxid in family_taxids:
                try:
                    seqs = fetch_sequences_for_taxid(
                        taxid=taxid,
                        email=email,
                        max_records=max_seq_count // 2,
                        api_key=api_key
                    )
                    family_exclusion_seqs.extend(seqs)
                    logger.info(f"Retrieved {len(seqs)} sequences for family exclusion taxid {taxid}")
                except Exception as e:
                    logger.warning(f"Error fetching sequences for {taxid}: {str(e)}")
            
            # Preprocess sequences
            processed_family_exclusion = preprocess_sequences(
                family_exclusion_seqs,
                max_sequences=max_seq_count,
                max_length=max_seq_length
            )
            
            # Start fresh with original inclusion sequences
            marker_result = find_discriminative_kmers(
                inclusion_sequences=processed_inclusion,
                exclusion_sequences=processed_family_exclusion,
                kmer_size=30,  # Try larger k-mer for more specificity
                min_conservation=min_conservation * 0.9,  # Slightly lower threshold
                max_sequences=max_seq_count,
                max_sequence_length=max_seq_length
            )
            
            if "error" not in marker_result and marker_result.get("marker_sequence"):
                logger.info(f"Found specific marker using family-level exclusion")
                if gene_selection_info:
                    marker_result['gene_selection_info'] = gene_selection_info
                if gene_name and auto_gene_selection:
                    marker_result['selected_gene'] = gene_name
                return marker_result
    
    # If all else fails, return the original conserved marker
    logger.info("Returning best available marker (may not be fully specific)")
    if processed_inclusion and len(processed_inclusion) > 0:
        sequence = str(processed_inclusion[0].seq)
        # Take a reasonable portion of the sequence (first 150bp or whole sequence if shorter)
        marker_length = min(150, len(sequence))
        conserved_marker = {
            "marker_sequence": sequence[:marker_length],
            "marker_length": marker_length,
            "conservation_score": 1.0, # Only one sequence used
            "specificity_score": 0.0, # Unknown specificity
            "description": "Fallback conserved marker (it's likely terrible, do your own validation please)"
        }
    else:
        # If there are no inclusion sequences, return an error
        conserved_marker = {
            "error": "No suitable marker found and no sequences available for fallback"
        }

    # Add gene selection info to fallback result
    if gene_selection_info:
        conserved_marker['gene_selection_info'] = gene_selection_info
    if gene_name and auto_gene_selection:
        conserved_marker['selected_gene'] = gene_name

    logger.info("Returning best available marker (may not be fully specific, please validate)")
    return conserved_marker

def process_with_sliding_windows(
    sequences: List[SeqRecord],
    window_size: int = 50000,
    step_size: int = 25000,
    max_windows: int = 5
) -> List[SeqRecord]:
    """
    Process sequences using a sliding window approach.
    
    Args:
        sequences: List of sequences to process
        window_size: Size of each window
        step_size: Step size between windows
        max_windows: Maximum windows per sequence
        
    Returns:
        List of windowed sequence segments
    """
    processed_sequences = []
    
    for seq in sequences:
        if len(seq.seq) <= window_size:
            # Keep sequences shorter than window_size as is
            processed_sequences.append(seq)
            continue
            
        windows_created = 0
        seq_len = len(seq.seq)
        
        # Create windows
        for start in range(0, seq_len - window_size + 1, step_size):
            if windows_created >= max_windows:
                break
                
            end = start + window_size
            
            # Create a new sequence record for this window
            from Bio.SeqRecord import SeqRecord
            from Bio.Seq import Seq
            
            window = SeqRecord(
                Seq(str(seq.seq[start:end])),
                id=f"{seq.id}_window_{start}_{end}",
                description=f"Window {start}-{end} of {seq.description}"
            )
            
            processed_sequences.append(window)
            windows_created += 1
    
    return processed_sequences
    
def find_specific_oligonucleotide_sites(
    inclusion_sequences: List[SeqRecord],
    exclusion_sequences: List[SeqRecord],
    max_amplicon_length: int = 200,
    min_oligo_length: int = 18,
    max_oligo_length: int = 30,
    min_amplicon_length: int = 50,
    min_conservation: float = 0.8
) -> Dict[str, Any]:
    """
    Find specific regions suitable for primers and probe.
    Instead of a single conserved region, identify multiple shorter regions.
    
    Args:
        inclusion_sequences: Sequences from target organisms
        exclusion_sequences: Sequences from non-target organisms
        min_oligo_length: Minimum length for oligos
        max_oligo_length: Maximum length for oligos
        min_amplicon_length: Minimum amplicon length
        max_amplicon_length: Maximum amplicon length
        min_conservation: Minimum conservation within inclusion sequences
        
    Returns:
        Dict with forward, probe, and reverse regions
    """
    # 1. Look for k-mers of appropriate oligo sizes
    conserved_regions = []
    
    # Try different k-mer sizes suitable for oligos
    for kmer_size in [25, 22, 20, 18]:
        # Find conserved k-mers within inclusion sequences
        conserved_kmers = kmer_analysis(
            inclusion_sequences, 
            exclusion_sequences,
            min_conservation=min_conservation,
            kmer_size=kmer_size,
            n_processes=parallel_processes
        ) if std_kmer else find_conserved_kmers(
            inclusion_sequences,
            min_conservation=min_conservation,
            kmer_size=kmer_size)
        
        # Map these k-mers to their positions in the sequences
        for kmer, conservation in conserved_kmers.items():
            positions = find_kmer_positions(kmer, inclusion_sequences)
            
            # Calculate specificity against exclusion sequences
            specificity = calculate_kmer_specificity(kmer, exclusion_sequences)
            
            conserved_regions.append({
                "kmer": kmer,
                "length": len(kmer),
                "conservation": conservation,
                "specificity": specificity,
                "positions": positions
            })
    
    # Sort regions by specificity
    conserved_regions.sort(key=lambda x: x["specificity"], reverse=True)
    
    # 2. Try to build an assay with three specific regions
    assay_candidates = design_assay_from_regions(
        conserved_regions,
        min_amplicon_length,
        max_amplicon_length
    )
    
    if assay_candidates:
        best_assay = assay_candidates[0]
        
        # Extract the full amplicon sequence
        template_seq = inclusion_sequences[0].seq
        amplicon_start = best_assay["forward"]["positions"][0][0]
        amplicon_end = best_assay["reverse"]["positions"][0][1]
        amplicon_seq = str(template_seq[amplicon_start:amplicon_end])
        
        return {
            "marker_sequence": amplicon_seq,
            "marker_length": len(amplicon_seq),
            "forward_primer_region": best_assay["forward"]["kmer"],
            "probe_region": best_assay["probe"]["kmer"],
            "reverse_primer_region": best_assay["reverse"]["kmer"],
            "conservation_score": min(
                best_assay["forward"]["conservation"],
                best_assay["probe"]["conservation"],
                best_assay["reverse"]["conservation"]
            ),
            "specificity_score": min(
                best_assay["forward"]["specificity"],
                best_assay["probe"]["specificity"],
                best_assay["reverse"]["specificity"]
            ),
            "description": "Amplicon with specific primer and probe binding sites"
        }
    else:
        # Fallback: Try to find the longest possible specific k-mer
        for kmer_size in range(50, min_oligo_length-1, -5):
            conserved_kmers = kmer_analysis(
                inclusion_sequences,
                exclusion_sequences, 
                min_conservation=min_conservation * 0.9,  # Lower threshold
                kmer_size=kmer_size
            )
            
            if conserved_kmers:
                best_kmer = max(conserved_kmers.items(), key=lambda x: x[1])[0]
                specificity = calculate_kmer_specificity(best_kmer, exclusion_sequences)
                
                return {
                    "marker_sequence": best_kmer,
                    "marker_length": len(best_kmer),
                    "conservation_score": conserved_kmers[best_kmer],
                    "specificity_score": specificity,
                    "description": "Single conserved specific region"
                }
    
    # If all else fails
    return {
        "error": "Could not identify suitable specific regions for assay design"
    }
    
def find_conserved_kmers(
    sequences: List[SeqRecord],
    min_conservation: float = 0.8,
    max_kmers: int = 1000,
    kmer_size: int = 20,
    use_optimized: bool = True
) -> Dict[str, float]:
    """
    Find k-mers conserved across sequences.
    
    Args:
        sequences: List of sequences to analyze
        min_conservation: Minimum conservation level (0-1)
        kmer_size: Size of k-mers to extract
        use_optimized: Whether to use the optimized implementation (do it)
        
    Returns:
        Dict mapping conserved k-mers to their conservation scores

    """
    if use_optimized:
        # Use optimized k-mer analysis
        return find_conserved_kmers_parallel(sequences, kmer_size, min_conservation)
        
    # Original k-mer analysis as a fallback
    if not sequences:
        return {}
        
    # Extract k-mers and count occurrences
    kmer_counts = {}
    
    for seq in sequences:
        seq_str = str(seq.seq).upper()
        # Use a set to only count each k-mer once per sequence
        seq_kmers = set()
        
        # Skip positions to improve speed (this checks every 3 positions so you're not sacrificing too much coverage)
        for i in range(0, len(seq_str) - kmer_size + 1,3): # Skip every 3 positions
            kmer = seq_str[i:i+kmer_size]
            if 'N' not in kmer and '-' not in kmer:
                seq_kmers.add(kmer)
                
        # Add each unique k-mer to the count
        for kmer in seq_kmers:
            kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
            
            # Early termination if we have too many k-mers
            if len(kmer_counts) > max_kmers:
                break
                
        # Calculation conservation scores and filter
        conserved_kmers = {}
        threshold = min_conservation * len(sequences)
        
        for kmer, count in kmer_counts.items():
            if count >= threshold:
                conserved_kmers[kmer] = count / len(sequences)
                
        return conserved_kmers
    
    
def find_kmer_positions(
    kmer: str,
    sequences: List[SeqRecord]
) -> List[tuple[int, int]]:
    """
    Find positions of a k-mer in sequences.
    
    Args:
        kmer: K-mer to locate
        sequences: Sequences to search in
        
    Returns:
        List of (start, end) positions in reference sequence
    """
    positions = []
    
    # Check if sequences is empty
    if not sequences:
        return positions
        
    # Use the first sequence as reference
    reference_sequence = sequences[0]
    
    # Check if reference_sequence is a SeqRecord or string
    if hasattr(reference_sequence, 'seq'):
        seq_str = str(reference_sequence.seq).upper()
    else:
        # It's already a string
        seq_str = str(reference_sequence).upper()
        
    # Find all occurrence
    pos = seq_str.find(kmer)
    while pos >= 0:
        end_pos = pos + len(kmer)
        positions.append((pos, end_pos))
        pos = seq_str.find(kmer, pos + 1)
        
    return positions
    
def calculate_kmer_specificity(
    kmer: str,
    exclusion_sequences: List[SeqRecord],
    max_mismatches: int = 3
) -> float:
    """
    Calculate specificity of a k-mer against exclusion sequences.
    
    Args:
        kmer: K-mer to evaluate
        exclusion_sequences: Sequences to check against
        max_mismatches: Maximum allowed mismatches
        
    Returns:
        Specificity score (0-1, higher is more specific)
    """
    if not exclusion_sequences:
        return 1.0  # Fully specific if no exclusion sequences
        
    else:
        for seq in exclusion_sequences:
            if not hasattr(seq, 'seq'):
                logger.warning(f"Found string instead of SeqRecord in exclusion sequences")
            
    
    # Count sequences with matches (allowing mismatches)
    match_count = 0
    
    for seq in exclusion_sequences:
        seq_str = str(seq.seq).upper()
        
        # Check for approximate matches
        for i in range(len(seq_str) - len(kmer) + 1):
            mismatches = 0
            for j in range(len(kmer)):
                if i+j < len(seq_str) and kmer[j] != seq_str[i+j]:
                    mismatches += 1
                    if mismatches > max_mismatches:
                        break
            
            if mismatches <= max_mismatches:
                match_count += 1
                break  # Found a match in this sequence, move to next
    
    # Calculate specificity as 1 - (match proportion)
    return 1.0 - (match_count / len(exclusion_sequences))

def design_assay_from_regions(
    regions: List[Dict],
    min_amplicon_length: int = 50,
    max_amplicon_length: int = 200
) -> List[Dict]:
    """
    Design assays using conserved specific regions.
    
    Args:
        regions: List of conserved and specific regions
        min_amplicon_length: Minimum amplicon length
        max_amplicon_length: Maximum amplicon length
        
    Returns:
        List of assay designs with forward, probe, and reverse regions
    """
    assay_candidates = []
    
    # Need at least 3 regions with positions
    valid_regions = [r for r in regions if r["positions"]]
    if len(valid_regions) < 3:
        return []
    
    # Try to find forward, probe, and reverse combinations
    for i, fwd in enumerate(valid_regions):
        for j, rev in enumerate(valid_regions):
            if i == j:
                continue
                
            # Check if these could form valid forward/reverse primers
            for fwd_pos in fwd["positions"]:
                for rev_pos in rev["positions"]:
                    # Check orientation and amplicon length
                    if fwd_pos[0] >= rev_pos[1]:  # Forward must be before reverse
                        continue
                        
                    amplicon_length = rev_pos[1] - fwd_pos[0]
                    if amplicon_length < min_amplicon_length or amplicon_length > max_amplicon_length:
                        continue
                    
                    # Look for a probe between them
                    for k, probe in enumerate(valid_regions):
                        if k == i or k == j:
                            continue
                            
                        for probe_pos in probe["positions"]:
                            # Check if probe is between primers
                            if probe_pos[0] > fwd_pos[1] and probe_pos[1] < rev_pos[0]:
                                # Found a valid assay design
                                assay = {
                                    "forward": fwd,
                                    "probe": probe,
                                    "reverse": rev,
                                    "amplicon_length": amplicon_length,
                                    "score": (fwd["specificity"] + probe["specificity"] + rev["specificity"]) / 3
                                }
                                assay_candidates.append(assay)
    
    # Sort by score
    assay_candidates.sort(key=lambda x: x["score"], reverse=True)
    
    return assay_candidates
    
def design_multi_region_assay(
    marker_info: Dict[str, Any],
    inclusion_sequences: List[SeqRecord],
    exclusion_sequences: List[SeqRecord],
    min_oligo_length: int = 18,
    max_oligo_length: int = 30,
    min_amplicon_length: int = 50,
    max_amplicon_length: int = 200
) -> Dict[str, Any]:
    """
    Design an assay using multiple specific regions rather than a single conserved marker.
    
    Args:
        marker_info (Dict[str, Any]): Information about marker region (may contain partial results)
        inclusion_sequences (List[SeqRecord]): Target sequences
        exclusion_sequences (List[SeqRecord]): Non-target sequences
        min_oligo_length (int): Minimum oligo length
        max_oligo_length (int): Maximum oligo length
        min_amplicon_length (int): Minimum amplicon length
        max_amplicon_length (int): Maximum amplicon length
        
    Returns:
        Dict[str, Any]: Assay information including primers, probe, and amplicon
    """
    # First, check if marker_info contains a usable sequence
    marker_sequence = marker_info.get("marker_sequence", "")
    
    if marker_sequence and len(marker_sequence) >= min_amplicon_length:
        # Try to design an assay from the existing marker
        assay_info = design_primers_and_probe(marker_info)
        if assay_info.get("status") == "success":
            return assay_info
        logger.info("Existing marker not suitable, searching for specific oligo regions")
    
    # We need to find specific regions for primers and probe
    logger.info("Searching for specific regions to use as primer and probe binding sites")
    
    # Step 1: Find k-mers that are specific to inclusion sequences
    forward_candidates = []
    probe_candidates = []
    reverse_candidates = []
    
    # Try different k-mer sizes
    for kmer_size in [25, 22, 20, 18]:
        if len(forward_candidates) >= 5 and len(probe_candidates) >= 5 and len(reverse_candidates) >= 5:
            break  # We have enough candidates
        
        # Find conserved k-mers within inclusion sequences
        conserved_kmers = find_conserved_kmers(inclusion_sequences, kmer_size=kmer_size)
        
        # Check each k-mer for specificity against exclusion sequences
        for kmer, conservation in conserved_kmers.items():
            # Calculate specificity score
            specificity = calculate_specificity(kmer, exclusion_sequences)
            
            # Calculate basic thermodynamic properties
            properties = calculate_oligo_properties(kmer)
            
            # Score this k-mer as an oligo
            oligo_score = calculate_oligo_score(properties)
            combined_score = (specificity * 0.6) + (oligo_score * 0.4)
            
            # Find positions in the reference sequence
            positions = find_positions(kmer, inclusion_sequences[0])
            
            if not positions:
                continue
                
            # Categorize based on position (beginning, middle, end)
            ref_seq_length = len(inclusion_sequences[0].seq)
            for pos in positions:
                relative_pos = pos / ref_seq_length
                
                if relative_pos < 0.3:  # Beginning of sequence
                    forward_candidates.append({
                        "sequence": kmer,
                        "length": len(kmer),
                        "position": pos,
                        "specificity": specificity,
                        "conservation": conservation,
                        "properties": properties,
                        "score": combined_score
                    })
                elif relative_pos > 0.7:  # End of sequence
                    # For reverse primer, we need the reverse complement
                    rev_comp = str(Seq(kmer).reverse_complement())
                    properties = calculate_oligo_properties(rev_comp)
                    
                    reverse_candidates.append({
                        "sequence": rev_comp,
                        "length": len(rev_comp),
                        "position": pos,
                        "specificity": specificity,
                        "conservation": conservation,
                        "properties": properties,
                        "score": combined_score
                    })
                else:  # Middle of sequence - potential probe
                    probe_candidates.append({
                        "sequence": kmer,
                        "length": len(kmer),
                        "position": pos,
                        "specificity": specificity,
                        "conservation": conservation,
                        "properties": properties,
                        "score": combined_score
                    })
    
    # Sort candidates by score
    forward_candidates.sort(key=lambda x: x["score"], reverse=True)
    probe_candidates.sort(key=lambda x: x["score"], reverse=True)
    reverse_candidates.sort(key=lambda x: x["score"], reverse=True)
    
    logger.info(f"Found {len(forward_candidates)} forward, {len(probe_candidates)} probe, and {len(reverse_candidates)} reverse candidates")
    
    # Step 2: Try to build an assay with appropriate spacing
    assay_designs = design_assays(
        forward_candidates,
        probe_candidates,
        reverse_candidates,
        min_amplicon_length,
        max_amplicon_length
    )
    
    if not assay_designs:
        logger.warning("Could not design a suitable assay with specific regions")
        return {
            "status": "failed",
            "message": "Could not find specific regions suitable for assay design",
            "primers": [],
            "probe": None,
            "amplicon": "",
            "amplicon_length": 0
        }
    
    # Choose the best assay design
    best_assay = assay_designs[0]
    
    # Step 3: Create the final assay components
    forward_primer = {
        "name": "Forward_Primer",
        "sequence": best_assay["forward"]["sequence"],
        "properties": best_assay["forward"]["properties"],
        "position": best_assay["forward"]["position"],
        "orientation": "forward",
        "specificity": best_assay["forward"]["specificity"]
    }
    
    reverse_primer = {
        "name": "Reverse_Primer",
        "sequence": best_assay["reverse"]["sequence"],
        "properties": best_assay["reverse"]["properties"],
        "position": best_assay["reverse"]["position"],
        "orientation": "reverse",
        "specificity": best_assay["reverse"]["specificity"]
    }
    
    probe = {
        "name": "Probe",
        "sequence": best_assay["probe"]["sequence"],
        "properties": best_assay["probe"]["properties"],
        "position": best_assay["probe"]["position"],
        "orientation": "forward",
        "specificity": best_assay["probe"]["specificity"]
    }
    
    # Extract the amplicon sequence from the reference sequence
    amplicon_start = best_assay["forward"]["position"]
    amplicon_end = best_assay["reverse"]["position"] + len(best_assay["reverse"]["sequence"])
    amplicon_sequence = str(inclusion_sequences[0].seq[amplicon_start:amplicon_end])
    
    return {
        "status": "success",
        "message": "Successfully designed assay using specific oligo regions",
        "primers": [forward_primer, reverse_primer],
        "probe": probe,
        "amplicon": amplicon_sequence,
        "amplicon_length": len(amplicon_sequence),
        "conservation_score": min(
            best_assay["forward"]["conservation"],
            best_assay["probe"]["conservation"],
            best_assay["reverse"]["conservation"]
        ),
        "specificity_score": min(
            best_assay["forward"]["specificity"],
            best_assay["probe"]["specificity"],
            best_assay["reverse"]["specificity"]
        )
    }

def find_conserved_kmers(
    sequences: List[SeqRecord],
    min_conservation: float = 0.8,
    kmer_size: int = 20
) -> Dict[str, float]:
    """
    Find k-mers conserved across sequences.
    
    Args:
        sequences: List of sequences to analyze
        min_conservation: Minimum conservation level (0-1)
        kmer_size: Size of k-mers to extract
        
    Returns:
        Dict mapping conserved k-mers to their conservation scores
    """
    if not sequences:
        return {}
    
    # Extract k-mers and count occurrences
    kmer_counts = {}
    
    for seq in sequences:
        seq_str = str(seq.seq).upper()
        # Use a set to only count each k-mer once per sequence
        seq_kmers = set()
        
        for i in range(len(seq_str) - kmer_size + 1):
            kmer = seq_str[i:i+kmer_size]
            if 'N' not in kmer and '-' not in kmer:
                seq_kmers.add(kmer)
        
        # Add each unique k-mer to the count
        for kmer in seq_kmers:
            kmer_counts[kmer] = kmer_counts.get(kmer, 0) + 1
    
    # Calculate conservation scores and filter
    conserved_kmers = {}
    threshold = min_conservation * len(sequences)
    
    for kmer, count in kmer_counts.items():
        if count >= threshold:
            conserved_kmers[kmer] = count / len(sequences)
    
    return conserved_kmers

def calculate_specificity(
    kmer: str,
    exclusion_sequences: List[SeqRecord],
    max_mismatches: int = 2
) -> float:
    """
    Calculate specificity of a k-mer against exclusion sequences.
    
    Args:
        kmer: K-mer to evaluate
        exclusion_sequences: Sequences to check against
        max_mismatches: Maximum allowed mismatches
        
    Returns:
        Specificity score (0-1, higher is more specific)
    """
    if not exclusion_sequences:
        return 1.0  # Fully specific if no exclusion sequences
    
    # Count sequences with matches (allowing mismatches)
    match_count = 0
    
    for seq in exclusion_sequences:
        seq_str = str(seq.seq).upper()
        
        # Check for approximate matches
        found_match = False
        for i in range(len(seq_str) - len(kmer) + 1):
            mismatches = 0
            for j in range(len(kmer)):
                if i+j < len(seq_str) and kmer[j] != seq_str[i+j]:
                    mismatches += 1
                    if mismatches > max_mismatches:
                        break
            
            if mismatches <= max_mismatches:
                found_match = True
                break
        
        if found_match:
            match_count += 1
    
    # Calculate specificity as 1 - (match proportion)
    return 1.0 - (match_count / len(exclusion_sequences))

def find_positions(kmer: str, reference_sequence: SeqRecord) -> List[int]:
    """
    Find positions of a k-mer in a reference sequence.
    
    Args:
        kmer: K-mer to locate
        reference_sequence: Reference sequence
        
    Returns:
        List of start positions
    """
    positions = []
    seq_str = str(reference_sequence.seq).upper()
    
    # Find all occurrences
    start = 0
    while True:
        start = seq_str.find(kmer, start)
        if start == -1:
            break
        positions.append(start)
        start += 1
    
    return positions

def design_assays(
    forward_candidates: List[Dict],
    probe_candidates: List[Dict],
    reverse_candidates: List[Dict],
    min_amplicon_length: int = 50,
    max_amplicon_length: int = 200
) -> List[Dict]:
    """
    Design assays by combining forward primer, probe, and reverse primer candidates.
    
    Args:
        forward_candidates: List of forward primer candidates
        probe_candidates: List of probe candidates
        reverse_candidates: List of reverse primer candidates
        min_amplicon_length: Minimum amplicon length
        max_amplicon_length: Maximum amplicon length
        
    Returns:
        List of assay designs
    """
    assay_designs = []
    
    # Limit to top candidates for efficiency
    forward_candidates = forward_candidates[:10]
    probe_candidates = probe_candidates[:15]
    reverse_candidates = reverse_candidates[:10]
    
    # Try all combinations
    for forward in forward_candidates:
        forward_pos = forward["position"]
        forward_len = forward["length"]
        
        for reverse in reverse_candidates:
            reverse_pos = reverse["position"]
            
            # Check amplicon length
            amplicon_length = reverse_pos - forward_pos + reverse["length"]
            if amplicon_length < min_amplicon_length or amplicon_length > max_amplicon_length:
                continue
            
            # Look for probes between primers
            valid_probes = []
            for probe in probe_candidates:
                probe_pos = probe["position"]
                
                # Check if probe is between primers
                if forward_pos + forward_len <= probe_pos and probe_pos + probe["length"] <= reverse_pos:
                    valid_probes.append(probe)
            
            if valid_probes:
                # Choose the best probe
                best_probe = max(valid_probes, key=lambda p: p["score"])
                
                # Calculate overall assay score
                assay_score = (
                    forward["score"] * 0.3 +
                    best_probe["score"] * 0.4 +  # Give higher weight to probe
                    reverse["score"] * 0.3
                )
                
                assay_designs.append({
                    "forward": forward,
                    "probe": best_probe,
                    "reverse": reverse,
                    "amplicon_length": amplicon_length,
                    "score": assay_score
                })
    
    # Sort designs by score
    assay_designs.sort(key=lambda x: x["score"], reverse=True)
    
    return assay_designs
    
def progressive_kmer_search(inclusion_seqs, exclusion_seqs, timeout_seconds=60):
    """Start with small regions and progressively expand only if needed."""
    start_time = time.time()
    
    # Try with small k-mers first (faster)
    for kmer_size in [20, 25, 30]:
        if time.time() - start_time > timeout_seconds:
            break
            
        result = find_discriminative_kmers(
            inclusion_seqs, 
            exclusion_seqs,
            kmer_size=kmer_size,
            max_sequence_length=10000  # Start with smaller sequences
        )
        
        # If we found good markers, return early
        if result.get("conservation_score", 0) > 0.8 and result.get("specificity_score", 0) > 0.7:
            return result
    
    # If we haven't found a good marker yet, try with a larger sequence window
    return find_discriminative_kmers(
        inclusion_seqs,
        exclusion_seqs,
        kmer_size=25,
        max_sequence_length=50000
    )
    
def sample_large_sequence(sequence, max_length=20000, samples=3):
    """Sample regions from a large sequence instead of just taking the beginning."""
    seq_length = len(sequence.seq)
    if seq_length <= max_length:
        return [sequence]
    
    samples_list = []
    # Take beginning, middle, and end segments
    segment_length = max_length // samples
    
    for i in range(samples):
        start = (seq_length * i) // samples
        end = start + segment_length
        
        segment = sequence[start:end]
        segment.id = f"{sequence.id}_segment_{start}-{end}"
        samples_list.append(segment)
    
    return samples_list
    
def find_discrete_binding_sites(
    inclusion_sequences: list[SeqRecord],
    exclusion_sequences: list[SeqRecord],
    primer_length_range: tuple[int, int] = (18, 25),
    probe_length_range: tuple[int, int] = (20, 30),
    min_amplicon_length: int = 75,
    max_amplicon_length: int = 200,
    min_conservation: float = 0.8,
    min_specificity: float = 0.8,
    std_kmer: bool = True,
    parallel_processes: Optional[int] = 8
) -> Dict[str, Any]:
    """
    Find discrete binding sites for primers and probe that meet specificity criteria
    even if the internal amplicon region is not specific.
    
    Args:
        inclusion_sequences: Target sequences
        exclusion_sequences: Non-target sequences
        primer_length_range: Min and max primer length
        probe_length_range: Min and max probe length
        min_amplicon_length: Minimum amplicon length
        max_amplicon_length: Maximum amplicon length
        min_conservation: Minimum conservation in inclusion sequences
        min_specificity: Minimum specificity compared to exclusion
        
    Returns:
        Dict with marker information and binding sites
    """
    logger.info("Searching for discrete specific binding sites for assay design")
    
    # 1. Find conserved k-mers in inclusion sequences
    min_kmer_size = min(primer_length_range[0], probe_length_range[0])
    conserved_kmers = {}
    
    for kmer_size in range(min_kmer_size, probe_length_range[1] + 1, 2):
        # Find conserved k-mers within inclusion sequences
        kmers = find_conserved_kmers(
            inclusion_sequences, 
            kmer_size=kmer_size, 
            min_conservation=min_conservation
        )
        conserved_kmers[kmer_size] = kmers
        logger.info(f"Found {len(kmers)} conserved {kmer_size}-mers")
        
        if len(kmers) > 50:
            break  # We have enough candidates
    
    # 2. Calculate specificity for each k-mer
    specific_sites = []
    
    for kmer_size, kmers in conserved_kmers.items():
        for kmer, conservation in kmers.items():
            # Calculate specificity against exclusion sequences
            specificity = calculate_kmer_specificity(kmer, exclusion_sequences)
            
            # Only keep k-mers that meet specificity threshold
            if specificity >= min_specificity:
                # Find positions in the reference sequence
                positions = find_kmer_positions(kmer, inclusion_sequences[0])
                
                if positions:
                    specific_sites.append({
                        "sequence": kmer,
                        "length": kmer_size,
                        "positions": positions,
                        "conservation": conservation,
                        "specificity": specificity
                    })
    
    if not specific_sites:
        logger.warning("No specific binding sites found")
        return {"error": "Could not find specific binding sites for assay design"}
    
    # 3. Sort sites by specificity
    specific_sites.sort(key=lambda x: x["specificity"], reverse=True)
    
    # 4. Find viable combinations for primers and probe
    viable_assays = []
    
    # Get a reference sequence for distance calculations
    ref_seq_length = len(inclusion_sequences[0].seq)
    
    # Try different combinations of sites as forward, probe, and reverse
    for fwd_idx, fwd_site in enumerate(specific_sites[:20]):  # Limit to top 20 for efficiency
        for rev_idx, rev_site in enumerate(specific_sites[:20]):
            # Skip if same site
            if fwd_idx == rev_idx:
                continue
                
            # Check each potential position pair
            for fwd_pos in fwd_site["positions"]:
                for rev_pos in rev_site["positions"]:
                    # Ensure correct orientation
                    if fwd_pos[0] >= rev_pos[0]:  # Forward must be before reverse
                        continue
                        
                    # Calculate amplicon length
                    amplicon_length = rev_pos[1] - fwd_pos[0]
                    if amplicon_length < min_amplicon_length or amplicon_length > max_amplicon_length:
                        continue
                    
                    # Find potential probes between primers
                    for probe_idx, probe_site in enumerate(specific_sites[:30]):
                        if probe_idx == fwd_idx or probe_idx == rev_idx:
                            continue
                            
                        for probe_pos in probe_site["positions"]:
                            # Check if probe is between primers
                            if probe_pos[0] > fwd_pos[1] and probe_pos[1] < rev_pos[0]:
                                # Calculate spacing
                                fwd_to_probe = probe_pos[0] - fwd_pos[1]
                                probe_to_rev = rev_pos[0] - probe_pos[1]
                                
                                # Ensure adequate spacing
                                if fwd_to_probe >= 5 and probe_to_rev >= 5:
                                    # Viable assay design found
                                    assay = {
                                        "forward_primer": fwd_site,
                                        "probe": probe_site,
                                        "reverse_primer": rev_site,
                                        "amplicon_length": amplicon_length,
                                        "forward_position": fwd_pos,
                                        "probe_position": probe_pos,
                                        "reverse_position": rev_pos,
                                        "score": (
                                            fwd_site["specificity"] +
                                            probe_site["specificity"] +
                                            rev_site["specificity"]
                                        ) / 3
                                    }
                                    viable_assays.append(assay)
    
    if not viable_assays:
        logger.warning("No viable assay designs found with specific binding sites")
        return {"error": "Could not find viable combination of specific binding sites"}
    
    # 5. Sort by score and return best assay
    viable_assays.sort(key=lambda x: x["score"], reverse=True)
    best_assay = viable_assays[0]
    
    # Extract full amplicon from reference sequence
    fwd_pos = best_assay["forward_position"][0]
    rev_pos = best_assay["reverse_position"][1]
    amplicon_seq = str(inclusion_sequences[0].seq[fwd_pos:rev_pos])
    
    return {
        "marker_sequence": amplicon_seq,
        "marker_length": len(amplicon_seq),
        "conservation_score": min(
            best_assay["forward_primer"]["conservation"],
            best_assay["probe"]["conservation"],
            best_assay["reverse_primer"]["conservation"]
        ),
        "specificity_score": best_assay["score"],
        "description": "Amplicon with specific primer and probe binding sites",
        "forward_primer_region": best_assay["forward_primer"]["sequence"],
        "probe_region": best_assay["probe"]["sequence"],
        "reverse_primer_region": best_assay["reverse_primer"]["sequence"]
    }