# assay_design/data_retrieval.py

import time
import logging
from typing import List, Optional, Dict, Any, Union
from Bio import Entrez, SeqIO
from Bio.SeqRecord import SeqRecord

from .cache_manager import SequenceCacheManager

# Initialize the cache manager (do this once at the module level)
cache_manager = SequenceCacheManager()

# Logger for the module
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

def fetch_sequence_by_accession(accession_id: str, email: str, db: str = "nucleotide") -> SeqIO.SeqRecord:
    """
    Fetch a single sequence by accession from NCBI.
    
    Args:
        accession_id (str): NCBI accession ID (e.g., 'NC_000913.3').
        email (str): User's email address (required by NCBI to track usage).
        db (str): NCBI database from which to fetch (default 'nucleotide').
        
    Returns:
        SeqIO.SeqRecord: A Biopython SeqRecord containing the requested sequence.
        
    Raises:
        ValueError: If no email is provided.
        RuntimeError: If retrieval fails or sequence is empty.
    """
    if not email:
        raise ValueError("You must provide a valid email address to comply with NCBI's usage policies.")
        
    # Set NCBI Entrez parameters
    Entrez.email = email
    
    try:
        logger.info(f"Fetching sequence for accession {accession_id} from NCBI {db} database.")
        handle = Entrez.efetch(db=db, id=accession_id, rettype="fasta", retmode="text")
        seq_record = SeqIO.read(handle, "fasta-pearson")
        handle.close()
        
        if not seq_record or len(seq_record.seq) == 0:
            raise RuntimeError(f"No sequence data returned for accession {accession_id}")
            
        logger.info(f"Successfully retrieved sequence: {seq_record.id}, length: {len(seq_record.seq)}")
        return seq_record
        
    except Exception as e:
        logger.error(f"Error fetching {accession_id} from NCBI: {str(e)}")
        raise

def fetch_sequences_for_taxid(
    taxid: str, 
    email: str, 
    query_term: Optional[str] = None, 
    max_records: int = 50,
    db: str = "nucleotide"
) -> List[SeqRecord]:        
    """
    Fetch sequences by NCBI TaxID with optional additional query terms.
    
    Args:
        taxid (str): NCBI TaxID, e.g., '562' for E. coli.
        email (str): User's email address (required by NCBI).
        db (str): NCBI database, default 'nucleotide'.
        max_records (int): Maximum number of sequences to fetch.
        query_term (Optional[str]): Optional specific NCBI query term.
            If None, will use "txid{taxid}[Organism]"
        
    Returns:
        List[SeqIO.SeqRecord]: A list of retrieved SeqRecords.
    """
    if not email:
        raise ValueError("Email address is required.")
        
    Entrez.email = email
    
    try:
        # Build the search query
        if query_term is None:
            query_term = f"txid{taxid}[Organism] AND biomol_genomic[PROP]"
            
        query_term += "NOT wgs[Filter]"
        
        logger.info(f"Searching for up to {max_records} sequences with query: {query_term}")
        search_handle = Entrez.esearch(db=db, term=query_term, retmax=max_records)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        ids = search_results.get("IdList", [])
        logger.info(f"Found {len(ids)} records. Fetching details...")
        
        records = []
        if ids:
            try:
                fetch_handle = Entrez.efetch(db=db, id=",".join(ids), rettype="gb", retmode="text")
                records = list(SeqIO.parse(fetch_handle, "genbank"))
                fetch_handle.close()
                
                # If we got no genbank records, try FASTA format as a fallback
                if not records:
                    fetch_handle = Entrez.efetch(db=db, id=",".join(ids), rettype="fasta", retmode="text")
                    records = list(SeqIO.parse(fetch_handle, "fasta"))
                    fetch_handle.close()
                    
            except Exception as inner_e:
                logger.warning(f"Error in first fetch attempt: {str(inner_e)}. Trying alternative format...")
                try:
                    fetch_handle = Entrez.efetch(db=db, id=",".join(ids), rettype="fasta", retmode="text")
                    records = list(SeqIO.parse(fetch_handle, "fasta"))
                    fetch_handle.close()
                except Exception as fallback_e:
                    logger.error(f"Both fetch attempts failed: {str(fallback_e)}")
                    
            logger.info(f"Fetched {len(records)} sequences.")
            
            # Validate sequences before returning
            records = ensure_sequences_are_defined(records, email)
            logger.info(f"Validated {len(records)} sequences with defined content.")
            
        return records
        
    except Exception as e:
        logger.error(f"Error fetching sequences: {str(e)}")
        raise
                
def fetch_gene_sequences(
    taxid: str,
    gene_name: str,
    email: str,
    max_records: int = 50,
    allow_wgs: bool = False
) -> List[SeqRecord]:
    """
    Fetch sequences for a specific gene from a specific taxon.
    
    Args:
        taxid (str): NCBI TaxID
        gene_name (str): Gene name to search for
        email (str): User's email address
        max_records (int): Maximum number of sequences to fetch
        
    Returns:
        List[SeqRecord]: List of gene sequences
    """
    # Define primary and fallback queries
    wgs_filter = "" if allow_wgs else "NOT wgs[Filter]"
    queries = [
        f"txid{taxid}[Organism] AND {gene_name}[Product] OR {gene_name}[Gene] OR {gene_name}[symbol] AND biomol_genomic[PROP] AND 200:20000[SLEN] {wgs_filter}",
        f"txid{taxid}[Organism] AND \"{gene_name}\"[Product] AND biomol_genomic[PROP] AND 200:20000[SLEN] {wgs_filter}",
        f"txid{taxid}[Organism] AND \"{gene_name}\" AND biomol_genomic[PROP] AND 200:20000[SLEN] {wgs_filter}",
        f"txid{taxid}[Organism] AND {gene_name} AND biomol_genomic[PROP] AND 200:20000[SLEN] {wgs_filter}"
    ]
    
    # Add alternative queries for common marker genes
    if "16S" in gene_name:
        queries.extend([
            f"txid{taxid}[Organism] AND 16S rRNA[Keyword] AND biomol_genomic[PROP] AND 200:20000[SLEN] {wgs_filter}",
            f"txid{taxid}[Organism] AND 16S ribosomal RNA[Keyword] AND biomol_genomic[PROP] AND 200:20000[SLEN] {wgs_filter}",
            f"txid{taxid}[Organism] AND 16S[Title] AND biomol_genomic[PROP] AND 200:20000[SLEN] {wgs_filter}",
            f"txid{taxid}[Organism] AND small subunit ribosomal RNA AND biomol_genomic[PROP] AND 200:20000[SLEN] {wgs_filter}"
        ])
    
    all_sequences = []
    
    for query in queries:
        try:
            logger.info(f"Trying query: {query}")
            sequences = fetch_sequences_for_taxid(
                taxid=taxid,
                email=email,
                query_term=query,
                max_records=max_records
            )
            
            if sequences:
                logger.info(f"Found {len(sequences)} sequences with query: {query}")
                all_sequences.extend(sequences)
                break  # Stop if we found sequences
                
        except Exception as e:
            logger.warning(f"Query failed: {str(e)}")
    
    # For 16S rRNA specifically, here are some alternative terms that might be used
    if not all_sequences and "16S" in gene_name:
        alternative_queries = [
            f"txid{taxid}[Organism] AND 16S rRNA[Keyword] AND biomol_genomic[PROP] AND 200:20000[SLEN] {wgs_filter}",
            f"txid{taxid}[Organism] AND 16S ribosomal RNA[Keyword] AND biomol_genomic[PROP] AND 200:20000[SLEN] {wgs_filter}",
            f"txid{taxid}[Organism] AND 16S[Title] AND biomol_genomic[PROP] AND 200:20000[SLEN] {wgs_filter}",
            f"txid{taxid}[Organism] AND small subunit ribosomal RNA AND biomol_genomic[PROP] AND 200:20000[SLEN] {wgs_filter}"
        ]
        
        for query in alternative_queries:
            try:
                logger.info(f"Trying alternative query: {query}")
                sequences = fetch_sequences_for_taxid(
                    taxid=taxid,
                    email=email,
                    query_term=query,
                    max_records=max_records
                ) 
                
                if sequences:
                    logger.info(f"Found {len(sequences)} sequences with alternative query: {query}")
                    all_sequences.extend(sequences)
                    break # Stop if we found sequences
                    
            except Exception as e:
                logger.warning(f"Alternative query failed: {str(e)}")
    
    return all_sequences    

def fetch_marker_genes(
    taxid: str,
    email: str,
    marker_gene: str = "16S ribosomal RNA",
    max_records: int = 50
) -> List[SeqRecord]:
    """
    Fetch common marker gene sequences for a specific taxon.
    
    Args:
        taxid (str): NCBI TaxID
        email (str): User's email address
        marker_gene (str): Marker gene to fetch (default: 16S rRNA)
        max_records (int): Maximum number of sequences to fetch
        
    Returns:
        List[SeqRecord]: List of marker gene sequences
    """
    # Create specific query for the marker gene
    query = f"txid{taxid}[Organism] AND {marker_gene}[Product]"
    return fetch_sequences_for_taxid(
        taxid=taxid,
        email=email,
        query_term=query,
        max_records=max_records
    )

def get_taxon_info(taxid: str, email: str) -> Dict[str, Any]:
    """
    Get information about a specific taxon.
    
    Args:
        taxid (str): NCBI TaxID
        email (str): User's email address
        
    Returns:
        Dict[str, Any]: Taxonomy information
    """
    if not email:
        raise ValueError("Email address is required.")
        
    Entrez.email = email
    
    try:
        logger.info(f"Fetching taxonomy information for taxid {taxid}")
        handle = Entrez.efetch(db="taxonomy", id=taxid, retmode="xml")
        records = Entrez.read(handle)
        handle.close()
        
        if not records:
            logger.error(f"No taxonomy information found for taxid {taxid}")
            return {}
        
        tax_info = records[0]
        
        # Extract relevant information
        result = {
            "taxid": taxid,
            "scientific_name": tax_info.get("ScientificName", ""),
            "rank": tax_info.get("Rank", ""),
            "division": tax_info.get("Division", ""),
            "lineage": tax_info.get("Lineage", ""),
            "genetic_code": tax_info.get("GeneticCode", {}).get("GCName", "")
        }
        
        # Extract lineage with ranks
        if "LineageEx" in tax_info:
            result["lineage_ex"] = [
                {
                    "taxid": entry.get("TaxId"),
                    "name": entry.get("ScientificName"),
                    "rank": entry.get("Rank")
                }
                for entry in tax_info["LineageEx"]
            ]
        
        return result
        
    except Exception as e:
        logger.error(f"Error fetching taxonomy info for {taxid}: {str(e)}")
        raise

def get_related_taxa(
    taxid: str, 
    email: str, 
    relationship: str = "sibling",
    rank: Optional[str] = None,
    max_results: int = 50
) -> List[Dict[str, Any]]:
    """
    Find related taxa for a given taxid.
    
    Args:
        taxid (str): NCBI TaxID
        email (str): User's email address
        relationship (str): Type of relationship ("sibling", "child", "parent")
        rank (Optional[str]): Taxonomic rank to consider
        max_results (int): Maximum number of results to return
        
    Returns:
        List[Dict[str, Any]]: List of related taxa info
    """
    if not email:
        raise ValueError("Email address is required.")
        
    Entrez.email = email
    
    try:
        # First get the taxonomy information for our target
        tax_info = get_taxon_info(taxid, email)
        if not tax_info:
            return []
        
        related_taxa = []
        
        if relationship == "sibling":
            # Find siblings by first getting the parent
            lineage = tax_info.get("lineage_ex", [])
            
            # If rank specified, find the parent at that rank
            parent_taxid = None
            if rank:
                for entry in lineage:
                    if entry.get("rank") == rank:
                        parent_taxid = entry.get("taxid")
                        break
            
            # If no parent at specified rank, use direct parent (usually genus)
            if not parent_taxid and lineage:
                # The direct parent is the last entry in the lineage
                parent_taxid = lineage[-1].get("taxid")
            
            # If we found a parent, get its children
            if parent_taxid:
                # Build a query to find all species in the same genus/family
                parent_info = get_taxon_info(parent_taxid, email)
                parent_name = parent_info.get("scientific_name", "")
                parent_rank = parent_info.get("rank", "")
                
                if parent_name:
                    # Search for all taxa under this parent
                    query = f"{parent_name}[{parent_rank}]"
                    search_handle = Entrez.esearch(db="taxonomy", term=query, retmax=max_results + 1)
                    search_results = Entrez.read(search_handle)
                    search_handle.close()
                    
                    child_ids = search_results.get("IdList", [])
                    
                    # Filter out the original taxid
                    child_ids = [tid for tid in child_ids if tid != taxid]
                    
                    # Get information for each child
                    for child_id in child_ids[:max_results]:
                        child_info = get_taxon_info(child_id, email)
                        if child_info:
                            related_taxa.append(child_info)
        
        return related_taxa
        
    except Exception as e:
        logger.error(f"Error finding related taxa for {taxid}: {str(e)}")
        return []

def get_related_taxa_weighted(
    taxid: str,
    email: str,
    max_results: int = 10,
    diversity_levels: Optional[List[str]] = None,
    min_sequence_count: int = 5
) -> List[Dict[str, Any]]:
    """
    Get phylogenetically weighted exclusion taxa with sequence availability.

    This enhanced version prioritizes taxa by:
    1. Phylogenetic proximity (closer relatives ranked higher)
    2. Sequence availability (taxa with more sequences preferred)
    3. Taxonomic diversity (spans multiple levels: genus, family, order)

    Args:
        taxid (str): NCBI TaxID for inclusion taxon
        email (str): User's email address
        max_results (int): Maximum number of related taxa to return
        diversity_levels (Optional[List[str]]): Taxonomic ranks to sample
            Default: ["genus", "family", "order"]
        min_sequence_count (int): Minimum sequences required (default: 5)

    Returns:
        List[Dict[str, Any]]: List of weighted taxa info with fields:
            - taxid: str
            - scientific_name: str
            - rank: str
            - phylogenetic_distance: int (lower = closer, in taxonomic steps)
            - sequence_count: int (estimated)
            - priority_score: float (higher = better candidate)
            - lineage: str
    """
    if not email:
        raise ValueError("Email address is required.")

    if diversity_levels is None:
        diversity_levels = ["genus", "family", "order"]

    Entrez.email = email

    try:
        # Get the taxonomy information for our target
        logger.info(f"Fetching weighted related taxa for taxid {taxid}")
        tax_info = get_taxon_info(taxid, email)

        if not tax_info:
            logger.error(f"Could not get taxonomy info for taxid {taxid}")
            return []

        inclusion_rank = tax_info.get("rank", "")
        lineage = tax_info.get("lineage_ex", [])

        logger.info(f"Inclusion taxon: {tax_info.get('scientific_name', 'Unknown')} (rank: {inclusion_rank})")

        # Collect candidate taxa from multiple taxonomic levels
        all_candidates = []

        for level in diversity_levels:
            logger.info(f"Searching for {level}-level relatives")

            # Get siblings at this taxonomic level
            level_siblings = get_related_taxa(
                taxid=taxid,
                email=email,
                relationship="sibling",
                rank=level,
                max_results=50  # Get more initially, will filter later
            )

            if not level_siblings:
                logger.info(f"No {level}-level siblings found")
                continue

            # Calculate phylogenetic distance for each candidate
            for sibling in level_siblings:
                sibling_taxid = sibling.get("taxid")
                sibling_name = sibling.get("scientific_name", "Unknown")
                sibling_rank = sibling.get("rank", "")

                # Calculate phylogenetic distance (number of taxonomic steps from inclusion)
                # Distance is determined by the level at which we found the sibling
                distance = _calculate_taxonomic_distance(level, inclusion_rank)

                # Estimate sequence availability
                try:
                    sequence_count = _estimate_sequence_count(sibling_taxid, email)
                except Exception as e:
                    logger.warning(f"Could not estimate sequence count for {sibling_name}: {e}")
                    sequence_count = 0

                # Skip taxa with insufficient sequences
                if sequence_count < min_sequence_count:
                    logger.debug(f"Skipping {sibling_name} (only {sequence_count} sequences)")
                    continue

                # Calculate priority score
                # Formula: priority = (1/distance) * log(seq_count + 1) * rank_weight
                import math
                rank_weights = {"genus": 1.0, "family": 0.8, "order": 0.6, "class": 0.4}
                rank_weight = rank_weights.get(level, 0.5)

                # Avoid division by zero
                distance_score = 1.0 / max(distance, 1.0)
                availability_score = math.log(sequence_count + 1)
                priority_score = distance_score * availability_score * rank_weight

                # Add to candidates
                all_candidates.append({
                    "taxid": sibling_taxid,
                    "scientific_name": sibling_name,
                    "rank": sibling_rank,
                    "phylogenetic_distance": distance,
                    "sequence_count": sequence_count,
                    "priority_score": priority_score,
                    "lineage": sibling.get("lineage", ""),
                    "diversity_level": level
                })

                logger.debug(f"  {sibling_name}: distance={distance}, seqs={sequence_count}, score={priority_score:.2f}")

        # Sort by priority score (descending)
        all_candidates.sort(key=lambda x: x["priority_score"], reverse=True)

        # Select top candidates, ensuring diversity across levels
        selected_taxa = _select_diverse_taxa(all_candidates, max_results, diversity_levels)

        logger.info(f"Selected {len(selected_taxa)} weighted exclusion taxa")
        for taxon in selected_taxa:
            logger.info(f"  {taxon['scientific_name']} ({taxon['diversity_level']}): "
                       f"score={taxon['priority_score']:.2f}, seqs={taxon['sequence_count']}")

        return selected_taxa

    except Exception as e:
        logger.error(f"Error finding weighted related taxa for {taxid}: {str(e)}")
        import traceback
        logger.debug(traceback.format_exc())
        return []

def _calculate_taxonomic_distance(level: str, inclusion_rank: str) -> int:
    """
    Calculate phylogenetic distance in taxonomic steps.

    Lower numbers indicate closer relationships.

    Args:
        level: Taxonomic level where sibling was found (genus, family, order)
        inclusion_rank: Rank of the inclusion taxon

    Returns:
        int: Distance in taxonomic steps
    """
    # Standard taxonomic hierarchy
    rank_hierarchy = {
        "subspecies": 1,
        "species": 2,
        "genus": 3,
        "family": 4,
        "order": 5,
        "class": 6,
        "phylum": 7,
        "kingdom": 8,
        "superkingdom": 9
    }

    # Get numeric ranks
    level_num = rank_hierarchy.get(level, 5)
    inclusion_num = rank_hierarchy.get(inclusion_rank, 2)

    # Distance is the absolute difference in rank levels
    distance = abs(level_num - inclusion_num)

    # If siblings are at the same level, distance is 1
    if distance == 0:
        distance = 1

    return distance

def _estimate_sequence_count(taxid: str, email: str, timeout: int = 5) -> int:
    """
    Estimate the number of sequences available for a taxon.

    Uses a quick NCBI ESearch query to count available sequences.

    Args:
        taxid: NCBI taxonomy ID
        email: User's email
        timeout: Maximum time to wait (seconds)

    Returns:
        int: Estimated sequence count
    """
    Entrez.email = email

    try:
        # Quick search to count sequences
        query = f"txid{taxid}[Organism] AND biomol_genomic[PROP] NOT wgs[Filter]"

        search_handle = Entrez.esearch(db="nucleotide", term=query, retmax=0)
        search_results = Entrez.read(search_handle)
        search_handle.close()

        count = int(search_results.get("Count", 0))
        return count

    except Exception as e:
        logger.debug(f"Error estimating sequence count for taxid {taxid}: {e}")
        return 0

def _select_diverse_taxa(
    candidates: List[Dict[str, Any]],
    max_results: int,
    diversity_levels: List[str]
) -> List[Dict[str, Any]]:
    """
    Select taxa ensuring diversity across taxonomic levels.

    Algorithm:
    - Allocate slots proportionally across diversity levels
    - Within each level, select highest-scoring taxa
    - If a level has fewer candidates than slots, redistribute

    Args:
        candidates: All candidate taxa (sorted by priority score)
        max_results: Maximum number to select
        diversity_levels: Taxonomic levels to sample

    Returns:
        List of selected taxa
    """
    if not candidates:
        return []

    # Group candidates by diversity level
    by_level = {level: [] for level in diversity_levels}
    for candidate in candidates:
        level = candidate.get("diversity_level")
        if level in by_level:
            by_level[level].append(candidate)

    # Calculate allocation per level
    # Prioritize closer relatives (genus > family > order)
    level_weights = {"genus": 0.5, "family": 0.3, "order": 0.2}

    # Adjust weights based on availability
    total_weight = 0
    available_weights = {}
    for level in diversity_levels:
        if by_level[level]:  # Only include levels with candidates
            weight = level_weights.get(level, 0.1)
            available_weights[level] = weight
            total_weight += weight

    # Normalize weights
    if total_weight > 0:
        for level in available_weights:
            available_weights[level] /= total_weight

    # Allocate slots
    selected = []
    remaining_slots = max_results

    for level in diversity_levels:
        if not by_level[level] or remaining_slots <= 0:
            continue

        # Calculate number of slots for this level
        target_slots = int(available_weights.get(level, 0) * max_results)
        actual_slots = min(target_slots, len(by_level[level]), remaining_slots)

        # Select top-scoring taxa from this level
        selected.extend(by_level[level][:actual_slots])
        remaining_slots -= actual_slots

    # If we still have remaining slots, add highest-scoring candidates regardless of level
    if remaining_slots > 0:
        already_selected_ids = {t["taxid"] for t in selected}
        for candidate in candidates:
            if candidate["taxid"] not in already_selected_ids:
                selected.append(candidate)
                remaining_slots -= 1
                if remaining_slots <= 0:
                    break

    return selected[:max_results]

def intelligent_exclusion_selection(
    inclusion_taxid: str,
    email: str,
    max_exclusion_taxa: int = 10,
    require_sequences: bool = True,
    min_sequence_count: int = 5,
    tier_allocation: Optional[Dict[str, float]] = None
) -> Dict[str, Any]:
    """
    Intelligently select exclusion taxa using tiered phylogenetic strategy.

    This function automatically builds a representative exclusion set by:
    1. Determining the inclusion taxon's rank
    2. Sampling taxa from multiple taxonomic levels (genus, family, order)
    3. Allocating slots proportionally by tier (default: 50% genus, 30% family, 20% order)
    4. Prioritizing by phylogenetic proximity and sequence availability

    Args:
        inclusion_taxid (str): NCBI TaxID for inclusion taxon
        email (str): User's email for NCBI queries
        max_exclusion_taxa (int): Maximum number of exclusion taxa (default: 10)
        require_sequences (bool): Only include taxa with sufficient sequences (default: True)
        min_sequence_count (int): Minimum sequences required per taxon (default: 5)
        tier_allocation (Optional[Dict[str, float]]): Custom tier weights
            Default: {"genus": 0.5, "family": 0.3, "order": 0.2}

    Returns:
        Dict[str, Any]: Exclusion selection results with fields:
            - exclusion_taxids: List[str] - Selected taxon IDs
            - taxa_info: List[Dict] - Full metadata for each taxon
            - selection_strategy: str - Description of selection approach
            - coverage_score: float - Phylogenetic breadth (0-1)
            - inclusion_info: Dict - Info about inclusion taxon
            - tier_summary: Dict - Count of taxa per tier

    Example:
        >>> result = intelligent_exclusion_selection(
        ...     inclusion_taxid="689",  # Vibrio mediterranei
        ...     email="user@example.com",
        ...     max_exclusion_taxa=10
        ... )
        >>> print(f"Selected {len(result['exclusion_taxids'])} exclusion taxa")
        >>> print(f"Coverage score: {result['coverage_score']:.2f}")
    """
    if not email:
        raise ValueError("Email address is required.")

    if tier_allocation is None:
        tier_allocation = {"genus": 0.5, "family": 0.3, "order": 0.2}

    Entrez.email = email

    try:
        # Get information about the inclusion taxon
        logger.info(f"Starting intelligent exclusion selection for taxid {inclusion_taxid}")
        inclusion_info = get_taxon_info(inclusion_taxid, email)

        if not inclusion_info:
            return {
                "error": f"Could not retrieve information for taxid {inclusion_taxid}",
                "exclusion_taxids": [],
                "taxa_info": []
            }

        inclusion_rank = inclusion_info.get("rank", "")
        inclusion_name = inclusion_info.get("scientific_name", "Unknown")

        logger.info(f"Inclusion taxon: {inclusion_name} (rank: {inclusion_rank})")

        # Determine appropriate diversity levels based on inclusion rank
        diversity_levels = _determine_diversity_levels(inclusion_rank)
        logger.info(f"Using diversity levels: {diversity_levels}")

        # Calculate target number of taxa per tier
        tier_targets = {}
        total_weight = sum(tier_allocation.get(level, 0) for level in diversity_levels)

        for level in diversity_levels:
            weight = tier_allocation.get(level, 0)
            if total_weight > 0:
                tier_targets[level] = max(1, int((weight / total_weight) * max_exclusion_taxa))
            else:
                tier_targets[level] = max(1, max_exclusion_taxa // len(diversity_levels))

        logger.info(f"Target allocation: {tier_targets}")

        # Get weighted related taxa spanning all diversity levels
        all_weighted_taxa = get_related_taxa_weighted(
            taxid=inclusion_taxid,
            email=email,
            max_results=max_exclusion_taxa * 3,  # Get more candidates than needed
            diversity_levels=diversity_levels,
            min_sequence_count=min_sequence_count if require_sequences else 0
        )

        if not all_weighted_taxa:
            logger.warning("No suitable exclusion taxa found with current criteria")

            # Try again with relaxed requirements
            if require_sequences and min_sequence_count > 1:
                logger.info("Retrying with relaxed sequence requirements...")
                all_weighted_taxa = get_related_taxa_weighted(
                    taxid=inclusion_taxid,
                    email=email,
                    max_results=max_exclusion_taxa * 2,
                    diversity_levels=diversity_levels,
                    min_sequence_count=1  # Very relaxed
                )

        if not all_weighted_taxa:
            # Last resort: try with even broader levels
            logger.warning("Still no taxa found, trying broader taxonomic levels...")
            broader_levels = diversity_levels + ["class", "phylum"]
            all_weighted_taxa = get_related_taxa_weighted(
                taxid=inclusion_taxid,
                email=email,
                max_results=max_exclusion_taxa,
                diversity_levels=broader_levels,
                min_sequence_count=0  # No minimum
            )

        if not all_weighted_taxa:
            return {
                "error": "Could not find any suitable exclusion taxa",
                "exclusion_taxids": [],
                "taxa_info": [],
                "inclusion_info": inclusion_info,
                "selection_strategy": "Failed - no suitable relatives found"
            }

        # Select top taxa (already sorted by priority score in get_related_taxa_weighted)
        selected_taxa = all_weighted_taxa[:max_exclusion_taxa]

        # Extract taxids
        exclusion_taxids = [taxon["taxid"] for taxon in selected_taxa]

        # Calculate coverage score based on phylogenetic diversity
        coverage_score = _calculate_coverage_score(selected_taxa, diversity_levels)

        # Summarize tier distribution
        tier_summary = {}
        for level in diversity_levels:
            count = sum(1 for t in selected_taxa if t.get("diversity_level") == level)
            tier_summary[level] = count

        # Build selection strategy description
        strategy_parts = []
        for level, count in tier_summary.items():
            if count > 0:
                strategy_parts.append(f"{count} {level}-level")

        selection_strategy = f"Tiered selection: {', '.join(strategy_parts)} relatives"

        logger.info(f"Selected {len(selected_taxa)} exclusion taxa")
        logger.info(f"Tier distribution: {tier_summary}")
        logger.info(f"Coverage score: {coverage_score:.2f}")

        return {
            "exclusion_taxids": exclusion_taxids,
            "taxa_info": selected_taxa,
            "selection_strategy": selection_strategy,
            "coverage_score": coverage_score,
            "inclusion_info": inclusion_info,
            "tier_summary": tier_summary,
            "diversity_levels": diversity_levels
        }

    except Exception as e:
        logger.error(f"Error in intelligent exclusion selection: {str(e)}")
        import traceback
        logger.debug(traceback.format_exc())

        return {
            "error": str(e),
            "exclusion_taxids": [],
            "taxa_info": []
        }

def _determine_diversity_levels(inclusion_rank: str) -> List[str]:
    """
    Determine appropriate diversity levels based on inclusion taxon rank.

    Args:
        inclusion_rank: Rank of the inclusion taxon

    Returns:
        List of taxonomic ranks to sample for exclusion

    Examples:
        - If inclusion is species → sample genus, family, order
        - If inclusion is genus → sample family, order, class
        - If inclusion is family → sample order, class, phylum
    """
    rank_hierarchy = [
        "subspecies",
        "species",
        "genus",
        "family",
        "order",
        "class",
        "phylum",
        "kingdom"
    ]

    # Find the index of the inclusion rank
    try:
        inclusion_idx = rank_hierarchy.index(inclusion_rank)
    except ValueError:
        # Default to species level if rank not recognized
        logger.warning(f"Unknown rank '{inclusion_rank}', defaulting to species-level sampling")
        inclusion_idx = 1

    # Select 3 levels above the inclusion rank
    diversity_levels = []
    for offset in [1, 2, 3]:
        target_idx = inclusion_idx + offset
        if target_idx < len(rank_hierarchy):
            diversity_levels.append(rank_hierarchy[target_idx])

    # Ensure we have at least some levels
    if not diversity_levels:
        diversity_levels = ["genus", "family", "order"]

    # Limit to commonly used levels
    valid_levels = ["genus", "family", "order", "class", "phylum"]
    diversity_levels = [level for level in diversity_levels if level in valid_levels]

    if not diversity_levels:
        diversity_levels = ["genus", "family", "order"]

    return diversity_levels

def _calculate_coverage_score(selected_taxa: List[Dict[str, Any]], diversity_levels: List[str]) -> float:
    """
    Calculate phylogenetic coverage score based on diversity.

    Score factors:
    1. Number of taxonomic levels represented (more = better)
    2. Number of taxa (more = better, up to diminishing returns)
    3. Spread of phylogenetic distances (wider = better)

    Args:
        selected_taxa: List of selected exclusion taxa with metadata
        diversity_levels: Target diversity levels

    Returns:
        float: Coverage score between 0 and 1
            - 1.0 = Excellent coverage across all levels
            - 0.7-0.9 = Good coverage
            - 0.5-0.7 = Fair coverage
            - <0.5 = Limited coverage
    """
    if not selected_taxa:
        return 0.0

    # Factor 1: Level diversity (0-0.4 points)
    levels_represented = set(t.get("diversity_level") for t in selected_taxa)
    level_score = len(levels_represented) / max(len(diversity_levels), 1)
    level_score = min(level_score, 1.0) * 0.4

    # Factor 2: Number of taxa (0-0.3 points, with diminishing returns)
    import math
    taxa_count = len(selected_taxa)
    # Use logarithmic scale to give diminishing returns after ~10 taxa
    taxa_score = min(math.log(taxa_count + 1) / math.log(15), 1.0) * 0.3

    # Factor 3: Distance spread (0-0.3 points)
    distances = [t.get("phylogenetic_distance", 1) for t in selected_taxa]
    if distances:
        min_dist = min(distances)
        max_dist = max(distances)
        distance_range = max_dist - min_dist

        # Ideal range is having both close (distance=1) and far (distance>=3) relatives
        ideal_range = 3
        distance_score = min(distance_range / ideal_range, 1.0) * 0.3
    else:
        distance_score = 0.0

    total_score = level_score + taxa_score + distance_score

    return round(total_score, 2)

def suggest_marker_genes(taxid: str, email: str) -> List[Dict[str, str]]:
    """
    Suggest appropriate marker genes for a given taxon.
    
    Args:
        taxid (str): NCBI TaxID
        email (str): User's email address
        
    Returns:
        List[Dict[str, str]]: List of suggested marker genes with descriptions
    """
    # Get taxonomy info to determine the organism type
    tax_info = get_taxon_info(taxid, email)
    
    if not tax_info:
        return []
    
    lineage = tax_info.get("lineage", "").lower()
    
    # Default markers for all organisms
    markers = [
        {"gene": "16S ribosomal RNA", "description": "Universal prokaryotic marker gene"},
    ]
    
    # Add organism-specific markers
    if "bacteria" in lineage:
        markers.extend([
            {"gene": "rpoB", "description": "RNA polymerase beta subunit"},
            {"gene": "gyrB", "description": "DNA gyrase subunit B"},
            {"gene": "recA", "description": "Recombination protein RecA"},
            {"gene": "groL", "description": "60 kDa chaperonin (cpn60; hsp60)"}
        ])
    elif "fungi" in lineage:
        markers.extend([
            {"gene": "ITS", "description": "Internal transcribed spacer"},
            {"gene": "18S ribosomal RNA", "description": "Small subunit ribosomal RNA"},
            {"gene": "28S ribosomal RNA", "description": "Large subunit ribosomal RNA"},
            {"gene": "RPB2", "description": "RNA polymerase II second largest subunit"}
        ])
    elif "virus" in lineage:
        markers.extend([
            {"gene": "polymerase", "description": "Viral polymerase gene"},
            {"gene": "capsid", "description": "Viral capsid protein gene"},
            {"gene": "envelope", "description": "Viral envelope protein gene"}
        ])
    
    return markers

def load_local_fasta(fasta_path: str) -> List[SeqRecord]:
    """
    Load sequences from a local FASTA file using Biopython.
    
    Args:
        fasta_path (str): Path to the local FASTA file.
        
    Returns:
        List[SeqRecord]: A list of sequences parsed from the FASTA file.
    """
    logger.info(f"Loading local FASTA file from {fasta_path}")
    try:
        seq_records = list(SeqIO.parse(fasta_path, "fasta"))
        logger.info(f"Loaded {len(seq_records)} sequences from {fasta_path}")
        return seq_records
    except Exception as e:
        logger.error(f"Error loading local FASTA file {fasta_path}: {str(e)}")
        raise

def fetch_cds_sequences(
    cds_name: str,
    email: str,
    max_seqs: int = 100,
    output_file: Optional[str] = None,
    database: str = "nuccore"
) -> List[SeqRecord]:
    """
    Fetch coding sequences (CDS) with specific gene annotation from NCBI RefSeq.
    Only downloads the annotated gene sequences, not entire genomes.
    
    Args:
        cds_name (str): Name of the CDS to search for (e.g., 'amoA')
        email (str): User's email address (required by NCBI)
        max_seqs (int): Maximum number of sequences to fetch
        output_file (Optional[str]): Output file path, if None will generate based on cds_name
        database (str): NCBI database to search (default: nuccore)
        
    Returns:
        List[SeqRecord]: List of retrieved CDS sequences
    """
    if not email:
        raise ValueError("Email address is required for NCBI queries")
        
    Entrez.email = email
    
    # Generate output filename if not provided
    if not output_file:
        output_file = f"{cds_name}_refseq_{max_seqs}.fna"
    
    logger.info(f"Searching for CDS '{cds_name}' in NCBI (max: {max_seqs} sequences)")
    
    # Construct search query for entries with specified CDS
    query = f"{cds_name}[Gene] OR {cds_name}[Product] OR {cds_name}[Function] AND CDS[Feature Key] NOT wgs[Filter]"
    
    try:
        # Search for matching entries
        search_handle = Entrez.esearch(db=database, term=query, retmax=max_seqs)
        search_results = Entrez.read(search_handle)
        search_handle.close()
        
        ids = search_results.get("IdList", [])
        found_count = len(ids)
        
        if not ids:
            logger.warning(f"No CDS entries found for '{cds_name}'")
            
            # Try with a more general search
            alt_query = f"{cds_name} AND genbank[Filter]"
            search_handle = Entrez.esearch(db=database, term=alt_query, retmax=max_seqs)
            search_results = Entrez.read(search_handle)
            search_handle.close()
            
            ids = search_results.get("IdList", [])
            if not ids:
                logger.warning(f"No entries found with alternative search for '{cds_name}'")
                return []
            
            logger.info(f"Found {len(ids)} entries with alternative search. Will extract CDS features.")
        
        logger.info(f"Found {len(ids)} entries. Retrieving CDS sequences...")
        
        # Fetch CDS sequences
        cds_records = []
        
        # Process in batches to avoid overwhelming NCBI
        batch_size = 50
        for i in range(0, len(ids), batch_size):
            batch_ids = ids[i:i+batch_size]
            logger.info(f"Fetching batch {i//batch_size + 1}/{(len(ids)-1)//batch_size + 1}...")
            
            # Fetch the GenBank records that contain the CDS features
            fetch_handle = Entrez.efetch(db=database, id=",".join(batch_ids), 
                                        rettype="gb", retmode="text")
            gb_records = list(SeqIO.parse(fetch_handle, "genbank"))
            fetch_handle.close()
            
            # Extract CDS features matching our target gene
            for gb_record in gb_records:
                for feature in gb_record.features:
                    if feature.type == "CDS":
                        try:
                            gene_qualifiers = feature.qualifiers.get("gene", [])
                            product_qualifiers = feature.qualifiers.get("product", [])
                            locus_tag = feature.qualifiers.get("locus_tag", [""])[0]
                        
                            # Check if this CDS matches our target gene
                            if any(cds_name.lower() in gene.lower() for gene in gene_qualifiers) or \
                               any(cds_name.lower() in product.lower() for product in product_qualifiers) or \
                               cds_name.lower() in locus_tag.lower():
                               
                                # Add error checking before extraction
                                if feature.location is None:
                                    logger.warning(f"Skipping CDS with undefined location in {gb_record.id}")
                                    continue
                                   
                                # Extract the CDS sequence only (not the entire genome)
                                try:
                                    cds_seq = feature.extract(gb_record.seq)
                                except Exception as extract_err:
                                    logger.warning(f"Could not extract CDS sequence from {gb_record.id}: {extract_err}")
                                    continue
                                    
                                # Verify sequence content is valid
                                if not cds_seq or len(cds_seq) == 0:
                                    logger.warning(f"Skipping empty CDS sequence in {gb_record.id}")
                                    continue
                            
                                # Get organism and taxonomy information
                                organism = gb_record.annotations.get("organism", "Unknown organism")
                                taxonomy = gb_record.annotations.get("taxonomy", [])
                                taxid = ""
                            
                                # Try to extract taxid from the source feature
                                for src_feature in gb_record.features:
                                    if src_feature.type == "source":
                                        db_xrefs = src_feature.qualifiers.get("db_xref", [])
                                        for xref in db_xrefs:
                                            if xref.startswith("taxon:"):
                                                taxid = xref.replace("taxon:", "")
                                                break
                            
                                # Create a descriptive ID with gene, organism, and accession
                                gene_name = gene_qualifiers[0] if gene_qualifiers else cds_name
                            
                                # Format description to include organism and taxid
                                description = f"{gene_name} gene from {organism}"
                                if taxid:
                                    description += f" (taxid:{taxid})"
                                
                                # Create a unique ID that includes accession and gene
                                record_id = f"{gb_record.id}_{gene_name}"
                            
                                # Add a check to ensure we're not getting overly short sequences
                                if len(cds_seq) < 30:
                                    logger.debug(f"Skipping very short CDS sequence: {len(cds_seq)} bp")
                                    continue
                                
                                # Create the SeqRecord
                                cds_record = SeqRecord(
                                    seq=cds_seq,
                                    id=record_id,
                                    name=gene_name,
                                    description=description
                                )
                            
                                # Add extra annotations
                                cds_record.annotations["organism"] = organism
                                if taxid:
                                    cds_record.annotations["taxid"] = taxid
                            
                                cds_records.append(cds_record)
                        except Exception as feature_err:
                            logger.warning(f"Error processing feature in {gb_record.id}: {feature_err}")
                            continue      
        
        # Log statistics about retrieved sequences
        if cds_records:
            # Filter out any problematic sequences
            valid_records = []
            for record in cds_records:
                try:
                    # Verify the sequence is valid
                    if record.seq and len(record.seq) > 0:
                        valid_records.append(record)
                    else:
                        logger.warning(f"Skipping record with invalid sequence: {record.id}")
                except Exception as e:
                    logger.warning(f"Error validating sequence {record.id}: {e}")
                    
            # Use the validated records
            cds_records = valid_records
            
            if not cds_records:
                logger.warning("No valid CDS sequences remain after validation")
                return []
                
                
            avg_length = sum(len(rec.seq) for rec in cds_records) / len(cds_records)
            logger.info(f"Retrieved {len(cds_records)} CDS sequences with average length {avg_length:.1f} bp")
            
            if len(cds_records) < max_seqs:
                logger.info(f"Note: Found fewer sequences ({len(cds_records)}) than requested ({max_seqs})")
                
            # Ensure some diversity in taxa if possible
            taxa_count = len(set(rec.annotations.get("taxid", "") for rec in cds_records if "taxid" in rec.annotations))
            logger.info(f"Retrieved sequences represent {taxa_count} different taxa")
            
            # Save to FASTA file
            try:
                SeqIO.write(cds_records, output_file, "fasta")
                logger.info(f"Saved {len(cds_records)} CDS sequences to {output_file}")
            except Exception as write_err:
                logger.error(f"Error writing sequences to file: {write_err}")
        
        return cds_records
        
    except Exception as e:
        logger.error(f"Error fetching CDS sequences: {str(e)}")
        raise    

def download_genome(
    accession_or_taxid: str,
    email: str,
    output_path: Optional[str] = None,
    is_taxid: bool = False
) -> str:
    """
    Download a complete genome by accession or taxid.
    
    Args:
        accession_or_taxid (str): NCBI accession or taxid
        email (str): User's email address
        output_path (Optional[str]): Path to save the genome
        is_taxid (bool): Whether the identifier is a taxid
        
    Returns:
        str: Path to the saved genome file
    """
    import os
    
    Entrez.email = email
    
    try:
        if is_taxid:
            # Search for the reference genome for this taxid
            logger.info(f"Searching for reference genome for taxid {accession_or_taxid}")
            query = f"txid{accession_or_taxid}[Organism] AND refseq[Filter] AND representative[Properties]"
            search_handle = Entrez.esearch(db="genome", term=query)
            search_results = Entrez.read(search_handle)
            search_handle.close()
            
            if not search_results.get("IdList"):
                # Try assembly database
                search_handle = Entrez.esearch(db="assembly", term=query)
                search_results = Entrez.read(search_handle)
                search_handle.close()
                
                if not search_results.get("IdList"):
                    raise ValueError(f"No reference genome found for taxid {accession_or_taxid}")
                
                # Get the assembly details
                assembly_id = search_results["IdList"][0]
                handle = Entrez.esummary(db="assembly", id=assembly_id)
                assembly_record = Entrez.read(handle)
                handle.close()
                
                # Extract the accession
                accession = assembly_record["DocumentSummarySet"]["DocumentSummary"][0]["AssemblyAccession"]
            else:
                # Get the genome record
                genome_id = search_results["IdList"][0]
                handle = Entrez.esummary(db="genome", id=genome_id)
                genome_record = Entrez.read(handle)
                handle.close()
                
                # Extract the accession
                accession = genome_record[0].get("AssemblyAccession")
                
                if not accession:
                    raise ValueError(f"Could not extract accession for taxid {accession_or_taxid}")
        else:
            accession = accession_or_taxid
        
        # Fetch the sequence
        logger.info(f"Fetching genome for accession {accession}")
        handle = Entrez.efetch(db="nucleotide", id=accession, rettype="fasta", retmode="text")
        record = SeqIO.read(handle, "fasta")
        handle.close()
        
        # Save to file
        if not output_path:
            output_path = f"{accession.replace('.', '_')}.fasta"
        
        SeqIO.write(record, output_path, "fasta")
        logger.info(f"Genome saved to {output_path}")
        
        return output_path
        
    except Exception as e:
        logger.error(f"Error downloading genome: {str(e)}")
        raise
        
def ensure_sequences_are_defined(sequences: List[SeqRecord], email: str) -> List[SeqRecord]:
    """
    Ensure all sequences have defined content before saving.
    
    Args:
        sequences: List of sequence records
        
    Returns:
        List of validated sequence records
    """
    from Bio.Seq import UndefinedSequenceError, Seq
    from Bio.SeqRecord import SeqRecord
    
    valid_sequences = []
    
    for i, record in enumerate(sequences):
        try:
            # Check if sequence is defined by attempting to convert to string
            seq_str = str(record.seq)
            valid_sequences.append(record)
        except UndefinedSequenceError:
            # Extract accession and try to fetch directly
            accession = record.id.split('.')[0]
            try:
                # Only try to fetch it if it looks like a sequence accession (and not a genome assembly)
                if not accession.endswith("000000000"):
                    logger.info(f"Trying to fetch sequence {accession} directly")
                    direct_record = fetch_sequence_by_accession(accession, email)
                    valid_sequences.append(direct_record)
                else:
                    logger.warning(f"Skipping genome assembly ID {accession}")
            except Exception as e:
                logger.warning(f"Failed to fetch direct sequence for {accession}: {e}")
            
    return valid_sequences
    
def load_local_fasta(fasta_path: str) -> List[SeqRecord]:
    """
    Load sequences from a local FASTA file using Biopython.
    
    Args:
        fasta_path (str): Path to the local FASTA file.
        
    Returns:
        List[SeqRecord]: A list of sequences parsed from the FASTA file.
    """
    logger.info(f"Loading local FASTA file from {fasta_path}")
    try:
        seq_records = list(SeqIO.parse(fasta_path, "fasta"))
        logger.info(f"Loaded {len(seq_records)} sequences from {fasta_path}")
        return seq_records
    except Exception as e:
        logger.error(f"Error loading local FASTA file {fasta_path}: {str(e)}")
        raise