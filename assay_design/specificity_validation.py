# assay_design/specificity_validation.py

import logging
import re
from typing import List, Dict, Any, Optional, Union
import subprocess
import tempfile
from Bio import SeqIO, Entrez
from Bio.Seq import Seq

logger = logging.getLogger(__name__)

def validate_primer_specificity(
    primers: List[Dict[str, Any]],
    inclusion_taxid: str,
    exclusion_taxids: List[str],
    email: str = None,
    max_blast_hits: int = 10
) -> Dict[str, Any]:
    """
    Validate primer specificity against inclusion and exclusion taxa.
    
    Args:
        primers (List[Dict[str, Any]]): List of designed primers
        inclusion_taxid (str): NCBI taxonomy ID for target organism
        exclusion_taxids (List[str]): NCBI taxonomy IDs for non-target organisms
        email (str): Email for NCBI queries
        max_blast_hits (int): Maximum number of BLAST hits to consider
        
    Returns:
        Dict[str, Any]: Validation results including specificity metrics
    """
    if not primers or len(primers) < 2:
        logger.error("At least two primers (forward and reverse) are required for validation")
        return {"error": "Insufficient primers provided"}
    
    # Extract primer sequences
    forward_primer = primers[0].get("sequence", "")
    reverse_primer = primers[1].get("sequence", "")
    
    if not forward_primer or not reverse_primer:
        logger.error("Invalid primer sequences provided")
        return {"error": "Invalid primer sequences"}
    
    # Placeholder for real validation results
    # In a production implementation, this would:
    # 1. Run in silico PCR using e-PCR or similar
    # 2. BLAST primers against NCBI database
    # 3. Check for potential cross-reactivity
    
    logger.info("Validating primer specificity")
    
    # Simulated validation results
    validation_results = {
        "inclusion_hits": simulate_inclusion_hits(forward_primer, reverse_primer, inclusion_taxid),
        "exclusion_hits": simulate_exclusion_hits(forward_primer, reverse_primer, exclusion_taxids),
        "specificity_score": 0.92,  # Placeholder value
        "potential_issues": []
    }
    
    # Calculate specificity metrics
    if validation_results["exclusion_hits"]:
        validation_results["potential_issues"].append(
            f"Potential cross-reactivity with {len(validation_results['exclusion_hits'])} non-target taxa"
        )
        validation_results["specificity_score"] = 0.5  # Reduce score for cross-reactivity
    
    return validation_results

def simulate_inclusion_hits(forward_primer: str, reverse_primer: str, inclusion_taxid: str) -> List[Dict[str, Any]]:
    """
    Simulate in silico PCR hits for inclusion taxid.
    
    Args:
        forward_primer (str): Forward primer sequence
        reverse_primer (str): Reverse primer sequence
        inclusion_taxid (str): Target taxid
        
    Returns:
        List[Dict[str, Any]]: Simulated PCR products from target taxon
    """
    # In a real implementation, this would query NCBI or run a local in silico PCR
    
    # Simulated positive hit
    return [
        {
            "accession": f"NZ_XXXX01000001.1",
            "organism": f"Target organism (taxid:{inclusion_taxid})",
            "product_size": 150,
            "sequence": "ATGC" * 50,  # Placeholder sequence
            "forward_position": 1000,
            "reverse_position": 1150
        }
    ]

def simulate_exclusion_hits(forward_primer: str, reverse_primer: str, exclusion_taxids: List[str]) -> List[Dict[str, Any]]:
    """
    Simulate in silico PCR hits for exclusion taxids.
    
    Args:
        forward_primer (str): Forward primer sequence
        reverse_primer (str): Reverse primer sequence
        exclusion_taxids (List[str]): Non-target taxids
        
    Returns:
        List[Dict[str, Any]]: Simulated PCR products from non-target taxa
    """
    # In a real implementation, this would query NCBI or run a local in silico PCR
    
    # For the prototype, we'll simulate no hits in exclusion group (good specificity)
    return []

def blast_primers(primer_sequence: str, email: str, taxid: Optional[str] = None) -> List[Dict[str, Any]]:
    """
    BLAST a primer sequence against NCBI database.
    
    Args:
        primer_sequence (str): Primer sequence to query
        email (str): Email for NCBI queries
        taxid (Optional[str]): Limit search to specific taxid
        
    Returns:
        List[Dict[str, Any]]: BLAST hits
    """
    from Bio.Blast import NCBIWWW, NCBIXML
    
    # Set NCBI email
    Entrez.email = email
    
    logger.info(f"BLASTing primer sequence of length {len(primer_sequence)}")
    
    taxon_limit = f"txid{taxid}" if taxid else ""
    
    # Perform BLAST search
    try:
        # This would be a real BLAST search in a production implementation
        # result_handle = NCBIWWW.qblast("blastn", "nt", primer_sequence, 
        #                               entrez_query=taxon_limit,
        #                               expect=1000, 
        #                               word_size=7)
        # blast_records = NCBIXML.parse(result_handle)
        # hits = process_blast_records(blast_records)
        
        # For the prototype, return simulated results
        hits = [
            {
                "accession": "NZ_XXXX01000001.1",
                "organism": "Simulated organism",
                "score": 40.1,
                "e_value": 0.001,
                "identity": 100.0,
                "alignment_length": len(primer_sequence)
            }
        ]
        
        return hits
        
    except Exception as e:
        logger.error(f"BLAST error: {str(e)}")
        return []

def run_in_silico_pcr(forward_primer: str, 
                     reverse_primer: str, 
                     template_sequences: List[Any],
                     max_product_size: int = 2000) -> List[Dict[str, Any]]:
    """
    Run in silico PCR on template sequences.
    
    Args:
        forward_primer (str): Forward primer sequence
        reverse_primer (str): Reverse primer sequence
        template_sequences (List[Any]): List of sequence templates (SeqRecord objects)
        max_product_size (int): Maximum allowed PCR product size
        
    Returns:
        List[Dict[str, Any]]: In silico PCR results
    """
    pcr_products = []
    
    # Reverse complement the reverse primer for accurate matching
    reverse_primer_rc = str(Seq(reverse_primer).reverse_complement())
    
    for template in template_sequences:
        template_seq = str(template.seq).upper()
        
        # Find forward primer binding sites
        forward_matches = find_approximate_matches(forward_primer, template_seq)
        
        # Find reverse primer binding sites
        reverse_matches = find_approximate_matches(reverse_primer_rc, template_seq)
        
        # Check for valid PCR products
        for fwd_pos in forward_matches:
            for rev_pos in reverse_matches:
                # Ensure correct orientation (forward -> reverse)
                if fwd_pos < rev_pos:
                    product_size = rev_pos - fwd_pos + len(reverse_primer)
                    
                    # Check if product size is reasonable
                    if product_size <= max_product_size:
                        product_seq = template_seq[fwd_pos:rev_pos + len(reverse_primer)]
                        
                        pcr_products.append({
                            "template_id": template.id,
                            "forward_position": fwd_pos,
                            "reverse_position": rev_pos + len(reverse_primer),
                            "product_size": product_size,
                            "product_sequence": product_seq
                        })
    
    return pcr_products

def find_approximate_matches(primer: str, sequence: str, max_mismatches: int = 1) -> List[int]:
    """
    Find approximate matches of a primer in a sequence, allowing mismatches.
    
    Args:
        primer (str): Primer sequence to search for
        sequence (str): Template sequence to search in
        max_mismatches (int): Maximum allowed mismatches
        
    Returns:
        List[int]: Positions of approximate matches
    """
    # For a real implementation, use regex with mismatches or a dedicated
    # sequence alignment algorithm. This is a simplified placeholder.
    
    primer = primer.upper()
    sequence = sequence.upper()
    primer_len = len(primer)
    positions = []
    
    # Naive implementation allowing mismatches
    for i in range(len(sequence) - primer_len + 1):
        mismatches = 0
        for j in range(primer_len):
            if sequence[i + j] != primer[j]:
                mismatches += 1
                if mismatches > max_mismatches:
                    break
        
        if mismatches <= max_mismatches:
            positions.append(i)
    
    return positions