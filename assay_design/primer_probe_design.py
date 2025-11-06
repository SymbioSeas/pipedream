# assay_design/primer_probe_design.py

import logging
import re
from typing import List, Dict, Any, Union, Tuple
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import MeltingTemp as mt
from Bio.Seq import Seq

logger = logging.getLogger(__name__)

# Quality control tresholds
MIN_PRIMER_LENGTH = 12
MAX_PRIMER_LENGTH = 30
MIN_PROBE_LENGTH = 15
MAX_PROBE_LENGTH = 30
MIN_AMPLICON_LENGTH = 50
MAX_AMPLICON_LENGTH = 200
MIN_PRIMER_TM = 50.0
MAX_PRIMER_TM = 65.0
MIN_PROBE_TM = 55.0
MAX_PROBE_TM = 55.0
MIN_GC_CONTENT = 40.0
MAX_GC_CONTENT = 60.0

def calculate_oligo_properties(sequence: str) -> Dict[str, Any]:
    """
    Calculate key properties of an oligonucleotide sequence (primer or probe).
    
    Args:
        sequence (str): Oligonucleotide sequence
        
    Returns:
        Dict[str, Any]: Dictionary of properties
    """
    seq = Seq(sequence)
    properties = {
        "length": len(sequence),
        "gc_content": calculate_gc_content(sequence),
        "tm": mt.Tm_Wallace(sequence),  # Basic melting temp calculation
        "self_complementarity": check_self_complementarity(sequence),
        "has_repeats": has_repeats(sequence),
        "end_stability": check_3prime_stability(sequence)
    }
    return properties

def calculate_gc_content(sequence: str) -> float:
    """Calculate GC content of a sequence."""
    sequence = sequence.upper()
    gc_count = sequence.count('G') + sequence.count('C')
    return round((gc_count / len(sequence)) * 100, 2) if len(sequence) > 0 else 0

def check_self_complementarity(sequence: str) -> float:
    """
    Check for self-complementarity in a sequence.
    Returns a score where higher values indicate more self-complementarity.
    """
    # Simple implementation - more sophisticated algorithms would be used in production
    seq = sequence.upper()
    rev_comp = str(Seq(seq).reverse_complement())
    
    max_complementary = 0
    for i in range(len(seq)):
        match_count = 0
        for j in range(min(len(seq) - i, len(rev_comp))):
            if seq[i + j] == rev_comp[j]:
                match_count += 1
            else:
                match_count = 0
            max_complementary = max(max_complementary, match_count)
    
    return max_complementary

def has_repeats(sequence: str, max_repeat: int = 4) -> bool:
    """Check for nucleotide repeats in a sequence."""
    for base in "ATGC":
        if base * max_repeat in sequence.upper():
            return True
    return False

def check_3prime_stability(sequence: str, window: int = 5) -> float:
    """
    Calculate 3' end stability.
    Lower values are more stable (better for specificity).
    """
    if len(sequence) < window:
        return 0
    
    # Get the 3' end subsequence
    three_prime = sequence[-window:].upper()
    
    # Count GC content at 3' end as a simple stability measure
    gc_count = three_prime.count('G') + three_prime.count('C')
    return gc_count / window

def calculate_oligo_score(properties: Dict[str, Any], oligo_type: str = "primer") -> float:
    """
    Calculate an overall score for an oligonucleotide based on its properties.
    
    Args:
        properties (Dict[str, Any]): Oligonucleotide properties
        oligo_type (str): Type of oligonucleotide ("primer" or "probe")
        
    Returns:
        float: Score between 0 and 1 (higher is better)
    """
    # Set optimal values based on oligo type
    if oligo_type == "probe":
        optimal_gc = 50
        optimal_tm = 68  # Probes typically have higher Tm
        optimal_length = 25  # Probes are often longer
    else:  # primer
        optimal_gc = 50
        optimal_tm = 60
        optimal_length = 22
    
    # Calculate component scores
    gc_score = 1.0 - abs(properties["gc_content"] - optimal_gc) / 40
    tm_score = 1.0 - abs(properties["tm"] - optimal_tm) / 15
    length_score = 1.0 - abs(properties["length"] - optimal_length) / 15
    complementarity_score = 1.0 - properties["self_complementarity"] / 8
    
    # Weight the components
    score = (
        0.3 * gc_score +
        0.3 * tm_score +
        0.2 * length_score +
        0.2 * complementarity_score
    )
    
    return max(0, min(1, score))

def design_primers_and_probe(marker_info: Dict[str, Any]) -> Dict[str, Any]:
    """
    Design forward primer, reverse primer, and probe for identified specific regions. Returns quality feedback if design doesn't meet criteria.
    """
    # Extract the marker sequence
    marker_sequence = marker_info.get("marker_sequence", "")
    
    if not marker_sequence:
        logger.error("No marker sequence provided for primer design")
        return {
            "status": "failed",
            "message": "No marker sequence available",
            "primers": [],
            "probe": None,
            "amplicon": "",
            "amplicon_length": 0
        }
    
     # Check if we have specific regions identified
    forward_region = marker_info.get("forward_primer_region", "")
    probe_region = marker_info.get("probe_region", "")
    reverse_region = marker_info.get("reverse_primer_region", "")
    
    # If specific regions are provided, use them
    if forward_region and probe_region and reverse_region:
        forward_primer = design_primer_from_region(forward_region, "forward")
        probe = design_probe_from_region(probe_region)
        reverse_primer = design_primer_from_region(reverse_region, "reverse")
    
        # Check if design was successful
        if forward_primer and probe and reverse_primer:
            # Find positions in marker sequence
            fwd_pos = marker_sequence.find(forward_region)
            probe_pos =marker_sequence.find(probe_region)
            rev_pos = marker_sequence.find(reverse_region)
            
            # Verify positions are valid and in correct order
            if fwd_pos >= 0 and probe_pos > fwd_pos + len(forward_region) and rev_pos > + len(probe_region):
                # Calculate amplicon
                reverse_primer_end = rev_pos + len(reverse_region)
                amplicon = marker_sequence[fwd_pos:reverse_primer_end]
                
                return {
                    "status": "success",
                    "message": "Successfully designed primers and probe from specific regions",
                    "primers": [forward_primer, reverse_primer],
                    "probe": probe,
                    "amplicon": amplicon,
                    "amplicon_length": len(amplicon)
                }
                
        # If specific regions weren't provided or weren't valid, design primers for the entire marker
        primers, amplicon_info = design_primers_for_sequence(marker_sequence)
        
        if primers and len(primers) >= 2:
            # Design a probe for the amplicon with position constraints
            probe = design_probe_for_amplicon(
                amplicon_info,
                min_distance_from_primers=5 # Ensures at least 5bp distance from primers
            )
            
            # Verify primer and probe positions
            if probe:
                # Get positions
                forward_pos = primers[0].get("position", 0)
                forward_len = len(primers[0].get("sequence", ""))
                probe_pos = probe.get("position", 0)
                probe_len = len(probe.get("sequence", ""))
                reverse_pos = primers[1].get("position", 0)
                
                # Check proper positioning
                if not (forward_pos + forward_len <= probe_pos and probe_pos + probe_len <= reverse_pos):
                    # If positioning isn't correct, redesign probe
                    logger.warning("Oops, initial probe positioning wasn't appropriate. Redesigning...")
                    
                    # Calculate the middle position for the probe
                    middle_pos = forward_pos + forward_len + (reverse_pos - (forward_pos + forward_len)) // 2
                    
                    # Adjust amplicon_info to force probe position
                    amplicon_info["preferred_probe_pos"] = middle_pos
                    
                    # Redesign probe
                    probe = design_probe_for_amplicon(
                        amplicon_info,
                        min_distance_from_primers=5,
                        force_position=True
                    )
                    
            # Final verification
            is_valid_design = False
            if probe:
                forward_pos = primers[0].get("position", 0)
                forward_len = len(primers[0].get("sequence", ""))
                probe_pos = probe.get("position", 0)
                probe_len = len(probe.get("sequence", ""))
                reverse_pos = primers[1].get("position", 0)
                
                is_valid_design = (forward_pos + forward_len <= probe_pos and probe_pos + probe_len <= reverse_pos)
                
            # Calculate the final amplicon
            forward_pos = primers[0].get("position", 0)
            reverse_primer_len = len(primers[1].get("sequence", ""))
            reverse_pos = primers[1].get("position", 0)
            
            # The amplicon should include both primers
            amplicon_sequence = marker_sequence[forward_pos:reverse_pos + reverse_primer_len]
            
            return {
                "status": "success" if is_valid_design else "partial",
                "message": "Woo, designed primers and probe with appropriate positioning!" if is_valid_design else "Designed primers but probe positioning may  not be optimal. You should really double check that.",
                "primers": primers,
                "probe": probe,
                "amplicon": amplicon_sequence,
                "amplicon_length": len(amplicon_sequence)
            }
            
        # If all else fails, return failure status
        return {
            "status": "failed",
            "message": "What a bummer. Could not design suitable primers for the marker sequence",
            "primers": [],
            "probe": None,
            "amplicon": "",
            "amplicon_length": 0
        }
                    
def design_probe_for_amplicon(
    amplicon: Dict[str, Any],
    min_length: int = MIN_PROBE_LENGTH,
    max_length: int = MAX_PROBE_LENGTH,
    min_tm: float = MIN_PROBE_TM,
    max_tm: float = MAX_PROBE_TM,
    min_gc: float = MIN_GC_CONTENT,
    max_gc: float = MAX_GC_CONTENT,
    min_distance_from_primers: int = 5,
    force_position: bool = False
) -> Dict[str, Any]:
    """
    Design a probe for the amplicon region with improved positioning.
    
    Args:
        amplicon (Dict[str, Any]): Amplicon information
        min_length (int): Minimum probe length
        max_length (int): Maximum probe length
        min_tm (float): Minimum melting temperature
        max_tm (float): Maximum melting temperature
        min_gc (float): Minimum GC content percentage
        max_gc (float): Maximum GC content percentage
        min_distance_from_primers (int): Minimum distance from primers
        force_position (bool): Force probe to be in preferred position
        
    Returns:
        Dict[str, Any]: Designed probe information or None if criteria can't be met
    """
    sequence = amplicon.get("sequence", "")
    
    if not sequence or len(sequence) < min_length:
        logger.warning("Amplicon too short for probe design")
        return None
        
    # Adjust parameters for short amplicons
    if len(sequence) < max_length * 1.5:
        max_length = len(sequence) // 2
        min_length = max(min_length, max_length - 5)
        logger.info(f"Adjusted probe length range to {min_length}-{max_length} for short amplicon")
        
    # Calculate safe region for probe (to make sure it doesn't overlap with primers)
    start_pos = min_distance_from_primers
    end_pos = len(sequence) - min_distance_from_primers
    
    # Ensure we have a valid region
    if end_pos - start_pos < min_length:
        logger.warning(f"Amplicon too short for probe with {min_distance_from_primers}bp buffer from primers")
        # Reduce the buffer if needed
        min_distance_from_primers = max(1, min_distance_from_primers // 2)
        start_pos = min_distance_from_primers
        end_pos = len(sequence) - min_distance_from_primers
        
        if end_pos - start_pos < min_length:
            logger.warning("Amplicon too short for non-overlapping probe")
            return None
            
    # If forcing probe position, use the preferred position if available
    if force_position and "preferred_probe_pos" in amplicon:
        preferred_pos = amplicon["preferred_probe_pos"]
        # Adjust to be within safe region
        start_pos = max(start_pos, preferred_pos - max_length)
        end_pos = min(end_pos, preferred_pos + max_length)
    else:
        # Default: aim for middle of amplicon
        middle = len(sequence) // 2
        # Adjust to ensure we have room for probe
        start_pos = max(start_pos, middle - (max_length // 2))
        end_pos = min(end_pos, middle + (max_length // 2))
        
    probe_region = sequence[start_pos:end_pos]
    
    # Generate probe candidates
    probe_candidates = []
    
    for length in range(min_length, min(max_length + 1, len(probe_region) +1)):
        for start in range(0, len(probe_region) - length + 1):
            probe_seq = probe_region[start:start + length]
            properties = calculate_oligo_properties(probe_seq)
            
            # Check if probe meets criteria
            if (min_gc <= properties["gc_content"] <= max_gc and
                min_tm <= properties["tm"] <= max_tm and
                not properties["has_repeats"] and
                properties["self_complementarity"] <= 5):
                
                probe_candidates.append({
                    "sequence": probe_seq,
                    "properties": properties,
                    "position": start_pos + start,
                    "orientation": "forward",
                    "score": calculate_oligo_score(properties, "probe")
                })
                
    # If no candidates found, try with relaxed criteria
    if not probe_candidates:
        logger.warning("No probe candidates found with initial criteria, relaxing constraints")
        for length in range(min_length, min(max_length + 1, len(sequence) + 1)):
            for start in range(start_pos, min(end_pos, len(sequence) - length + 1)):
                probe_seq = sequence[start:start + length]
                properties = calculate_oligo_properties(probe_seq)
                
                # Relaxed criteria
                if (min_gc - 5 <= properties["gc_content"] <= max_gc + 5 and
                    min_tm - 3 <= properties["tm"] <= max_tm + 3):
                    
                    probe_candidates.append({
                        "sequence": probe_seq,
                        "properties": properties,
                        "position": start,
                        "orientation": "forward",
                        "score": calculate_oligo_score(properties, "probe") * 0.8 # Penalty for using relaxed criteria
                    })
                    
    # If still no probe candidates....create a basic probe
    if not probe_candidates:
        logger.warning("No suitable probes found, creating a basic probe. Proceed with caution.")
        middle = len(sequence) // 2
        probe_start = max(start_pos, middle - 10)
        probe_end = min(end_pos, probe_start + 20)
        probe_seq = sequence[probe_start:probe_end]
        
        return {
            "name": "Probe",
            "sequence": probe_seq,
            "properties": calculate_oligo_properties(probe_seq),
            "position": probe_start,
            "orientation": "forward",
            "score": 0.5
        }
        
    # Sort candidates by score
    probe_candidates.sort(key=lambda x: x["score"], reverse=True)
    
    # Select the best probe
    best_probe = probe_candidates[0]
    best_probe["name"] = "Probe"
    
    return best_probe
                 
def design_primers_for_sequence(
    sequence: str,
    min_length: int = MIN_PRIMER_LENGTH, 
    max_length: int = MAX_PRIMER_LENGTH,
    min_tm: float = MIN_PRIMER_TM,
    max_tm: float = MAX_PRIMER_TM,
    min_gc: float = MIN_GC_CONTENT,
    max_gc: float = MAX_GC_CONTENT,
    target_amplicon_size: int = 150
) -> Tuple[List[Dict[str, Any]], Dict[str, Any]]:
    """
    Design primers for a given sequence and return the resulting amplicon.
    
    Args:
        sequence (str): The marker sequence
        min_length (int): Minimum primer length
        max_length (int): Maximum primer length
        min_tm (float): Minimum melting temperature
        max_tm (float): Maximum melting temperature
        min_gc (float): Minimum GC content percentage
        max_gc (float): Maximum GC content percentage
        target_amplicon_size (int): Target amplicon size
        
    Returns:
        Tuple[List[Dict[str, Any]], Dict[str, Any]]: 
            1. List of designed primers
            2. Dictionary with amplicon information
    """
    logger.info(f"Designing primers for sequence of length {len(sequence)}")
    
    # Check if sequence is too short
    if len(sequence) < MIN_AMPLICON_LENGTH:
        logger.warning(f"Sequence too short ({len(sequence)} bp) for meaningful primer design")
        return [], {"sequence": "", "start": 0, "end": 0, "length": 0}
        
    # Generate forward primer candidates
    forward_candidates = []
    forward_region = sequence[:min(100, len(sequence) // 3)]
    
    for length in range(min_length, min(max_length + 1, len(forward_region) + 1)):
        for start in range(0, min(20, len(forward_region) - length + 1)):
            primer_seq = forward_region[start:start + length]
            properties = calculate_oligo_properties(primer_seq)
            
            # Check if primer meets criteria
            if (min_gc <= properties["gc_content"] <= max_gc and
                min_tm <= properties["tm"] <= max_tm and
                not properties["has_repeats"] and
                properties["self_complementarity"] <= 4):
                
                forward_candidates.append({
                    "sequence": primer_seq,
                    "properties": properties,
                    "position": start,
                    "orientation": "forward",
                    "score": calculate_oligo_score(properties, "primer")
                })
    
    # Generate reverse primer candidates
    reverse_candidates = []
    reverse_region = sequence[max(0, len(sequence) - 100):]
    reverse_comp = str(Seq(reverse_region).reverse_complement())
    
    for length in range(min_length, min(max_length + 1, len(reverse_region) + 1)):
        for start in range(0, min(20, len(reverse_region) - length + 1)):
            primer_seq = reverse_comp[start:start + length]
            properties = calculate_oligo_properties(primer_seq)
            
            # Check if primer meets criteria
            if (min_gc <= properties["gc_content"] <= max_gc and
                min_tm <= properties["tm"] <= max_tm and
                not properties["has_repeats"] and
                properties["self_complementarity"] <= 4):
                
                reverse_candidates.append({
                    "sequence": primer_seq,
                    "properties": properties,
                    "position": len(sequence) - (len(reverse_region) - start) + length,
                    "orientation": "reverse",
                    "score": calculate_oligo_score(properties, "primer")
                })
    
    # Sort candidates by score
    forward_candidates.sort(key=lambda x: x["score"], reverse=True)
    reverse_candidates.sort(key=lambda x: x["score"], reverse=True)
    
    logger.info(f"Generated {len(forward_candidates)} forward and {len(reverse_candidates)} reverse primer candidates")
    
    # If we couldn't find primers that meet all criteria, create basic ones
    if not forward_candidates:
        logger.warning("No good forward primers found, creating basic primer")
        primer_seq = sequence[:min(20, len(sequence) // 3)]
        forward_candidates = [{
            "sequence": primer_seq,
            "properties": calculate_oligo_properties(primer_seq),
            "position": 0,
            "orientation": "forward",
            "score": 0.5
        }]
    
    if not reverse_candidates:
        logger.warning("No good reverse primers found, creating basic primer")
        rev_seq = str(Seq(sequence[-min(20, len(sequence) // 3):]).reverse_complement())
        reverse_candidates = [{
            "sequence": rev_seq,
            "properties": calculate_oligo_properties(rev_seq),
            "position": len(sequence) - min(20, len(sequence) // 3),
            "orientation": "reverse",
            "score": 0.5
        }]
    
    # Find the best primer pair
    primer_pair = optimize_primer_pair(
        forward_candidates, 
        reverse_candidates, 
        target_amplicon_size
    )
    
    # Extract the amplicon region
    forward_primer = primer_pair["forward_primer"]
    reverse_primer = primer_pair["reverse_primer"]
    
    # Get the amplicon sequence
    amplicon_start = forward_primer["position"]
    amplicon_end = reverse_primer["position"] + 1
    
    if amplicon_start < amplicon_end and amplicon_end <= len(sequence):
        amplicon_sequence = sequence[amplicon_start:amplicon_end]
    else:
        # Fallback if positions are invalid
        logger.warning("Invalid amplicon coordinates, using placeholder")
        amplicon_sequence = sequence[:min(150, len(sequence))]
    
    # Create amplicon dictionary
    amplicon = {
        "sequence": amplicon_sequence,
        "start": amplicon_start,
        "end": amplicon_end,
        "length": len(amplicon_sequence)
    }
    
    # Format the final primers
    final_primers = []
    
    forward = forward_primer.copy()
    forward["name"] = "Forward_Primer"
    final_primers.append(forward)
    
    reverse = reverse_primer.copy()
    reverse["name"] = "Reverse_Primer"
    final_primers.append(reverse)
    
    return final_primers, amplicon

def optimize_primer_pair(
    forward_candidates: List[Dict[str, Any]],
    reverse_candidates: List[Dict[str, Any]],
    target_amplicon_size: int = 150
) -> Dict[str, Any]:
    """
    Find the optimal primer pair for amplification.
    
    Args:
        forward_candidates (List[Dict[str, Any]]): Forward primer candidates
        reverse_candidates (List[Dict[str, Any]]): Reverse primer candidates
        target_amplicon_size (int): Target amplicon size
        
    Returns:
        Dict[str, Any]: Optimal primer pair information
    """
    if not forward_candidates or not reverse_candidates:
        logger.warning("Not enough primers to form a pair")
        return {
            "forward_primer": forward_candidates[0] if forward_candidates else None,
            "reverse_primer": reverse_candidates[0] if reverse_candidates else None,
            "product_size": target_amplicon_size,  # Placeholder
            "compatibility_score": 0.5  # Placeholder
        }
    
    # Find the best pair
    best_pair = None
    best_score = -1
    
    for forward in forward_candidates[:5]:  # Limit to top 5 candidates for efficiency
        for reverse in reverse_candidates[:5]:
            # Calculate product size
            if forward.get("position") is not None and reverse.get("position") is not None:
                product_size = reverse["position"] - forward["position"]
                
                # Only consider pairs with reasonable product size
                if product_size <= 50 or product_size > 1000:
                    continue
                
                # Calculate compatibility score
                tm_diff = abs(forward["properties"]["tm"] - reverse["properties"]["tm"])
                tm_compatibility = max(0, 1.0 - (tm_diff / 5))
                
                # Adjust score based on product size proximity to target
                size_diff = abs(product_size - target_amplicon_size)
                size_compatibility = max(0, 1.0 - (size_diff / target_amplicon_size))
                
                # Calculate overall compatibility
                compatibility_score = 0.7 * tm_compatibility + 0.3 * size_compatibility
                
                # Calculate overall score considering individual primer scores
                overall_score = (
                    0.4 * forward["score"] +
                    0.4 * reverse["score"] +
                    0.2 * compatibility_score
                )
                
                if overall_score > best_score:
                    best_score = overall_score
                    best_pair = {
                        "forward_primer": forward,
                        "reverse_primer": reverse,
                        "product_size": product_size,
                        "compatibility_score": compatibility_score
                    }
    
    # If no good pair found, just take the best primers
    if not best_pair:
        logger.warning("Could not find optimal primer pair, using best individual primers")
        return {
            "forward_primer": forward_candidates[0],
            "reverse_primer": reverse_candidates[0],
            "product_size": target_amplicon_size,  # Placeholder
            "compatibility_score": 0.5  # Placeholder
        }
    
    return best_pair

def design_probe_for_amplicon(
    amplicon: Dict[str, Any],
    min_length: int = MIN_PROBE_LENGTH,
    max_length: int = MAX_PROBE_LENGTH,
    min_tm: float = MIN_PROBE_TM,
    max_tm: float = MAX_PROBE_TM,
    min_gc: float = MIN_GC_CONTENT,
    max_gc: float = MAX_GC_CONTENT
) -> Dict[str, Any]:
    """
    Design a probe for the amplicon region.
    
    Args:
        amplicon (Dict[str, Any]): Amplicon information
        min_length (int): Minimum probe length
        max_length (int): Maximum probe length
        min_tm (float): Minimum melting temperature
        max_tm (float): Maximum melting temperature
        min_gc (float): Minimum GC content percentage
        max_gc (float): Maximum GC content percentage
        
    Returns:
        Dict[str, Any]: Designed probe information or None if criteria can't be met
    """
    sequence = amplicon.get("sequence", "")
    
    if not sequence or len(sequence) < min_length:
        logger.warning("Amplicon too short for probe design")
        return None
    
    # Adjust parameters for short amplicons
    if len(sequence) < max_length * 1.5:
        max_length = len(sequence) // 2
        min_length = max(min_length, max_length - 5)
        logger.info(f"Adjusted probe length range to {min_length}-{max_length} for short amplicon")
    
    # Try to position the probe in the middle third of the amplicon
    start_pos = max(0, len(sequence) // 3)
    end_pos = min(len(sequence), 2 * len(sequence) // 3)
    
    # Adjust if amplicon is very short
    if end_pos - start_pos < min_length:
        start_pos = 0
        end_pos = len(sequence)
    
    probe_region = sequence[start_pos:end_pos]
    
    # Generate probe candidates
    probe_candidates = []
    
    for length in range(min_length, min(max_length + 1, len(probe_region) + 1)):
        for start in range(0, len(probe_region) - length + 1):
            probe_seq = probe_region[start:start + length]
            properties = calculate_oligo_properties(probe_seq)
            
            # Check if probe meets criteria
            if (min_gc <= properties["gc_content"] <= max_gc and
                min_tm <= properties["tm"] <= max_tm and
                not properties["has_repeats"] and
                properties["self_complementarity"] <= 5):
                
                probe_candidates.append({
                    "sequence": probe_seq,
                    "properties": properties,
                    "position": start_pos + start,
                    "orientation": "forward",
                    "score": calculate_oligo_score(properties, "probe")
                })
    
    # If no candidates found, try with relaxed criteria
    if not probe_candidates:
        logger.warning("No probe candidates found with initial criteria, relaxing constraints")
        for length in range(min_length, min(max_length + 1, len(sequence) + 1)):
            for start in range(0, len(sequence) - length + 1):
                probe_seq = sequence[start:start + length]
                properties = calculate_oligo_properties(probe_seq)
                
                # Relaxed criteria
                if (min_gc - 5 <= properties["gc_content"] <= max_gc + 5 and
                    min_tm - 3 <= properties["tm"] <= max_tm + 3):
                    
                    probe_candidates.append({
                        "sequence": probe_seq,
                        "properties": properties,
                        "position": start,
                        "orientation": "forward",
                        "score": calculate_oligo_score(properties, "probe") * 0.8  # Penalty for relaxed criteria
                    })
    
    # If still no candidates, create a basic probe
    if not probe_candidates:
        logger.warning("No suitable probe found, creating basic probe")
        middle = len(sequence) // 2
        probe_start = max(0, middle - 10)
        probe_end = min(len(sequence), probe_start + 20)
        probe_seq = sequence[probe_start:probe_end]
        
        return {
            "name": "Probe",
            "sequence": probe_seq,
            "properties": calculate_oligo_properties(probe_seq),
            "position": probe_start,
            "orientation": "forward",
            "score": 0.5
        }
    
    # Sort candidates by score
    probe_candidates.sort(key=lambda x: x["score"], reverse=True)
    
    # Select the best probe
    best_probe = probe_candidates[0]
    best_probe["name"] = "Probe"
    
    return best_probe
    
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
        # Check if assay_info is not None before checking its status
        if assay_info is not None and assay_info.get("status") == "success":
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