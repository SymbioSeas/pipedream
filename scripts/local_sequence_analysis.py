import os
from Bio import SeqIO, AlignIO
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter
import tempfile
import subprocess
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import time

def analyze_fasta(fasta_file):
    """Analyze a FASTA file and return basic statistics"""
    print(f"Analyzing {fasta_file}...")
    
    # Load sequences
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    seq_count = len(sequences)
    
    if seq_count == 0:
        print("No sequences found in the file.")
        return
    
    # Address duplicate IDs by adding a unique suffix
    seen_ids = {}
    for i, seq in enumerate(sequences):
        if seq.id in seen_ids:
            seen_ids[seq.id] += 1
            seq.id = f"{seq.id}_{seen_ids[seq.id]}"
        else:
            seen_ids[seq.id] = 0
    
    # Calculate length statistics
    lengths = [len(seq) for seq in sequences]
    min_len = min(lengths)
    max_len = max(lengths)
    avg_len = sum(lengths) / len(lengths)
    
    print(f"Found {seq_count} sequences")
    print(f"Length range: {min_len} - {max_len} bp (average: {avg_len:.1f} bp)")
    
    # Plot length distribution
    plt.figure(figsize=(10, 6))
    plt.hist(lengths, bins=20)
    plt.title('Sequence Length Distribution')
    plt.xlabel('Length (bp)')
    plt.ylabel('Count')
    plt.savefig('sequence_length_distribution.png')
    print("Saved sequence length distribution plot to 'sequence_length_distribution.png'")
    
    # Are the sequences similar in length? (important for alignment)
    length_variation = np.std(lengths) / avg_len
    print(f"Length variation coefficient: {length_variation:.3f}")
    
    # Take a sample of sequences for alignment if there are too many
    max_for_alignment = 100  # Limit for alignment visualization
    if seq_count > max_for_alignment:
        print(f"Sampling {max_for_alignment} sequences for alignment visualization")
        import random
        random.seed(42)  # For reproducibility
        sample_seqs = random.sample(sequences, max_for_alignment)
    else:
        sample_seqs = sequences
    
    # Create a temporary file for alignment
    with tempfile.NamedTemporaryFile(suffix='.fasta', delete=False) as temp:
        sample_file = temp.name
        SeqIO.write(sample_seqs, sample_file, "fasta")
    
    # Align sequences with MUSCLE 5.2
    try:
        print("Aligning sequences with MUSCLE...")
        output_file = "full_alignment.fasta"
        
        # Updated command for MUSCLE 5.2
        muscle_cmd = f"muscle -align {sample_file} -output {output_file}"
        subprocess.run(muscle_cmd, shell=True, check=True)
        
        # Read the alignment
        alignment = AlignIO.read(output_file, "fasta")
        
        # Create a phylogenetic tree from the full alignment
        print("Creating phylogenetic tree from full alignment...")
        create_phylogenetic_tree(alignment, "full_alignment_tree.png", "Full Sequence Phylogenetic Tree")
        
        # Calculate conservation per position
        alignment_length = alignment.get_alignment_length()
        conservation = []
        
        for i in range(alignment_length):
            column = alignment[:, i]
            # Count non-gap characters
            non_gaps = [c for c in column if c != '-']
            if len(non_gaps) == 0:
                conservation.append(0)
                continue
                
            # Count occurrences of each base
            base_counts = Counter(non_gaps)
            most_common = base_counts.most_common(1)[0][1]
            cons_score = most_common / len(non_gaps)
            conservation.append(cons_score)
        
        # Plot conservation
        plt.figure(figsize=(12, 6))
        plt.plot(range(alignment_length), conservation)
        plt.axhline(y=0.9, color='r', linestyle='--')  # 90% conservation threshold
        plt.title('Sequence Conservation by Position')
        plt.xlabel('Position')
        plt.ylabel('Conservation Score')
        plt.savefig('sequence_conservation.png')
        print("Saved conservation plot to 'sequence_conservation.png'")
        
        # Find highly conserved regions (potential primer/probe sites)
        window_size = 18  # Typical primer length
        region_scores = []
        
        for i in range(alignment_length - window_size + 1):
            window_score = sum(conservation[i:i+window_size]) / window_size
            region_scores.append((i, window_score))
        
        # Sort regions by conservation score
        region_scores.sort(key=lambda x: x[1], reverse=True)
        
        # Print top conserved regions
        print("\nTop conserved regions (potential primer sites):")
        for i, (pos, score) in enumerate(region_scores[:5]):
            if score < 0.8:  # Skip low conservation regions
                continue
            print(f"Region {i+1}: Position {pos}-{pos+window_size}, Score: {score:.3f}")
            # Extract consensus sequence for this region
            consensus = ""
            for j in range(pos, pos+window_size):
                column = alignment[:, j]
                # Remove gaps
                bases = [b for b in column if b != '-']
                if not bases:
                    consensus += 'N'
                    continue
                # Get most common base
                most_common = Counter(bases).most_common(1)[0][0]
                consensus += most_common
            print(f"Consensus: {consensus}")
        
        # Find all high conserved regions
        min_conservation = 0.9
        max_amplicon = 300
        high_conserved_regions = [(pos, score) for pos, score in region_scores if score >= min_conservation]
        
        # Start timer for optimized search
        start_time = time.time()
        print(f"\nSearching for primer-probe combinations using optimized approach...")
        
        # Use optimized search for assay combinations
        assay_candidates = find_assay_combinations(
            high_conserved_regions, 
            alignment, 
            window_size,
        )
        
        end_time = time.time()
        print(f"Assay search completed in {end_time - start_time:.2f} seconds")
        
        if assay_candidates:
            # Print top assay candidates
            print(f"\nFound {len(assay_candidates)} potential assays. Here are the top candidates:")
            
            for i, assay in enumerate(assay_candidates[:5]):  # Show top 5
                # Extract consensus sequences
                fwd_consensus = "".join([consensus_base(alignment, p) for p in range(assay['forward']['position'], 
                                                                                    assay['forward']['position']+window_size)])
                probe_consensus = "".join([consensus_base(alignment, p) for p in range(assay['probe']['position'], 
                                                                                      assay['probe']['position']+window_size)])
                rev_consensus = "".join([consensus_base(alignment, p) for p in range(assay['reverse']['position'], 
                                                                                    assay['reverse']['position']+window_size)])
                rev_consensus_rc = str(Seq(rev_consensus).reverse_complement())
                
                # Calculate GC content
                fwd_gc = calculate_gc_content(fwd_consensus)
                probe_gc = calculate_gc_content(probe_consensus)
                rev_gc = calculate_gc_content(rev_consensus_rc)
                
                print(f"\nAssay {i+1} (Quality Score: {assay['quality_score']:.3f}):")
                print(f"  Forward: Pos {assay['forward']['position']}-{assay['forward']['position']+window_size}, " +
                      f"Score: {assay['forward']['score']:.3f}, GC: {fwd_gc:.1f}%")
                print(f"  Primer: {fwd_consensus}")
                
                print(f"  Probe: Pos {assay['probe']['position']}-{assay['probe']['position']+window_size}, " +
                      f"Score: {assay['probe']['score']:.3f}, GC: {probe_gc:.1f}%")
                print(f"  Probe: {probe_consensus}")
                print(f"  Spacing: {assay['probe']['spacing_to_fwd']} bp after forward, " +
                      f"{assay['probe']['spacing_to_rev']} bp before reverse")
                
                print(f"  Reverse: Pos {assay['reverse']['position']}-{assay['reverse']['position']+window_size}, " +
                      f"Score: {assay['reverse']['score']:.3f}, GC: {rev_gc:.1f}%")
                print(f"  Primer (RC): {rev_consensus_rc}")
                
                print(f"  Amplicon size: {assay['amplicon_size']} bp")
                
                # Store sequence information in the assay
                assay['forward']['consensus'] = fwd_consensus
                assay['forward']['gc'] = fwd_gc
                assay['probe']['consensus'] = probe_consensus
                assay['probe']['gc'] = probe_gc
                assay['reverse']['consensus'] = rev_consensus
                assay['reverse']['consensus_rc'] = rev_consensus_rc
                assay['reverse']['gc'] = rev_gc
                
                # Extract the in silico amplicon for this assay
                assay['amplicon_alignment'] = alignment[:, assay['forward']['position']:assay['reverse']['position']+window_size]
            
            # Best assay is the first one (they're sorted by quality score)
            best_assay = assay_candidates[0]
            print("\nBest assay selected:")
            print(f"  Forward: {best_assay['forward']['consensus']}")
            print(f"  Probe: {best_assay['probe']['consensus']}")
            print(f"  Reverse (RC): {best_assay['reverse']['consensus_rc']}")
            print(f"  Amplicon size: {best_assay['amplicon_size']} bp")
            print(f"  Quality score: {best_assay['quality_score']:.3f}")
            
            # Save the amplicon alignment to a file
            amplicon_file = "amplicon_alignment.fasta"
            with open(amplicon_file, "w") as f:
                for i, seq in enumerate(alignment):
                    name = seq.id
                    # Get the sequence region and ensure it's handled as a string
                    start = best_assay['forward']['position']
                    end = best_assay['reverse']['position'] + window_size
                    amplicon_seq = str(seq.seq[start:end])
                    
                    # Make sure to remove any gaps that might cause length differences
                    amplicon_seq = amplicon_seq.replace('-', '')
                    
                    # Make sure the sequence isn't empty
                    if len(amplicon_seq) > 0:
                        f.write(f">{name}\n{amplicon_seq}\n")
                    else:
                        print(f"Warning: Empty amplicon sequences for {name}")
            
            # Then create a new alignment from these ungapped sequences
            print("Creating alignment from in silico amplicon sequences...")
            amplicon_output = "amplicon_alignment_processed.fasta"
            muscle_cmd = f"muscle -align {amplicon_file} -output {amplicon_output}"
            subprocess.run(muscle_cmd, shell=True, check=True)
            
            # Use this new alignment for the tree
            print("Creating phylogenetic tree from in silico amplicon...")
            amplicon_alignment = AlignIO.read(amplicon_output, "fasta")
            create_phylogenetic_tree(amplicon_alignment, "amplicon_tree.png", "In Silico Amplicon Phylogenetic Tree")
                                    
            # Generate assay design parameters for assay-design command
            print("\nRecommended assay-design command parameters:")
            print(f"--max-amplicon-length {max(best_assay['amplicon_size'] + 50, 200)} " +
                  f"--min-conservation {min(0.9, min(best_assay['forward']['score'], best_assay['probe']['score'], best_assay['reverse']['score']) - 0.05)} " +
                  "--conserved-mode")
        else:
            print("No suitable assay combinations found.")
            print("Try adjusting the parameters:")
            print("  1. Increase maximum amplicon length")
            print("  2. Decrease minimum conservation threshold")
            print("  3. Try a different primer/probe length")
        
        # Clean up
        os.remove(sample_file)
        
        return alignment, conservation, region_scores
    
    except Exception as e:
        print(f"Error during alignment and analysis: {e}")
        import traceback
        traceback.print_exc()
        # Clean up
        try:
            os.remove(sample_file)
        except:
            pass
        return None, None, None

def find_assay_combinations(high_conserved_regions, alignment, window_size):
    """Optimized search for primer-probe combinations"""
    # Configuration
    min_amplicon = 70
    max_amplicon = 300
    max_assays_to_find = 10
    alignment_length = alignment.get_alignment_length()
    
    # 1. Pre-filtering candidate regions by position
    forward_candidates = [(pos, score) for pos, score in high_conserved_regions 
                         if pos < alignment_length * 0.7]  # First 70% of sequence
    
    reverse_candidates = [(pos, score) for pos, score in high_conserved_regions 
                         if pos > alignment_length * 0.3]  # Last 70% of sequence
                         
    # Check if we have enough candidates with position filtering
    if len(forward_candidates) < 3 or len(reverse_candidates) < 3:
        print("Not enough position-based candidates, using all conserved regions")
        # Fall back: use all conserved regions and rely on amplicon size constraints
        all_regions = high_conserved_regions.copy()
        # Sort by position for better organization
        all_regions.sort(key=lambda x: x[0])
        
        # For forward, prefer earlier positions
        forward_candidates = all_regions.copy()
        # For reverse, prefer later positions
        reverse_candidates = list(reversed(all_regions.copy()))
    
    # 2. Sort by conservation score
    forward_candidates.sort(key=lambda x: x[1], reverse=True)
    reverse_candidates.sort(key=lambda x: x[1], reverse=True)
    
    # 3. Take top candidates for an initial pool
    forward_candidates = forward_candidates[:min(20, len(forward_candidates))]
    reverse_candidates = reverse_candidates[:min(20, len(reverse_candidates))]
    
    print(f"Evaluating {len(forward_candidates)} forward and {len(reverse_candidates)} reverse candidates")
    
    # 4. Find promising primer pairs
    promising_pairs = []
    for fwd_pos, fwd_score in forward_candidates:
        for rev_pos, rev_score in reverse_candidates:
            if fwd_pos >= rev_pos:
                continue
                
            amplicon_size = rev_pos + window_size - fwd_pos
            if amplicon_size < min_amplicon or amplicon_size > max_amplicon:
                continue
            
            # Calculate preliminary pair score
            size_score = 1.0 - abs(amplicon_size - 150) / 250
            pair_score = (fwd_score + rev_score) / 2 * 0.7 + size_score * 0.3
            
            promising_pairs.append({
                "fwd_pos": fwd_pos, 
                "fwd_score": fwd_score,
                "rev_pos": rev_pos, 
                "rev_score": rev_score,
                "amplicon_size": amplicon_size,
                "pair_score": pair_score
            })
    
    # 5. Sort promising pairs and take top candidates
    promising_pairs.sort(key=lambda x: x["pair_score"], reverse=True)
    promising_pairs = promising_pairs[:min(50, len(promising_pairs))]
    
    print(f"Found {len(promising_pairs)} promising primer pairs")


    if len(promising_pairs) == 0:
        # Diagnose why no promising pairs were found
        # Check amplicon size constraint failures
        amplicon_size_failures = 0
        orientation_failures = 0
        
        for fwd_pos, fwd_score in forward_candidates:
            for rev_pos, rev_score in reverse_candidates:
                if fwd_pos >= rev_pos:
                    orientation_failures += 1
                    continue
                    
                amplicon_size = rev_pos + window_size - fwd_pos
                if amplicon_size < min_amplicon or amplicon_size > max_amplicon:
                    amplicon_size_failures += 1
        
        print(f"Diagnostic information:")
        print(f"  - Orientation failures: {orientation_failures} (forward position >= reverse position)")
        print(f"  - Amplicon size failures: {amplicon_size_failures} (outside {min_amplicon}-{max_amplicon} bp range)")
        
        # Check conservation distribution
        conservation_counts = {}
        for threshold in [0.9, 0.8, 0.7, 0.6, 0.5]:
            count = sum(1 for _, score in high_conserved_regions if score >= threshold)
            conservation_counts[threshold] = count
        
        print("Conservation distribution:")
        for threshold, count in conservation_counts.items():
            print(f"  - {threshold*100:.0f}% conservation: {count} regions")
        
        # Show position distribution of conserved regions (at lower threshold)
        lower_threshold = 0.7  # Try 70% conservation
        lower_conserved = [(pos, score) for pos, score in region_scores if score >= lower_threshold]
        
        # Group by position in 10% increments of sequence length
        position_bins = [0] * 10
        for pos, _ in lower_conserved:
            bin_idx = min(9, int(pos * 10 / alignment_length))
            position_bins[bin_idx] += 1
        
        print("Position distribution of regions with ≥70% conservation:")
        for i, count in enumerate(position_bins):
            print(f"  - {i*10}-{(i+1)*10}% of sequence: {count} regions")
    
    # 6. For each promising pair, find the best probe
    assay_candidates = []
    
    for pair in promising_pairs:
        fwd_pos = pair["fwd_pos"]
        rev_pos = pair["rev_pos"]
        
        # Find potential probe regions between primers
        probe_candidates = []
        for probe_pos, probe_score in high_conserved_regions:
            if probe_pos <= fwd_pos or probe_pos + window_size >= rev_pos:
                continue
                
            # Calculate spacing
            fwd_to_probe = probe_pos - (fwd_pos + window_size)
            probe_to_rev = rev_pos - (probe_pos + window_size)
            
            if fwd_to_probe < 5 or probe_to_rev < 5:
                continue  # Too close to primers
                
            # Calculate spacing score
            spacing_score = min(1.0, fwd_to_probe / 20) * min(1.0, probe_to_rev / 20)
            
            # Combined probe score
            combined_score = probe_score * 0.7 + spacing_score * 0.3
            
            probe_candidates.append({
                "position": probe_pos,
                "score": probe_score,
                "spacing_to_fwd": fwd_to_probe,
                "spacing_to_rev": probe_to_rev,
                "combined_score": combined_score
            })
        
        if probe_candidates:
            # Take the best probe
            probe_candidates.sort(key=lambda x: x["combined_score"], reverse=True)
            best_probe = probe_candidates[0]
            
            # Calculate full assay quality score
            quality_score = (
                (pair["fwd_score"] + pair["rev_score"] + best_probe["score"]) / 3 * 0.6 +
                (1.0 - abs(pair["amplicon_size"] - 150) / 250) * 0.2 +
                (min(1.0, best_probe["spacing_to_fwd"] / 15) * 
                 min(1.0, best_probe["spacing_to_rev"] / 15)) * 0.2
            )
            
            # Add to assay candidates
            assay_candidates.append({
                "forward": {
                    "position": fwd_pos,
                    "score": pair["fwd_score"]
                },
                "probe": best_probe,
                "reverse": {
                    "position": rev_pos,
                    "score": pair["rev_score"]
                },
                "amplicon_size": pair["amplicon_size"],
                "quality_score": quality_score
            })
            
            # Early termination if we have enough good candidates
            if len(assay_candidates) >= max_assays_to_find:
                break
    
    # Sort final candidates by quality score
    assay_candidates.sort(key=lambda x: x["quality_score"], reverse=True)
    
    return assay_candidates

def consensus_base(alignment, position):
    """Get the consensus base at a given position in an alignment"""
    column = alignment[:, position]
    # Remove gaps
    bases = [b for b in column if b != '-']
    if not bases:
        return 'N'
    # Get most common base
    most_common = Counter(bases).most_common(1)[0][0]
    return most_common

def calculate_gc_content(sequence):
    """Calculate GC content as a percentage"""
    if not sequence:
        return 0
    gc_count = sequence.count('G') + sequence.count('C')
    return (gc_count / len(sequence)) * 100

def create_phylogenetic_tree(alignment, output_file, title):
    """Create and save a phylogenetic tree from an alignment"""
    # Calculate distance matrix
    calculator = DistanceCalculator('identity')
    dm = calculator.get_distance(alignment)
    
    # Build tree
    constructor = DistanceTreeConstructor()
    tree = constructor.upgma(dm)
    
    # Remove inner node labels
    for clade in tree.get_nonterminals():
        clade.name = ""
    
    # Draw and save tree
    plt.figure(figsize=(14, 16))
    Phylo.draw(tree, do_show=False, label_func=lambda x: x.name)
    plt.title(title)
    plt.savefig(output_file, dpi=300, bbox_inches='tight')
    plt.close()
    
    # Also save tree in Newick format
    newick_file = output_file.replace('.png', '.newick')
    Phylo.write(tree, newick_file, 'newick')
    
    print(f"Phylogenetic tree saved to {output_file} and {newick_file}")
    
    return tree

def save_assay_results(assay, filename="assay_design_results.txt"):
    """Save assay design results to a text file"""
    with open(filename, "w") as f:
        f.write("==== Assay Design Results ====\n\n")
        f.write(f"Forward Primer: {assay['forward']['consensus']}\n")
        f.write(f"Probe: {assay['probe']['consensus']}\n")
        f.write(f"Reverse Primer (RC): {assay['reverse']['consensus_rc']}\n")
        f.write(f"Amplicon Size: {assay['amplicon_size']} bp\n\n")
        
        f.write("Primer Details:\n")
        f.write(f"Forward: Position {assay['forward']['position']}, Conservation {assay['forward']['score']:.3f}, GC {assay['forward']['gc']:.1f}%\n")
        f.write(f"Probe: Position {assay['probe']['position']}, Conservation {assay['probe']['score']:.3f}, GC {assay['probe']['gc']:.1f}%\n")
        f.write(f"Reverse: Position {assay['reverse']['position']}, Conservation {assay['reverse']['score']:.3f}, GC {assay['reverse']['gc']:.1f}%\n\n")
        
        f.write(f"Spacing: {assay['probe']['spacing_to_fwd']} bp between Forward and Probe\n")
        f.write(f"         {assay['probe']['spacing_to_rev']} bp between Probe and Reverse\n\n")
        
        f.write(f"Overall Quality Score: {assay['quality_score']:.3f}\n")
    
    print(f"Assay design results saved to {filename}")

# If run as a script, analyze the provided file
if __name__ == "__main__":
    import sys
    
    if len(sys.argv) > 1:
        fasta_file = sys.argv[1]
    else:
        # Default file path if none provided
        fasta_file = input("Enter path to FASTA file: ")
    
    # Run the analysis
    alignment, conservation, region_scores = analyze_fasta(fasta_file)