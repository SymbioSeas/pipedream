# assay_design/cli.py

import argparse
import logging
import sys
import time
import os
import json
import datetime
from typing import List, Optional, Dict, Any

from .data_retrieval import (
    fetch_sequences_for_taxid,
    fetch_gene_sequences,
    get_related_taxa,
    suggest_marker_genes,
    get_taxon_info,
    fetch_cds_sequences,
    load_local_fasta,
    intelligent_exclusion_selection
)
#Import hierarchical search function
from .hierarchical_search import hierarchical_marker_search
from .target_identification import (
    find_optimal_marker,
    find_conserved_without_exclusion,
    find_discriminative_kmers,
    extend_kmer
)
# Import the primer and probe design
from .primer_probe_design import (
    design_primers_and_probe, 
    design_multi_region_assay,
    find_conserved_kmers
)
from .specificity_validation import validate_primer_specificity
from .kmer_analysis import kmer_analysis
# Import visualization functions
from .visualization import create_visualizations

# Set up logging
logging.basicConfig(
    level=logging.INFO,
    format='%(asctime)s - %(name)s - %(levelname)s - %(message)s'
)
logger = logging.getLogger(__name__)

def create_parser():
    """Create an argument parser for the CLI."""
    parser = argparse.ArgumentParser(
        description='Design specific molecular assays for target organisms',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Design an assay with automatic gene and exclusion selection (recommended)
  assay-design --inclusion 689 --email user@example.com --auto-gene --exclusion-strategy intelligent --output-dir ./vibrio_auto

  # Design a quantification assay (single-copy gene, multi-level exclusion)
  assay-design --inclusion 689 --email user@example.com --auto-gene --gene-use-case quantification --exclusion-strategy intelligent --max-exclusion-taxa 15 --output-dir ./vibrio_quant

  # Design an assay for Vibrio mediterranei using automatic exclusion taxa
  assay-design --inclusion 689 --email user@example.com --output-dir ./vibrio_assay

  # Design an assay for a specific gene in Vibrio mediterranei
  assay-design --inclusion 689 --email user@example.com --gene "16S ribosomal RNA" --output-dir ./vibrio_16S_assay

  # Get suggestions for marker genes for a taxon
  assay-design --inclusion 689 --email user@example.com --suggest-genes --output-dir ./vibrio_marker_assay

  # Design with manual exclusion taxa
  assay-design --inclusion 689 --exclusion 717,670,672 --email user@example.com --output-dir ./vibrio_specific_assay

  # Design with intelligent multi-level exclusion
  assay-design --inclusion 689 --email user@example.com --exclusion-strategy intelligent --max-exclusion-taxa 10 --output-dir ./vibrio_intelligent

  # Design an assay using local FASTA files
  assay-design --inclusion-fasta nifH_refseq_100.fna --email user@example.com --output-dir ./nifH_assay

  # Design an assay using local FASTA files for both inclusion and exclusion
  assay-design --inclusion-fasta nifH_refseq_100.fna --exclusion-fasta non_target_sequences.fasta --email user@example.com --output-dir ./custom_assay
"""
    )
    
    # Required arguments
    required = parser.add_argument_group('required arguments')
    required.add_argument(
        '--email', 
        type=str, 
        required=True,
        help='Your email for NCBI queries (required by NCBI)'
    )
    required.add_argument(
        '--output-dir',
        type=str,
        required=True,
        help='Output directory for assay design results'
    )
    
    # Make inclusion group conditionally required (required if --cds isn't specified)
    inclusion_group = parser.add_argument_group('target specification (one is required)')
    inclusion_group = inclusion_group.add_mutually_exclusive_group(required=True)
    inclusion_group.add_argument(
        '--inclusion',
        type=str,
        help='NCBI taxid for the target organism'
    )
    inclusion_group.add_argument(
        '--cds',
        type=str,
        help='Fetch CDS sequences for the specified gene (e.g., "amoA") and use for assay design'
    )
    inclusion_group.add_argument(
        '--inclusion-fasta',
        type=str,
        help='Path to a local multi-FASTA file to use for inclusion sequences'
    )
    
    
    # Optional arguments
    parser.add_argument(
        '--conserved-mode',
        action='store_true',
        help='Optimize for conservation across diverse taxa without exclusions (recommended with --cds)')
    parser.add_argument(
        '--exclusion', 
        type=str,
        help='NCBI taxids for non-target organisms (comma-separated)'
    )
    parser.add_argument(
        '--gene',
        type=str,
        help='Specific gene target (e.g., "16S ribosomal RNA")'
    )
    parser.add_argument(
        '--auto-gene',
        action='store_true',
        help='Automatically select optimal gene for assay design (overrides --gene)'
    )
    parser.add_argument(
        '--gene-use-case',
        type=str,
        choices=['quantification', 'phylogeny', 'detection'],
        default='quantification',
        help='Use case for gene selection: quantification (single-copy, HGT-resistant), phylogeny (conserved, informative), or detection (highly specific) (default: quantification)'
    )
    parser.add_argument(
        '--max-genes-to-try',
        type=int,
        default=3,
        help='Maximum number of genes to try in fallback strategy (default: 3)'
    )
    parser.add_argument(
        '--min-gene-score',
        type=float,
        default=0.4,
        help='Minimum acceptable gene suitability score 0-1 (default: 0.4)'
    )
    parser.add_argument(
        '--auto-exclusion',
        action='store_true',
        help='Automatically select related taxa for exclusion'
    )
    parser.add_argument(
        '--exclusion-strategy',
        type=str,
        choices=['intelligent', 'manual', 'siblings', 'family', 'order'],
        default='siblings',
        help='Strategy for exclusion taxa selection: intelligent (multi-level phylogenetic), manual (user-specified only), siblings (genus-level), family, or order (default: siblings)'
    )
    parser.add_argument(
        '--max-exclusion-taxa',
        type=int,
        default=10,
        help='Maximum number of exclusion taxa to select (default: 10, only used with intelligent strategy)'
    )
    parser.add_argument(
        '--auto-exclusion-rank', 
        type=str,
        choices=['species', 'genus', 'family', 'order'],
        default='genus',
        help='Taxonomic rank for auto-exclusion (default: genus)'
    )
    parser.add_argument(
        '--auto-exclusion-count', 
        type=int,
        default=3,
        help='Number of taxa to auto-select for exclusion (default: 3)'
    )
    parser.add_argument(
        '--suggest-genes', 
        action='store_true',
        help='Suggest appropriate marker genes for the target'
    )
    parser.add_argument(
        '--max-runtime', 
        type=int,
        default=120,
        help='Maximum runtime in seconds for marker identification (default: 120)'
    )
    parser.add_argument(
        '--max-amplicon-length',
        type=int,
        default=200,
        help='Maximum length of the amplicon (default: 200 bp)'
    )
    parser.add_argument(
        '--max-seq-length', 
        type=int,
        default=100000,
        help='Maximum sequence length to analyze (default: 100,000 bp)'
    )
    parser.add_argument(
        '--max-seq-count', 
        type=int,
        default=20,
        help='Maximum number of sequences to analyze per group (default: 5)'
    )
    parser.add_argument(
        '--verbose', 
        action='store_true',
        help='Enable verbose output'
    )
    parser.add_argument(
        '--no-visualization',
        action='store_true',
        help='Disable visualization generation'
    )
    parser.add_argument(
        '--processing-strategy',
        type=str,
        choices=['window', 'adaptive', 'gene-specific', 'progressive'],
        default='window',
        help='Strategy for processing large sequences'
    )
    parser.add_argument(
        '--compute-intensity',
        type=str,
        choices=['low', 'medium', 'high'],
        default='medium',
        help='Compute intensity level (affects search thoroughness)'
    )
    parser.add_argument(
        '--search-strategy',
        type=str,
        choices=['basic', 'hierarchical'],
        default='hierarchical',
        help='Strategy for marker identification (default: hierarchical)'
    )
    parser.add_argument(
        '--min-conservation',
        type=float,
        default=0.9,
        help='Minimum conservation within inclusion sequences (default: 0.9)'
    )
    parser.add_argument(
        '--min-specificity',
        type=float,
        default=0.8,
        help='Minimum specificity compared to exclusion sequences (default: 0.8)'
    )
    parser.add_argument(
        '--single-region-design',
        action='store_true',
        help='Use traditional single-region design instead of multi-region design'
    )
    parser.add_argument(
        '--no-lsh',
        action='store_true',
        help='Disable Locality-Sensitive Hashing (LSH) optimization (using traditional k-mer comparison instead)'
    )
    parser.add_argument(
        '--lsh-timeout',
        type=int,
        default=60,
        help='Timeout in seconds for LSH marker identification (default: 60, only relevant when LSH is enabled)'
    )
    parser.add_argument(
        '--whole-amplicon-specificity',
        action='store_true',
        help='Requires the entire amplicon to be specific (not just the primer and probe binding sites)'
    )
    parser.add_argument(
        '--std-kmer',
        action='store_true',
        help='Uses standard k-mer analysis instead of optimized approach (less effecient but requires fewer dependencies)'
    )
    parser.add_argument(
        '--parallel',
        type=int,
        default=8,
        help='Number of parallel processes to use for k-mer analysis (default: CPU count)'
    )
    parser.add_argument(
        '--max-seqs',
        type=int,
        default=100,
        help='Maximum number of sequences to retrieve (default: 100)'
    )
    parser.add_argument(
        '--output-fasta',
        type=str,
        help='Output FASTA file path (default: {cds}_refseq_{max_seqs}.fna)'
    )
    parser.add_argument(
        '--exclusion-fasta',
        type=str,
        help='Path to a local multi-FASTA file to use for exclusion sequences'
    )
    
    return parser

def create_output_directory(output_dir: str) -> str:
    """
    Create output directory structure.
    
    Args:
        output_dir (str): Base output directory
        
    Returns:
        str: Path to created directory
    """
    # Create main directory
    try:
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)
            logger.info(f"Created output directory: {output_dir}")
        
        # Create subdirectories
        subdirs = ["sequences", "results", "visualizations", "logs"]
        for subdir in subdirs:
            subdir_path = os.path.join(output_dir, subdir)
            if not os.path.exists(subdir_path):
                os.makedirs(subdir_path)
                logger.info(f"Created subdirectory: {subdir_path}")
        
        return output_dir
    
    except Exception as e:
        logger.error(f"Error creating output directory: {e}")
        sys.exit(1)

def save_sequences(output_dir: str, inclusion_sequences: List, exclusion_sequences: List):
    """
    Save inclusion and exclusion sequences to files.
    
    Args:
        output_dir (str): Output directory
        inclusion_sequences (List): Inclusion sequences
        exclusion_sequences (List): Exclusion sequences
    """
    from Bio import SeqIO
    from Bio.Seq import UndefinedSequenceError
    
    seq_dir = os.path.join(output_dir, "sequences")
    
    # Save inclusion sequences
    if inclusion_sequences:
        # Filter out sequences with undefined content
        valid_inclusion = []
        for seq in inclusion_sequences:
            try:
                # Test if sequence is defined
                _ = str(seq.seq)
                valid_inclusion.append(seq)
            except UndefinedSequenceError:
                logger.warning(f"Skipping sequence {seq.id} with undefined content")
                
        if valid_inclusion:        
            inclusion_path = os.path.join(seq_dir, "inclusion_sequences.fasta")
            SeqIO.write(inclusion_sequences, inclusion_path, "fasta")
            logger.info(f"Saved {len(inclusion_sequences)} inclusion sequences to {inclusion_path}")
        else:
            logger.warning("No valid inclusion sequences to save")
    
    # Save exclusion sequences
    if exclusion_sequences:
        valid_exclusion = []
        for seq in exclusion_sequences:
            try:
                _ = str(seq.seq)
                valid_exclusion.append(seq)
            except UndefinedSequenceError:
                logger.warning(f"Skipping sequence {seq.id} with undefined content")
                
        if valid_exclusion:        
            exclusion_path = os.path.join(seq_dir, "exclusion_sequences.fasta")
            SeqIO.write(exclusion_sequences, exclusion_path, "fasta")
            logger.info(f"Saved {len(exclusion_sequences)} exclusion sequences to {exclusion_path}")
        else:
            logger.warning("No valid exclusion sequences to save")
    
def save_marker_sequence(output_dir: str, marker_info: Dict[str, Any]):
    """
    Save marker sequence to file.
    
    Args:
        output_dir (str): Output directory
        marker_info (Dict[str, Any]): Marker information
    """
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio import SeqIO
    
    marker_sequence = marker_info.get("marker_sequence", "")
    if marker_sequence:
        seq_record = SeqRecord(
            Seq(marker_sequence),
            id="marker",
            description=f"Marker sequence: {marker_info.get('description', 'Unknown')}"
        )
        
        marker_path = os.path.join(output_dir, "sequences", "marker_sequence.fasta")
        SeqIO.write([seq_record], marker_path, "fasta")
        logger.info(f"Saved marker sequence ({len(marker_sequence)} bp) to {marker_path}")

def save_primers_and_probe(output_dir: str, assay_info: Dict[str, Any]):
    """
    Save primers and probe to files.
    
    Args:
        output_dir (str): Output directory
        assay_info (Dict[str, Any]): Assay information including primers and probe
    """
    from Bio.SeqRecord import SeqRecord
    from Bio.Seq import Seq
    from Bio import SeqIO
    
    results_dir = os.path.join(output_dir, "results")
    
    # Save primers
    primers = assay_info.get("primers", [])
    if primers:
        primer_records = []
        for primer in primers:
            record = SeqRecord(
                Seq(primer.get("sequence", "")),
                id=primer.get("name", "primer"),
                description=f"Tm={primer.get('properties', {}).get('tm', 0):.1f}, GC={primer.get('properties', {}).get('gc_content', 0):.1f}%"
            )
            primer_records.append(record)
        
        primers_path = os.path.join(results_dir, "primers.fasta")
        SeqIO.write(primer_records, primers_path, "fasta")
        logger.info(f"Saved {len(primers)} primers to {primers_path}")
    
    # Save probe
    probe = assay_info.get("probe")
    if probe:
        probe_record = SeqRecord(
            Seq(probe.get("sequence", "")),
            id=probe.get("name", "probe"),
            description=f"Tm={probe.get('properties', {}).get('tm', 0):.1f}, GC={probe.get('properties', {}).get('gc_content', 0):.1f}%"
        )
        
        probe_path = os.path.join(results_dir, "probe.fasta")
        SeqIO.write([probe_record], probe_path, "fasta")
        logger.info(f"Saved probe to {probe_path}")
    
    # Save amplicon
    amplicon = assay_info.get("amplicon", "")
    if amplicon:
        amplicon_record = SeqRecord(
            Seq(amplicon),
            id="amplicon",
            description=f"Length: {len(amplicon)} bp"
        )
        
        amplicon_path = os.path.join(results_dir, "amplicon.fasta")
        SeqIO.write([amplicon_record], amplicon_path, "fasta")
        logger.info(f"Saved amplicon sequence to {amplicon_path}")

def save_metadata(output_dir: str, args: argparse.Namespace, runtime: float, taxa_info: Dict[str, Any]):
    """
    Save metadata about the assay design run.
    
    Args:
        output_dir (str): Output directory
        args (argparse.Namespace): Command line arguments
        runtime (float): Total runtime in seconds
        taxa_info (Dict[str, Any]): Information about taxa used
    """
    metadata = {
        "timestamp": datetime.datetime.now().isoformat(),
        "runtime_seconds": runtime,
        "arguments": vars(args),
        "taxa_info": taxa_info
    }
    
    metadata_path = os.path.join(output_dir, "results", "metadata.json")
    with open(metadata_path, 'w') as f:
        json.dump(metadata, f, indent=2)
    
    logger.info(f"Saved metadata to {metadata_path}")

def save_results_text(output_dir: str, results_text: str):
    """
    Save results as text file.
    
    Args:
        output_dir (str): Output directory
        results_text (str): Formatted results text
    """
    results_path = os.path.join(output_dir, "results", "assay_results.txt")
    with open(results_path, 'w') as f:
        f.write(results_text)
    
    logger.info(f"Saved results to {results_path}")

def main():
    """Main entry point for the CLI."""
    parser = create_parser()
    args = parser.parse_args()
    
    start_time = time.time()
    
    # Set verbosity
    if args.verbose:
        logging.getLogger().setLevel(logging.DEBUG)
    
    # Create output directory structure
    output_dir = create_output_directory(args.output_dir)
    
    # Configure file logging
    log_file = os.path.join(output_dir, "logs", "assay_design.log")
    file_handler = logging.FileHandler(log_file)
    file_handler.setFormatter(logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s'))
    logging.getLogger().addHandler(file_handler)
    
    # Initialize variables for sequences
    inclusion_sequences = []
    exclusion_sequences = []
    target_info = {}
    exclusion_info = []
    exclusion_taxids = []
    selected_gene = None  # Track the selected gene for reporting
    
    # Check if local inclusion FASTA is provided
    if args.inclusion_fasta:
        try:
            logger.info(f"Loading inclusion sequences from local FASTA file: {args.inclusion_fasta}")
            inclusion_sequences = load_local_fasta(args.inclusion_fasta)
            
            if not inclusion_sequences:
                logger.error("No valid sequences found in inclusion FASTA file")
                print("\nNo valid sequences found in inclusion FASTA file")
                sys.exit(1)
                
            # Create a placeholder target info
            fasta_filename = os.path.basename(args.inclusion_fasta)
            target_name = f"Local FASTA: {fasta_filename}"
            target_info = {
                "taxid": "local",
                "scientific_name": target_name
            }
            
            logger.info(f"Loaded {len(inclusion_sequences)} sequences from inclusion FASTA file")
            print(f"Loaded {len(inclusion_sequences)} sequences from inclusion FASTA file")
            
        except Exception as e:
            logger.error(f"Error loading inclusion FASTA file: {e}")
            if args.verbose:
                import traceback
                traceback.print_exc()
            sys.exit(1)
            
    # Check if local exclusion FASTA is provided
    if args.exclusion_fasta:
        try:
            logger.info(f"Loading exclusion sequences from local FASTA file: {args.exclusion_fasta}")
            exclusion_sequences = load_local_fasta(args.exclusion_fasta)
            
            if not exclusion_sequences:
                logger.warning("No valid sequences found in exclusion FASTA file")
                print("\nWarning: No valid sequences found in exclusion FASTA file")
            else:
                # Create placeholder exclusion info
                fasta_filename = os.path.basename(args.exclusion_fasta)
                exclusion_info.append({
                    "taxid": "local",
                    "scientific_name": f"Local FASTA: {fasta_filename}"
                })
                
                logger.info(f"Loaded {len(exclusion_sequences)} sequences from exclusion FASTA file")
                print(f"Loaded {len(exclusion_sequences)} sequences from exclusion FASTA file")
                
        except Exception as e:
            logger.error(f"Error loading exclusion FASTA file: {e}")
            if args.verbose:
                import traceback
                traceback.print_exc()
                
    # Set up CDS workflow if specified and no local FASTA files provided
    if args.cds and not args.inclusion_fasta:
        logger.info(f"Fetching CDS sequences for '{args.cds}'")
        try:
            output_file = getattr(args, 'output_fasta', None)
            if not output_file:
                output_file = os.path.join(output_dir, "sequences", f"{args.cds}_refseq_{args.max_seqs}.fna")
            
            cds_sequences = []
            try:
                cds_sequences = fetch_cds_sequences(
                    cds_name=args.cds,
                    email=args.email,
                    max_seqs=getattr(args, 'max_seqs', 100),
                    output_file=output_file
                )
            except Exception as e:
                logger.error(f"Error fetching CDS sequences: {e}")
                # Continue and check if we got any sequences despite the error
            
            if cds_sequences:
                logger.info(f"Retrieved {len(cds_sequences)} CDS sequences for '{args.cds}'")
                print(f"Retrieved {len(cds_sequences)} CDS sequences for '{args.cds}'")
                print(f"Sequences saved to: {output_file}")
                
                from Bio.Seq import UndefinedSequenceError
                valid_sequences = []
                for seq in cds_sequences:
                    try:
                        # Test if sequence is defined
                        _ = str(seq.seq)
                        valid_sequences.append(seq)
                    except UndefinedSequenceError:
                        logger.warning(f"Skipping sequence {seq.id} with undefined content")
                        
                # Update with only valid sequences
                cds_sequences = valid_sequences
                
                # Check if we still have sequences after filtering
                if not cds_sequences:
                    logger.error("No valid CDS sequences found after filtering")
                    print("\nNo valid CDS sequences found after filtering")
                    sys.exit(1)
                
                # Use CDS sequences as inclusion sequences
                inclusion_sequences = cds_sequences
                
                # Create a placeholder target info
                target_name = f"{args.cds} gene (multiple taxa)"
                target_info = {
                    "taxid": "multiple",
                    "scientific_name": target_name
                }
                
                # Track selected gene
                selected_gene = args.cds
                
                # Force conserved mode when using CDS for diverse taxa
                if not args.conserved_mode:
                    logger.info("Enabling conserved mode for CDS-based assay design")
                    args.conserved_mode = True
            else:
                logger.error(f"No CDS sequences found for '{args.cds}'")
                print(f"\nNo CDS sequences found for '{args.cds}'")
                sys.exit(1)
                
        except Exception as e:
            logger.error(f"Error in CDS processing workflow: {e}")
            if args.verbose:
                import traceback
                traceback.print_exc()
            sys.exit(1)
        
    # Process inclusion taxid if neither CDS nor local FASTA is specified
    elif args.inclusion and not args.inclusion_fasta:
        inclusion_taxid = args.inclusion
        logger.info(f"Target organism: taxid {inclusion_taxid}")
        
        # Get taxonomy information for the target
        try:
            target_info = get_taxon_info(inclusion_taxid, args.email)
            target_name = target_info.get("scientific_name", f"Taxid:{inclusion_taxid}")
            logger.info(f"Target organism {target_name} ({inclusion_taxid})")
        except Exception as e:
            logger.error(f"Error fetching taxonomy information: {e}")
            target_info = {"taxid": inclusion_taxid}
            target_name = f"Taxid:{inclusion_taxid}"
        
        # Process exclusion taxa if not in conserved mode
        if not args.conserved_mode:
            if args.exclusion:
                exclusion_taxids = args.exclusion.split(',')
                logger.info(f"Manually specified exclusion taxa: {exclusion_taxids}")

            # Handle intelligent exclusion strategy
            if args.exclusion_strategy == 'intelligent':
                logger.info(f"Using intelligent multi-level exclusion strategy (max: {args.max_exclusion_taxa} taxa)")
                try:
                    exclusion_result = intelligent_exclusion_selection(
                        inclusion_taxid=inclusion_taxid,
                        email=args.email,
                        max_exclusion_taxa=args.max_exclusion_taxa,
                        require_sequences=True,
                        min_sequence_count=5
                    )

                    if 'error' not in exclusion_result:
                        new_taxids = exclusion_result.get('exclusion_taxids', [])
                        taxa_info_dict = exclusion_result.get('taxa_info', {})

                        for taxid in new_taxids:
                            if taxid not in exclusion_taxids:
                                exclusion_taxids.append(taxid)
                                taxon_info = taxa_info_dict.get(taxid, {'taxid': taxid})
                                exclusion_info.append(taxon_info)
                                logger.info(f"Intelligent exclusion selected: {taxon_info.get('scientific_name', 'Unknown')} "
                                          f"(taxid: {taxid}, level: {taxon_info.get('diversity_level', 'unknown')})")

                        coverage_score = exclusion_result.get('coverage_score', 0)
                        logger.info(f"Phylogenetic coverage score: {coverage_score:.2f}")
                        logger.info(f"Tier summary: {exclusion_result.get('tier_summary', {})}")
                    else:
                        logger.warning(f"Intelligent exclusion failed: {exclusion_result.get('error')}")
                        logger.info("Falling back to default exclusion strategy")

                except Exception as e:
                    logger.error(f"Error in intelligent exclusion selection: {e}")
                    logger.info("Falling back to default exclusion strategy")

            # Auto-generate exclusion taxa with specific rank if requested
            elif args.auto_exclusion or args.exclusion_strategy in ['siblings', 'family', 'order']:
                # Determine rank from strategy
                if args.exclusion_strategy in ['siblings', 'family', 'order']:
                    rank = 'genus' if args.exclusion_strategy == 'siblings' else args.exclusion_strategy
                else:
                    rank = args.auto_exclusion_rank

                logger.info(f"Automatically selecting related taxa at rank '{rank}'")
                related_taxa = get_related_taxa(
                    taxid=inclusion_taxid,
                    email=args.email,
                    relationship="sibling",
                    rank=rank,
                    max_results=args.auto_exclusion_count
                )

                for taxon in related_taxa:
                    taxid = taxon.get("taxid")
                    if taxid and taxid not in exclusion_taxids:
                        exclusion_taxids.append(taxid)
                        exclusion_info.append(taxon)
                        logger.info(f"Auto-selected exclusion taxon: {taxon.get('scientific_name')} (taxid: {taxid})")

            # Get exclusion taxa info if not already obtained
            if not exclusion_info and exclusion_taxids:
                for taxid in exclusion_taxids:
                    try:
                        taxon_info = get_taxon_info(taxid, args.email)
                        exclusion_info.append(taxon_info)
                        logger.info(f"Exclusion taxon: {taxon_info.get('scientific_name')} ({taxid})")
                    except Exception as e:
                        logger.error(f"Error fetching exclusion taxon info: {e}")
                        exclusion_info.append({"taxid": taxid})
                
        # Determine gene selection approach
        if args.auto_gene and not args.gene:
            logger.info(f"Auto-gene selection enabled (use case: {args.gene_use_case})")
            # Gene will be selected by hierarchical_marker_search
            selected_gene = None
        elif args.gene:
            selected_gene = args.gene
            logger.info(f"Using manually specified gene: {selected_gene}")
        else:
            selected_gene = None

        # Fetch sequences for inclusion taxid only if not using auto-gene
        # (auto-gene will fetch sequences internally)
        if not args.auto_gene:
            logger.info(f"Fetching sequences for inclusion taxid {inclusion_taxid}")
            try:
                if selected_gene:
                    inclusion_sequences = fetch_gene_sequences(
                        taxid=inclusion_taxid,
                        gene_name=selected_gene,
                        email=args.email,
                        max_records=args.max_seq_count * 2
                    )
                else:
                    inclusion_sequences = fetch_sequences_for_taxid(
                        taxid=inclusion_taxid,
                        email=args.email,
                        max_records=args.max_seq_count * 2
                    )

                logger.info(f"Retrieved {len(inclusion_sequences)} sequences for inclusion taxid")

                if not inclusion_sequences:
                    logger.error(f"No sequences found for taxid {inclusion_taxid}" +
                                (f" with gene {selected_gene}" if selected_gene else ""))
                    sys.exit(1)

            except Exception as e:
                logger.error(f"Error fetching inclusion sequences: {e}")
                sys.exit(1)

    # Save sequences to files
    save_sequences(output_dir, inclusion_sequences, exclusion_sequences)
    
    # Fetch sequences for exclusion taxids if not in conserved mode
    if not args.conserved_mode and exclusion_taxids:
        logger.info(f"Fetching sequences for {len(exclusion_taxids)} exclusion taxids")
        
        for taxid in exclusion_taxids:
            try:
                sequences = []
                if args.gene:
                    sequences = fetch_gene_sequences(
                        taxid=taxid,
                        gene_name=args.gene,
                        email=args.email,
                        max_records=args.max_seq_count
                    )
                else:
                    sequences = fetch_sequences_for_taxid(
                        taxid=taxid,
                        email=args.email,
                        max_records=args.max_seq_count
                    )
                    
                exclusion_sequences.extend(sequences)
                logger.info(f"Retrieved {len(sequences)} sequences for exclusion taxid {taxid}")
            except Exception as e:
                logger.error(f"Error fetching exclusion sequences for taxid {taxid}: {e}")
    
    # Find conserved markers using the specified approach
    logger.info("Identifying conserved marker regions")
    
    try:
        if args.conserved_mode:
            # Use optimized approach for conserved gene mode (no exclusions)
            logger.info("Using conserved gene mode (optimizing for conservation across diverse taxa)")
            marker_info = find_conserved_without_exclusion(
                inclusion_sequences=inclusion_sequences,
                min_conservation=args.min_conservation,
                max_amplicon_length=args.max_amplicon_length,
                timeout_seconds=args.max_runtime
            )
        elif args.search_strategy == 'basic':
            # Use basic search approach
            marker_info = find_optimal_marker(
                inclusion_sequences,
                exclusion_sequences,
                timeout_seconds=args.max_runtime
            )
        else:
            # Use hierarchical search (default)
            if args.inclusion_fasta:
                # For local FASTA files, use "local" as the taxid
                inclusion_id = "local"
            else:
                inclusion_id = args.inclusion if args.inclusion else "multiple"
                
            marker_info = hierarchical_marker_search(
                inclusion_taxid=inclusion_id,
                email=args.email,
                gene_name=selected_gene,
                exclusion_taxids=exclusion_taxids if exclusion_taxids else None,
                max_amplicon_length=args.max_amplicon_length,
                max_seq_count=args.max_seq_count,
                max_seq_length=args.max_seq_length,
                min_conservation=args.min_conservation,
                min_specificity=args.min_specificity,
                timeout_seconds=args.max_runtime,
                use_lsh=not args.no_lsh,
                binding_sites_only=not args.whole_amplicon_specificity,
                inclusion_sequences=inclusion_sequences if (args.inclusion_fasta or not args.auto_gene) else None,
                exclusion_sequences=exclusion_sequences if args.exclusion_fasta else None,
                auto_gene_selection=args.auto_gene,
                gene_use_case=args.gene_use_case,
                max_genes_to_try=args.max_genes_to_try,
                min_gene_score=args.min_gene_score
            )

            # Extract gene selection info if auto-gene was used
            if args.auto_gene and 'gene_selection_info' in marker_info:
                gene_info = marker_info['gene_selection_info']
                if 'selected_gene' in gene_info:
                    selected_gene = gene_info['selected_gene']
                    logger.info(f"Auto-selected gene: {selected_gene} (score: {gene_info.get('gene_score', 0):.2f})")
                elif 'error' in gene_info:
                    logger.warning(f"Auto-gene selection encountered error: {gene_info['error']}")
    except Exception as e:
        logger.error(f"Error in marker identification: {e}")
        marker_info = {
            "error": f"Marker identification failed: {str(e)}",
            "marker_sequence": "",
            "marker_length": 0
        }
    
    # Save marker sequence
    save_marker_sequence(output_dir, marker_info)
    
    # Design primers and probe
    logger.info("Designing primers and probe for the identified marker region")

    if args.single_region_design:
        # Use single region approach
        assay_info = design_primers_and_probe(marker_info)
    else:
        # Use multi-region approach (default)
        assay_info = design_multi_region_assay(
            marker_info=marker_info,
            inclusion_sequences=inclusion_sequences,
            exclusion_sequences=exclusion_sequences
        )
    
    # Check if design was successful
    if assay_info.get("status") == "failed":
        logger.error(f"Assay design failed: {assay_info.get('message')}")
    
        # Save primers and probe
        save_primers_and_probe(output_dir, assay_info)
        
        # Create a simple result text
        gene_info = f"Target gene: {selected_gene}" if selected_gene else ""
        target_id = args.inclusion if args.inclusion else "multiple taxa"
        
        failure_text = f"""=== Assay Design Results ===
        
Target: {target_name} (taxid: {target_id})
{gene_info}

ERROR: Could not design suitable assay

Reason: {assay_info.get('message')}

Please consider the following options:
1. Try a different target gene
2. Adjust parameters (e.g., increase maximum amplicon length)
3. Use manual exclusion instead of auto-exclusion
"""
        # Save results text
        save_results_text(output_dir, failure_text)
    
        # Print results to console
        print("\n" + failure_text)
        print(f"\nResults saved to: {output_dir}")
    
        logger.info(f"Assay design completed with errors in {time.time() - start_time:.2f} seconds")
        return 1
    else:
        # Success path
        # Validate specificity
        logger.info("Validating primer specificity")
        validation_results = validate_primer_specificity(
            assay_info["primers"], 
            args.inclusion if args.inclusion else "multiple",
            exclusion_taxids,
            email=args.email
        )
    
        # Create a composite result with all information
        complete_results = {
            "target": {
                "taxid": args.inclusion if args.inclusion else "multiple",
                "name": target_name,
                "info": target_info
            },
            "exclusion": [
                {"taxid": info.get("taxid"), "name": info.get("scientific_name", "Unknown")}
                for info in exclusion_info
            ],
            "gene": selected_gene,
            "marker": marker_info,
            "assay": assay_info,
            "validation": validation_results,
            "runtime": time.time() - start_time
        }
    
        # Save complete results as JSON
        results_json_path = os.path.join(output_dir, "results", "complete_results.json")
        with open(results_json_path, 'w') as f:
            json.dump(complete_results, f, indent=2)
        
        logger.info(f"Saved complete results to {results_json_path}")
        
        # Format and output human-readable results
        target_id = args.inclusion if args.inclusion else "multiple"
        results_text = format_results(
            inclusion_taxid=target_id,
            target_name=target_name,
            gene=selected_gene,
            exclusion_taxids=exclusion_taxids,
            exclusion_info=exclusion_info,
            marker_info=marker_info,
            assay_info=assay_info,
            validation_results=validation_results,
            runtime=time.time() - start_time
        )
        
        # Save results text
        save_results_text(output_dir, results_text)
        
        # Create visualizations
        if not args.no_visualization:
            logger.info("Creating visualizations")
            try:
                create_visualizations(output_dir, complete_results)
            except Exception as e:
                logger.error(f"Error creating visualizations: {e}")
        
        # Save metadata
        taxa_info = {
            "inclusion": target_info,
            "exclusion": exclusion_info
        }
        save_metadata(output_dir, args, time.time() - start_time, taxa_info)
        
        # Print results to console
        print("\n" + results_text)
        print(f"\nResults saved to: {output_dir}")
        
        logger.info(f"Assay design completed in {time.time() - start_time:.2f} seconds")
        return 0
        
def format_results(
    inclusion_taxid: str,
    target_name: str,
    gene: Optional[str],
    exclusion_taxids: List[str],
    exclusion_info: List[Dict[str, Any]],
    marker_info: Dict[str, Any],
    assay_info: Dict[str, Any],
    validation_results: Dict[str, Any],
    runtime: float
) -> str:
    """Format the results as a string."""
    lines = ["=== Assay Design Results ===", ""]
    
    lines.append(f"Target: {target_name}")
    # Only show taxid if it's not a local file
    if inclusion_taxid != "local":
        lines.append(f"Taxid: {inclusion_taxid}")
        
    if gene:
        lines.append(f"Target gene: {gene}")
    
    # Format exclusion taxa
    exclusion_names = []
    for info in exclusion_info:
        if info.get("taxid") == "local":
            exclusion_names.append(info.get("scientific_name", "Local FASTA file"))
        else:
            name = info.get("scientific_name", f"Taxid:{info.get('taxid')}")
            exclusion_names.append(f"{name} ({info.get('taxid')})")
    
    if exclusion_names:
        lines.append(f"Exclusion taxa: {', '.join(exclusion_names)}")
    lines.append(f"Total runtime: {runtime:.2f} seconds")
    
    # Marker region information
    lines.append("\nIdentified marker region:")
    lines.append(f"Length: {marker_info.get('marker_length')} bp")
    lines.append(f"Conservation score: {marker_info.get('conservation_score', 0):.2f}")
    if 'specificity_score' in marker_info:
        lines.append(f"Specificity score: {marker_info.get('specificity_score', 0):.2f}")
    lines.append(f"Description: {marker_info.get('description')}")
    
    # Only show part of the marker sequence if it's long
    marker_seq = marker_info.get('marker_sequence', '')
    if len(marker_seq) > 100:
        lines.append(f"Marker sequence (first 50 bp): {marker_seq[:50]}...")
        lines.append(f"Marker sequence (last 50 bp): ...{marker_seq[-50:]}")
    else:
        lines.append(f"Marker sequence: {marker_seq}")
    
    # Primers
    lines.append("\nDesigned primers:")
    for primer in assay_info.get("primers", []):
        lines.append(f"{primer.get('name')} ({primer.get('orientation', 'unknown')}):")
        lines.append(f"  Sequence: {primer.get('sequence', '')}")
        
        properties = primer.get('properties', {})
        if properties:
            lines.append(f"  Length: {properties.get('length', 0)} bp")
            lines.append(f"  GC Content: {properties.get('gc_content', 0):.1f}%")
            lines.append(f"  Tm: {properties.get('tm', 0):.1f}°C")
    
    # Probe
    probe = assay_info.get("probe")
    if probe:
        lines.append("\nDesigned probe:")
        lines.append(f"  Sequence: {probe.get('sequence', '')}")
        
        properties = probe.get('properties', {})
        if properties:
            lines.append(f"  Length: {properties.get('length', 0)} bp")
            lines.append(f"  GC Content: {properties.get('gc_content', 0):.1f}%")
            lines.append(f"  Tm: {properties.get('tm', 0):.1f}°C")
    
    # Amplicon
    amplicon = assay_info.get("amplicon", "")
    if amplicon:
        lines.append("\nAmplicon:")
        lines.append(f"  Length: {len(amplicon)} bp")
        if len(amplicon) > 100:
            lines.append(f"  Sequence (first 50 bp): {amplicon[:50]}...")
            lines.append(f"  Sequence (last 50 bp): ...{amplicon[-50:]}")
        else:
            lines.append(f"  Sequence: {amplicon}")
    
    # Specificity validation
    lines.append("\nSpecificity validation:")
    if "error" in validation_results:
        lines.append(f"Error: {validation_results.get('error')}")
    else:
        lines.append(f"Specificity score: {validation_results.get('specificity_score', 0):.2f}")
        
        if validation_results.get('inclusion_hits'):
            lines.append("\nExpected amplification in target:")
            for hit in validation_results.get('inclusion_hits', []):
                lines.append(f"  {hit.get('organism')} - Product size: {hit.get('product_size')} bp")
        
        if validation_results.get('exclusion_hits'):
            lines.append("\nPotential cross-reactivity:")
            for hit in validation_results.get('exclusion_hits', []):
                lines.append(f"  {hit.get('organism')} - Product size: {hit.get('product_size')} bp")
        
        if validation_results.get('potential_issues'):
            lines.append("\nPotential issues:")
            for issue in validation_results.get('potential_issues', []):
                lines.append(f"  - {issue}")
    
    return "\n".join(lines)

if __name__ == "__main__":
    sys.exit(main())