# assay_design/__init__.py

"""
Assay Design Package

A Python package for designing specific molecular assays targeting particular
taxa while excluding others.
"""

__version__ = '0.1.0'

# Import primary functions to expose at package level
from .data_retrieval import (
    fetch_sequence_by_accession,
    fetch_sequences_for_taxid,
    fetch_gene_sequences,
    get_taxon_info,
    get_related_taxa,
    suggest_marker_genes,
    load_local_fasta,
    fetch_cds_sequences
)

# Import the optimized target identification function
from .target_identification import (
    find_optimal_marker,
    find_conserved_without_exclusion,
    find_discriminative_kmers,
    extend_kmer
)

# Import primer design functions
from .primer_probe_design import (
    design_primers_and_probe,
    calculate_oligo_properties,
    find_conserved_kmers
)

# Import specificity validation
from .specificity_validation import validate_primer_specificity

# Import visualization functions
from .visualization import (
    create_visualizations,
    create_assay_schematic,
    create_conservation_plot,
    create_summary_figure,
    create_html_report
)

# Import CLI
from .cli import main

__all__ = [
    # Data retrieval
    'fetch_sequence_by_accession',
    'fetch_sequences_for_taxid',
    'fetch_gene_sequences',
    'get_taxon_info',
    'get_related_taxa',
    'suggest_marker_genes',
    'load_local_fasta',
    'fetch_cds_sequences',
    
    # Target identification
    'find_optimal_marker',
    'find_discriminative_kmers',
    'extend_kmer',
    'find_conserved_without_exclusion',
    
    # Primer and probe design
    'design_primers_and_probe',
    'calculate_oligo_properties',
    'find_conserved_kmers',
    
    # Specificity validation
    'validate_primer_specificity',
    
    # Visualization
    'create_visualizations',
    'create_assay_schematic',
    'create_conservation_plot',
    'create_temperature_plot',
    'create_summary_figure',
    'create_html_report',
    
    # CLI
    'main'
]