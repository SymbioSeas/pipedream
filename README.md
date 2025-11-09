# Assay Design Package

A Python package for designing specific molecular assays targeting particular taxa while excluding others.

## Features

- **Automated Gene Selection**: Intelligently selects the most suitable gene markers based on sequence availability, copy number, and evolutionary properties
- **Phylogenetic-Aware Exclusion**: Automatically identify related taxa for exclusion using phylogenetic distance weighting
- **Flexible Exclusion Strategies**: Choose from genus-level, family-level, order-level, or custom exclusion approaches
- **Multi-Domain Support**: Automatic detection and handling of Bacteria, Archaea, Eukaryota, and Viruses
- **Gene Suitability Scoring**: Evaluates candidate genes based on sequence availability, copy number stability, HGT resistance, and optimal length
- **Fetch sequence data from NCBI**: Retrieve genomic and gene sequences for any taxonomic ID
- **Find conserved marker regions**: Identify taxa-specific conserved regions
- **Design and validate primers**: Automated primer design for marker regions
- **Command-line interface**: Easy-to-use CLI with extensive customization options

## Installation

### From PyPI

```bash
pip install assay_design
```

### From source

```bash
git clone https://github.com/yourusername/assay_design.git
cd assay_design
pip install -e .
```

## Usage

### Command Line Interface

Design an assay for a specific organism (using taxid 689 for Vibrio mediterranei as an example):

#### Fully Automated Workflow (Recommended)

```bash
# Fully automated: auto-select gene AND exclusion taxa
assay-design --inclusion 689 --email your.email@example.com --auto-gene --auto-exclusion

# Automated with phylogenetic distance weighting (prioritizes phylogenetically closer taxa)
assay-design --inclusion 689 --email your.email@example.com --auto-gene --auto-exclusion --exclusion-phylo-distance 0.8

# Use family-level exclusion strategy (broader exclusion scope)
assay-design --inclusion 689 --email your.email@example.com --auto-gene --exclusion-strategy family
```

#### Automated Gene Selection

```bash
# Let the system choose the best gene automatically
assay-design --inclusion 689 --email your.email@example.com --auto-gene --auto-exclusion

# Customize gene selection criteria (JSON format)
assay-design --inclusion 689 --email your.email@example.com --auto-gene \
  --gene-criteria '{"min_sequence_count": 20, "ideal_length_range": [800, 2000], "single_copy_preferred": true}'

# View available genes and their suitability scores
assay-design --inclusion 689 --email your.email@example.com --suggest-genes
```

#### Exclusion Strategy Options

```bash
# Genus-level exclusion (default) - excludes sibling species in the same genus
assay-design --inclusion 689 --email your.email@example.com --exclusion-strategy genus

# Family-level exclusion - excludes related genera in the same family
assay-design --inclusion 689 --email your.email@example.com --exclusion-strategy family

# Order-level exclusion - excludes related families in the same order
assay-design --inclusion 689 --email your.email@example.com --exclusion-strategy order

# Custom exclusion - manually specify taxa to exclude
assay-design --inclusion 689 --exclusion 717,670,672 --email your.email@example.com
```

#### Manual Gene Selection

```bash
# Specify a particular gene
assay-design --inclusion 689 --gene "rpoB" --email your.email@example.com --auto-exclusion

# Specify gene with custom exclusion
assay-design --inclusion 689 --gene "16S ribosomal RNA" --exclusion 717,670,672 --email your.email@example.com
```

#### Advanced Options

```bash
# Disable LSH optimization for traditional k-mer comparison
assay-design --inclusion 689 --email your.email@example.com --auto-gene --auto-exclusion --no-lsh

# Combine multiple options for fine-tuned control
assay-design --inclusion 689 --email your.email@example.com \
  --auto-gene \
  --gene-criteria '{"min_sequence_count": 15, "ideal_length_range": [1000, 1800]}' \
  --exclusion-strategy family \
  --exclusion-phylo-distance 0.75
```

### Python API

#### Fully Automated Workflow (Recommended)

```python
from assay_design.gene_selection import (
    auto_select_gene_for_taxon,
    GeneSelectionCriteria
)
from assay_design.hierarchical_search import intelligent_exclusion_selection
from assay_design.data_retrieval import fetch_gene_sequences
from assay_design.target_identification import find_conserved_marker
from assay_design.primer_design import design_primers

# Step 1: Automatically select the best gene
gene_result = auto_select_gene_for_taxon(
    taxid="689",  # Vibrio mediterranei
    email="your.email@example.com"
)
selected_gene = gene_result['gene']
print(f"Selected gene: {selected_gene} (score: {gene_result['total_score']:.2f})")

# Step 2: Intelligently select exclusion taxa using phylogenetic distance weighting
exclusion_taxa = intelligent_exclusion_selection(
    taxid="689",
    email="your.email@example.com",
    strategy="genus",  # or "family", "order"
    phylo_distance_weight=0.8,
    max_exclusion_taxa=5
)
exclusion_taxids = [t['taxid'] for t in exclusion_taxa]

# Step 3: Fetch sequences for inclusion and exclusion
inclusion_sequences = fetch_gene_sequences(
    taxid="689",
    gene_name=selected_gene,
    email="your.email@example.com",
    max_records=10
)

exclusion_sequences = []
for taxid in exclusion_taxids:
    sequences = fetch_gene_sequences(
        taxid=taxid,
        gene_name=selected_gene,
        email="your.email@example.com",
        max_records=5
    )
    exclusion_sequences.extend(sequences)

# Step 4: Find conserved markers
marker_info = find_conserved_marker(
    inclusion_sequences=inclusion_sequences,
    exclusion_sequences=exclusion_sequences
)

# Step 5: Design primers
primers = design_primers(marker_info)

print(f"Designed primers: {primers}")
```

#### Manual Gene Selection with Custom Criteria

```python
from assay_design.gene_selection import (
    rank_candidate_genes,
    evaluate_gene_suitability,
    GeneSelectionCriteria,
    BACTERIA_GENES
)

# Define custom gene selection criteria
custom_criteria = GeneSelectionCriteria(
    min_sequence_count=20,
    ideal_length_range=(1000, 1800),
    single_copy_preferred=True,
    hgt_resistant_preferred=True
)

# Evaluate all candidate genes
gene_scores = []
for gene_name, gene_info in BACTERIA_GENES.items():
    result = evaluate_gene_suitability(
        taxid="689",
        gene_name=gene_name,
        email="your.email@example.com",
        criteria=custom_criteria
    )
    gene_scores.append(result)

# Rank genes by suitability
ranked_genes = rank_candidate_genes(gene_scores)

# Use the top-ranked gene
best_gene = ranked_genes[0]
print(f"Best gene: {best_gene['gene']}")
print(f"Score: {best_gene['total_score']:.2f}")
print(f"Sequences available: {best_gene['sequence_count']}")
```

#### Traditional Workflow with get_related_taxa

```python
from assay_design.data_retrieval import fetch_sequences_for_taxid, get_related_taxa
from assay_design.target_identification import find_conserved_marker
from assay_design.primer_design import design_primers

# Fetch sequences for inclusion taxid
inclusion_sequences = fetch_sequences_for_taxid(
    taxid="689",  # Vibrio mediterranei
    email="your.email@example.com",
    max_records=10
)

# Get related taxa for exclusion with phylogenetic distance weighting
related_taxa = get_related_taxa(
    taxid="689",
    email="your.email@example.com",
    relationship="sibling",
    max_results=5,
    phylo_distance_weight=0.8  # Prioritize phylogenetically closer taxa
)
exclusion_taxids = [taxon["taxid"] for taxon in related_taxa]

# Fetch sequences for exclusion taxids
exclusion_sequences = []
for taxid in exclusion_taxids:
    sequences = fetch_sequences_for_taxid(
        taxid=taxid,
        email="your.email@example.com",
        max_records=5
    )
    exclusion_sequences.extend(sequences)

# Find conserved markers
marker_info = find_conserved_marker(
    inclusion_sequences=inclusion_sequences,
    exclusion_sequences=exclusion_sequences
)

# Design primers
primers = design_primers(marker_info)

print(f"Designed primers: {primers}")
```

## Gene Selection Algorithm

### How Automated Gene Selection Works

The automated gene selection feature evaluates candidate genes based on multiple criteria to identify the most suitable marker for your target taxon:

1. **Sequence Availability**: Evaluates how many sequences are available in NCBI databases
   - Optimal: 20+ sequences available
   - Score decreases for sparse sequence data

2. **Copy Number Stability**: Prefers single-copy genes to avoid PCR bias
   - Single-copy genes: Full score
   - Multi-copy genes (e.g., rRNA): Partial score

3. **HGT (Horizontal Gene Transfer) Resistance**: Prioritizes genes with low HGT rates
   - Core housekeeping genes: Higher score
   - Mobile or frequently transferred genes: Lower score

4. **Sequence Length**: Targets optimal length range for primer design
   - Ideal: 800-1500 bp
   - Penalties for sequences too short (<500 bp) or too long (>3000 bp)

### Available Gene Databases

The package includes curated gene databases for different taxonomic domains:

**Bacteria** (30 genes):
- rpoB, rpoD, gyrB, recA, dnaK, groEL, atpD, fusA
- 16S rRNA, 23S rRNA, infB, tuf, rplB, rpsB
- And 16 additional housekeeping genes

**Archaea** (20 genes):
- rpoB, rpoA1, rpoA2, EF-2, SecY, RecA
- 16S rRNA, 23S rRNA, and 12 additional marker genes

**Eukaryota** (20 genes):
- 18S rRNA, 28S rRNA, ITS, COI, CytB
- actin, tubulin, EF-1α, and 12 additional markers

**Viruses** (15 genes):
- RdRp, capsid, polymerase, and 12 viral-specific markers

### Customizing Gene Selection Criteria

You can customize the gene selection criteria using the `GeneSelectionCriteria` class:

```python
from assay_design.gene_selection import GeneSelectionCriteria

custom_criteria = GeneSelectionCriteria(
    min_sequence_count=15,           # Minimum sequences required
    ideal_length_range=(1000, 1800), # Optimal length range (bp)
    single_copy_preferred=True,      # Prefer single-copy genes
    hgt_resistant_preferred=True     # Prefer HGT-resistant genes
)
```

## Intelligent Exclusion Selection

### Phylogenetic Distance Weighting

The `get_related_taxa()` function now supports phylogenetic distance weighting, which prioritizes taxa that are phylogenetically closer to your target:

- **Weight = 0.0**: No distance weighting (all taxa equally likely)
- **Weight = 0.5**: Moderate preference for closer taxa
- **Weight = 0.8**: Strong preference for closer taxa (recommended)
- **Weight = 1.0**: Maximum preference for closest taxa

This helps design more specific assays by focusing on the most relevant exclusion taxa.

### Exclusion Strategies

Three hierarchical exclusion strategies are available:

1. **Genus-level** (default): Excludes sibling species within the same genus
   - Best for species-specific assays
   - Example: Targeting *Vibrio mediterranei*, excludes other *Vibrio* species

2. **Family-level**: Excludes related genera within the same family
   - Best for genus-specific assays
   - Example: Targeting *Vibrio* genus, excludes other Vibrionaceae genera

3. **Order-level**: Excludes related families within the same order
   - Best for family-specific assays
   - Example: Targeting Vibrionaceae, excludes other Vibrionales families

## Performance Optimization

By default, the package uses Locality-Sensitive Hashing (LSH) for efficient sequence comparison when identifying specific marker regions. This significantly improves performance for large datasets.

If you need to use the traditional k-mer comparison method (e.g., for debugging or validation), you can disable LSH with the '--no-lsh' flag:

```bash
assay-design --inclusion 689 --email your.email@example.com --no-lsh
```

## Requirements

- Python 3.7+
- Biopython
- External tools (optional but recommended):
  - MAFFT, MUSCLE, or ClustalW for multiple sequence alignment
  - Primer3 for advanced primer design (the package includes a basic implementation)

## License

This project is licensed under the MIT License - see the LICENSE file for details.