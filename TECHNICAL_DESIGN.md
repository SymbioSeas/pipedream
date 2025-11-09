# Technical Design Document: Automated Gene Selection and Intelligent Exclusion

This document provides technical details about the implementation of automated gene selection and intelligent exclusion features in the assay_design package.

## Table of Contents

1. [Architecture Overview](#architecture-overview)
2. [Gene Selection Module](#gene-selection-module)
3. [Intelligent Exclusion Selection](#intelligent-exclusion-selection)
4. [Data Structures](#data-structures)
5. [Algorithms](#algorithms)
6. [Integration Points](#integration-points)
7. [Testing Strategy](#testing-strategy)
8. [Performance Considerations](#performance-considerations)

---

## Architecture Overview

### Module Structure

```
assay_design/
├── gene_selection.py          # Gene evaluation and selection logic
├── hierarchical_search.py     # Intelligent exclusion selection
├── data_retrieval.py          # Enhanced with phylo-distance weighting
├── cli.py                     # Updated CLI with new flags
└── target_identification.py   # Conserved marker finding (existing)
```

### Data Flow

```
User Input (taxid)
    ↓
[Domain Detection] → Auto-select gene database
    ↓
[Gene Evaluation] → Evaluate all candidate genes
    ↓
[Gene Ranking] → Rank by composite score
    ↓
[Exclusion Selection] → Select exclusion taxa with phylo-weighting
    ↓
[Sequence Retrieval] → Fetch sequences for inclusion + exclusion
    ↓
[Marker Finding] → Identify conserved, specific regions
    ↓
[Primer Design] → Design primers for markers
    ↓
Results
```

---

## Gene Selection Module

### File: `assay_design/gene_selection.py`

#### Core Classes

##### `GeneSelectionCriteria`

```python
@dataclass
class GeneSelectionCriteria:
    """Criteria for evaluating gene suitability."""
    min_sequence_count: int = 10
    ideal_length_range: Tuple[int, int] = (800, 1500)
    single_copy_preferred: bool = True
    hgt_resistant_preferred: bool = True
```

**Design rationale:**
- Immutable dataclass for type safety
- Sensible defaults based on molecular biology best practices
- Easily customizable via CLI JSON or Python API

#### Core Functions

##### `evaluate_gene_suitability()`

**Signature:**
```python
def evaluate_gene_suitability(
    taxid: str,
    gene_name: str,
    email: str,
    criteria: Optional[GeneSelectionCriteria] = None
) -> Dict[str, Any]
```

**Algorithm:**
1. Fetch sequences from NCBI using `fetch_gene_sequences()`
2. Calculate sequence availability score (0-1 scale)
3. Calculate copy number score (1.0 for single-copy, 0.5 for multi-copy)
4. Calculate HGT resistance score (1.0 for resistant, 0.7 for susceptible)
5. Calculate length score based on ideal range
6. Compute weighted composite score

**Scoring weights:**
```python
WEIGHTS = {
    'sequence_availability': 0.35,
    'copy_number': 0.25,
    'hgt_resistance': 0.20,
    'length': 0.20
}
```

**Output format:**
```python
{
    'gene': 'rpoB',
    'sequence_count': 147,
    'sequence_availability_score': 1.0,
    'copy_number_score': 1.0,
    'hgt_score': 1.0,
    'length_score': 0.95,
    'total_score': 9.2,
    'avg_length': 3900,
    'domain': 'Bacteria'
}
```

##### `rank_candidate_genes()`

**Signature:**
```python
def rank_candidate_genes(
    gene_evaluations: List[Dict[str, Any]]
) -> List[Dict[str, Any]]
```

**Algorithm:**
1. Filter out genes with 0 sequences
2. Sort by total_score (descending)
3. Return sorted list

**Time complexity:** O(n log n) where n = number of genes

##### `auto_select_gene_for_taxon()`

**Signature:**
```python
def auto_select_gene_for_taxon(
    taxid: str,
    email: str,
    criteria: Optional[GeneSelectionCriteria] = None
) -> Dict[str, Any]
```

**Algorithm:**
1. Detect taxonomic domain using `get_taxon_info()`
2. Select appropriate gene database (BACTERIA_GENES, ARCHAEA_GENES, etc.)
3. Evaluate all genes in parallel (potential optimization point)
4. Rank genes by composite score
5. Return top-ranked gene

**Error handling:**
- Raises `ValueError` if no suitable genes found
- Logs warnings if best score < 6.0

---

## Intelligent Exclusion Selection

### File: `assay_design/hierarchical_search.py`

#### Core Function

##### `intelligent_exclusion_selection()`

**Signature:**
```python
def intelligent_exclusion_selection(
    taxid: str,
    email: str,
    strategy: str = "genus",
    phylo_distance_weight: float = 0.8,
    max_exclusion_taxa: int = 5
) -> List[Dict[str, Any]]
```

**Strategy mapping:**
```python
STRATEGY_MAP = {
    'genus': 'sibling',      # Species in same genus
    'family': 'cousin',      # Genera in same family
    'order': 'second_cousin' # Families in same order
}
```

**Algorithm:**
1. Map strategy to NCBI relationship type
2. Call `get_related_taxa()` with phylo_distance_weight
3. Limit results to max_exclusion_taxa
4. Return list of taxa dictionaries

**Example output:**
```python
[
    {
        'taxid': '717',
        'scientific_name': 'Vibrio cholerae',
        'rank': 'species',
        'phylo_distance': 0.15,
        'weighted_score': 0.85
    },
    # ... more taxa
]
```

### Enhanced `get_related_taxa()`

**File:** `assay_design/data_retrieval.py`

**New parameter:** `phylo_distance_weight: float = 0.0`

**Phylogenetic distance weighting algorithm:**

```python
def calculate_weighted_score(taxon, phylo_distance_weight):
    """
    Calculate weighted score favoring phylogenetically closer taxa.

    Score = (1 - phylo_distance) * phylo_distance_weight +
            random_factor * (1 - phylo_distance_weight)
    """
    phylo_distance = taxon.get('phylo_distance', 0.5)  # 0 = very close, 1 = distant
    random_factor = random.random()

    score = (1 - phylo_distance) * phylo_distance_weight + \
            random_factor * (1 - phylo_distance_weight)

    return score
```

**Behavior by weight:**
- `weight = 0.0`: Pure random selection (original behavior)
- `weight = 0.5`: Balanced random + phylogenetic preference
- `weight = 0.8`: Strong phylogenetic preference (recommended)
- `weight = 1.0`: Purely deterministic phylogenetic ordering

**Implementation details:**
1. Fetch related taxa from NCBI Taxonomy
2. Calculate phylogenetic distances (based on taxonomic hierarchy)
3. Apply weighting formula to each taxon
4. Sort by weighted score (descending)
5. Return top N taxa

---

## Data Structures

### Gene Databases

Each gene database is a dictionary mapping gene names to metadata:

```python
BACTERIA_GENES = {
    'rpoB': {
        'full_name': 'RNA polymerase beta subunit',
        'typical_length': 3900,
        'copy_number': 'single',
        'hgt_susceptibility': 'low',
        'phylogenetic_utility': 'high'
    },
    'gyrB': {
        'full_name': 'DNA gyrase subunit B',
        'typical_length': 2100,
        'copy_number': 'single',
        'hgt_susceptibility': 'low',
        'phylogenetic_utility': 'high'
    },
    # ... 28 more genes
}
```

**Gene categories:**
1. **Core housekeeping genes**: rpoB, rpoD, gyrB, recA, dnaK
2. **Ribosomal genes**: 16S rRNA, 23S rRNA, rplB, rpsB
3. **Translation factors**: tuf, infB, fusA
4. **Metabolic genes**: atpD, pgk, pyrH

**Selection rationale:**
- Genes selected based on literature review
- Focus on phylogenetically informative markers
- Balance between sequence availability and specificity
- Include both single-copy and multi-copy options

---

## Algorithms

### Sequence Availability Scoring

```python
def score_sequence_availability(count: int, min_count: int) -> float:
    """
    Score sequence availability on 0-1 scale.

    Sigmoid function for smooth scoring:
    - count >= 20: score = 1.0
    - count = 10: score ≈ 0.7
    - count = 5: score ≈ 0.5
    - count < min_count: score = 0.0
    """
    if count < min_count:
        return 0.0

    # Sigmoid scaling
    optimal_count = 20
    if count >= optimal_count:
        return 1.0

    # Smooth curve between min and optimal
    x = (count - min_count) / (optimal_count - min_count)
    return x ** 0.5  # Square root for gentle scaling
```

### Length Scoring

```python
def score_length(length: float, ideal_range: Tuple[int, int]) -> float:
    """
    Score sequence length based on ideal range.

    Penalties:
    - Too short (<500 bp): Heavy penalty
    - Too long (>3000 bp): Moderate penalty
    - Within ideal range: Maximum score
    """
    min_ideal, max_ideal = ideal_range

    if min_ideal <= length <= max_ideal:
        return 1.0

    if length < 500:
        return 0.3  # Heavy penalty

    if length > 3000:
        return 0.6  # Moderate penalty

    if length < min_ideal:
        # Penalty proportional to distance from min_ideal
        return 0.7 + 0.3 * (length / min_ideal)

    if length > max_ideal:
        # Penalty proportional to distance from max_ideal
        penalty = min(0.4, (length - max_ideal) / 1000 * 0.2)
        return 1.0 - penalty

    return 0.8  # Default fallback
```

### Phylogenetic Distance Calculation

**Simplified taxonomic distance:**

```python
def calculate_phylo_distance(taxon1_lineage, taxon2_lineage) -> float:
    """
    Calculate normalized phylogenetic distance (0-1 scale).

    Distance based on taxonomic rank divergence:
    - Same species: 0.0
    - Same genus: 0.1
    - Same family: 0.3
    - Same order: 0.5
    - Same class: 0.7
    - Same phylum: 0.9
    - Different phyla: 1.0
    """
    ranks = ['species', 'genus', 'family', 'order', 'class', 'phylum']
    distances = [0.0, 0.1, 0.3, 0.5, 0.7, 0.9]

    for rank, distance in zip(ranks, distances):
        if taxon1_lineage[rank] == taxon2_lineage[rank]:
            return distance

    return 1.0  # Different phyla
```

---

## Integration Points

### CLI Integration

**File:** `assay_design/cli.py`

**New flags:**

```python
@click.option('--auto-gene', is_flag=True,
              help='Automatically select the best gene')

@click.option('--gene-criteria', type=str,
              help='JSON string with gene selection criteria')

@click.option('--exclusion-strategy',
              type=click.Choice(['genus', 'family', 'order', 'custom']),
              default='genus',
              help='Hierarchical exclusion strategy')

@click.option('--exclusion-phylo-distance', type=float, default=0.8,
              help='Phylogenetic distance weighting (0-1)')
```

**Integration logic:**

```python
if auto_gene:
    # Parse criteria if provided
    criteria = None
    if gene_criteria:
        criteria_dict = json.loads(gene_criteria)
        criteria = GeneSelectionCriteria(**criteria_dict)

    # Auto-select gene
    gene_result = auto_select_gene_for_taxon(inclusion, email, criteria)
    gene = gene_result['gene']
    click.echo(f"Auto-selected gene: {gene} (score: {gene_result['total_score']:.2f})")

if auto_exclusion:
    # Intelligent exclusion selection
    exclusion_taxa = intelligent_exclusion_selection(
        taxid=inclusion,
        email=email,
        strategy=exclusion_strategy,
        phylo_distance_weight=exclusion_phylo_distance,
        max_exclusion_taxa=5
    )
    exclusion_taxids = [t['taxid'] for t in exclusion_taxa]
```

### Python API Integration

Users can import and use functions directly:

```python
from assay_design.gene_selection import (
    auto_select_gene_for_taxon,
    evaluate_gene_suitability,
    rank_candidate_genes,
    GeneSelectionCriteria,
    BACTERIA_GENES
)

from assay_design.hierarchical_search import intelligent_exclusion_selection
```

---

## Testing Strategy

### Unit Tests

**File:** `tests/test_gene_selection.py` (Phase 1 tests - 55 tests)

**Coverage:**
- `GeneSelectionCriteria` dataclass initialization
- Gene database completeness (all 4 domains)
- Domain detection logic
- `get_related_taxa()` with phylogenetic weighting
- `intelligent_exclusion_selection()` for all strategies

**File:** `tests/test_gene_selection_phase2.py` (Phase 2 tests - 20 tests)

**Coverage:**
- `evaluate_gene_suitability()` with various sequence counts
- Scoring algorithms (availability, copy number, HGT, length)
- `rank_candidate_genes()` sorting logic
- `auto_select_gene_for_taxon()` domain detection
- Integration with hierarchical search

### Mocking Strategy

**Mock NCBI API calls:**

```python
@pytest.fixture
def mock_sequences_abundant():
    """Mock abundant sequences (good availability)."""
    sequences = []
    for i in range(25):
        seq = SeqRecord(
            Seq("ATCG" * 400),  # 1600 bp
            id=f"seq_{i}",
            description=f"Mock sequence {i}"
        )
        sequences.append(seq)
    return sequences

@patch('assay_design.data_retrieval.fetch_gene_sequences')
def test_abundant_sequences(mock_fetch, mock_sequences_abundant):
    mock_fetch.return_value = mock_sequences_abundant
    result = evaluate_gene_suitability("562", "rpoB", "test@example.com")
    assert result['sequence_availability_score'] > 0.9
```

**Benefits:**
- Fast test execution (no network calls)
- Deterministic results
- Can test edge cases (no sequences, sparse data, etc.)

### Integration Tests

**Test hierarchical search integration:**

```python
def test_hierarchical_search_integration():
    """Test full workflow from gene selection to exclusion selection."""
    with patch('assay_design.data_retrieval.fetch_gene_sequences'), \
         patch('assay_design.data_retrieval.get_taxon_info'), \
         patch('assay_design.data_retrieval.get_related_taxa'):

        # Step 1: Auto-select gene
        gene_result = auto_select_gene_for_taxon("562", "test@example.com")

        # Step 2: Intelligent exclusion
        exclusion_taxa = intelligent_exclusion_selection(
            "562", "test@example.com", strategy="genus"
        )

        # Verify integration
        assert gene_result['gene'] in BACTERIA_GENES
        assert len(exclusion_taxa) > 0
```

---

## Performance Considerations

### Optimization Opportunities

#### 1. Parallel Gene Evaluation

**Current:** Sequential evaluation of genes
**Proposed:** Parallel evaluation using `concurrent.futures`

```python
from concurrent.futures import ThreadPoolExecutor

def auto_select_gene_for_taxon_parallel(taxid, email, criteria=None):
    """Parallel version of auto_select_gene_for_taxon."""

    with ThreadPoolExecutor(max_workers=5) as executor:
        futures = {
            executor.submit(evaluate_gene_suitability, taxid, gene, email, criteria): gene
            for gene in gene_database.keys()
        }

        results = []
        for future in futures:
            try:
                result = future.result(timeout=30)
                results.append(result)
            except Exception as e:
                logger.warning(f"Failed to evaluate {futures[future]}: {e}")

    return rank_candidate_genes(results)[0]
```

**Expected speedup:** 3-5x for 30 genes (network-bound)

#### 2. Caching NCBI Queries

**Current:** Fresh queries for every run
**Proposed:** Cache results with TTL

```python
from functools import lru_cache
import hashlib

@lru_cache(maxsize=1000)
def fetch_gene_sequences_cached(taxid, gene_name, email, max_records):
    """Cached version of fetch_gene_sequences."""
    return fetch_gene_sequences(taxid, gene_name, email, max_records)
```

**Benefits:**
- Faster repeated queries
- Reduced NCBI API load
- Better user experience during experimentation

#### 3. Early Termination

**Current:** Evaluate all genes even if excellent candidate found
**Proposed:** Stop if score > 9.0

```python
def auto_select_gene_smart(taxid, email, criteria=None):
    """Smart version with early termination."""

    results = []
    for gene in gene_database.keys():
        result = evaluate_gene_suitability(taxid, gene, email, criteria)
        results.append(result)

        # Early termination if excellent gene found
        if result['total_score'] > 9.0:
            logger.info(f"Excellent gene found: {gene} (score: {result['total_score']:.2f})")
            return result

    return rank_candidate_genes(results)[0]
```

### Memory Considerations

**Sequence storage:**
- Each SeqRecord ≈ 2-10 KB
- For 30 genes × 20 sequences = 600 sequences ≈ 1.2-6 MB
- Acceptable memory footprint

**Gene databases:**
- 4 databases × 30 genes × 200 bytes ≈ 24 KB
- Negligible memory usage

### Network Considerations

**NCBI API rate limits:**
- 3 requests/second without API key
- 10 requests/second with API key
- Recommend users provide API key for better performance

**Retry logic:**
```python
from Bio import Entrez
import time

def fetch_with_retry(query_func, max_retries=3):
    """Retry NCBI queries with exponential backoff."""
    for attempt in range(max_retries):
        try:
            return query_func()
        except Exception as e:
            if attempt == max_retries - 1:
                raise
            wait_time = 2 ** attempt
            logger.warning(f"Retry {attempt+1}/{max_retries} after {wait_time}s")
            time.sleep(wait_time)
```

---

## Future Enhancements

### 1. Machine Learning Gene Selection

Use historical success rates to train a model:
```python
from sklearn.ensemble import RandomForestClassifier

def ml_gene_selection(taxid, email):
    """ML-based gene selection using historical data."""
    # Features: sequence_count, avg_length, GC_content, etc.
    # Labels: success/failure from past assays
    # Train RandomForest to predict best gene
    pass
```

### 2. Multi-Objective Optimization

Balance multiple objectives (specificity, sensitivity, cost):
```python
from scipy.optimize import minimize

def optimize_gene_and_exclusion(taxid, email, objectives):
    """Multi-objective optimization for gene and exclusion selection."""
    # Objectives: [specificity, sensitivity, assay_cost, turnaround_time]
    # Use Pareto optimization to find optimal balance
    pass
```

### 3. Dynamic Gene Database Updates

Periodically update gene databases from NCBI:
```python
def update_gene_database(domain='Bacteria'):
    """Update gene database with latest NCBI data."""
    # Query NCBI for gene statistics
    # Update typical_length, sequence_counts
    # Re-rank genes by utility
    pass
```

### 4. Phylogenomic Integration

Use real phylogenetic trees instead of taxonomy:
```python
from Bio import Phylo

def phylogenomic_distance(taxon1, taxon2, tree):
    """Calculate true phylogenetic distance from tree."""
    # Use real phylogenetic tree (e.g., from PhyloT or TimeTree)
    # Calculate patristic distance
    # More accurate than taxonomy-based distance
    pass
```

---

## References

1. **Phylogenetic markers:** Yarza P, et al. (2014) "Uniting the classification of cultured and uncultured bacteria and archaea using 16S rRNA gene sequences" Nat Rev Microbiol.

2. **Gene selection for molecular diagnostics:** Rossen JW, et al. (2015) "Practical issues in implementing whole-genome-sequencing in routine diagnostic microbiology" Clin Microbiol Infect.

3. **HGT and gene stability:** Puigbo P, et al. (2014) "Genomes in turmoil: quantification of genome dynamics in prokaryote supergenomes" BMC Biol.

4. **NCBI Taxonomy:** Schoch CL, et al. (2020) "NCBI Taxonomy: a comprehensive update on curation, resources and tools" Database.

5. **Primer design principles:** Untergasser A, et al. (2012) "Primer3—new capabilities and interfaces" Nucleic Acids Res.

---

## Version History

- **v1.0.0** (2024): Initial implementation
  - Phase 1: Intelligent exclusion with phylogenetic weighting
  - Phase 2: Automated gene selection with multi-criteria scoring
  - Phase 3: Comprehensive documentation

---

## Contributors

For questions or contributions, please open an issue on GitHub.
