# Pipedream Enhancement Implementation Plan

**Version:** 1.0
**Date:** 2025-11-09
**Objective:** Implement intelligent exclusion selection and automated gene target selection for taxa-specific dPCR assay design

---

## Overview

This plan outlines the implementation of three major enhancement phases to make pipedream fully autonomous in designing taxa-specific primers and probes. The enhancements address two critical user challenges:

1. **Determining appropriate exclusion taxa** - Users know what to include (e.g., "Bennettella") but not what to exclude
2. **Selecting optimal gene targets** - Automated selection of single-copy, HGT-resistant genes suitable for dPCR assays

---

## Phase 1: Core Enhancements (Intelligent Exclusion & Gene Criteria)

**Goal:** Establish foundational infrastructure for intelligent taxonomic filtering and gene selection criteria

### 1.1 Enhanced `get_related_taxa()` with Phylogenetic Distance Weighting

**File:** `assay_design/data_retrieval.py`

**Current Limitation:**
- Returns siblings at a single taxonomic rank (genus OR family OR order)
- No prioritization based on phylogenetic distance
- No sequence availability checking

**Enhancement Specifications:**

```python
def get_related_taxa_weighted(
    taxid: str,
    email: str,
    max_results: int = 10,
    diversity_levels: List[str] = ["genus", "family", "order"]
) -> List[Dict[str, Any]]:
    """
    Get phylogenetically weighted exclusion taxa.

    Returns taxa sorted by:
    1. Phylogenetic proximity (closer relatives ranked higher)
    2. Sequence availability (taxa with >10 sequences preferred)
    3. Taxonomic diversity (span multiple levels)

    Returns:
        List of dicts with:
        - taxid: str
        - scientific_name: str
        - rank: str
        - phylogenetic_distance: int (lower = closer)
        - sequence_count: int
        - priority_score: float
    """
```

**Implementation Steps:**
1. Fetch full lineage for inclusion taxid using `get_taxon_info()`
2. For each diversity level (genus → family → order):
   - Get siblings at that rank
   - Query NCBI to estimate sequence availability
   - Calculate phylogenetic distance score
3. Compute priority score: `priority = (1/distance) * log(seq_count + 1) * rank_weight`
4. Return top N taxa spanning all diversity levels

**Testing:**
- Test with Vibrio mediterranei (taxid: 689) → should return other Vibrio species + related genera
- Verify sequence count estimates match NCBI
- Ensure diversity across taxonomic levels

---

### 1.2 Intelligent Exclusion Selection with Tiered Approach

**File:** `assay_design/data_retrieval.py`

**Purpose:** Automatically build a representative exclusion set without user input

**Enhancement Specifications:**

```python
def intelligent_exclusion_selection(
    inclusion_taxid: str,
    email: str,
    max_exclusion_taxa: int = 10,
    require_sequences: bool = True
) -> Dict[str, Any]:
    """
    Intelligently select exclusion taxa using tiered strategy.

    Tier 1 (50%): Closest relatives (genus-level siblings)
    Tier 2 (30%): Family-level relatives
    Tier 3 (20%): Order-level representatives

    Returns:
        {
            'exclusion_taxids': List[str],
            'taxa_info': List[Dict],  # Full metadata
            'selection_strategy': str,  # Description of selection
            'coverage_score': float  # 0-1, phylogenetic breadth
        }
    """
```

**Algorithm:**
1. Get inclusion taxon lineage and rank
2. Allocate exclusion slots by tier:
   - If inclusion is species → 5 genus siblings, 3 family relatives, 2 order relatives
   - If inclusion is genus → 7 family relatives, 3 order relatives
3. For each tier, call `get_related_taxa_weighted()`
4. Filter taxa with insufficient sequences (if `require_sequences=True`)
5. Calculate coverage score based on phylogenetic span
6. Return selected taxids with metadata

**Edge Cases:**
- Handle cases where insufficient taxa exist at a level (e.g., monotypic genera)
- Default to higher-level taxa if lower levels unavailable
- Minimum 3 exclusion taxa, maximum 20

**Testing:**
- Test with diverse taxa: bacteria (Vibrio), fungi (Candida), viruses (Influenza)
- Verify tier allocation percentages
- Test edge cases (monotypic genera, poorly sequenced clades)

---

### 1.3 Gene Selection Module (`gene_selection.py`)

**File:** `assay_design/gene_selection.py` (NEW)

**Purpose:** Centralized gene selection criteria and metadata

**Module Structure:**

```python
# assay_design/gene_selection.py

import logging
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass

logger = logging.getLogger(__name__)

@dataclass
class GeneMetadata:
    """Metadata for a candidate gene"""
    gene_name: str
    copy_number: str  # 'single', 'low', 'multi'
    hgt_propensity: str  # 'low', 'medium', 'high'
    evolutionary_rate: str  # 'slow', 'moderate', 'fast'
    typical_length_bp: Tuple[int, int]  # (min, max)
    functional_category: str  # 'housekeeping', 'marker', 'metabolic'
    description: str

class GeneSelectionCriteria:
    """
    Curated database of gene suitability for dPCR assays.

    Criteria for ideal dPCR target genes:
    1. Single or low copy number (avoid multi-copy bias)
    2. Low HGT propensity (phylogenetically informative)
    3. Moderate evolutionary rate (conserved but informative)
    4. Appropriate length (500-2000 bp for primer design flexibility)
    5. Housekeeping or marker gene (universally present)
    """

    # Bacteria single-copy housekeeping genes
    BACTERIA_SINGLE_COPY = {
        'rpoB': GeneMetadata(
            gene_name='rpoB',
            copy_number='single',
            hgt_propensity='low',
            evolutionary_rate='moderate',
            typical_length_bp=(3500, 4200),
            functional_category='housekeeping',
            description='RNA polymerase beta subunit - highly conserved, single copy'
        ),
        'gyrB': GeneMetadata(
            gene_name='gyrB',
            copy_number='single',
            hgt_propensity='low',
            evolutionary_rate='moderate',
            typical_length_bp=(1800, 2400),
            functional_category='housekeeping',
            description='DNA gyrase subunit B - single copy, HGT-resistant'
        ),
        'recA': GeneMetadata(
            gene_name='recA',
            copy_number='single',
            hgt_propensity='low',
            evolutionary_rate='moderate',
            typical_length_bp=(1000, 1200),
            functional_category='housekeeping',
            description='Recombinase A - universal bacterial gene'
        ),
        'dnaK': GeneMetadata(
            gene_name='dnaK',
            copy_number='single',
            hgt_propensity='low',
            evolutionary_rate='moderate',
            typical_length_bp=(1800, 2000),
            functional_category='housekeeping',
            description='Heat shock protein 70 - universally conserved'
        ),
        'groEL': GeneMetadata(
            gene_name='groEL',
            copy_number='single',
            hgt_propensity='low',
            evolutionary_rate='moderate',
            typical_length_bp=(1600, 1650),
            functional_category='housekeeping',
            description='60 kDa chaperonin - single copy chaperone'
        ),
    }

    # Bacterial marker genes (may be multi-copy but widely used)
    BACTERIA_MARKER = {
        '16S rRNA': GeneMetadata(
            gene_name='16S ribosomal RNA',
            copy_number='multi',  # Typically 1-15 copies
            hgt_propensity='low',
            evolutionary_rate='slow',
            typical_length_bp=(1400, 1600),
            functional_category='marker',
            description='Small subunit ribosomal RNA - universal phylogenetic marker'
        ),
        'rpoD': GeneMetadata(
            gene_name='rpoD',
            copy_number='single',
            hgt_propensity='low',
            evolutionary_rate='moderate',
            typical_length_bp=(1900, 2100),
            functional_category='housekeeping',
            description='RNA polymerase sigma factor - single copy'
        ),
    }

    # Archaea housekeeping genes
    ARCHAEA_SINGLE_COPY = {
        'rps3': GeneMetadata(
            gene_name='rps3',
            copy_number='single',
            hgt_propensity='low',
            evolutionary_rate='slow',
            typical_length_bp=(600, 700),
            functional_category='housekeeping',
            description='Ribosomal protein S3 - conserved archaeal marker'
        ),
        'rpl2': GeneMetadata(
            gene_name='rpl2',
            copy_number='single',
            hgt_propensity='low',
            evolutionary_rate='slow',
            typical_length_bp=(700, 900),
            functional_category='housekeeping',
            description='Ribosomal protein L2 - conserved archaeal marker'
        ),
    }

    # Eukaryote (fungi) housekeeping genes
    FUNGI_SINGLE_COPY = {
        'RPB1': GeneMetadata(
            gene_name='RPB1',
            copy_number='single',
            hgt_propensity='low',
            evolutionary_rate='moderate',
            typical_length_bp=(3500, 4500),
            functional_category='housekeeping',
            description='RNA polymerase II largest subunit - single copy'
        ),
        'RPB2': GeneMetadata(
            gene_name='RPB2',
            copy_number='single',
            hgt_propensity='low',
            evolutionary_rate='moderate',
            typical_length_bp=(3000, 4000),
            functional_category='housekeeping',
            description='RNA polymerase II second largest subunit - single copy'
        ),
        'TEF1': GeneMetadata(
            gene_name='TEF1',
            copy_number='single',
            hgt_propensity='low',
            evolutionary_rate='moderate',
            typical_length_bp=(1200, 1500),
            functional_category='housekeeping',
            description='Translation elongation factor 1-alpha - single copy'
        ),
    }

    # Fungal marker genes
    FUNGI_MARKER = {
        'ITS': GeneMetadata(
            gene_name='ITS',
            copy_number='multi',
            hgt_propensity='low',
            evolutionary_rate='fast',
            typical_length_bp=(500, 700),
            functional_category='marker',
            description='Internal transcribed spacer - standard fungal barcode'
        ),
        '18S rRNA': GeneMetadata(
            gene_name='18S ribosomal RNA',
            copy_number='multi',
            hgt_propensity='low',
            evolutionary_rate='slow',
            typical_length_bp=(1700, 1900),
            functional_category='marker',
            description='Small subunit ribosomal RNA - eukaryotic marker'
        ),
    }

    @classmethod
    def get_genes_for_domain(cls, domain: str, prefer_single_copy: bool = True) -> Dict[str, GeneMetadata]:
        """
        Get recommended genes for a taxonomic domain.

        Args:
            domain: 'bacteria', 'archaea', 'fungi', 'eukaryote'
            prefer_single_copy: If True, prioritize single-copy genes

        Returns:
            Dict of gene_name -> GeneMetadata
        """
        domain = domain.lower()

        if domain == 'bacteria':
            if prefer_single_copy:
                return cls.BACTERIA_SINGLE_COPY
            else:
                return {**cls.BACTERIA_SINGLE_COPY, **cls.BACTERIA_MARKER}
        elif domain == 'archaea':
            return cls.ARCHAEA_SINGLE_COPY
        elif domain in ['fungi', 'eukaryote']:
            if prefer_single_copy:
                return cls.FUNGI_SINGLE_COPY
            else:
                return {**cls.FUNGI_SINGLE_COPY, **cls.FUNGI_MARKER}
        else:
            logger.warning(f"Unknown domain '{domain}', defaulting to bacteria")
            return cls.BACTERIA_SINGLE_COPY

    @classmethod
    def rank_genes_by_criteria(cls, genes: Dict[str, GeneMetadata], criteria_weights: Dict[str, float] = None) -> List[Tuple[str, float]]:
        """
        Rank genes by suitability for dPCR assays.

        Args:
            genes: Dict of gene_name -> GeneMetadata
            criteria_weights: Custom weights for scoring criteria

        Returns:
            List of (gene_name, score) tuples, sorted by score (descending)
        """
        if criteria_weights is None:
            criteria_weights = {
                'single_copy': 0.4,
                'hgt_resistance': 0.3,
                'evolutionary_rate': 0.2,
                'length': 0.1
            }

        scored_genes = []

        for gene_name, metadata in genes.items():
            score = 0.0

            # Copy number score
            copy_score = {'single': 1.0, 'low': 0.7, 'multi': 0.3}.get(metadata.copy_number, 0.5)
            score += copy_score * criteria_weights['single_copy']

            # HGT resistance score
            hgt_score = {'low': 1.0, 'medium': 0.5, 'high': 0.1}.get(metadata.hgt_propensity, 0.5)
            score += hgt_score * criteria_weights['hgt_resistance']

            # Evolutionary rate score (moderate is best for primers)
            evo_score = {'slow': 0.8, 'moderate': 1.0, 'fast': 0.5}.get(metadata.evolutionary_rate, 0.5)
            score += evo_score * criteria_weights['evolutionary_rate']

            # Length score (prefer 1000-2000 bp)
            min_len, max_len = metadata.typical_length_bp
            ideal_length = 1500
            avg_length = (min_len + max_len) / 2
            length_score = 1.0 - min(abs(avg_length - ideal_length) / ideal_length, 1.0)
            score += length_score * criteria_weights['length']

            scored_genes.append((gene_name, score))

        return sorted(scored_genes, key=lambda x: x[1], reverse=True)
```

**Testing:**
- Verify gene metadata completeness
- Test ranking algorithm with different weight configurations
- Validate domain-specific gene retrieval

---

### 1.4 Unit Tests for Phase 1

**File:** `tests/test_gene_selection.py` (NEW)

**Test Coverage:**
- `test_get_related_taxa_weighted()` - verify phylogenetic distance calculation
- `test_intelligent_exclusion_selection()` - test tier allocation
- `test_gene_metadata_completeness()` - ensure all genes have required fields
- `test_gene_ranking()` - verify scoring algorithm
- `test_domain_detection()` - test automatic domain inference from taxid

---

## Phase 2: Automated Gene Selection (Ranking & Integration)

**Goal:** Enable fully automated gene target selection based on sequence availability and suitability

### 2.1 Evaluate Gene Suitability Function

**File:** `assay_design/gene_selection.py`

**Purpose:** Dynamically evaluate gene suitability by querying NCBI for sequence availability

**Enhancement Specifications:**

```python
def evaluate_gene_suitability(
    taxid: str,
    gene_name: str,
    email: str,
    metadata: Optional[GeneMetadata] = None
) -> Dict[str, Any]:
    """
    Evaluate a specific gene for assay design suitability.

    Scoring criteria:
    1. Sequence availability (more sequences = more robust design)
    2. Copy number (single copy preferred)
    3. HGT resistance (low HGT preferred)
    4. Sequence length distribution (ideal: 1000-2000 bp)

    Args:
        taxid: NCBI taxonomy ID
        gene_name: Gene to evaluate
        email: NCBI email
        metadata: Pre-loaded gene metadata (if available)

    Returns:
        {
            'gene': str,
            'sequence_count': int,
            'avg_sequence_length': float,
            'copy_number_score': float,  # 0-1
            'hgt_resistance_score': float,  # 0-1
            'sequence_availability_score': float,  # 0-1
            'length_appropriateness_score': float,  # 0-1
            'overall_score': float,  # weighted average
            'recommendation': str  # 'excellent', 'good', 'fair', 'poor'
        }
    """
```

**Implementation:**
1. Use `fetch_gene_sequences()` to query NCBI for the gene in the given taxon
2. Count sequences and calculate average length
3. Retrieve gene metadata from `GeneSelectionCriteria` if available
4. Calculate component scores:
   - `sequence_availability_score = min(1.0, sequence_count / 20)`
   - `copy_number_score`, `hgt_resistance_score` from metadata
   - `length_score = 1.0 - abs(avg_length - 1500) / 1500`
5. Compute overall score (weighted average)
6. Return recommendation tier

**Testing:**
- Test with well-sequenced genes (16S rRNA, rpoB)
- Test with poorly sequenced genes
- Verify score calculations

---

### 2.2 Rank Candidate Genes Function

**File:** `assay_design/gene_selection.py`

**Purpose:** Automatically rank multiple candidate genes for a given taxon

**Enhancement Specifications:**

```python
def rank_candidate_genes(
    taxid: str,
    email: str,
    domain: Optional[str] = None,
    prefer_single_copy: bool = True,
    max_genes_to_test: int = 5,
    timeout_per_gene: int = 10
) -> List[Tuple[str, Dict[str, Any]]]:
    """
    Rank all candidate genes for a taxon.

    Args:
        taxid: NCBI taxonomy ID
        email: NCBI email
        domain: Taxonomic domain (auto-detected if None)
        prefer_single_copy: Prioritize single-copy genes
        max_genes_to_test: Maximum genes to evaluate (top-ranked by criteria)
        timeout_per_gene: Max seconds to spend querying each gene

    Returns:
        List of (gene_name, evaluation_results) sorted by overall_score
    """
```

**Algorithm:**
1. Auto-detect domain from taxid lineage if not provided
2. Get candidate genes from `GeneSelectionCriteria.get_genes_for_domain()`
3. Pre-rank genes using `GeneSelectionCriteria.rank_genes_by_criteria()`
4. For top N genes, call `evaluate_gene_suitability()`
5. Re-rank based on empirical sequence availability
6. Return sorted list

**Testing:**
- Test with Vibrio mediterranei (bacteria)
- Test with Candida albicans (fungi)
- Verify timeout enforcement

---

### 2.3 Integration into Hierarchical Search

**File:** `assay_design/hierarchical_search.py`

**Purpose:** Add automatic gene selection to the hierarchical search workflow

**Modifications:**

```python
def hierarchical_marker_search(
    inclusion_taxid: str,
    email: str,
    gene_name: Optional[str] = None,
    auto_select_gene: bool = False,  # NEW
    gene_selection_criteria: str = 'balanced',  # NEW: 'single_copy', 'balanced', 'marker'
    # ... existing parameters
) -> Dict[str, Any]:
    """
    Enhanced hierarchical search with automatic gene selection.

    If auto_select_gene=True and gene_name=None:
        1. Detect domain from inclusion_taxid
        2. Rank candidate genes using rank_candidate_genes()
        3. Try top 3 genes in order until markers found
        4. Fall back to next gene if first fails
    """
```

**Algorithm:**
1. Check if `auto_select_gene=True` and `gene_name=None`
2. If yes, call `rank_candidate_genes()`
3. Iterate through top 3 genes:
   - Fetch sequences for gene
   - Run marker identification
   - If markers found with sufficient quality, break
   - Otherwise, try next gene
4. Log gene selection results

**Testing:**
- Test with taxa having poor 16S coverage
- Verify fallback to alternative genes
- Test performance impact

---

### 2.4 CLI Integration

**File:** `assay_design/cli.py`

**Purpose:** Add command-line flags for automated features

**New Arguments:**

```python
parser.add_argument(
    '--auto-gene',
    action='store_true',
    help='Automatically select optimal gene target for the taxon'
)

parser.add_argument(
    '--gene-selection-criteria',
    type=str,
    choices=['single_copy', 'balanced', 'marker'],
    default='balanced',
    help='Gene selection strategy: single_copy (prefer single-copy genes), '
         'balanced (balance copy number and availability), '
         'marker (include marker genes like 16S rRNA)'
)

parser.add_argument(
    '--exclusion-strategy',
    type=str,
    choices=['intelligent', 'manual', 'siblings', 'family'],
    default='intelligent',
    help='Exclusion taxa selection strategy'
)

parser.add_argument(
    '--max-exclusion-taxa',
    type=int,
    default=10,
    help='Maximum number of exclusion taxa to use (default: 10)'
)
```

**Main Function Updates:**

Integrate new functions into the main CLI workflow:

```python
def main():
    # ... existing argument parsing

    # Automatic exclusion selection
    if args.exclusion_strategy == 'intelligent':
        exclusion_result = intelligent_exclusion_selection(
            inclusion_taxid=args.inclusion,
            email=args.email,
            max_exclusion_taxa=args.max_exclusion_taxa
        )
        exclusion_taxids = exclusion_result['exclusion_taxids']
        logger.info(f"Intelligently selected {len(exclusion_taxids)} exclusion taxa")
        logger.info(f"Coverage score: {exclusion_result['coverage_score']:.2f}")

    # Automatic gene selection
    if args.auto_gene and not args.gene:
        ranked_genes = rank_candidate_genes(
            taxid=args.inclusion,
            email=args.email,
            prefer_single_copy=(args.gene_selection_criteria == 'single_copy')
        )

        if ranked_genes:
            selected_gene = ranked_genes[0][0]
            logger.info(f"Auto-selected gene: {selected_gene}")
            logger.info(f"Gene score: {ranked_genes[0][1]['overall_score']:.2f}")
            args.gene = selected_gene
```

**Testing:**
- Test full CLI workflow with `--auto-gene --exclusion-strategy intelligent`
- Verify output directory structure
- Test error handling for edge cases

---

### 2.5 Unit Tests for Phase 2

**File:** `tests/test_automated_gene_selection.py` (NEW)

**Test Coverage:**
- `test_evaluate_gene_suitability()` - verify scoring
- `test_rank_candidate_genes()` - test ranking algorithm
- `test_hierarchical_search_auto_gene()` - integration test
- `test_cli_auto_gene_flag()` - CLI integration
- `test_multi_gene_fallback()` - verify fallback behavior

---

## Phase 3: Advanced Optimizations

**Goal:** Maximize computational efficiency and provide rich user experience

### 3.1 Pre-computed Gene Database

**File:** `assay_design/gene_database.py` (NEW)

**Purpose:** Cache frequently accessed gene metadata to avoid redundant NCBI queries

**Enhancement Specifications:**

```python
class GeneDatabase:
    """
    Persistent cache of gene metadata across taxa.

    Stores:
    - Gene names and aliases
    - Typical sequence lengths
    - Copy number estimates
    - Sequence availability statistics
    - Last updated timestamps
    """

    def __init__(self, cache_file: str = "~/.assay_design/gene_cache.json"):
        """Initialize database with optional cache file"""

    def get_gene_info(self, gene_name: str, taxid: Optional[str] = None) -> Optional[Dict]:
        """Retrieve cached gene information"""

    def update_gene_info(self, gene_name: str, taxid: str, info: Dict):
        """Update cache with new gene information"""

    def get_common_aliases(self, gene_name: str) -> List[str]:
        """Get common aliases for a gene (e.g., '16S rRNA' -> ['16S', 'SSU rRNA'])"""
```

**Data Structure:**

```json
{
  "genes": {
    "rpoB": {
      "full_name": "RNA polymerase beta subunit",
      "aliases": ["rpoB", "beta subunit"],
      "copy_number": "single",
      "typical_length": [3500, 4200],
      "taxa_stats": {
        "562": {  // E. coli taxid
          "sequence_count": 1547,
          "avg_length": 3900,
          "last_updated": "2025-11-09"
        }
      }
    }
  }
}
```

**Testing:**
- Test cache read/write
- Test cache expiration (30 days)
- Verify thread safety for parallel queries

---

### 3.2 Adaptive Sequence Sampling

**File:** `assay_design/target_identification.py`

**Purpose:** Intelligently sample sequences to maximize phylogenetic diversity while minimizing computation

**Enhancement Specifications:**

```python
def adaptive_sequence_sampling(
    sequences: List[SeqRecord],
    target_count: int = 20,
    diversity_weight: float = 0.6,
    quality_weight: float = 0.4
) -> List[SeqRecord]:
    """
    Sample sequences optimizing for diversity and quality.

    Algorithm:
    1. Cluster sequences by k-mer similarity (fast clustering)
    2. Select representatives from each cluster
    3. Weight by sequence quality (length, N-content)
    4. Return balanced sample

    Args:
        sequences: All available sequences
        target_count: Desired number of sequences
        diversity_weight: Weight for phylogenetic diversity (0-1)
        quality_weight: Weight for sequence quality (0-1)

    Returns:
        Sampled sequences maximizing diversity and quality
    """
```

**Implementation:**
1. Use MinHash LSH to cluster sequences rapidly
2. Calculate cluster sizes
3. Sample proportionally from each cluster
4. Within clusters, prioritize by sequence quality metrics
5. Return balanced sample

**Testing:**
- Test with highly redundant sequence sets
- Verify diversity improvement over random sampling
- Benchmark performance vs. full dataset

---

### 3.3 Enhanced Early Termination

**File:** `assay_design/hierarchical_search.py`

**Purpose:** Stop search early when high-quality markers are found

**Enhancement Specifications:**

```python
def assess_marker_quality(marker_info: Dict, primers_info: Dict) -> str:
    """
    Assess overall marker quality.

    Returns: 'excellent', 'good', 'fair', 'poor'

    Criteria for 'excellent':
    - Conservation ≥ 0.95
    - Specificity ≥ 0.99
    - Primer Tm difference ≤ 2°C
    - No primer dimers
    - GC content 40-60%
    """
```

Add early termination logic:

```python
# In hierarchical_marker_search()
if assess_marker_quality(marker_info, primers_info) == 'excellent':
    logger.info("Excellent marker found, terminating search early")
    return result
```

**Testing:**
- Verify termination with high-quality markers
- Ensure continuation with poor-quality markers
- Test timeout behavior

---

### 3.4 Documentation and Examples

**Files:**
- `docs/user_guide.md` - Comprehensive user guide
- `docs/api_reference.md` - API documentation
- `examples/auto_assay_design.py` - Example scripts

**Content:**

Create detailed documentation covering:
1. Automatic exclusion selection
2. Automatic gene selection
3. Customizing gene selection criteria
4. Interpreting quality scores
5. Troubleshooting common issues

Example usage:

```bash
# Fully automated assay design for Bennettella
assay-design \
    --inclusion 123456 \
    --email user@example.com \
    --output-dir ./bennettella_assay \
    --auto-gene \
    --gene-selection-criteria single_copy \
    --exclusion-strategy intelligent \
    --max-exclusion-taxa 10
```

---

### 3.5 Integration Testing and Benchmarking

**File:** `tests/integration/test_full_workflow.py` (NEW)

**Test Cases:**
1. End-to-end test with Vibrio mediterranei
2. End-to-end test with poorly sequenced organism
3. Performance benchmarking (time and memory)
4. Stress testing with large exclusion sets

**Benchmarks:**
- Full workflow should complete in <5 minutes for typical cases
- Memory usage should stay <4GB
- Cache hit rate should be >70% after first run

---

## Success Criteria

### Phase 1 Success Metrics
- [ ] `intelligent_exclusion_selection()` returns 5-10 diverse taxa
- [ ] Phylogenetic coverage score >0.7 for typical cases
- [ ] Gene metadata database includes ≥20 genes across domains
- [ ] All unit tests pass

### Phase 2 Success Metrics
- [ ] `rank_candidate_genes()` successfully ranks ≥3 genes per taxon
- [ ] Auto-gene selection works for 90% of well-sequenced taxa
- [ ] Multi-gene fallback successfully finds markers when first gene fails
- [ ] CLI integration works seamlessly
- [ ] All unit tests pass

### Phase 3 Success Metrics
- [ ] Gene database cache reduces NCBI queries by ≥50%
- [ ] Adaptive sampling reduces computation time by ≥30%
- [ ] Early termination triggers for high-quality markers
- [ ] Complete documentation published
- [ ] All integration tests pass
- [ ] Performance benchmarks meet targets

---

## Implementation Timeline

### Phase 1: 3-4 days
- Day 1: Enhanced `get_related_taxa()` and `intelligent_exclusion_selection()`
- Day 2: `gene_selection.py` module and gene metadata database
- Day 3: Unit tests and debugging
- Day 4: Documentation and code review

### Phase 2: 3-4 days
- Day 1: `evaluate_gene_suitability()` and `rank_candidate_genes()`
- Day 2: Hierarchical search integration
- Day 3: CLI integration
- Day 4: Unit tests and debugging

### Phase 3: 2-3 days
- Day 1: Gene database and adaptive sampling
- Day 2: Enhanced termination and documentation
- Day 3: Integration testing and benchmarking

**Total Estimated Time:** 8-11 days

---

## Risk Mitigation

### Risk 1: NCBI API Rate Limiting
**Mitigation:** Implement exponential backoff, respect NCBI rate limits (3 requests/second), cache aggressively

### Risk 2: Insufficient Sequence Data
**Mitigation:** Multi-gene fallback strategy, graceful degradation, clear user warnings

### Risk 3: Performance Degradation
**Mitigation:** Parallel processing, LSH optimization, early termination, adaptive sampling

### Risk 4: Edge Cases (Monotypic Taxa)
**Mitigation:** Comprehensive error handling, fallback to higher taxonomic levels, user notifications

---

## Dependencies

### Required Libraries (already installed)
- biopython ≥ 1.79
- numpy ≥ 1.20.0
- click ≥ 8.0

### Optional Libraries (for enhanced features)
- datasketch - MinHash LSH (already installed)
- pybloom-live - Bloom filters (already installed)

### External Dependencies
- NCBI Entrez API (requires email)
- Internet connection for NCBI queries

---

## Maintenance and Future Enhancements

### Ongoing Maintenance
- Update gene metadata database quarterly
- Monitor NCBI API changes
- Refresh cached data periodically

### Future Enhancement Ideas
- Support for viral assays
- Integration with Primer3 for advanced primer design
- Web interface for non-command-line users
- Multi-species assays (detect multiple taxa)
- Machine learning for gene selection optimization

---

## Appendix: File Structure After Implementation

```
pipedream/
├── assay_design/
│   ├── __init__.py
│   ├── cli.py                      # Enhanced with new CLI flags
│   ├── data_retrieval.py           # Enhanced with intelligent_exclusion_selection()
│   ├── gene_selection.py           # NEW - Gene selection module
│   ├── gene_database.py            # NEW - Pre-computed gene cache
│   ├── target_identification.py    # Enhanced with adaptive_sequence_sampling()
│   ├── hierarchical_search.py      # Enhanced with auto-gene selection
│   ├── primer_probe_design.py      # (unchanged)
│   ├── specificity_validation.py   # (unchanged)
│   ├── kmer_analysis.py            # (unchanged)
│   ├── lsh_sequence_search.py      # (unchanged)
│   ├── cache_manager.py            # (unchanged)
│   ├── visualization.py            # (unchanged)
│   └── io_utils.py                 # (unchanged)
├── tests/
│   ├── test_gene_selection.py      # NEW - Phase 1 tests
│   ├── test_automated_gene_selection.py  # NEW - Phase 2 tests
│   ├── integration/
│   │   └── test_full_workflow.py   # NEW - Integration tests
│   └── (existing test files)
├── docs/
│   ├── IMPLEMENTATION_PLAN.md      # This document
│   ├── user_guide.md               # NEW - User guide
│   ├── api_reference.md            # NEW - API docs
│   └── index.md
├── examples/
│   └── auto_assay_design.py        # NEW - Example scripts
└── (existing files)
```

---

**End of Implementation Plan**
