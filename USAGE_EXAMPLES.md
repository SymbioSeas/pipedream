# Assay Design Usage Examples

This document provides comprehensive, real-world examples of using the assay_design package for various scenarios.

## Table of Contents

1. [Quick Start Examples](#quick-start-examples)
2. [Automated Gene Selection Examples](#automated-gene-selection-examples)
3. [Intelligent Exclusion Selection Examples](#intelligent-exclusion-selection-examples)
4. [Custom Criteria Examples](#custom-criteria-examples)
5. [Multi-Domain Examples](#multi-domain-examples)
6. [Python API Examples](#python-api-examples)
7. [Troubleshooting Common Issues](#troubleshooting-common-issues)

---

## Quick Start Examples

### Example 1: Fully Automated Species-Specific Assay

Design a species-specific assay for *Vibrio mediterranei* with everything automated:

```bash
assay-design --inclusion 689 --email your.email@example.com --auto-gene --auto-exclusion
```

**What this does:**
- Automatically selects the best gene marker (likely rpoB or gyrB)
- Automatically identifies sibling Vibrio species for exclusion
- Uses phylogenetic distance weighting to prioritize closest relatives
- Designs primers specific to V. mediterranei

**Expected output:**
```
Selected gene: rpoB (score: 8.75/10)
Exclusion taxa: 5 sibling species in genus Vibrio
Primers designed: Forward: ATCG..., Reverse: GCTA...
```

### Example 2: View Gene Suggestions Before Deciding

Check which genes are available for your target organism:

```bash
assay-design --inclusion 562 --email your.email@example.com --suggest-genes
```

**What this does:**
- Evaluates all candidate genes for E. coli (taxid 562)
- Shows suitability scores for each gene
- Displays sequence availability and other metrics

**Expected output:**
```
Gene Suitability Report for taxid 562 (Escherichia coli):

Rank  Gene        Score  Sequences  Copy  HGT-Res  Avg Length
----  ----        -----  ---------  ----  -------  ----------
1     rpoB        9.2    147        Yes   Yes      3900 bp
2     gyrB        8.8    132        Yes   Yes      2100 bp
3     recA        8.5    98         Yes   Yes      1050 bp
4     16S rRNA    7.2    856        No    Yes      1550 bp
5     dnaK        7.0    76         Yes   Yes      1800 bp
```

---

## Automated Gene Selection Examples

### Example 3: Auto-Select Gene with Custom Criteria

Automatically select a gene, but use custom criteria requiring at least 30 sequences:

```bash
assay-design --inclusion 689 --email your.email@example.com --auto-gene \
  --gene-criteria '{"min_sequence_count": 30, "ideal_length_range": [1000, 2000]}' \
  --auto-exclusion
```

**When to use this:**
- You need high sequence availability for robust alignment
- You want longer amplicons (1000-2000 bp)
- You're designing assays for well-characterized organisms

### Example 4: Auto-Select Single-Copy Genes Only

Force selection of only single-copy genes (exclude rRNA genes):

```bash
assay-design --inclusion 562 --email your.email@example.com --auto-gene \
  --gene-criteria '{"single_copy_preferred": true, "min_sequence_count": 20}' \
  --auto-exclusion
```

**When to use this:**
- Quantitative PCR applications where copy number matters
- Avoiding amplification bias from multi-copy genes
- Phylogenetic studies requiring single-copy markers

### Example 5: Manual Gene Selection

If you already know which gene to use:

```bash
assay-design --inclusion 689 --gene "rpoB" --email your.email@example.com --auto-exclusion
```

**When to use this:**
- You have prior knowledge of the best marker for your organism
- You're following a published protocol
- You're designing primers for a specific gene of interest

---

## Intelligent Exclusion Selection Examples

### Example 6: Genus-Level Exclusion (Species-Specific Assay)

Design a species-specific assay that excludes other species in the same genus:

```bash
assay-design --inclusion 689 --email your.email@example.com \
  --auto-gene \
  --exclusion-strategy genus \
  --exclusion-phylo-distance 0.8
```

**Use case:** Detecting *Vibrio mediterranei* but NOT other Vibrio species
**Exclusion taxa:** V. cholerae, V. parahaemolyticus, V. vulnificus, etc.

### Example 7: Family-Level Exclusion (Genus-Specific Assay)

Design a genus-specific assay that excludes other genera in the same family:

```bash
assay-design --inclusion 662 --email your.email@example.com \
  --auto-gene \
  --exclusion-strategy family \
  --exclusion-phylo-distance 0.75
```

**Use case:** Detecting any *Vibrio* species but NOT other Vibrionaceae genera
**Exclusion taxa:** Photobacterium, Aliivibrio, Grimontia, etc.

### Example 8: Order-Level Exclusion (Family-Specific Assay)

Design a family-specific assay with broader exclusion:

```bash
assay-design --inclusion 641 --email your.email@example.com \
  --auto-gene \
  --exclusion-strategy order \
  --exclusion-phylo-distance 0.7
```

**Use case:** Detecting Vibrionaceae family but NOT other families in Vibrionales
**Exclusion taxa:** Other families in the Vibrionales order

### Example 9: Manual Exclusion (Custom Taxa)

Manually specify which taxa to exclude:

```bash
assay-design --inclusion 689 \
  --exclusion 717,670,672 \
  --gene "rpoB" \
  --email your.email@example.com
```

**When to use this:**
- You have specific knowledge of problematic taxa
- You're targeting an unusual phylogenetic relationship
- Published literature suggests specific exclusion taxa

---

## Custom Criteria Examples

### Example 10: High Sequence Availability Required

For well-studied organisms with abundant sequence data:

```bash
assay-design --inclusion 562 --email your.email@example.com --auto-gene \
  --gene-criteria '{"min_sequence_count": 50, "ideal_length_range": [800, 1500]}' \
  --auto-exclusion
```

**Typical organisms:** E. coli, Salmonella, Staphylococcus aureus

### Example 11: Relaxed Criteria for Rare Organisms

For poorly characterized organisms with limited sequence data:

```bash
assay-design --inclusion 12345 --email your.email@example.com --auto-gene \
  --gene-criteria '{"min_sequence_count": 5, "ideal_length_range": [500, 3000]}' \
  --auto-exclusion
```

**When to use this:**
- Novel or rare organisms
- Environmental isolates with limited genomic data
- Preliminary assay design before more data becomes available

### Example 12: Long Amplicons for Sequencing

When you need longer amplicons for Sanger sequencing or other applications:

```bash
assay-design --inclusion 689 --email your.email@example.com --auto-gene \
  --gene-criteria '{"ideal_length_range": [1500, 2500], "min_sequence_count": 15}' \
  --auto-exclusion
```

**Use case:** Sanger sequencing, MinION sequencing, or phylogenetic analysis

---

## Multi-Domain Examples

### Example 13: Bacterial Assay (Most Common)

```bash
# E. coli detection
assay-design --inclusion 562 --email your.email@example.com --auto-gene --auto-exclusion

# Salmonella detection
assay-design --inclusion 590 --email your.email@example.com --auto-gene --auto-exclusion
```

**Automatic gene database:** BACTERIA_GENES (30 genes including rpoB, gyrB, recA, 16S rRNA)

### Example 14: Archaeal Assay

```bash
# Methanococcus maripaludis
assay-design --inclusion 267377 --email your.email@example.com --auto-gene --auto-exclusion
```

**Automatic gene database:** ARCHAEA_GENES (20 genes including rpoB, EF-2, 16S rRNA)

### Example 15: Eukaryotic Assay

```bash
# Saccharomyces cerevisiae
assay-design --inclusion 4932 --email your.email@example.com --auto-gene --auto-exclusion
```

**Automatic gene database:** EUKARYOTA_GENES (20 genes including 18S rRNA, COI, ITS)

### Example 16: Viral Assay

```bash
# SARS-CoV-2
assay-design --inclusion 2697049 --email your.email@example.com --auto-gene --auto-exclusion
```

**Automatic gene database:** VIRUS_GENES (15 genes including RdRp, capsid, polymerase)

---

## Python API Examples

### Example 17: Full Automated Workflow in Python

```python
from assay_design.gene_selection import auto_select_gene_for_taxon
from assay_design.hierarchical_search import intelligent_exclusion_selection
from assay_design.data_retrieval import fetch_gene_sequences
from assay_design.target_identification import find_conserved_marker
from assay_design.primer_design import design_primers

def design_automated_assay(taxid, email):
    """Fully automated assay design workflow."""

    # Step 1: Auto-select best gene
    print(f"Step 1: Selecting best gene for taxid {taxid}...")
    gene_result = auto_select_gene_for_taxon(taxid, email)
    gene = gene_result['gene']
    print(f"  Selected: {gene} (score: {gene_result['total_score']:.2f})")

    # Step 2: Intelligent exclusion selection
    print(f"Step 2: Selecting exclusion taxa...")
    exclusion_taxa = intelligent_exclusion_selection(
        taxid=taxid,
        email=email,
        strategy="genus",
        phylo_distance_weight=0.8,
        max_exclusion_taxa=5
    )
    exclusion_taxids = [t['taxid'] for t in exclusion_taxa]
    print(f"  Found {len(exclusion_taxids)} exclusion taxa")

    # Step 3: Fetch sequences
    print(f"Step 3: Fetching sequences...")
    inclusion_seqs = fetch_gene_sequences(taxid, gene, email, max_records=10)
    exclusion_seqs = []
    for ex_taxid in exclusion_taxids:
        seqs = fetch_gene_sequences(ex_taxid, gene, email, max_records=5)
        exclusion_seqs.extend(seqs)
    print(f"  Inclusion: {len(inclusion_seqs)} sequences")
    print(f"  Exclusion: {len(exclusion_seqs)} sequences")

    # Step 4: Find conserved markers
    print(f"Step 4: Finding conserved markers...")
    marker_info = find_conserved_marker(inclusion_seqs, exclusion_seqs)

    # Step 5: Design primers
    print(f"Step 5: Designing primers...")
    primers = design_primers(marker_info)

    return {
        'gene': gene,
        'exclusion_taxa': exclusion_taxa,
        'primers': primers
    }

# Usage
result = design_automated_assay(taxid="689", email="your.email@example.com")
print(f"\nFinal Result:")
print(f"Gene: {result['gene']}")
print(f"Primers: {result['primers']}")
```

### Example 18: Evaluate Multiple Genes in Python

```python
from assay_design.gene_selection import (
    evaluate_gene_suitability,
    rank_candidate_genes,
    BACTERIA_GENES
)

def compare_all_genes(taxid, email):
    """Evaluate and compare all candidate genes."""

    results = []
    for gene_name in BACTERIA_GENES.keys():
        try:
            result = evaluate_gene_suitability(taxid, gene_name, email)
            results.append(result)
            print(f"{gene_name}: {result['total_score']:.2f} "
                  f"({result['sequence_count']} sequences)")
        except Exception as e:
            print(f"{gene_name}: Failed ({str(e)})")

    # Rank genes
    ranked = rank_candidate_genes(results)

    print("\n=== Top 5 Genes ===")
    for i, gene in enumerate(ranked[:5], 1):
        print(f"{i}. {gene['gene']}: {gene['total_score']:.2f}")
        print(f"   Sequences: {gene['sequence_count']}")
        print(f"   Avg Length: {gene['avg_length']:.0f} bp")

    return ranked

# Usage
ranked_genes = compare_all_genes(taxid="562", email="your.email@example.com")
```

### Example 19: Custom Exclusion Strategy in Python

```python
from assay_design.data_retrieval import get_related_taxa, get_taxon_info

def smart_exclusion_selection(taxid, email):
    """Custom exclusion strategy combining multiple approaches."""

    # Get taxon info
    taxon_info = get_taxon_info(taxid, email)

    # Strategy 1: Get phylogenetically close siblings
    close_siblings = get_related_taxa(
        taxid, email,
        relationship="sibling",
        max_results=3,
        phylo_distance_weight=0.9  # Very close relatives
    )

    # Strategy 2: Get some distant relatives
    distant_siblings = get_related_taxa(
        taxid, email,
        relationship="sibling",
        max_results=2,
        phylo_distance_weight=0.3  # More distant relatives
    )

    # Combine results
    all_exclusion = close_siblings + distant_siblings

    print(f"Close exclusion taxa ({len(close_siblings)}):")
    for taxon in close_siblings:
        print(f"  - {taxon['scientific_name']} (taxid: {taxon['taxid']})")

    print(f"\nDistant exclusion taxa ({len(distant_siblings)}):")
    for taxon in distant_siblings:
        print(f"  - {taxon['scientific_name']} (taxid: {taxon['taxid']})")

    return all_exclusion

# Usage
exclusion_taxa = smart_exclusion_selection(taxid="689", email="your.email@example.com")
```

---

## Troubleshooting Common Issues

### Issue 1: No Sequences Found for Gene

**Problem:**
```
Error: No sequences found for gene 'gyrB' in taxid 12345
```

**Solutions:**
1. Try auto-gene selection instead of manual:
   ```bash
   assay-design --inclusion 12345 --email your.email@example.com --auto-gene --auto-exclusion
   ```

2. Use --suggest-genes to see what's available:
   ```bash
   assay-design --inclusion 12345 --email your.email@example.com --suggest-genes
   ```

3. Relax gene criteria:
   ```bash
   assay-design --inclusion 12345 --email your.email@example.com --auto-gene \
     --gene-criteria '{"min_sequence_count": 3}'
   ```

### Issue 2: Too Few Exclusion Taxa

**Problem:**
```
Warning: Only 1 exclusion taxon found
```

**Solutions:**
1. Use a broader exclusion strategy:
   ```bash
   assay-design --inclusion 689 --email your.email@example.com \
     --auto-gene --exclusion-strategy family
   ```

2. Reduce phylo-distance weighting to include more distant taxa:
   ```bash
   assay-design --inclusion 689 --email your.email@example.com \
     --auto-gene --auto-exclusion --exclusion-phylo-distance 0.3
   ```

### Issue 3: Gene Selection Scores Are Low

**Problem:**
```
Warning: Best gene score is only 5.2/10
```

**Meaning:** Your target organism has limited sequence data available

**Solutions:**
1. Proceed with the best available gene:
   ```bash
   assay-design --inclusion 12345 --email your.email@example.com --auto-gene --auto-exclusion
   ```

2. Use 16S rRNA as a fallback (usually has more sequences):
   ```bash
   assay-design --inclusion 12345 --gene "16S ribosomal RNA" --email your.email@example.com
   ```

### Issue 4: Primers Not Specific Enough

**Problem:** Designed primers cross-react with exclusion taxa

**Solutions:**
1. Add more exclusion taxa:
   ```bash
   assay-design --inclusion 689 --email your.email@example.com \
     --auto-gene --exclusion-strategy family  # Broader exclusion
   ```

2. Try a different gene:
   ```bash
   assay-design --inclusion 689 --email your.email@example.com --suggest-genes
   # Then select a different top-ranked gene
   ```

3. Manually add problematic taxa to exclusion:
   ```bash
   assay-design --inclusion 689 --exclusion 717,670,672,999 --gene "rpoB" \
     --email your.email@example.com
   ```

---

## Best Practices

### 1. Start with Full Automation
Begin with the fully automated workflow:
```bash
assay-design --inclusion <TAXID> --email <EMAIL> --auto-gene --auto-exclusion
```

### 2. Check Gene Suggestions
If results aren't satisfactory, review available genes:
```bash
assay-design --inclusion <TAXID> --email <EMAIL> --suggest-genes
```

### 3. Adjust Exclusion Strategy
Fine-tune specificity by adjusting the exclusion strategy:
- Species-specific → `--exclusion-strategy genus`
- Genus-specific → `--exclusion-strategy family`
- Family-specific → `--exclusion-strategy order`

### 4. Customize for Your Organism
- Well-characterized organisms: Increase `min_sequence_count` (30-50)
- Rare organisms: Decrease `min_sequence_count` (5-10)
- Phylogenetic studies: Use `single_copy_preferred=true`
- Multi-copy acceptable: Allow rRNA genes for better sequence availability

### 5. Validate Results
Always validate your primers using the built-in specificity validation or external tools like Primer-BLAST.

---

## Additional Resources

- **NCBI Taxonomy Browser**: https://www.ncbi.nlm.nih.gov/Taxonomy/
- **Primer3 Documentation**: https://primer3.org/
- **Biopython Documentation**: https://biopython.org/

For more help, see the main [README.md](README.md) or open an issue on GitHub.
