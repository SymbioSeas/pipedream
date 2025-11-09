# assay_design/gene_selection.py

"""
Gene selection module for taxa-specific dPCR assay design.

This module provides curated databases of gene suitability for dPCR assays,
including single-copy genes, HGT-resistant genes, and phylogenetic markers
across bacterial, archaeal, and fungal domains.

Key criteria for dPCR target genes:
1. Single or low copy number (avoid multi-copy bias)
2. Low horizontal gene transfer (HGT) propensity
3. Moderate evolutionary rate (conserved but informative)
4. Appropriate length (500-2000 bp for primer design flexibility)
5. Housekeeping or marker gene (universally present)
"""

import logging
from typing import Dict, List, Optional, Tuple
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass
class GeneMetadata:
    """
    Metadata for a candidate gene for dPCR assay design.

    Attributes:
        gene_name: Official gene name
        copy_number: 'single', 'low' (2-3 copies), or 'multi' (>3 copies)
        hgt_propensity: 'low', 'medium', or 'high' likelihood of horizontal transfer
        evolutionary_rate: 'slow', 'moderate', or 'fast' rate of sequence evolution
        typical_length_bp: (min, max) typical gene length range in base pairs
        functional_category: 'housekeeping', 'marker', or 'metabolic'
        description: Brief description of gene function and suitability
    """
    gene_name: str
    copy_number: str  # 'single', 'low', 'multi'
    hgt_propensity: str  # 'low', 'medium', 'high'
    evolutionary_rate: str  # 'slow', 'moderate', 'fast'
    typical_length_bp: Tuple[int, int]  # (min, max)
    functional_category: str  # 'housekeeping', 'marker', 'metabolic'
    description: str


class GeneSelectionCriteria:
    """
    Curated database of gene suitability for dPCR assay design.

    This class provides pre-vetted gene candidates organized by domain
    (bacteria, archaea, fungi) and suitability criteria.

    Gene databases include:
    - Single-copy housekeeping genes (ideal for quantification)
    - HGT-resistant genes (phylogenetically informative)
    - Marker genes (widely used, may be multi-copy)
    """

    # ==================== BACTERIA ====================

    BACTERIA_SINGLE_COPY = {
        'rpoB': GeneMetadata(
            gene_name='rpoB',
            copy_number='single',
            hgt_propensity='low',
            evolutionary_rate='moderate',
            typical_length_bp=(3500, 4200),
            functional_category='housekeeping',
            description='RNA polymerase beta subunit - highly conserved, single copy, '
                       'universally present in bacteria. Excellent for species ID and quantification.'
        ),
        'gyrB': GeneMetadata(
            gene_name='gyrB',
            copy_number='single',
            hgt_propensity='low',
            evolutionary_rate='moderate',
            typical_length_bp=(1800, 2400),
            functional_category='housekeeping',
            description='DNA gyrase subunit B - single copy, HGT-resistant, '
                       'good phylogenetic resolution at species/genus level.'
        ),
        'recA': GeneMetadata(
            gene_name='recA',
            copy_number='single',
            hgt_propensity='low',
            evolutionary_rate='moderate',
            typical_length_bp=(1000, 1200),
            functional_category='housekeeping',
            description='Recombinase A - universal bacterial gene, single copy, '
                       'widely used for phylogenetic studies.'
        ),
        'dnaK': GeneMetadata(
            gene_name='dnaK',
            copy_number='single',
            hgt_propensity='low',
            evolutionary_rate='moderate',
            typical_length_bp=(1800, 2000),
            functional_category='housekeeping',
            description='Heat shock protein 70 (HSP70) - universally conserved chaperone, '
                       'single copy in most bacteria.'
        ),
        'groEL': GeneMetadata(
            gene_name='groEL',
            copy_number='single',
            hgt_propensity='low',
            evolutionary_rate='moderate',
            typical_length_bp=(1600, 1650),
            functional_category='housekeeping',
            description='60 kDa chaperonin (cpn60, hsp60) - single copy essential chaperone, '
                       'highly conserved across bacteria.'
        ),
        'rpoD': GeneMetadata(
            gene_name='rpoD',
            copy_number='single',
            hgt_propensity='low',
            evolutionary_rate='moderate',
            typical_length_bp=(1900, 2100),
            functional_category='housekeeping',
            description='RNA polymerase sigma factor 70 - primary sigma factor, '
                       'single copy, essential gene.'
        ),
        'gyrA': GeneMetadata(
            gene_name='gyrA',
            copy_number='single',
            hgt_propensity='low',
            evolutionary_rate='moderate',
            typical_length_bp=(2400, 2700),
            functional_category='housekeeping',
            description='DNA gyrase subunit A - single copy, pairs with gyrB, '
                       'good for species discrimination.'
        ),
        'fusA': GeneMetadata(
            gene_name='fusA',
            copy_number='single',
            hgt_propensity='low',
            evolutionary_rate='moderate',
            typical_length_bp=(1900, 2100),
            functional_category='housekeeping',
            description='Elongation factor G (EF-G) - single copy translation factor, '
                       'highly conserved.'
        ),
        'atpD': GeneMetadata(
            gene_name='atpD',
            copy_number='single',
            hgt_propensity='low',
            evolutionary_rate='moderate',
            typical_length_bp=(1400, 1600),
            functional_category='housekeeping',
            description='ATP synthase F1 beta subunit - single copy essential metabolic gene, '
                       'good phylogenetic marker.'
        ),
    }

    BACTERIA_MARKER = {
        '16S ribosomal RNA': GeneMetadata(
            gene_name='16S ribosomal RNA',
            copy_number='multi',  # Typically 1-15 copies depending on species
            hgt_propensity='low',
            evolutionary_rate='slow',
            typical_length_bp=(1400, 1600),
            functional_category='marker',
            description='Small subunit ribosomal RNA - universal phylogenetic marker, '
                       'most widely sequenced gene. Multi-copy may cause quantification bias.'
        ),
        '23S ribosomal RNA': GeneMetadata(
            gene_name='23S ribosomal RNA',
            copy_number='multi',
            hgt_propensity='low',
            evolutionary_rate='slow',
            typical_length_bp=(2800, 3000),
            functional_category='marker',
            description='Large subunit ribosomal RNA - alternative to 16S, more conserved, '
                       'useful for higher-level phylogeny.'
        ),
    }

    # ==================== ARCHAEA ====================

    ARCHAEA_SINGLE_COPY = {
        'rps3': GeneMetadata(
            gene_name='rps3',
            copy_number='single',
            hgt_propensity='low',
            evolutionary_rate='slow',
            typical_length_bp=(600, 700),
            functional_category='housekeeping',
            description='Ribosomal protein S3 - conserved archaeal marker, '
                       'single copy, good for phylogenetic analysis.'
        ),
        'rpl2': GeneMetadata(
            gene_name='rpl2',
            copy_number='single',
            hgt_propensity='low',
            evolutionary_rate='slow',
            typical_length_bp=(700, 900),
            functional_category='housekeeping',
            description='Ribosomal protein L2 - conserved archaeal marker, '
                       'single copy, commonly used in archaeal phylogenomics.'
        ),
        'ef2': GeneMetadata(
            gene_name='EF-2',
            copy_number='single',
            hgt_propensity='low',
            evolutionary_rate='moderate',
            typical_length_bp=(2200, 2500),
            functional_category='housekeeping',
            description='Elongation factor 2 - single copy translation factor, '
                       'good resolution for archaeal species.'
        ),
    }

    ARCHAEA_MARKER = {
        '16S ribosomal RNA': GeneMetadata(
            gene_name='16S ribosomal RNA',
            copy_number='multi',
            hgt_propensity='low',
            evolutionary_rate='slow',
            typical_length_bp=(1400, 1500),
            functional_category='marker',
            description='Small subunit ribosomal RNA - widely used archaeal marker, '
                       'similar to bacterial 16S.'
        ),
    }

    # ==================== FUNGI/EUKARYOTES ====================

    FUNGI_SINGLE_COPY = {
        'RPB1': GeneMetadata(
            gene_name='RPB1',
            copy_number='single',
            hgt_propensity='low',
            evolutionary_rate='moderate',
            typical_length_bp=(3500, 4500),
            functional_category='housekeeping',
            description='RNA polymerase II largest subunit - single copy essential gene, '
                       'excellent for fungal phylogenetics.'
        ),
        'RPB2': GeneMetadata(
            gene_name='RPB2',
            copy_number='single',
            hgt_propensity='low',
            evolutionary_rate='moderate',
            typical_length_bp=(3000, 4000),
            functional_category='housekeeping',
            description='RNA polymerase II second largest subunit - single copy, '
                       'widely used fungal phylogenetic marker.'
        ),
        'TEF1': GeneMetadata(
            gene_name='TEF1',
            copy_number='single',
            hgt_propensity='low',
            evolutionary_rate='moderate',
            typical_length_bp=(1200, 1500),
            functional_category='housekeeping',
            description='Translation elongation factor 1-alpha - single copy, '
                       'good species-level resolution in fungi.'
        ),
        'TUB2': GeneMetadata(
            gene_name='beta-tubulin',
            copy_number='single',
            hgt_propensity='low',
            evolutionary_rate='moderate',
            typical_length_bp=(1200, 1400),
            functional_category='housekeeping',
            description='Beta-tubulin - single copy cytoskeletal protein, '
                       'commonly used for fungal species identification.'
        ),
        'ACT': GeneMetadata(
            gene_name='actin',
            copy_number='low',  # Usually 1-2 copies
            hgt_propensity='low',
            evolutionary_rate='slow',
            typical_length_bp=(1100, 1300),
            functional_category='housekeeping',
            description='Actin - low copy cytoskeletal protein, '
                       'conserved housekeeping gene.'
        ),
    }

    FUNGI_MARKER = {
        'ITS': GeneMetadata(
            gene_name='ITS',
            copy_number='multi',
            hgt_propensity='low',
            evolutionary_rate='fast',
            typical_length_bp=(500, 700),
            functional_category='marker',
            description='Internal transcribed spacer (ITS1-5.8S-ITS2) - standard fungal DNA barcode, '
                       'high variability for species discrimination. Multi-copy.'
        ),
        '18S ribosomal RNA': GeneMetadata(
            gene_name='18S ribosomal RNA',
            copy_number='multi',
            hgt_propensity='low',
            evolutionary_rate='slow',
            typical_length_bp=(1700, 1900),
            functional_category='marker',
            description='Small subunit ribosomal RNA - eukaryotic universal marker, '
                       'conserved, useful for higher-level phylogeny.'
        ),
        '28S ribosomal RNA': GeneMetadata(
            gene_name='28S ribosomal RNA',
            copy_number='multi',
            hgt_propensity='low',
            evolutionary_rate='slow',
            typical_length_bp=(3300, 3600),
            functional_category='marker',
            description='Large subunit ribosomal RNA - alternative eukaryotic marker, '
                       'contains variable D1/D2 regions.'
        ),
    }

    # ==================== CLASS METHODS ====================

    @classmethod
    def get_genes_for_domain(
        cls,
        domain: str,
        prefer_single_copy: bool = True
    ) -> Dict[str, GeneMetadata]:
        """
        Get recommended genes for a taxonomic domain.

        Args:
            domain: Taxonomic domain - 'bacteria', 'archaea', 'fungi', or 'eukaryote'
            prefer_single_copy: If True, only return single-copy genes;
                              if False, include marker genes (may be multi-copy)

        Returns:
            Dict mapping gene_name -> GeneMetadata

        Examples:
            >>> genes = GeneSelectionCriteria.get_genes_for_domain('bacteria')
            >>> print(len(genes))  # Single-copy genes only
            9
            >>> all_genes = GeneSelectionCriteria.get_genes_for_domain('bacteria', False)
            >>> print(len(all_genes))  # Includes 16S, 23S
            11
        """
        domain = domain.lower()

        if domain == 'bacteria':
            if prefer_single_copy:
                return cls.BACTERIA_SINGLE_COPY.copy()
            else:
                return {**cls.BACTERIA_SINGLE_COPY, **cls.BACTERIA_MARKER}

        elif domain == 'archaea':
            if prefer_single_copy:
                return cls.ARCHAEA_SINGLE_COPY.copy()
            else:
                return {**cls.ARCHAEA_SINGLE_COPY, **cls.ARCHAEA_MARKER}

        elif domain in ['fungi', 'eukaryote']:
            if prefer_single_copy:
                return cls.FUNGI_SINGLE_COPY.copy()
            else:
                return {**cls.FUNGI_SINGLE_COPY, **cls.FUNGI_MARKER}

        else:
            logger.warning(f"Unknown domain '{domain}', defaulting to bacteria")
            return cls.BACTERIA_SINGLE_COPY.copy()

    @classmethod
    def rank_genes_by_criteria(
        cls,
        genes: Dict[str, GeneMetadata],
        criteria_weights: Optional[Dict[str, float]] = None
    ) -> List[Tuple[str, float]]:
        """
        Rank genes by suitability for dPCR assays.

        Scoring criteria:
        1. Single copy status (40% weight) - single > low > multi
        2. HGT resistance (30% weight) - low > medium > high
        3. Evolutionary rate (20% weight) - moderate > slow > fast
        4. Length appropriateness (10% weight) - ideal ~1500 bp

        Args:
            genes: Dict of gene_name -> GeneMetadata
            criteria_weights: Optional custom weights for each criterion
                Default: {'single_copy': 0.4, 'hgt_resistance': 0.3,
                         'evolutionary_rate': 0.2, 'length': 0.1}

        Returns:
            List of (gene_name, score) tuples, sorted by score (descending)
            Score range: 0.0 to 1.0

        Examples:
            >>> genes = GeneSelectionCriteria.get_genes_for_domain('bacteria')
            >>> ranked = GeneSelectionCriteria.rank_genes_by_criteria(genes)
            >>> print(ranked[0])  # Top gene
            ('rpoB', 0.95)
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

            # 1. Copy number score (higher is better)
            copy_scores = {'single': 1.0, 'low': 0.7, 'multi': 0.3}
            copy_score = copy_scores.get(metadata.copy_number, 0.5)
            score += copy_score * criteria_weights['single_copy']

            # 2. HGT resistance score (low HGT is better)
            hgt_scores = {'low': 1.0, 'medium': 0.5, 'high': 0.1}
            hgt_score = hgt_scores.get(metadata.hgt_propensity, 0.5)
            score += hgt_score * criteria_weights['hgt_resistance']

            # 3. Evolutionary rate score (moderate is best for primer design)
            evo_scores = {'slow': 0.8, 'moderate': 1.0, 'fast': 0.5}
            evo_score = evo_scores.get(metadata.evolutionary_rate, 0.5)
            score += evo_score * criteria_weights['evolutionary_rate']

            # 4. Length score (prefer genes around 1500 bp for primer flexibility)
            min_len, max_len = metadata.typical_length_bp
            avg_length = (min_len + max_len) / 2
            ideal_length = 1500

            # Calculate distance from ideal
            length_diff = abs(avg_length - ideal_length)
            # Normalize to 0-1 (perfect score at ideal length, 0 at >3000 bp difference)
            length_score = max(1.0 - (length_diff / 3000), 0.0)
            score += length_score * criteria_weights['length']

            scored_genes.append((gene_name, round(score, 3)))

        # Sort by score (descending)
        return sorted(scored_genes, key=lambda x: x[1], reverse=True)

    @classmethod
    def get_gene_metadata(cls, gene_name: str, domain: Optional[str] = None) -> Optional[GeneMetadata]:
        """
        Retrieve metadata for a specific gene.

        Args:
            gene_name: Name of the gene
            domain: Optional domain to search in (bacteria, archaea, fungi)
                   If None, searches all domains

        Returns:
            GeneMetadata if found, None otherwise

        Examples:
            >>> metadata = GeneSelectionCriteria.get_gene_metadata('rpoB')
            >>> print(metadata.copy_number)
            single
            >>> print(metadata.description)
            RNA polymerase beta subunit - highly conserved...
        """
        # Search in specified domain
        if domain:
            genes = cls.get_genes_for_domain(domain, prefer_single_copy=False)
            return genes.get(gene_name)

        # Search all domains
        all_databases = [
            cls.BACTERIA_SINGLE_COPY,
            cls.BACTERIA_MARKER,
            cls.ARCHAEA_SINGLE_COPY,
            cls.ARCHAEA_MARKER,
            cls.FUNGI_SINGLE_COPY,
            cls.FUNGI_MARKER,
        ]

        for db in all_databases:
            if gene_name in db:
                return db[gene_name]

        return None

    @classmethod
    def list_all_genes(cls, domain: Optional[str] = None) -> List[str]:
        """
        List all available genes, optionally filtered by domain.

        Args:
            domain: Optional domain filter (bacteria, archaea, fungi)

        Returns:
            List of gene names

        Examples:
            >>> all_genes = GeneSelectionCriteria.list_all_genes()
            >>> bacterial_genes = GeneSelectionCriteria.list_all_genes('bacteria')
        """
        if domain:
            genes = cls.get_genes_for_domain(domain, prefer_single_copy=False)
            return sorted(genes.keys())

        # Return all genes from all domains
        all_genes = set()

        for db in [cls.BACTERIA_SINGLE_COPY, cls.BACTERIA_MARKER,
                   cls.ARCHAEA_SINGLE_COPY, cls.ARCHAEA_MARKER,
                   cls.FUNGI_SINGLE_COPY, cls.FUNGI_MARKER]:
            all_genes.update(db.keys())

        return sorted(all_genes)

    @classmethod
    def get_recommended_genes(
        cls,
        domain: str,
        use_case: str = 'quantification',
        max_genes: int = 5
    ) -> List[Tuple[str, GeneMetadata, float]]:
        """
        Get top recommended genes for a specific use case.

        Use cases:
        - 'quantification': Prioritize single-copy genes (for absolute quantification)
        - 'phylogeny': Balance between conservation and variation
        - 'detection': Include marker genes (higher sensitivity)

        Args:
            domain: Taxonomic domain
            use_case: 'quantification', 'phylogeny', or 'detection'
            max_genes: Maximum number of recommendations

        Returns:
            List of (gene_name, metadata, score) tuples

        Examples:
            >>> recs = GeneSelectionCriteria.get_recommended_genes('bacteria', 'quantification')
            >>> print(f"Top gene: {recs[0][0]} (score: {recs[0][2]})")
        """
        # Adjust weights based on use case
        if use_case == 'quantification':
            weights = {
                'single_copy': 0.6,  # Very important
                'hgt_resistance': 0.2,
                'evolutionary_rate': 0.1,
                'length': 0.1
            }
            prefer_single_copy = True

        elif use_case == 'phylogeny':
            weights = {
                'single_copy': 0.3,
                'hgt_resistance': 0.4,  # Very important
                'evolutionary_rate': 0.2,
                'length': 0.1
            }
            prefer_single_copy = True

        elif use_case == 'detection':
            weights = {
                'single_copy': 0.2,  # Less important
                'hgt_resistance': 0.2,
                'evolutionary_rate': 0.3,
                'length': 0.3
            }
            prefer_single_copy = False  # Include marker genes

        else:
            logger.warning(f"Unknown use case '{use_case}', using default weights")
            weights = None
            prefer_single_copy = True

        # Get genes for domain
        genes = cls.get_genes_for_domain(domain, prefer_single_copy=prefer_single_copy)

        # Rank genes
        ranked = cls.rank_genes_by_criteria(genes, criteria_weights=weights)

        # Return top genes with full metadata
        results = []
        for gene_name, score in ranked[:max_genes]:
            metadata = genes[gene_name]
            results.append((gene_name, metadata, score))

        return results
