"""
Criterion 1: Pathway Mechanism Alignment  [0,1]

Sub-metric A: Pathway Overlap (Jaccard)
  - Disease pathway set from Reactome
  - Drug target pathway set from Reactome (using expanded STRING neighborhood)
  - Score = Jaccard similarity

Sub-metric B: Perturbation Signal Similarity
  - Enrichr GO-BP enrichment comparison (disease gene set vs drug target gene set)
  - Uses expanded target set to capture downstream effectors

Composite: 0.5 * A + 0.5 * B

Target expansion: When a drug has < 5 direct MOA targets, expands to include
1-hop STRING neighbors (confidence ≥ 700) to capture downstream effectors
(e.g. MTOR, RPTOR for FKBP1A-binding mTOR inhibitors). This expansion is used
only for pathway scoring; network proximity uses direct MOA targets only.
"""

import logging
from src.loaders.reactome_loader import get_combined_pathways_for_genes
from src.loaders.enrichr_loader import enrich_gene_list, get_pathway_score_vector, cosine_similarity_pathway_vectors
from src.loaders.opentargets_loader import get_drug_target_genes as ot_get_drug_targets
from src.loaders.chembl_loader import get_targets_for_drug, get_gene_symbols_for_target
from src.loaders.string_loader import get_interaction_partners

logger = logging.getLogger(__name__)

# Genes to exclude from STRING expansion — not meaningful for pathway scoring
_EXPANSION_BLOCKLIST = {
    "RYR1", "RYR2", "RYR3",  # Ryanodine receptors — calcium, not signaling
    "TRDN", "CASQ1", "CASQ2",  # SR proteins
}


def score_pathway_alignment(
    drug_name: str,
    drug_chembl_id: str,
    disease_module_genes: list[str],
    disease_pathway_vector: dict[str, float] | None = None,
) -> dict:
    """
    Compute Pathway Mechanism Alignment score for a drug vs. a disease module.

    Parameters:
        drug_name: human-readable drug name
        drug_chembl_id: ChEMBL compound ID
        disease_module_genes: gene symbols in the disease module
        disease_pathway_vector: precomputed pathway score vector for disease (or None to compute)

    Returns dict with:
        score: float [0,1]
        jaccard: float  — Sub-metric A
        signal_similarity: float  — Sub-metric B
        drug_target_genes: list[str]  — direct MOA targets (for network proximity)
        drug_target_genes_expanded: list[str]  — MOA + STRING neighbors (for pathway scoring)
        shared_pathways: list[str]  — top shared pathways
    """
    # Resolve direct MOA targets
    direct_targets = _get_direct_targets(drug_chembl_id)

    if not direct_targets:
        return _zero_result("No high-confidence targets found for drug")

    # Expand to STRING neighbors for richer pathway signature
    expanded_targets = _expand_with_string_neighbors(direct_targets)

    # Sub-metric A: Pathway Overlap (Jaccard via Reactome — use expanded set)
    jaccard, shared_pathways = _compute_jaccard(expanded_targets, disease_module_genes)

    # Sub-metric B: Perturbation Signal Similarity (Enrichr — use expanded set)
    signal_sim = _compute_signal_similarity(expanded_targets, disease_module_genes, disease_pathway_vector)

    # Composite
    score = round(0.5 * jaccard + 0.5 * signal_sim, 3)

    return {
        "score": score,
        "jaccard": round(jaccard, 3),
        "signal_similarity": round(signal_sim, 3),
        "drug_target_genes": direct_targets,          # Used by network proximity
        "drug_target_genes_expanded": expanded_targets,
        "shared_pathways": shared_pathways[:5],
        "n_drug_targets": len(direct_targets),
    }


def _get_direct_targets(chembl_id: str) -> list[str]:
    """
    Resolve ChEMBL ID to direct drug target gene symbols.
    Primary: OpenTargets drug MOA.
    Fallback: ChEMBL activity data (pChEMBL >= 7).
    """
    genes = ot_get_drug_targets(chembl_id)
    if genes:
        return genes

    try:
        targets = get_targets_for_drug(chembl_id, min_pchembl=7.0)
        chembl_genes = []
        for target in targets:
            syms = get_gene_symbols_for_target(target["target_chembl_id"])
            chembl_genes.extend(syms)
        return list(set(chembl_genes))
    except Exception:
        return []


def _expand_with_string_neighbors(direct_targets: list[str]) -> list[str]:
    """
    Expand drug target gene set with top STRING 1-hop neighbors (confidence ≥ 700).
    Only expands when fewer than 5 direct targets found.
    Excludes genes in blocklist (structural proteins, ion channels unrelated to signaling).
    """
    if len(direct_targets) >= 5:
        return direct_targets

    try:
        interactions = get_interaction_partners(direct_targets, confidence=700, limit=30)
        neighbors: dict[str, float] = {}
        for e in interactions:
            for gene in direct_targets:
                neighbor = None
                if e["gene_a"] == gene and e["gene_b"] not in direct_targets:
                    neighbor = e["gene_b"]
                elif e["gene_b"] == gene and e["gene_a"] not in direct_targets:
                    neighbor = e["gene_a"]
                if neighbor and neighbor not in _EXPANSION_BLOCKLIST:
                    if neighbor not in neighbors or e["score"] > neighbors[neighbor]:
                        neighbors[neighbor] = e["score"]

        # Take top 10 neighbors by STRING score
        top_neighbors = [g for g, _ in sorted(neighbors.items(), key=lambda x: x[1], reverse=True)[:10]]
        expanded = list(set(direct_targets) | set(top_neighbors))
        logger.debug("Expanded %s -> %d genes via STRING", direct_targets, len(expanded))
        return expanded
    except Exception as e:
        logger.debug("STRING expansion failed: %s", e)
        return direct_targets


def _compute_jaccard(drug_genes: list[str], disease_genes: list[str]) -> tuple[float, list[str]]:
    """Compute Jaccard similarity of Reactome pathway sets."""
    drug_pathways = get_combined_pathways_for_genes(drug_genes[:15])
    disease_pathways = get_combined_pathways_for_genes(disease_genes[:15])

    if not drug_pathways or not disease_pathways:
        return 0.0, []

    intersection = drug_pathways & disease_pathways
    union = drug_pathways | disease_pathways

    jaccard = len(intersection) / len(union) if union else 0.0
    shared = sorted(intersection)

    return jaccard, shared


def _compute_signal_similarity(
    drug_genes: list[str],
    disease_genes: list[str],
    disease_pathway_vector: dict[str, float] | None = None,
) -> float:
    """
    Compute cosine similarity of Enrichr GO-BP pathway score vectors.
    Uses the disease module gene set as the disease reference.
    """
    drug_enrichment = enrich_gene_list(drug_genes[:50])
    drug_vector = get_pathway_score_vector(drug_enrichment)

    if not drug_vector:
        return 0.0

    if disease_pathway_vector is None:
        disease_enrichment = enrich_gene_list(disease_genes[:50])
        disease_pathway_vector = get_pathway_score_vector(disease_enrichment)

    if not disease_pathway_vector:
        return 0.0

    return cosine_similarity_pathway_vectors(drug_vector, disease_pathway_vector)


def _zero_result(reason: str) -> dict:
    return {
        "score": 0.0,
        "jaccard": 0.0,
        "signal_similarity": 0.0,
        "drug_target_genes": [],
        "drug_target_genes_expanded": [],
        "shared_pathways": [],
        "n_drug_targets": 0,
        "note": reason,
    }
