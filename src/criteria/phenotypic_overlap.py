"""
Criterion 3: Phenotypic Overlap  [0,1]

Pipeline:
  Step 1 — HPO Semantic Similarity (Resnik BMA)     → hpo_similarity
  Step 2 — HPO → Gene → Cosine similarity           → gene_similarity
  Step 3 — Gene → Pathway Enrichment → Cosine       → pathway_similarity
  Step 4 — Composite: 0.4*hpo + 0.3*gene + 0.3*pathway

Sources: HPO files (hp.obo, phenotype_to_genes.txt), Enrichr API
"""

import logging
import numpy as np
from src.loaders.hpo_loader import HPOSimilarity, get_disease_hpo_terms
from src.loaders.enrichr_loader import enrich_gene_list, get_pathway_score_vector, cosine_similarity_pathway_vectors

logger = logging.getLogger(__name__)

# Singleton HPO similarity instance (loaded once per pipeline run)
_hpo_sim: HPOSimilarity | None = None


def get_hpo_similarity() -> HPOSimilarity:
    global _hpo_sim
    if _hpo_sim is None:
        _hpo_sim = HPOSimilarity()
    return _hpo_sim


def score_phenotypic_overlap(
    disease_a_omim_ids: list[str],
    disease_b_omim_ids: list[str],
    disease_a_genes: list[str],
    disease_b_genes: list[str],
    disease_a_pathway_vector: dict[str, float] | None = None,
) -> dict:
    """
    Compute phenotypic overlap between two diseases.

    Can be used for:
    - Disease A = root disease (e.g. TSC), Disease B = drug target disease
    - Disease A = root disease, Disease B = drug's mechanism disease

    Returns dict with:
        score: float [0,1]
        hpo_similarity: float
        gene_similarity: float
        pathway_similarity: float
        confidence: str — High/Med/Low
        shared_hpo_terms: list[dict]
        shared_genes: list[str]
        shared_pathways: list[str]
    """
    hpo_sim = get_hpo_similarity()

    # Get HPO terms for both diseases
    terms_a = _get_hpo_terms(disease_a_omim_ids, disease_a_genes, hpo_sim)
    terms_b = _get_hpo_terms(disease_b_omim_ids, disease_b_genes, hpo_sim)

    # Step 1: HPO semantic similarity
    if terms_a and terms_b:
        hpo_score = hpo_sim.bma_similarity(terms_a, terms_b)
        shared_hpo = hpo_sim.top_shared_hpo_terms(terms_a, terms_b, top_k=5)
    else:
        hpo_score = 0.0
        shared_hpo = []

    # Step 2: Gene similarity (cosine of IC-weighted gene vectors)
    gene_score, shared_genes = _compute_gene_similarity(terms_a, terms_b, hpo_sim)

    # Step 3: Pathway similarity (Enrichr enrichment cosine)
    pathway_score, shared_pathways = _compute_pathway_similarity(
        disease_a_genes, disease_b_genes, disease_a_pathway_vector
    )

    # Composite
    score = round(0.4 * hpo_score + 0.3 * gene_score + 0.3 * pathway_score, 3)
    confidence = _confidence_rating(len(terms_a), len(terms_b), len(shared_genes))

    return {
        "score": score,
        "hpo_similarity": round(hpo_score, 3),
        "gene_similarity": round(gene_score, 3),
        "pathway_similarity": round(pathway_score, 3),
        "confidence": confidence,
        "shared_hpo_terms": shared_hpo,
        "shared_genes": shared_genes[:10],
        "shared_pathways": shared_pathways[:5],
        "n_hpo_terms_a": len(terms_a),
        "n_hpo_terms_b": len(terms_b),
    }


def _get_hpo_terms(omim_ids: list[str], genes: list[str], hpo_sim: HPOSimilarity) -> list[str]:
    """
    Get HPO terms from OMIM disease annotations.
    Fallback: when no OMIM IDs, look up HPO terms annotated to the given genes
    via the reverse of the phenotype_to_genes mapping.
    """
    terms: set[str] = set()
    for omim_id in (omim_ids or []):
        terms.update(get_disease_hpo_terms(omim_id))
    # Gene-based fallback: find HPO terms that mention any of the given genes
    if not terms and genes:
        gene_set = set(genes)
        for term_id, term_genes in hpo_sim.gene_map.items():
            if gene_set & set(term_genes):
                terms.add(term_id)
    return list(terms)


def _compute_gene_similarity(
    terms_a: list[str],
    terms_b: list[str],
    hpo_sim: HPOSimilarity,
) -> tuple[float, list[str]]:
    """Compute cosine similarity of IC-weighted gene vectors from HPO term sets."""
    if not terms_a or not terms_b:
        return 0.0, []

    genes_a = hpo_sim.get_genes_for_terms(terms_a)
    genes_b = hpo_sim.get_genes_for_terms(terms_b)

    if not genes_a or not genes_b:
        return 0.0, []

    all_genes = list(set(genes_a.keys()) | set(genes_b.keys()))
    vec_a = np.array([genes_a.get(g, 0.0) for g in all_genes])
    vec_b = np.array([genes_b.get(g, 0.0) for g in all_genes])

    norm_a = np.linalg.norm(vec_a)
    norm_b = np.linalg.norm(vec_b)
    if norm_a == 0 or norm_b == 0:
        return 0.0, []

    sim = float(np.dot(vec_a, vec_b) / (norm_a * norm_b))

    # Find shared genes (present in both with non-zero score)
    shared = sorted(
        [g for g in all_genes if genes_a.get(g, 0) > 0 and genes_b.get(g, 0) > 0],
        key=lambda g: genes_a.get(g, 0) + genes_b.get(g, 0),
        reverse=True,
    )

    return round(sim, 3), shared[:10]


def _compute_pathway_similarity(
    genes_a: list[str],
    genes_b: list[str],
    disease_a_pathway_vector: dict[str, float] | None = None,
) -> tuple[float, list[str]]:
    """Compute cosine similarity of Enrichr pathway score vectors."""
    if not genes_a or not genes_b:
        return 0.0, []

    # Enrich disease A genes (or use precomputed)
    if disease_a_pathway_vector is None:
        enrichment_a = enrich_gene_list(genes_a[:50])
        vec_a = get_pathway_score_vector(enrichment_a)
    else:
        vec_a = disease_a_pathway_vector

    # Enrich disease B genes
    enrichment_b = enrich_gene_list(genes_b[:50])
    vec_b = get_pathway_score_vector(enrichment_b)

    if not vec_a or not vec_b:
        return 0.0, []

    sim = cosine_similarity_pathway_vectors(vec_a, vec_b)

    # Find shared top pathways
    shared = sorted(
        set(vec_a.keys()) & set(vec_b.keys()),
        key=lambda p: vec_a.get(p, 0) + vec_b.get(p, 0),
        reverse=True,
    )[:5]

    return round(sim, 3), shared


def _confidence_rating(n_terms_a: int, n_terms_b: int, n_shared_genes: int) -> str:
    if n_terms_a >= 10 and n_terms_b >= 10 and n_shared_genes >= 5:
        return "High"
    elif n_terms_a >= 5 or n_terms_b >= 5:
        return "Med"
    else:
        return "Low"
