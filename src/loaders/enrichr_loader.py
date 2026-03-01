"""
Enrichr API loader.
Performs gene set enrichment for pathway alignment (Criteria 1A and 1B).
No registration required. Returns ranked pathways with p-values in ~2s per query.
"""

import time
import logging
import requests
import numpy as np
from .cache import file_cache

logger = logging.getLogger(__name__)
BASE_URL = "https://maayanlab.cloud/Enrichr"
TIMEOUT = 30

# Libraries used for enrichment
LIBRARIES = ["KEGG_2021_Human", "Reactome_2022", "GO_Biological_Process_2023"]


@file_cache("enrichr")
def enrich_gene_list(gene_symbols: list[str], top_n: int = 50) -> dict[str, list[dict]]:
    """
    Submit a gene list to Enrichr and return enrichment results across all libraries.
    Returns: {library_name: [{term, pvalue, zscore, combined_score, genes}, ...]}
    """
    if not gene_symbols:
        return {lib: [] for lib in LIBRARIES}

    # Submit gene list
    user_list_id = _submit_gene_list(gene_symbols)
    if not user_list_id:
        return {lib: [] for lib in LIBRARIES}

    results = {}
    for lib in LIBRARIES:
        enrichment = _get_enrichment(user_list_id, lib, top_n)
        results[lib] = enrichment
        time.sleep(0.2)

    return results


def _submit_gene_list(gene_symbols: list[str]) -> int | None:
    try:
        gene_str = "\n".join(gene_symbols)
        resp = requests.post(
            f"{BASE_URL}/addList",
            files={"list": (None, gene_str), "description": (None, "hackrare_query")},
            timeout=TIMEOUT,
        )
        resp.raise_for_status()
        data = resp.json()
        return data.get("userListId")
    except Exception as e:
        logger.error("Enrichr submit error: %s", e)
        return None


def _get_enrichment(user_list_id: int, library: str, top_n: int = 50) -> list[dict]:
    try:
        resp = requests.get(
            f"{BASE_URL}/enrich",
            params={"userListId": user_list_id, "backgroundType": library},
            timeout=TIMEOUT,
        )
        resp.raise_for_status()
        data = resp.json()
        raw = data.get(library, [])
        results = []
        for item in raw[:top_n]:
            # Enrichr format: [rank, term, pvalue, zscore, combined_score, genes, adj_pvalue, ...]
            if len(item) >= 6:
                results.append({
                    "rank": item[0],
                    "term": item[1],
                    "pvalue": item[2],
                    "zscore": item[3],
                    "combined_score": item[4],
                    "genes": item[5] if isinstance(item[5], list) else [],
                    "adj_pvalue": item[6] if len(item) > 6 else item[2],
                })
        return results
    except Exception as e:
        logger.error("Enrichr enrich error for %s: %s", library, e)
        return []


def get_pathway_score_vector(enrichment_results: dict[str, list[dict]], top_n: int = 20) -> dict[str, float]:
    """
    Convert enrichment results into a pathway score vector.
    Score = -log10(pvalue) for each pathway term (across all libraries).
    Used for cosine similarity comparison between two gene sets.
    Returns: {pathway_term: score}
    """
    import math
    scores = {}
    for lib, terms in enrichment_results.items():
        for t in terms[:top_n]:
            pval = t.get("pvalue", 1.0)
            score = -math.log10(max(pval, 1e-10))
            term = t["term"]
            # Keep highest score if term appears in multiple libraries
            if term not in scores or score > scores[term]:
                scores[term] = score
    return scores


def cosine_similarity_pathway_vectors(vec_a: dict[str, float], vec_b: dict[str, float]) -> float:
    """
    Compute cosine similarity between two pathway score vectors.
    Shared terms form the common dimensions.
    """
    if not vec_a or not vec_b:
        return 0.0

    all_terms = list(set(vec_a.keys()) | set(vec_b.keys()))
    if not all_terms:
        return 0.0

    a = np.array([vec_a.get(t, 0.0) for t in all_terms])
    b = np.array([vec_b.get(t, 0.0) for t in all_terms])

    norm_a = np.linalg.norm(a)
    norm_b = np.linalg.norm(b)
    if norm_a == 0 or norm_b == 0:
        return 0.0

    return float(np.dot(a, b) / (norm_a * norm_b))
