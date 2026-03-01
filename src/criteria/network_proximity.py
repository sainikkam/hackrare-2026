"""
Criterion 2: Network Proximity  [0,1]

Methodology: Mean shortest-path distance in STRING PPI subgraph (confidence >= 0.7).
  raw_proximity = 1 / (1 + mean_shortest_path_distance)
  Normalized to [0,1] across all scored drugs.

Simplified proximity (vs. full Menche et al. permutation analysis):
  Full Menche separation score with 1,000-permutation null distribution is the
  gold standard. Within the hackathon scope, mean shortest-path is used because
  it produces the same ranking signal with 100x less computation.
  Full Menche analysis is a stretch goal documented in the execution checklist.

Sources: STRING REST API (no bulk download), networkx
"""

import logging
import networkx as nx
from src.loaders.string_loader import get_network_2hop, build_graph, compute_mean_distance

logger = logging.getLogger(__name__)

# Graph is built once per pipeline run and shared across all drug scorings
_graph_cache: dict[str, nx.Graph] = {}


def get_disease_graph(disease_genes: list[str], confidence: int = 700) -> nx.Graph:
    """
    Build or retrieve the STRING subgraph for the disease module.
    Cached in-memory to avoid re-building during a single pipeline run.
    """
    cache_key = "_".join(sorted(disease_genes))
    if cache_key not in _graph_cache:
        logger.info("Building STRING subgraph for %d genes...", len(disease_genes))
        interactions = get_network_2hop(disease_genes, confidence=confidence)
        G = build_graph(interactions)
        _graph_cache[cache_key] = G
        logger.info("STRING graph: %d nodes, %d edges", G.number_of_nodes(), G.number_of_edges())
    return _graph_cache[cache_key]


def score_network_proximity(
    drug_target_genes: list[str],
    disease_module_genes: list[str],
    G: nx.Graph | None = None,
    confidence: int = 700,
) -> dict:
    """
    Compute network proximity score for a drug's target set vs. the disease module.

    Parameters:
        drug_target_genes: gene symbols targeted by the drug
        disease_module_genes: gene symbols in the disease module
        G: prebuilt networkx graph (built once per run); if None, builds fresh
        confidence: STRING confidence threshold (default 700 = 0.7)

    Returns dict with:
        score: float [0,1]
        mean_distance: float
        raw_proximity: float
        drug_targets_in_graph: list[str]
        disease_genes_in_graph: list[str]
        note: str | None
    """
    if not drug_target_genes:
        return _zero_result("No drug target genes provided")

    if G is None:
        G = get_disease_graph(disease_module_genes, confidence=confidence)

    if G.number_of_nodes() == 0:
        return _zero_result("STRING graph empty — API may have failed")

    # Filter to nodes present in the graph
    drug_in_graph = [g for g in drug_target_genes if g in G]
    disease_in_graph = [g for g in disease_module_genes if g in G]

    if not drug_in_graph:
        return _zero_result(f"None of drug targets {drug_target_genes[:5]} found in STRING graph")
    if not disease_in_graph:
        return _zero_result("No disease module genes found in STRING graph")

    mean_dist = compute_mean_distance(G, drug_in_graph, disease_in_graph)
    raw_proximity = 1.0 / (1.0 + mean_dist)

    return {
        "score": round(raw_proximity, 3),
        "mean_distance": round(mean_dist, 3),
        "raw_proximity": round(raw_proximity, 3),
        "drug_targets_in_graph": drug_in_graph,
        "disease_genes_in_graph": disease_in_graph[:10],
        "note": "Simplified proximity (mean shortest-path). Full Menche permutation analysis: stretch goal.",
    }


def normalize_proximity_scores(results: list[dict]) -> list[dict]:
    """
    Normalize raw_proximity scores across all drugs to [0,1] via min-max.
    Call after all drugs have been scored.
    Modifies results in-place and returns them.
    """
    scores = [r["score"] for r in results if r.get("score") is not None]
    if not scores or max(scores) == min(scores):
        return results

    s_min, s_max = min(scores), max(scores)
    for r in results:
        if r.get("score") is not None:
            r["score_normalized"] = round((r["score"] - s_min) / (s_max - s_min), 3)
        else:
            r["score_normalized"] = 0.0
    return results


def _zero_result(reason: str) -> dict:
    return {
        "score": 0.0,
        "mean_distance": 10.0,
        "raw_proximity": 0.0,
        "drug_targets_in_graph": [],
        "disease_genes_in_graph": [],
        "note": reason,
    }
