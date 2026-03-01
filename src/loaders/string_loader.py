"""
STRING REST API loader.
Builds a PPI subgraph for the disease module (2-hop neighborhood).
No registration required. No bulk download needed (~500-800 proteins via API).
"""

import logging
import requests
import networkx as nx
from .cache import file_cache

logger = logging.getLogger(__name__)
BASE_URL = "https://string-db.org/api/json"
TIMEOUT = 60
SPECIES_HUMAN = 9606


@file_cache("string")
def get_interaction_partners(gene_symbols: list[str], confidence: int = 700, limit: int = 500) -> list[dict]:
    """
    Fetch STRING interaction partners for a list of gene symbols.
    confidence: 0-1000 (700 = high confidence >= 0.7)
    Returns list of {gene_a, gene_b, score}.
    """
    if not gene_symbols:
        return []
    try:
        identifiers = "\r".join(gene_symbols)
        resp = requests.post(
            f"{BASE_URL}/interaction_partners",
            data={
                "identifiers": identifiers,
                "species": SPECIES_HUMAN,
                "required_score": confidence,
                "limit": limit,
            },
            timeout=TIMEOUT,
        )
        resp.raise_for_status()
        data = resp.json()
        interactions = []
        for item in data:
            interactions.append({
                "gene_a": item.get("preferredName_A", ""),
                "gene_b": item.get("preferredName_B", ""),
                "score": item.get("score", 0),
                "string_id_a": item.get("stringId_A", ""),
                "string_id_b": item.get("stringId_B", ""),
            })
        return interactions
    except Exception as e:
        logger.error("STRING API error for %s: %s", gene_symbols, e)
        return []


@file_cache("string")
def get_network_2hop(seed_genes: list[str], confidence: int = 700) -> list[dict]:
    """
    Build 2-hop neighborhood: first fetch partners of seeds, then partners of those partners.
    Returns all interactions as edge list.
    """
    logger.info("STRING: fetching 1-hop neighbors for %d seed genes", len(seed_genes))
    hop1 = get_interaction_partners(seed_genes, confidence=confidence, limit=300)

    hop1_genes = set()
    for e in hop1:
        hop1_genes.add(e["gene_a"])
        hop1_genes.add(e["gene_b"])
    hop1_genes -= set(seed_genes)

    logger.info("STRING: found %d 1-hop partners, fetching 2-hop", len(hop1_genes))

    # Limit 2-hop expansion to keep graph manageable
    hop1_sample = list(hop1_genes)[:100]
    hop2 = get_interaction_partners(hop1_sample, confidence=confidence, limit=100) if hop1_sample else []

    all_interactions = hop1 + hop2
    # Deduplicate
    seen = set()
    unique = []
    for e in all_interactions:
        key = tuple(sorted([e["gene_a"], e["gene_b"]]))
        if key not in seen:
            seen.add(key)
            unique.append(e)

    logger.info("STRING: total unique interactions: %d", len(unique))
    return unique


def build_graph(interactions: list[dict]) -> nx.Graph:
    """Build a networkx graph from STRING interaction list."""
    G = nx.Graph()
    for e in interactions:
        if e["gene_a"] and e["gene_b"]:
            G.add_edge(e["gene_a"], e["gene_b"], weight=e["score"] / 1000.0)
    return G


def compute_mean_distance(G: nx.Graph, source_genes: list[str], target_genes: list[str]) -> float:
    """
    Compute mean shortest-path distance from source_genes to target_genes in G.
    Nodes not in the graph are skipped. Returns large value if no path found.
    """
    sources = [g for g in source_genes if g in G]
    targets = [g for g in target_genes if g in G]

    if not sources or not targets:
        return 10.0  # Max distance — not connected

    distances = []
    for s in sources:
        for t in targets:
            if s == t:
                distances.append(0)
                continue
            try:
                d = nx.shortest_path_length(G, s, t)
                distances.append(d)
            except nx.NetworkXNoPath:
                distances.append(10)

    return sum(distances) / len(distances) if distances else 10.0
