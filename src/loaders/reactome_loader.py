"""
Reactome REST API loader.
Fetches pathway memberships for gene sets.
No registration required.
"""

import logging
import requests
from .cache import file_cache

logger = logging.getLogger(__name__)
BASE_URL = "https://reactome.org/ContentService"
TIMEOUT = 30


@file_cache("reactome")
def get_pathways_for_gene(gene_symbol: str, species: str = "Homo sapiens") -> list[dict]:
    """
    Fetch Reactome pathways containing a given gene symbol.
    Returns list of {pathway_id, pathway_name, top_level}.
    """
    try:
        # First resolve gene symbol to Reactome entity
        url = f"{BASE_URL}/data/mapping/UniProt/{gene_symbol}/pathways"
        resp = requests.get(url, params={"species": species}, timeout=TIMEOUT)
        if resp.status_code == 404:
            # Try by gene name directly
            url2 = f"{BASE_URL}/data/query/enhanced/{gene_symbol}"
            resp = requests.get(url2, timeout=TIMEOUT)
            if resp.status_code == 404:
                return []
        resp.raise_for_status()
        data = resp.json()
        if isinstance(data, list):
            return [
                {
                    "pathway_id": p.get("stId", ""),
                    "pathway_name": p.get("displayName", ""),
                    "top_level": p.get("isInferred", False),
                }
                for p in data
            ]
        return []
    except Exception as e:
        logger.error("Reactome pathway error for %s: %s", gene_symbol, e)
        return []


@file_cache("reactome")
def get_pathways_for_gene_list(gene_symbols: list[str]) -> dict[str, list[str]]:
    """
    Batch fetch pathways for multiple genes.
    Returns dict: {gene_symbol: [pathway_name, ...]}
    """
    result = {}
    for gene in gene_symbols:
        pathways = get_pathways_for_gene(gene)
        result[gene] = [p["pathway_name"] for p in pathways]
    return result


@file_cache("reactome")
def get_pathway_genes(pathway_id: str) -> list[str]:
    """Fetch all genes in a Reactome pathway."""
    try:
        url = f"{BASE_URL}/data/pathway/{pathway_id}/participatingPhysicalEntities"
        resp = requests.get(url, timeout=TIMEOUT)
        resp.raise_for_status()
        data = resp.json()
        genes = []
        for entity in data:
            if "gene" in entity.get("className", "").lower():
                genes.append(entity.get("displayName", ""))
        return genes
    except Exception as e:
        logger.error("Reactome pathway genes error for %s: %s", pathway_id, e)
        return []


def get_combined_pathways_for_genes(gene_symbols: list[str]) -> set[str]:
    """Return the union of all Reactome pathway names for a list of genes."""
    all_pathways = set()
    for gene in gene_symbols:
        pathways = get_pathways_for_gene(gene)
        all_pathways.update(p["pathway_name"] for p in pathways)
    return all_pathways
