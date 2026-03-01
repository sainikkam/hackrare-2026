"""
DisGeNET REST API loader.
Fetches gene-disease associations to augment disease module construction.
Requires DISGENET_API_KEY environment variable.
Registration: https://www.disgenet.org/signup (free, <5 min)
"""

import os
import logging
import requests
from .cache import file_cache

logger = logging.getLogger(__name__)
BASE_URL = "https://www.disgenet.org/api"
TIMEOUT = 30


def _get_headers() -> dict:
    key = os.environ.get("DISGENET_API_KEY", "")
    if not key:
        logger.warning("DISGENET_API_KEY not set — DisGeNET calls will be skipped")
    return {"Authorization": f"Bearer {key}", "accept": "application/json"}


@file_cache("disgenet")
def get_genes_for_disease(disease_id: str, min_score: float = 0.3) -> list[dict]:
    """
    Fetch gene-disease associations for a DisGeNET disease concept ID (e.g. 'C0041341').
    Returns list of {gene_symbol, gene_id, score}.
    """
    key = os.environ.get("DISGENET_API_KEY", "")
    if not key:
        return []

    url = f"{BASE_URL}/gda/disease/{disease_id}"
    params = {"format": "json", "min_score": min_score}
    try:
        resp = requests.get(url, headers=_get_headers(), params=params, timeout=TIMEOUT)
        if resp.status_code == 401:
            logger.error("DisGeNET: invalid API key")
            return []
        resp.raise_for_status()
        # DisGeNET REST API sometimes returns HTML (API restructuring).
        # Detect and fall back gracefully.
        ct = resp.headers.get("content-type", "")
        if "html" in ct or resp.text.strip().startswith("<!"):
            logger.warning("DisGeNET returned HTML (API may have changed) — skipping")
            return []
        data = resp.json()
        return [
            {
                "gene_symbol": item.get("gene_symbol", ""),
                "gene_id": item.get("geneid", ""),
                "score": item.get("score", 0),
                "disease_name": item.get("disease_name", ""),
            }
            for item in (data if isinstance(data, list) else [])
            if item.get("score", 0) >= min_score
        ]
    except Exception as e:
        logger.error("DisGeNET error for %s: %s", disease_id, e)
        return []


@file_cache("disgenet")
def get_diseases_for_gene(entrez_id: int, min_score: float = 0.3) -> list[dict]:
    """
    Fetch disease associations for a given Entrez gene ID.
    Returns list of {disease_id, disease_name, score}.
    """
    key = os.environ.get("DISGENET_API_KEY", "")
    if not key:
        return []

    url = f"{BASE_URL}/gda/gene/{entrez_id}"
    params = {"format": "json", "min_score": min_score}
    try:
        resp = requests.get(url, headers=_get_headers(), params=params, timeout=TIMEOUT)
        resp.raise_for_status()
        ct = resp.headers.get("content-type", "")
        if "html" in ct or resp.text.strip().startswith("<!"):
            logger.warning("DisGeNET returned HTML — skipping gene %s", entrez_id)
            return []
        data = resp.json()
        return [
            {
                "disease_id": item.get("diseaseid", ""),
                "disease_name": item.get("disease_name", ""),
                "score": item.get("score", 0),
            }
            for item in (data if isinstance(data, list) else [])
            if item.get("score", 0) >= min_score
        ]
    except Exception as e:
        logger.error("DisGeNET error for gene %s: %s", entrez_id, e)
        return []
