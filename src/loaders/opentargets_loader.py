"""
OpenTargets GraphQL API loader.
Fetches: disease module genes, drugs by disease, disease search by name.
No registration required.
"""

import logging
import requests
from .cache import file_cache

logger = logging.getLogger(__name__)
OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"
TIMEOUT = 30


def _gql(query: str, variables: dict) -> dict:
    resp = requests.post(OT_URL, json={"query": query, "variables": variables}, timeout=TIMEOUT)
    resp.raise_for_status()
    return resp.json().get("data", {})


@file_cache("opentargets")
def search_disease(query_string: str) -> list[dict]:
    """Search for a disease by name; returns list of {id, name} hits."""
    q = """
    query SearchDisease($q: String!) {
      search(queryString: $q, entityNames: ["disease"], page: {index: 0, size: 5}) {
        hits { id name entity }
      }
    }
    """
    data = _gql(q, {"q": query_string})
    hits = data.get("search", {}).get("hits", [])
    return [{"id": h["id"], "name": h["name"]} for h in hits if h.get("entity") == "disease"]


@file_cache("opentargets")
def get_disease_associated_targets(disease_id: str, min_score: float = 0.3, size: int = 100) -> list[dict]:
    """Returns gene symbols associated with a disease above min_score."""
    q = """
    query DiseaseTargets($id: String!, $size: Int!) {
      disease(efoId: $id) {
        associatedTargets(page: {index: 0, size: $size}) {
          rows {
            target { id approvedSymbol }
            score
          }
        }
      }
    }
    """
    data = _gql(q, {"id": disease_id, "size": size})
    rows = data.get("disease", {}).get("associatedTargets", {}).get("rows", [])
    return [
        {"ensembl_id": r["target"]["id"], "gene_symbol": r["target"]["approvedSymbol"], "score": r["score"]}
        for r in rows if r["score"] >= min_score
    ]


@file_cache("opentargets")
def get_drugs_for_disease(disease_id: str) -> list[dict]:
    """Returns known drugs for a disease with their targets."""
    q = """
    query DiseaseKnownDrugs($id: String!) {
      disease(efoId: $id) {
        knownDrugs(size: 100) {
          rows {
            drug { id name }
            phase
            mechanismOfAction
            approvedName
          }
        }
      }
    }
    """
    data = _gql(q, {"id": disease_id})
    rows = data.get("disease", {}).get("knownDrugs", {}).get("rows", [])
    drugs = []
    seen = set()
    for r in rows:
        drug = r.get("drug", {})
        cid = drug.get("id")
        if cid and cid not in seen:
            seen.add(cid)
            drugs.append({
                "chembl_id": cid,
                "name": drug.get("name", ""),
                "phase": r.get("phase", 0),
                "moa": r.get("mechanismOfAction", ""),
                "target_name": r.get("approvedName", ""),
                "max_phase": r.get("phase", 0),
            })
    return drugs


@file_cache("opentargets")
def get_target_associated_diseases(ensembl_id: str, size: int = 20) -> list[dict]:
    """Returns diseases associated with a target gene — used for Layer 3 expansion."""
    q = """
    query TargetDiseases($id: String!, $size: Int!) {
      target(ensemblId: $id) {
        associatedDiseases(page: {index: 0, size: $size}) {
          rows {
            disease { id name }
            score
          }
        }
      }
    }
    """
    data = _gql(q, {"id": ensembl_id, "size": size})
    rows = data.get("target", {}).get("associatedDiseases", {}).get("rows", [])
    return [{"disease_id": r["disease"]["id"], "disease_name": r["disease"]["name"], "score": r["score"]} for r in rows]


@file_cache("opentargets")
def get_drug_target_genes(chembl_id: str) -> list[str]:
    """
    Fetch gene symbols targeted by a drug from OpenTargets drug MOA.
    More reliable than ChEMBL activity endpoint for hackathon use.
    Returns list of gene symbols.
    """
    q = """
    query DrugTargets($id: String!) {
      drug(chemblId: $id) {
        name
        mechanismsOfAction {
          rows {
            targets { id approvedSymbol }
            mechanismOfAction
          }
        }
      }
    }
    """
    try:
        data = _gql(q, {"id": chembl_id})
        drug = data.get("drug") or {}
        rows = (drug.get("mechanismsOfAction") or {}).get("rows", [])
        genes = []
        for row in rows:
            for target in row.get("targets", []):
                sym = target.get("approvedSymbol")
                if sym:
                    genes.append(sym)
        return list(set(genes))
    except Exception as e:
        logger.error("OpenTargets drug targets error for %s: %s", chembl_id, e)
        return []


def resolve_disease_id(search_query: str) -> str | None:
    """Resolve a disease name to its OpenTargets EFO/MONDO ID. Returns first hit ID."""
    hits = search_disease(search_query)
    if hits:
        logger.info("Resolved '%s' -> %s (%s)", search_query, hits[0]["id"], hits[0]["name"])
        return hits[0]["id"]
    logger.warning("Could not resolve disease: %s", search_query)
    return None
