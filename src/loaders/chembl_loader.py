"""
ChEMBL REST API loader.
Fetches: drug targets (high-confidence pChEMBL >= threshold), drug properties,
drug indications, and target-to-pathway mappings.
No registration required. Rate limit: ~10 req/s.
"""

import time
import logging
import requests
from .cache import file_cache

logger = logging.getLogger(__name__)
BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
TIMEOUT = 30
_last_call = 0.0


def _get(endpoint: str, params: dict) -> dict:
    global _last_call
    elapsed = time.time() - _last_call
    if elapsed < 0.12:
        time.sleep(0.12 - elapsed)
    resp = requests.get(f"{BASE_URL}/{endpoint}.json", params=params, timeout=TIMEOUT)
    _last_call = time.time()
    resp.raise_for_status()
    return resp.json()


@file_cache("chembl")
def get_drug_info(chembl_id: str) -> dict | None:
    """Fetch basic drug info: name, molecular properties, max_phase."""
    try:
        data = _get(f"molecule/{chembl_id}", {})
        props = data.get("molecule_properties") or {}
        return {
            "chembl_id": chembl_id,
            "name": data.get("pref_name", chembl_id),
            "max_phase": data.get("max_phase", 0),
            "mw": props.get("mw_freebase"),
            "xlogp": props.get("cx_logp"),
            "tpsa": props.get("psa"),
            "hbd": props.get("hbd"),
            "hba": props.get("hba"),
            "alogp": props.get("alogp"),
            "ro5_violations": props.get("num_ro5_violations"),
            "smiles": data.get("molecule_structures", {}).get("canonical_smiles", ""),
        }
    except Exception as e:
        logger.error("ChEMBL drug info error for %s: %s", chembl_id, e)
        return None


@file_cache("chembl")
def get_targets_for_drug(chembl_id: str, min_pchembl: float = 7.0) -> list[dict]:
    """Fetch protein targets for a drug with pChEMBL >= min_pchembl."""
    results = []
    try:
        # Get all activities for this compound
        data = _get("activity", {
            "molecule_chembl_id": chembl_id,
            "pchembl_value__gte": min_pchembl,
            "limit": 100,
            "assay_type": "B",  # Binding assays only
        })
        for item in data.get("activities", []):
            target_id = item.get("target_chembl_id")
            if target_id and item.get("pchembl_value"):
                results.append({
                    "target_chembl_id": target_id,
                    "pchembl_value": float(item["pchembl_value"]),
                    "assay_type": item.get("assay_type", ""),
                })
        # Deduplicate, keep highest pChEMBL per target
        best = {}
        for r in results:
            tid = r["target_chembl_id"]
            if tid not in best or r["pchembl_value"] > best[tid]["pchembl_value"]:
                best[tid] = r
        return list(best.values())
    except Exception as e:
        logger.error("ChEMBL targets error for %s: %s", chembl_id, e)
        return []


@file_cache("chembl")
def get_gene_symbols_for_target(target_chembl_id: str) -> list[str]:
    """Resolve a ChEMBL target ID to gene symbols."""
    try:
        data = _get(f"target/{target_chembl_id}", {})
        components = data.get("target_components", [])
        symbols = []
        for comp in components:
            for xref in comp.get("target_component_xrefs", []):
                if xref.get("xref_src_db") == "UniProt":
                    pass  # Use gene_name from component directly
            gene = comp.get("component_description", "")
            # Extract gene symbol from synonyms
            for syn in comp.get("target_component_synonyms", []):
                if syn.get("syn_type") in ("GENE_SYMBOL", "GENE_NAME"):
                    symbols.append(syn["component_synonym"])
                    break
        return list(set(symbols))
    except Exception as e:
        logger.error("ChEMBL target resolve error for %s: %s", target_chembl_id, e)
        return []


@file_cache("chembl")
def get_drugs_for_indication(disease_name: str, min_phase: int = 2) -> list[dict]:
    """Search ChEMBL drug indications by disease name."""
    try:
        data = _get("drug_indication", {
            "mesh_heading__icontains": disease_name,
            "max_phase_for_ind__gte": min_phase,
            "limit": 100,
        })
        drugs = []
        seen = set()
        for item in data.get("drug_indications", []):
            cid = item.get("molecule_chembl_id")
            if cid and cid not in seen:
                seen.add(cid)
                drugs.append({
                    "chembl_id": cid,
                    "indication": item.get("mesh_heading", ""),
                    "phase": item.get("max_phase_for_ind", 0),
                    "efo_id": item.get("efo_id", ""),
                })
        return drugs
    except Exception as e:
        logger.error("ChEMBL indication error for %s: %s", disease_name, e)
        return []


@file_cache("chembl")
def search_drug_by_name(name: str) -> dict | None:
    """Search for a drug by preferred name or synonym."""
    try:
        data = _get("molecule", {"pref_name__iexact": name, "limit": 1})
        mols = data.get("molecules", [])
        if mols:
            return {"chembl_id": mols[0]["molecule_chembl_id"], "name": mols[0].get("pref_name", name)}
        # Try synonym search
        data2 = _get("molecule", {"molecule_synonyms__molecule_synonym__iexact": name, "limit": 1})
        mols2 = data2.get("molecules", [])
        if mols2:
            return {"chembl_id": mols2[0]["molecule_chembl_id"], "name": mols2[0].get("pref_name", name)}
        return None
    except Exception as e:
        logger.error("ChEMBL search error for %s: %s", name, e)
        return None
