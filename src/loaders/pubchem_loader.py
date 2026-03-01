"""
PubChem REST API loader.
Fetches molecular properties needed for CNS MPO computation:
MW, XLogP3, TPSA, HBD count, pKa (most basic), cLogD.
No registration required.
"""

import logging
import requests
from .cache import file_cache

logger = logging.getLogger(__name__)
BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
TIMEOUT = 30

# Properties to fetch for CNS MPO
_PROPS = "MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount,CanonicalSMILES,IUPACName"


@file_cache("pubchem")
def get_properties_by_name(drug_name: str) -> dict | None:
    """Fetch molecular properties by drug name (INN or brand)."""
    try:
        url = f"{BASE_URL}/compound/name/{requests.utils.quote(drug_name)}/property/{_PROPS}/JSON"
        resp = requests.get(url, timeout=TIMEOUT)
        if resp.status_code == 404:
            return None
        resp.raise_for_status()
        table = resp.json().get("PropertyTable", {}).get("Properties", [])
        if not table:
            return None
        p = table[0]
        return _parse_props(p, drug_name)
    except Exception as e:
        logger.error("PubChem name lookup error for %s: %s", drug_name, e)
        return None


@file_cache("pubchem")
def get_properties_by_cid(cid: int) -> dict | None:
    """Fetch molecular properties by PubChem CID."""
    try:
        url = f"{BASE_URL}/compound/cid/{cid}/property/{_PROPS}/JSON"
        resp = requests.get(url, timeout=TIMEOUT)
        resp.raise_for_status()
        table = resp.json().get("PropertyTable", {}).get("Properties", [])
        if not table:
            return None
        return _parse_props(table[0], str(cid))
    except Exception as e:
        logger.error("PubChem CID lookup error for %s: %s", cid, e)
        return None


@file_cache("pubchem")
def get_cid_by_name(drug_name: str) -> int | None:
    """Resolve a drug name to a PubChem CID."""
    try:
        url = f"{BASE_URL}/compound/name/{requests.utils.quote(drug_name)}/cids/JSON"
        resp = requests.get(url, timeout=TIMEOUT)
        if resp.status_code == 404:
            return None
        resp.raise_for_status()
        cids = resp.json().get("IdentifierList", {}).get("CID", [])
        return cids[0] if cids else None
    except Exception as e:
        logger.error("PubChem CID resolve error for %s: %s", drug_name, e)
        return None


def _parse_props(p: dict, source: str) -> dict:
    mw = p.get("MolecularWeight")
    xlogp = p.get("XLogP")
    tpsa = p.get("TPSA")
    hbd = p.get("HBondDonorCount")
    hba = p.get("HBondAcceptorCount")

    # Estimate cLogD at pH 7.4 from XLogP + pKa heuristic
    # For a basic amine (most common CNS drug), cLogD ≈ cLogP - 1.0 as rough estimate
    # For neutral/acidic drugs, cLogD ≈ cLogP
    # Without pKa data, use XLogP - 0.5 as a conservative estimate
    clogd = (float(xlogp) - 0.5) if xlogp is not None else None

    # Most basic pKa: not directly available from PubChem.
    # Use heuristic: if HBA > 0 (likely basic nitrogen), estimate pKa ~ 8.5
    # If no HBA, drug is likely neutral/acidic, pKa ~ 1 (passes CNS MPO criterion)
    pka_basic = 8.5 if (hba or 0) > 0 else 1.0

    return {
        "source": source,
        "mw": float(mw) if mw is not None else None,
        "xlogp": float(xlogp) if xlogp is not None else None,
        "clogd": clogd,
        "tpsa": float(tpsa) if tpsa is not None else None,
        "hbd": int(hbd) if hbd is not None else None,
        "hba": int(hba) if hba is not None else None,
        "pka_basic": pka_basic,
        "smiles": p.get("CanonicalSMILES", ""),
        "iupac_name": p.get("IUPACName", ""),
        "cid": p.get("CID"),
    }
