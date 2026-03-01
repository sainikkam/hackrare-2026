"""
openFDA adverse event API loader.
Fetches pediatric serious adverse event counts for Gate 2 Component C.
No registration required. Rate limit: 240 req/min.
"""

import logging
import requests
from .cache import file_cache

logger = logging.getLogger(__name__)
BASE_URL = "https://api.fda.gov/drug/event.json"
TIMEOUT = 30


@file_cache("openfda")
def get_pediatric_ae_stats(drug_name: str) -> dict:
    """
    Fetch pediatric (age < 18) serious adverse event statistics for a drug.
    Returns: {total_reports, pediatric_serious, pediatric_ae_rate, ae_score}
    """
    # Total reports for the drug
    total = _count_drug_events(drug_name, serious_only=False, pediatric_only=False)
    if total == 0:
        return _default_ae_stats()

    # Pediatric serious reports (age groups 1-4: neonate, infant, child, adolescent)
    ped_serious = _count_drug_events(drug_name, serious_only=True, pediatric_only=True)

    # Compute normalized rate
    rate = ped_serious / total if total > 0 else 0.0

    # AE score: inverted so lower rate = higher score
    ae_score = max(0.0, 1.0 - min(rate * 10, 1.0))  # rate of 10% -> score 0; rate 0% -> score 1

    return {
        "total_reports": total,
        "pediatric_serious": ped_serious,
        "pediatric_ae_rate": round(rate, 4),
        "ae_score": round(ae_score, 3),
    }


def _count_drug_events(drug_name: str, serious_only: bool, pediatric_only: bool) -> int:
    try:
        search_parts = [f'patient.drug.openfda.generic_name:"{drug_name}"']
        if serious_only:
            search_parts.append("serious:1")
        if pediatric_only:
            # Age groups: 1=Neonate, 2=Infant, 3=Child, 4=Adolescent
            search_parts.append("patient.patientagegroup:(1+2+3+4)")

        params = {
            "search": " AND ".join(search_parts),
            "limit": 1,
        }
        resp = requests.get(BASE_URL, params=params, timeout=TIMEOUT)
        if resp.status_code == 404:
            return 0
        resp.raise_for_status()
        data = resp.json()
        return data.get("meta", {}).get("results", {}).get("total", 0)
    except Exception as e:
        logger.error("openFDA AE count error for %s: %s", drug_name, e)
        return 0


def _default_ae_stats() -> dict:
    """Default stats when no data is available — assume moderate safety (neutral)."""
    return {
        "total_reports": 0,
        "pediatric_serious": 0,
        "pediatric_ae_rate": 0.0,
        "ae_score": 0.5,  # Neutral when no data
    }


@file_cache("openfda")
def get_pediatric_trial_status(drug_name: str) -> dict:
    """
    Query ClinicalTrials.gov v2 API for pediatric trial status.
    Returns: {has_pediatric_trial, trial_phase, nct_ids}
    """
    try:
        url = "https://clinicaltrials.gov/api/v2/studies"
        params = {
            "query.term": drug_name,
            "query.patient": "CHILD",
            "filter.overallStatus": "RECRUITING,ACTIVE_NOT_RECRUITING,COMPLETED",
            "pageSize": 10,
        }
        resp = requests.get(url, params=params, timeout=TIMEOUT)
        resp.raise_for_status()
        data = resp.json()
        studies = data.get("studies", [])
        if not studies:
            return {"has_pediatric_trial": False, "trial_count": 0, "nct_ids": []}
        nct_ids = [s.get("protocolSection", {}).get("identificationModule", {}).get("nctId", "") for s in studies[:5]]
        return {
            "has_pediatric_trial": True,
            "trial_count": len(studies),
            "nct_ids": [n for n in nct_ids if n],
        }
    except Exception as e:
        logger.error("ClinicalTrials lookup error for %s: %s", drug_name, e)
        return {"has_pediatric_trial": False, "trial_count": 0, "nct_ids": []}
