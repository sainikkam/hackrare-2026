"""
Gate 2: Pediatric Safety Compatibility
Pass threshold: score >= 0.6

Scoring:
  Hard disqualifiers (auto-fail, score = 0):
    - Black box warning referencing pediatric patients
    - Explicit pediatric contraindication

  Component A: FDA Pediatric Approval Status   (weight: 0.40)
  Component B: Minimum Approved Age            (weight: 0.30)
  Component C: Serious Pediatric AE Rate       (weight: 0.30)

  Composite: 0.40*A + 0.30*B + 0.30*C

Sources: DailyMed SPL XML, openFDA adverse events, ClinicalTrials.gov v2
"""

import logging
from src.loaders.dailymed_loader import get_label_text
from src.loaders.openfda_loader import get_pediatric_ae_stats, get_pediatric_trial_status

logger = logging.getLogger(__name__)
SAFETY_THRESHOLD = 0.6


def score_pediatric_safety(drug_name: str) -> dict:
    """
    Compute pediatric safety gate score for a drug.

    Returns dict with:
        score: float [0,1]
        passes: bool
        caution: bool  — True if 0.60 <= score < 0.80
        hard_fail: bool
        hard_fail_reason: str | None
        component_a: float
        component_b: float
        component_c: float
        min_age_years: float | None
        pediatric_approved: bool
        ae_stats: dict
        trial_status: dict
        label_data: dict
    """
    # Fetch label data
    label = get_label_text(drug_name)

    # --- Hard disqualifier check ---
    if label.get("has_pediatric_bb_warning"):
        return _fail_result(drug_name, "Black box warning references pediatric patients", label)
    if label.get("has_pediatric_contraindication"):
        return _fail_result(drug_name, "Explicit pediatric contraindication in label", label)

    # --- Component A: FDA Pediatric Approval Status ---
    comp_a = _score_component_a(label, drug_name)

    # --- Component B: Minimum Approved Age ---
    comp_b = _score_component_b(label)

    # --- Component C: Serious Pediatric AE Rate ---
    ae_stats = get_pediatric_ae_stats(drug_name)
    comp_c = ae_stats.get("ae_score", 0.5)

    # --- Composite ---
    score = round(0.40 * comp_a + 0.30 * comp_b + 0.30 * comp_c, 3)
    passes = score >= SAFETY_THRESHOLD
    caution = passes and score < 0.80

    trial_status = get_pediatric_trial_status(drug_name)

    return {
        "score": score,
        "passes": passes,
        "caution": caution,
        "caution_label": "Caution: limited pediatric data" if caution else None,
        "hard_fail": False,
        "hard_fail_reason": None,
        "component_a": round(comp_a, 3),
        "component_b": round(comp_b, 3),
        "component_c": round(comp_c, 3),
        "min_age_years": label.get("min_age_years"),
        "pediatric_approved": label.get("pediatric_approved", False),
        "ae_stats": ae_stats,
        "trial_status": trial_status,
        "threshold": SAFETY_THRESHOLD,
        "label_data": {
            "boxed_warning_excerpt": label.get("boxed_warning", "")[:200],
            "pediatric_use_excerpt": label.get("pediatric_use", "")[:400],
        },
    }


def _score_component_a(label: dict, drug_name: str) -> float:
    """
    Component A: FDA Pediatric Approval Status
    1.0 = approved >= 1 pediatric group (not just adolescents)
    0.6 = adolescent only (>= 12 years)
    0.4 = active pediatric trial
    0.3 = no approval, no contraindication
    """
    if label.get("pediatric_approved"):
        min_age = label.get("min_age_years")
        if min_age is not None and min_age >= 12:
            return 0.6  # Adolescent only
        return 1.0

    # Check for active pediatric trial via ClinicalTrials
    trial = get_pediatric_trial_status(drug_name)
    if trial.get("has_pediatric_trial"):
        return 0.4

    return 0.3


def _score_component_b(label: dict) -> float:
    """
    Component B: Minimum Approved Age (inverse — younger = higher score)
    """
    min_age = label.get("min_age_years")
    if min_age is None:
        return 0.2  # No age data — conservative

    if min_age <= 2:
        return 1.0
    elif min_age <= 5:
        return 0.8
    elif min_age <= 11:
        return 0.6
    elif min_age <= 17:
        return 0.4
    else:
        return 0.1  # Adults only


def _fail_result(drug_name: str, reason: str, label: dict) -> dict:
    logger.debug("Gate 2 FAIL (hard disqualifier): %s — %s", drug_name, reason)
    return {
        "score": 0.0,
        "passes": False,
        "caution": False,
        "caution_label": None,
        "hard_fail": True,
        "hard_fail_reason": reason,
        "component_a": 0.0,
        "component_b": 0.0,
        "component_c": 0.0,
        "min_age_years": None,
        "pediatric_approved": False,
        "ae_stats": {},
        "trial_status": {},
        "threshold": SAFETY_THRESHOLD,
        "label_data": {
            "boxed_warning_excerpt": label.get("boxed_warning", "")[:200],
            "pediatric_use_excerpt": "",
        },
    }
