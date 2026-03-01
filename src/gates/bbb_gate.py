"""
Gate 1: Blood-Brain Barrier Permeability
Pass threshold: score > 0.5 (composite of CNS MPO + clinical/experimental CNS evidence)

Methodology: Wager et al. CNS MPO (ACS Chem. Neurosci. 2010) + Pgp efflux adjustment.
For drugs with documented CNS efficacy (FDA neurological approval or preclinical CNS data),
a clinical evidence component is incorporated.

NOTE ON THRESHOLD: The original methodology specifies > 0.8 for physicochemical-only scores.
This implementation uses > 0.5 for the composite score, which incorporates clinical CNS
evidence alongside CNS MPO. This is a deliberate design choice for rare neurological disease
repurposing, where many relevant drugs (mTOR inhibitors, MEK inhibitors) achieve CNS efficacy
via mechanisms not captured by CNS MPO alone.
"""

import json
import logging
from pathlib import Path

from src.loaders.pubchem_loader import get_properties_by_name

logger = logging.getLogger(__name__)

BBB_THRESHOLD = 0.5
_EXP_BBB_PATH = Path(__file__).resolve().parents[2] / "data" / "experimental_bbb.json"
_exp_bbb_cache: dict | None = None


def _load_exp_bbb() -> dict:
    global _exp_bbb_cache
    if _exp_bbb_cache is None:
        if _EXP_BBB_PATH.exists():
            _exp_bbb_cache = json.loads(_EXP_BBB_PATH.read_text())
        else:
            _exp_bbb_cache = {}
    return _exp_bbb_cache


# Pgp substrate classification (rule-based + curated overrides)
_PGP_KNOWN_SUBSTRATES = {
    "everolimus", "sirolimus", "rapamycin", "tacrolimus", "cyclosporine",
    "imatinib", "dasatinib", "nilotinib", "lapatinib", "sunitinib",
    "gefitinib", "erlotinib", "paclitaxel", "docetaxel", "vincristine",
}
_PGP_KNOWN_NON_SUBSTRATES = {
    "metformin", "selumetinib", "trametinib", "binimetinib", "cobimetinib",
    "valproic acid", "levetiracetam", "memantine", "lithium",
    "vigabatrin", "ketamine", "minocycline",
}


def _pgp_adjustment(drug_name: str, props: dict | None) -> tuple[float, str]:
    """
    Return (adjustment_factor, classification_source).
    1.0 = not a substrate; 0.8 = predicted; 0.6 = confirmed substrate.
    """
    name_lower = drug_name.lower()
    if name_lower in _PGP_KNOWN_SUBSTRATES:
        return 0.6, "confirmed_substrate"
    if name_lower in _PGP_KNOWN_NON_SUBSTRATES:
        return 1.0, "confirmed_non_substrate"

    # Rule-based heuristic: large MW + high PSA suggests Pgp substrate
    if props:
        mw = props.get("mw", 0) or 0
        tpsa = props.get("tpsa", 0) or 0
        if mw > 400 and tpsa > 90:
            return 0.8, "predicted_substrate_rule"
        if mw > 500:
            return 0.8, "predicted_substrate_mw"

    return 1.0, "assumed_non_substrate"


def _cns_mpo_desirability(props: dict) -> tuple[float, dict]:
    """
    Compute CNS MPO score from Wager et al. 2010.
    Returns (normalized_score, component_scores).
    """
    def linear_desirability(value, optimal_max, hard_max):
        """1.0 at <= optimal_max, linearly decreasing to 0 at > hard_max."""
        if value is None:
            return 0.5  # Neutral when data missing
        if value <= optimal_max:
            return 1.0
        if value >= hard_max:
            return 0.0
        return (hard_max - value) / (hard_max - optimal_max)

    mw = props.get("mw") or props.get("xlogp")  # fallback handled below
    mw = props.get("mw")
    xlogp = props.get("xlogp")
    clogd = props.get("clogd")
    tpsa = props.get("tpsa")
    hbd = props.get("hbd")
    pka = props.get("pka_basic")

    components = {
        "mw":    linear_desirability(mw,    360.0, 500.0),
        "clogp": linear_desirability(xlogp,   3.0,   5.0),
        "clogd": linear_desirability(clogd,   2.0,   4.0),
        "tpsa":  linear_desirability(tpsa,   60.0,  90.0),
        "hbd":   linear_desirability(hbd,     0.0,   3.0),
        "pka":   linear_desirability(pka,     8.0,  10.0),
    }

    score = sum(components.values()) / 6.0
    return round(score, 3), components


def score_bbb(drug_name: str, props: dict | None = None) -> dict:
    """
    Compute BBB gate score for a drug.

    Returns dict with:
        score: float [0,1] — final BBB score
        passes: bool — True if score > BBB_THRESHOLD
        cns_mpo: float — raw CNS MPO score
        pgp_adjustment: float — Pgp multiplier applied
        pgp_source: str — classification source
        experimental_score: float | None — clinical/experimental CNS evidence score
        method: str — 'physicochemical' | 'composite'
        components: dict — individual CNS MPO desirability scores
    """
    # Fetch molecular properties if not provided
    if props is None:
        props = get_properties_by_name(drug_name)

    # CNS MPO component
    if props:
        mpo_score, components = _cns_mpo_desirability(props)
    else:
        mpo_score, components = 0.3, {}  # Unknown — assume poor without data

    # Pgp adjustment
    pgp_adj, pgp_source = _pgp_adjustment(drug_name, props)
    physicochemical_score = round(mpo_score * pgp_adj, 3)

    # Check experimental / clinical CNS evidence
    exp_bbb = _load_exp_bbb()
    exp_entry = exp_bbb.get(drug_name.lower())

    if exp_entry:
        exp_score = exp_entry["score"]
        # Composite: weight 0.35 physicochemical + 0.65 experimental
        final_score = round(0.35 * physicochemical_score + 0.65 * exp_score, 3)
        method = "composite"
    else:
        final_score = physicochemical_score
        exp_score = None
        method = "physicochemical"

    passes = final_score > BBB_THRESHOLD

    result = {
        "score": final_score,
        "passes": passes,
        "threshold": BBB_THRESHOLD,
        "cns_mpo": mpo_score,
        "pgp_adjustment": pgp_adj,
        "pgp_source": pgp_source,
        "physicochemical_score": physicochemical_score,
        "experimental_score": exp_score,
        "experimental_evidence": exp_entry.get("evidence", "") if exp_entry else None,
        "method": method,
        "components": components,
    }

    if not passes:
        logger.debug("Gate 1 FAIL: %s (score=%.3f)", drug_name, final_score)

    return result
