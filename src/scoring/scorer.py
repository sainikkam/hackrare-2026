"""
Main scoring engine.

For each drug/disease pair:
  1. Gate 1: BBB Permeability
  2. Gate 2: Pediatric Safety
  3. (If both gates pass) Criteria 1, 2, 3
  4. Final Score = (c1 + c2 + c3) / 3  [arithmetic mean = cos(theta)*magnitude]

Evidence ledger: full breakdown of all scores, shared biology, safety flags.
"""

import logging
import networkx as nx
from src.gates.bbb_gate import score_bbb
from src.gates.pediatric_safety_gate import score_pediatric_safety
from src.criteria.pathway_alignment import score_pathway_alignment
from src.criteria.network_proximity import score_network_proximity
from src.criteria.phenotypic_overlap import score_phenotypic_overlap
from src.loaders.pubchem_loader import get_properties_by_name

logger = logging.getLogger(__name__)


def score_drug(
    drug: dict,
    disease_config: dict,
    disease_module_genes: list[str],
    string_graph: nx.Graph | None = None,
    disease_pathway_vector: dict[str, float] | None = None,
    disease_hpo_terms: list[str] | None = None,
) -> dict:
    """
    Full scoring pipeline for one drug against the root disease.

    drug: {chembl_id, name, phase, moa, ...}
    disease_config: loaded YAML config for the root disease
    disease_module_genes: list of gene symbols in the disease module
    string_graph: prebuilt STRING networkx graph (shared across all drugs for speed)
    disease_pathway_vector: precomputed Enrichr pathway vector for disease genes
    disease_hpo_terms: precomputed HPO terms for root disease

    Returns full evidence ledger dict.
    """
    drug_name = drug.get("name", "")
    chembl_id = drug.get("chembl_id", "")
    logger.info("Scoring %s (%s)", drug_name, chembl_id)

    # Fetch molecular properties once
    mol_props = get_properties_by_name(drug_name)

    # ── Gate 1: BBB ──────────────────────────────────────────────────────────
    bbb = score_bbb(drug_name, props=mol_props)
    if not bbb["passes"]:
        return _build_result(drug, bbb, None, None, None, None, "gate1_fail")

    # ── Gate 2: Pediatric Safety ─────────────────────────────────────────────
    safety = score_pediatric_safety(drug_name)
    if not safety["passes"]:
        return _build_result(drug, bbb, safety, None, None, None, "gate2_fail")

    # ── Criterion 1: Pathway Alignment ───────────────────────────────────────
    c1 = score_pathway_alignment(
        drug_name=drug_name,
        drug_chembl_id=chembl_id,
        disease_module_genes=disease_module_genes,
        disease_pathway_vector=disease_pathway_vector,
    )

    # ── Criterion 2: Network Proximity ───────────────────────────────────────
    # Use direct MOA targets (not expanded) for mechanistically accurate distance
    c2 = score_network_proximity(
        drug_target_genes=c1.get("drug_target_genes", []),
        disease_module_genes=disease_module_genes,
        G=string_graph,
    )

    # ── Criterion 3: Phenotypic Overlap ──────────────────────────────────────
    # Use expanded targets for richer pathway coverage in phenotypic comparison
    root_omim_ids = disease_config.get("omim_ids", [])
    drug_disease_omim = _infer_drug_disease_omim(drug)
    c3 = score_phenotypic_overlap(
        disease_a_omim_ids=root_omim_ids,
        disease_b_omim_ids=drug_disease_omim,
        disease_a_genes=disease_module_genes,
        disease_b_genes=c1.get("drug_target_genes_expanded") or c1.get("drug_target_genes", []),
        disease_a_pathway_vector=disease_pathway_vector,
    )

    return _build_result(drug, bbb, safety, c1, c2, c3, "scored")


def _build_result(drug, bbb, safety, c1, c2, c3, status) -> dict:
    drug_name = drug.get("name", "")
    chembl_id = drug.get("chembl_id", "")

    if status == "scored":
        final_score = round((c1["score"] + c2["score"] + c3["score"]) / 3.0, 3)
    else:
        final_score = None

    result = {
        "drug_name": drug_name,
        "chembl_id": chembl_id,
        "drug_phase": drug.get("phase", 0) or drug.get("max_phase", 0),
        "moa": drug.get("moa", ""),
        "status": status,
        "final_score": final_score,
        "gate1_bbb": bbb,
        "gate2_safety": safety,
        "criterion1_pathway": c1,
        "criterion2_network": c2,
        "criterion3_phenotypic": c3,
        # Convenience fields for UI
        "passes_gate1": bbb.get("passes", False),
        "passes_gate2": safety.get("passes", False) if safety else False,
        "caution": safety.get("caution", False) if safety else False,
        "caution_label": safety.get("caution_label") if safety else None,
    }

    if status == "scored":
        result["evidence_summary"] = {
            "final_score": final_score,
            "bbb_score": bbb.get("score"),
            "safety_score": safety.get("score"),
            "pathway_alignment": c1.get("score"),
            "network_proximity": c2.get("score"),
            "phenotypic_overlap": c3.get("score"),
            "top_shared_pathways": c1.get("shared_pathways", []),
            "top_shared_hpo_terms": c3.get("shared_hpo_terms", []),
            "shared_genes": c3.get("shared_genes", []),
            "drug_targets": c1.get("drug_target_genes_expanded") or c1.get("drug_target_genes", []),
            "mean_network_distance": c2.get("mean_distance"),
            "hpo_confidence": c3.get("confidence"),
        }

    return result


def _infer_drug_disease_omim(drug: dict) -> list[str]:
    """Best-effort: return empty list (no OMIM override — use gene-based overlap)."""
    return []


def run_sensitivity_analysis(results: list[dict], weight_perturbation: float = 0.1) -> dict:
    """
    Perturb criterion weights by ±perturbation and check top-3 ranking stability.
    Returns: {weight_set_name: [top3_drug_names]}
    """
    scored = [r for r in results if r.get("status") == "scored"]
    if not scored:
        return {}

    weight_sets = {
        "baseline":           (1/3, 1/3, 1/3),
        "pathway_emphasized": (1/3 + weight_perturbation, 1/3, 1/3 - weight_perturbation),
        "network_emphasized": (1/3, 1/3 + weight_perturbation, 1/3 - weight_perturbation),
        "phenotypic_emphasized": (1/3 - weight_perturbation, 1/3, 1/3 + weight_perturbation),
    }

    stability = {}
    for label, (w1, w2, w3) in weight_sets.items():
        ranked = sorted(
            scored,
            key=lambda r: (
                w1 * (r["criterion1_pathway"] or {}).get("score", 0) +
                w2 * (r["criterion2_network"] or {}).get("score", 0) +
                w3 * (r["criterion3_phenotypic"] or {}).get("score", 0)
            ),
            reverse=True,
        )
        stability[label] = [r["drug_name"] for r in ranked[:3]]

    # Check stability: are the same 3 drugs in top-3 across all weight sets?
    top3_sets = [set(names) for names in stability.values()]
    all_same = all(s == top3_sets[0] for s in top3_sets)
    stability["_stable"] = all_same
    stability["_baseline_top3"] = stability["baseline"]

    return stability
