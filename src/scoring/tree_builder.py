"""
Tree builder: constructs the nested probability tree JSON from scored results.

Output structure:
  Layer 1: root disease node
  Layer 2: related disease nodes, each with ranked drug list
    - drugs are sorted by final_score (descending)
    - each drug carries its full evidence ledger

Layer 3 is deferred (documented as roadmap item).
"""

import logging

logger = logging.getLogger(__name__)


def build_tree(
    root_config: dict,
    layer2_results: dict[str, list[dict]],
    layer2_overlap_scores: dict[str, dict] | None = None,
) -> dict:
    """
    Build the probability tree JSON.

    root_config: disease YAML config (dict)
    layer2_results: {disease_name: [scored_drug_results, ...]}
    layer2_overlap_scores: {disease_name: phenotypic_overlap_score_dict}

    Returns nested dict (JSON-serializable).
    """
    root_disease = root_config.get("name", "Unknown Disease")

    layer2_nodes = []
    all_scored_drugs = []

    for disease_name, drug_results in layer2_results.items():
        # Separate scored vs. failed drugs
        scored = [r for r in drug_results if r.get("status") == "scored"]
        gate1_fail = [r for r in drug_results if r.get("status") == "gate1_fail"]
        gate2_fail = [r for r in drug_results if r.get("status") == "gate2_fail"]

        # Sort by final score
        scored_sorted = sorted(scored, key=lambda r: r.get("final_score", 0), reverse=True)
        all_scored_drugs.extend(scored_sorted)

        # Layer 2 disease node
        overlap = (layer2_overlap_scores or {}).get(disease_name, {})
        node = {
            "disease": disease_name,
            "overlap_score_to_root": overlap.get("score"),
            "hpo_similarity": overlap.get("hpo_similarity"),
            "pathway_similarity": overlap.get("pathway_similarity"),
            "n_candidates_evaluated": len(drug_results),
            "n_gate1_pass": len(drug_results) - len(gate1_fail),
            "n_gate2_pass": len(drug_results) - len(gate1_fail) - len(gate2_fail),
            "n_scored": len(scored),
            "drugs": [_format_drug_node(r) for r in scored_sorted],
            "gate1_failures": [{"drug_name": r["drug_name"], "bbb_score": r["gate1_bbb"].get("score")} for r in gate1_fail],
            "gate2_failures": [{"drug_name": r["drug_name"], "safety_score": r["gate2_safety"].get("score") if r.get("gate2_safety") else None} for r in gate2_fail],
        }
        layer2_nodes.append(node)

    # Sort Layer 2 nodes by overlap score to root (if available)
    layer2_nodes.sort(
        key=lambda n: n.get("overlap_score_to_root") or 0,
        reverse=True,
    )

    # Global top drugs (across all Layer 2 diseases)
    global_top = sorted(
        all_scored_drugs,
        key=lambda r: r.get("final_score", 0),
        reverse=True,
    )[:10]

    tree = {
        "root": {
            "disease": root_disease,
            "omim_ids": root_config.get("omim_ids", []),
            "genes": root_config.get("genes", []),
            "layer2": layer2_nodes,
        },
        "summary": {
            "root_disease": root_disease,
            "n_layer2_diseases": len(layer2_nodes),
            "n_total_candidates": sum(n["n_candidates_evaluated"] for n in layer2_nodes),
            "n_gate1_pass": sum(n["n_gate1_pass"] for n in layer2_nodes),
            "n_gate2_pass": sum(n["n_gate2_pass"] for n in layer2_nodes),
            "n_scored": sum(n["n_scored"] for n in layer2_nodes),
            "global_top_drugs": [_format_drug_node(r) for r in global_top],
        },
    }

    return tree


def _format_drug_node(result: dict) -> dict:
    """Flatten a scored drug result into a clean tree node."""
    ev = result.get("evidence_summary", {})
    safety = result.get("gate2_safety", {}) or {}
    bbb = result.get("gate1_bbb", {}) or {}
    c2 = result.get("criterion2_network", {}) or {}

    return {
        "drug_name": result.get("drug_name", ""),
        "chembl_id": result.get("chembl_id", ""),
        "drug_phase": result.get("drug_phase", 0),
        "moa": result.get("moa", ""),
        "final_score": result.get("final_score"),
        "caution": result.get("caution", False),
        "caution_label": result.get("caution_label"),
        "scores": {
            "bbb": bbb.get("score"),
            "safety": safety.get("score"),
            "pathway_alignment": ev.get("pathway_alignment"),
            "network_proximity": ev.get("network_proximity"),
            "phenotypic_overlap": ev.get("phenotypic_overlap"),
        },
        "evidence": {
            "top_shared_pathways": ev.get("top_shared_pathways", []),
            "top_shared_hpo_terms": ev.get("top_shared_hpo_terms", []),
            "shared_genes": ev.get("shared_genes", []),
            "drug_targets": ev.get("drug_targets", []),
            "mean_network_distance": c2.get("mean_distance"),
            "hpo_confidence": ev.get("hpo_confidence"),
            "scoring_method": ev.get("scoring_method"),
            "bbb_method": bbb.get("method"),
            "bbb_experimental_evidence": bbb.get("experimental_evidence"),
            "min_age_years": safety.get("min_age_years"),
            "pediatric_approved": safety.get("pediatric_approved"),
            "pediatric_ae_rate": (safety.get("ae_stats") or {}).get("pediatric_ae_rate"),
        },
    }


def get_validation_summary(tree: dict, positive_controls: list[dict], negative_controls: list[dict]) -> dict:
    """
    Check positive and negative controls against the final tree.
    Returns validation pass/fail for each control.
    """
    # Collect all scored drug names (lowercase)
    all_scored = {}
    for node in tree.get("root", {}).get("layer2", []):
        for drug in node.get("drugs", []):
            name = drug.get("drug_name", "").lower()
            if name not in all_scored or (drug.get("final_score") or 0) > (all_scored[name].get("final_score") or 0):
                all_scored[name] = drug

    # Global top 3
    top3 = sorted(all_scored.values(), key=lambda d: d.get("final_score") or 0, reverse=True)[:3]
    top3_names = {d["drug_name"].lower() for d in top3}

    validation = {"positive_controls": [], "negative_controls": [], "overall_pass": True}

    for ctrl in positive_controls:
        name = ctrl["drug_name"].lower()
        in_tree = name in all_scored
        in_top3 = name in top3_names
        score = all_scored.get(name, {}).get("final_score")
        passes = in_tree and score and score >= 0.4
        if not passes:
            validation["overall_pass"] = False
        validation["positive_controls"].append({
            "drug": ctrl["drug_name"],
            "in_tree": in_tree,
            "in_top3": in_top3,
            "final_score": score,
            "passes": passes,
            "expected": "Top 3 candidate",
        })

    for ctrl in negative_controls:
        name = ctrl["drug_name"].lower()
        in_tree = name in all_scored
        score = all_scored.get(name, {}).get("final_score")
        passes = not in_tree or (score is not None and score < 0.3)
        if not passes:
            validation["overall_pass"] = False
        validation["negative_controls"].append({
            "drug": ctrl["drug_name"],
            "in_tree": in_tree,
            "final_score": score,
            "passes": passes,
            "expected": "Excluded from tree or low score",
        })

    return validation
