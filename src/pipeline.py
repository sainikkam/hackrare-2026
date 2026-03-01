"""
Main pipeline orchestrator.

run_pipeline(config_path) -> tree dict

Steps:
  1. Load config
  2. Build disease module (OpenTargets + DisGeNET)
  3. Build STRING graph for disease module (shared across all drug scorings)
  4. Precompute disease pathway vector (Enrichr) — shared across Criteria 1 & 3
  5. Collect candidate drugs from Layer 2 diseases
  6. Score each candidate (gates + 3 criteria)
  7. Build tree JSON
  8. Validate against positive/negative controls
"""

import json
import logging
import os
from pathlib import Path

import yaml
from dotenv import load_dotenv

# Load .env before anything else
load_dotenv(Path(__file__).resolve().parents[1] / ".env")

from src.loaders.opentargets_loader import (
    resolve_disease_id,
    get_disease_associated_targets,
    get_drugs_for_disease,
)
from src.loaders.disgenet_loader import get_diseases_for_gene
from src.loaders.enrichr_loader import enrich_gene_list, get_pathway_score_vector
from src.loaders.hpo_loader import get_disease_hpo_terms
from src.criteria.network_proximity import get_disease_graph
from src.criteria.phenotypic_overlap import score_phenotypic_overlap, get_hpo_similarity
from src.scoring.scorer import score_drug, run_sensitivity_analysis
from src.scoring.tree_builder import build_tree, get_validation_summary

logger = logging.getLogger(__name__)

OUTPUT_DIR = Path(__file__).resolve().parents[1] / "data"


def run_pipeline(config_path: str, max_candidates: int | None = None) -> dict:
    """
    Run the full pipeline for a disease config file.

    Parameters:
        config_path: path to YAML config (e.g. 'config/tsc.yaml')
        max_candidates: override max_candidates from config (for testing)

    Returns:
        tree dict (JSON-serializable)
    """
    config = _load_config(config_path)
    disease_name = config["disease"]["name"]
    logger.info("=" * 60)
    logger.info("Starting pipeline for: %s", disease_name)
    logger.info("=" * 60)

    filters = config.get("filters", {})
    n_max = max_candidates or filters.get("max_candidates", 40)
    min_score = filters.get("opentargets_min_genetic_score", 0.3)
    string_conf = filters.get("string_confidence", 700)

    # ── Step 1: Build disease module ─────────────────────────────────────────
    logger.info("[Step 1] Building disease module for %s", disease_name)
    disease_module_genes = _build_disease_module(config, min_score)
    logger.info("Disease module: %d genes", len(disease_module_genes))

    if not disease_module_genes:
        logger.error("Empty disease module — check config and API keys")
        return {}

    # ── Step 2: Build STRING graph ───────────────────────────────────────────
    logger.info("[Step 2] Building STRING subgraph (confidence=%d)...", string_conf)
    string_graph = get_disease_graph(disease_module_genes, confidence=string_conf)
    logger.info("STRING graph: %d nodes, %d edges",
                string_graph.number_of_nodes(), string_graph.number_of_edges())

    # ── Step 3: Precompute disease pathway vector ────────────────────────────
    logger.info("[Step 3] Computing disease pathway vector (Enrichr)...")
    disease_enrichment = enrich_gene_list(disease_module_genes[:50])
    disease_pathway_vector = get_pathway_score_vector(disease_enrichment)
    logger.info("Disease pathway vector: %d terms", len(disease_pathway_vector))

    # ── Step 4: Collect candidate drugs from Layer 2 diseases ───────────────
    logger.info("[Step 4] Collecting candidate drugs from %d Layer 2 diseases...",
                len(config.get("layer2_diseases", [])))
    all_candidates, layer2_disease_drugs = _collect_candidates(config, n_max)
    logger.info("Total unique candidate drugs: %d", len(all_candidates))

    # ── Step 5: Score candidates ─────────────────────────────────────────────
    logger.info("[Step 5] Scoring %d candidates...", len(all_candidates))
    layer2_results: dict[str, list[dict]] = {}

    for disease_entry in config.get("layer2_diseases", []):
        dname = disease_entry["name"]
        drugs_for_disease = layer2_disease_drugs.get(dname, [])
        logger.info("  Scoring %d drugs for %s", len(drugs_for_disease), dname)

        results = []
        for drug in drugs_for_disease:
            result = score_drug(
                drug=drug,
                disease_config=config["disease"],
                disease_module_genes=disease_module_genes,
                string_graph=string_graph,
                disease_pathway_vector=disease_pathway_vector,
            )
            results.append(result)

        layer2_results[dname] = results

    # ── Step 6: Compute Layer 2 disease-to-root overlap scores ───────────────
    logger.info("[Step 6] Computing Layer 2 disease overlap scores...")
    layer2_overlap_scores = _compute_layer2_overlaps(config, disease_pathway_vector)

    # ── Step 7: Build tree ───────────────────────────────────────────────────
    logger.info("[Step 7] Building probability tree...")
    tree = build_tree(
        root_config=config["disease"],
        layer2_results=layer2_results,
        layer2_overlap_scores=layer2_overlap_scores,
    )

    # ── Step 8: Validate ─────────────────────────────────────────────────────
    logger.info("[Step 8] Validating against controls...")
    validation = get_validation_summary(
        tree,
        positive_controls=config.get("positive_controls", []),
        negative_controls=config.get("negative_controls", []),
    )
    tree["validation"] = validation
    _log_validation(validation)

    # Sensitivity analysis
    all_results = [r for results in layer2_results.values() for r in results]
    sensitivity = run_sensitivity_analysis(all_results)
    tree["sensitivity_analysis"] = sensitivity

    # ── Save output ──────────────────────────────────────────────────────────
    out_path = OUTPUT_DIR / f"tree_{_slugify(disease_name)}.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(tree, indent=2, default=str))
    logger.info("Tree saved to: %s", out_path)

    # Print summary
    _print_summary(tree)
    return tree


def _load_config(config_path: str) -> dict:
    path = Path(config_path)
    if not path.exists():
        path = Path(__file__).resolve().parents[1] / config_path
    with open(path) as f:
        return yaml.safe_load(f)


def _build_disease_module(config: dict, min_score: float) -> list[str]:
    """Collect disease module genes from OpenTargets + config-specified genes."""
    genes = set(config["disease"].get("genes", []))

    # OpenTargets
    search_query = config["disease"].get("opentargets_search", config["disease"]["name"])
    disease_id = resolve_disease_id(search_query)
    if disease_id:
        ot_targets = get_disease_associated_targets(disease_id, min_score=min_score, size=80)
        genes.update(t["gene_symbol"] for t in ot_targets)
        logger.info("OpenTargets: added %d genes (disease_id=%s)", len(ot_targets), disease_id)

    # DisGeNET (gene-centric — fetch diseases associated with seed genes)
    for entrez_id in config["disease"].get("entrez_ids", []):
        try:
            disgenet_results = get_diseases_for_gene(entrez_id, min_score=0.3)
            logger.info("DisGeNET: %d disease associations for gene %d", len(disgenet_results), entrez_id)
        except Exception as e:
            logger.warning("DisGeNET error for gene %d: %s", entrez_id, e)

    return sorted(genes)


def _collect_candidates(config: dict, max_per_disease: int) -> tuple[dict[str, dict], dict[str, list[dict]]]:
    """
    Collect candidate drugs for each Layer 2 disease.
    Returns:
        all_candidates: {drug_name: drug_dict} — deduplicated
        layer2_disease_drugs: {disease_name: [drug_dict, ...]}
    """
    all_candidates: dict[str, dict] = {}
    layer2_disease_drugs: dict[str, list[dict]] = {}

    for disease_entry in config.get("layer2_diseases", []):
        dname = disease_entry["name"]
        search_query = disease_entry.get("opentargets_search", dname)
        disease_id = resolve_disease_id(search_query)

        drugs = []
        if disease_id:
            ot_drugs = get_drugs_for_disease(disease_id)
            drugs.extend(ot_drugs)
            logger.info("  %s: %d drugs from OpenTargets", dname, len(ot_drugs))

        # Deduplicate per disease and cap at max_per_disease
        seen_in_disease = set()
        deduped = []
        for d in drugs:
            cid = d.get("chembl_id", "")
            if cid and cid not in seen_in_disease:
                seen_in_disease.add(cid)
                deduped.append(d)

        # Prioritize by clinical phase (higher phase = more validated)
        deduped.sort(key=lambda d: d.get("phase", 0) or d.get("max_phase", 0), reverse=True)
        deduped = deduped[:max_per_disease]

        layer2_disease_drugs[dname] = deduped

        for d in deduped:
            name = d.get("name", "").lower()
            if name and name not in all_candidates:
                all_candidates[name] = d

    # Add explicit controls from config (ensure they are included)
    for ctrl_list in [config.get("positive_controls", []), config.get("negative_controls", [])]:
        for ctrl in ctrl_list:
            name = ctrl["drug_name"].lower()
            if name not in all_candidates:
                all_candidates[name] = {"name": ctrl["drug_name"], "chembl_id": ctrl.get("chembl_id", ""), "phase": 4}
                # Add to first Layer 2 disease for scoring
                if layer2_disease_drugs:
                    first_disease = list(layer2_disease_drugs.keys())[0]
                    layer2_disease_drugs[first_disease].append(all_candidates[name])

    return all_candidates, layer2_disease_drugs


def _compute_layer2_overlaps(config: dict, disease_a_pathway_vector: dict) -> dict[str, dict]:
    """Compute phenotypic overlap score between root disease and each Layer 2 disease."""
    overlaps = {}
    root_omim_ids = config["disease"].get("omim_ids", [])
    root_genes = config["disease"].get("genes", [])

    for disease_entry in config.get("layer2_diseases", []):
        dname = disease_entry["name"]
        b_omim = [str(eid) for eid in disease_entry.get("entrez_ids", [])]  # Approximation
        b_genes = disease_entry.get("genes", [])
        b_omim_ids = []  # We don't have OMIM IDs for layer2 in config — use gene-based overlap

        try:
            overlap = score_phenotypic_overlap(
                disease_a_omim_ids=root_omim_ids,
                disease_b_omim_ids=b_omim_ids,
                disease_a_genes=root_genes,
                disease_b_genes=b_genes,
                disease_a_pathway_vector=disease_a_pathway_vector,
            )
            overlaps[dname] = overlap
        except Exception as e:
            logger.warning("Layer 2 overlap error for %s: %s", dname, e)
            overlaps[dname] = {"score": None}

    return overlaps


def _log_validation(validation: dict):
    logger.info("── Validation ──────────────────────────────────")
    for ctrl in validation.get("positive_controls", []):
        status = "✓ PASS" if ctrl["passes"] else "✗ FAIL"
        logger.info("  [+] %s: %s (score=%.3f, top3=%s)",
                    ctrl["drug"], status,
                    ctrl.get("final_score") or 0,
                    ctrl.get("in_top3", False))
    for ctrl in validation.get("negative_controls", []):
        status = "✓ PASS" if ctrl["passes"] else "✗ FAIL"
        logger.info("  [-] %s: %s (in_tree=%s, score=%s)",
                    ctrl["drug"], status, ctrl.get("in_tree"),
                    ctrl.get("final_score"))
    overall = "✓ ALL CONTROLS PASS" if validation.get("overall_pass") else "✗ SOME CONTROLS FAIL"
    logger.info("  Overall: %s", overall)


def _print_summary(tree: dict):
    summary = tree.get("summary", {})
    print("\n" + "=" * 60)
    print(f"PIPELINE COMPLETE: {summary.get('root_disease')}")
    print("=" * 60)
    print(f"  Candidates evaluated:  {summary.get('n_total_candidates', 0)}")
    print(f"  Passed Gate 1 (BBB):   {summary.get('n_gate1_pass', 0)}")
    print(f"  Passed Gate 2 (Safety):{summary.get('n_gate2_pass', 0)}")
    print(f"  Fully scored:          {summary.get('n_scored', 0)}")
    print("\nTop Candidates:")
    for i, drug in enumerate(summary.get("global_top_drugs", [])[:5], 1):
        caution = " ⚠ CAUTION" if drug.get("caution") else ""
        print(f"  {i}. {drug['drug_name']:<20} score={drug.get('final_score', 'N/A'):.3f}{caution}")
    val = tree.get("validation", {})
    print(f"\nValidation: {'✓ PASS' if val.get('overall_pass') else '✗ FAIL'}")
    print("=" * 60 + "\n")


def _slugify(s: str) -> str:
    return s.lower().replace(" ", "_").replace("/", "_").replace("-", "_")[:40]
