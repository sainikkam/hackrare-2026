#!/usr/bin/env python3
"""
webapp.py — Web interface for the Rare Disease Drug Repurposing Navigator.

Run:
    pip install streamlit pyvis requests networkx
    streamlit run webapp.py

Uses the same algorithm as repurpose.py — no API keys required.
"""

import json
import sys
import tempfile
import os
from pathlib import Path

import streamlit as st
import pandas as pd

# Import all algorithm functions from repurpose.py (same directory)
from repurpose import (
    ot_search_disease,
    ot_resolve_disease,
    ot_disease_genes,
    ot_disease_drugs,
    build_ppi_graph,
    enrichr_pathway_vector,
    score_drug,
    magnitude_aware_cosine,
    DEFAULT_PANEL,
    BBB_THRESHOLD,
    SAFETY_THRESHOLD,
    CACHE_DIR,
)

st.set_page_config(
    page_title="HackRare 2026 — Drug Repurposing Navigator",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)

# ──────────────────────────────────────────────────────────────────────────────
# Helpers
# ──────────────────────────────────────────────────────────────────────────────

def fmt(v) -> str:
    return f"{float(v):.3f}" if v is not None else "—"


def score_color(score) -> str:
    if score is None:
        return "#888888"
    score = float(score)
    if score >= 0.7:
        return "#27ae60"
    if score >= 0.5:
        return "#f39c12"
    if score >= 0.3:
        return "#e67e22"
    return "#e74c3c"


def _slugify(s: str) -> str:
    return s.lower().replace(" ", "_").replace("/", "_").replace("-", "_")[:40]


# ──────────────────────────────────────────────────────────────────────────────
# Pipeline runner (wraps repurpose.py functions for the web UI)
# ──────────────────────────────────────────────────────────────────────────────

def run_pipeline_web(disease_name: str, progress_callback=None) -> dict:
    """Run the full scoring pipeline, returning a structured result dict for the UI."""

    def update(step, total, msg):
        if progress_callback:
            progress_callback(step / total, msg)

    update(0, 6, "Searching OpenTargets...")
    hits = ot_search_disease(disease_name)
    if hits:
        disease_id = hits[0]["id"]
        matched_name = hits[0]["name"]
    else:
        disease_id = None
        matched_name = disease_name

    genes = []
    if disease_id:
        genes = ot_disease_genes(disease_id)

    update(1, 6, f"Building PPI network ({len(genes)} genes)...")
    import networkx as nx
    G = build_ppi_graph(genes, confidence=700) if genes else nx.Graph()

    update(2, 6, "Computing disease pathway signature...")
    disease_pathway_vec = enrichr_pathway_vector(genes[:50])

    update(3, 6, "Collecting candidate drugs from comparison panel...")
    all_drugs: dict = {}
    panel_names = []
    panel_drug_map: dict[str, list] = {}  # panel_disease -> [drug_dict, ...]
    for entry in DEFAULT_PANEL:
        pname = entry["name"]
        panel_names.append(pname)
        did = ot_resolve_disease(entry["opentargets_search"])
        if not did:
            panel_drug_map[pname] = []
            continue
        drugs_for_panel = []
        for d in ot_disease_drugs(did):
            key = d["name"].lower()
            if key and key not in all_drugs:
                all_drugs[key] = d
                drugs_for_panel.append(d)
            elif key in all_drugs:
                drugs_for_panel.append(all_drugs[key])
        panel_drug_map[pname] = drugs_for_panel

    update(4, 6, f"Scoring {len(all_drugs)} candidates (2 gates + 3 criteria)...")
    results = []
    for drug in all_drugs.values():
        results.append(score_drug(drug, genes, G, disease_pathway_vec))

    scored = [r for r in results if r["status"] == "scored"]
    g1_fail = [r for r in results if r["status"] == "gate1_fail"]
    g2_fail = [r for r in results if r["status"] == "gate2_fail"]
    scored.sort(key=lambda r: r["final_score"], reverse=True)

    update(5, 6, "Building probability tree...")

    # Build tree structure for visualization
    layer2_nodes = []
    for pname in panel_names:
        drugs_in_panel = panel_drug_map.get(pname, [])
        drug_names_in_panel = {d["name"].lower() for d in drugs_in_panel}

        scored_drugs = []
        for r in scored:
            if r["drug_name"].lower() in drug_names_in_panel:
                scored_drugs.append({
                    "drug_name": r["drug_name"],
                    "final_score": r["final_score"],
                    "scores": {
                        "bbb": r.get("bbb_score"),
                        "safety": r.get("safety_score"),
                        "pathway_alignment": r.get("c1_pathway"),
                        "network_proximity": r.get("c2_network"),
                        "phenotypic_overlap": r.get("c3_phenotypic"),
                    },
                    "evidence": {
                        "top_shared_pathways": r.get("shared_pathways", []),
                        "drug_targets": r.get("drug_targets", []),
                        "shared_genes": r.get("shared_genes", []),
                        "mean_network_distance": r.get("mean_network_distance"),
                        "bbb_experimental_evidence": (r.get("gate1_bbb") or {}).get("experimental_evidence"),
                    },
                    "caution": r.get("caution", False),
                    "caution_label": r.get("caution_label"),
                    "drug_phase": r.get("phase", 0),
                    "moa": r.get("moa", ""),
                })
        layer2_nodes.append({
            "disease": pname,
            "n_candidates": len(drugs_in_panel),
            "n_scored": len(scored_drugs),
            "overlap_score_to_root": 0.5,  # placeholder
            "drugs": scored_drugs,
        })

    tree = {
        "root": {
            "disease": matched_name,
            "disease_module_size": len(genes),
            "layer2": layer2_nodes,
        },
        "summary": {
            "root_disease": matched_name,
            "n_total_candidates": len(results),
            "n_gate1_pass": len(results) - len(g1_fail),
            "n_gate2_pass": len(results) - len(g1_fail) - len(g2_fail),
            "n_scored": len(scored),
            "n_layer2_diseases": len(panel_names),
            "global_top_drugs": [
                {"drug_name": r["drug_name"], "final_score": r["final_score"],
                 "caution": r.get("caution", False)}
                for r in scored[:10]
            ],
        },
        "all_results": results,
        "scored": scored,
        "scoring_method": "magnitude_aware_cosine",
    }

    update(6, 6, "Done!")
    return tree


# ──────────────────────────────────────────────────────────────────────────────
# Rendering
# ──────────────────────────────────────────────────────────────────────────────

def render_tree_viz(tree: dict, min_score: float):
    try:
        from pyvis.network import Network
        import streamlit.components.v1 as components

        net = Network(
            height="580px", width="100%", directed=True,
            bgcolor="#0e1117", font_color="white",
        )
        net.set_options(json.dumps({
            "physics": {"enabled": True, "stabilization": {"iterations": 100}},
            "interaction": {"hover": True, "tooltipDelay": 200},
        }))

        root_disease = tree.get("root", {}).get("disease", "Root Disease")
        net.add_node(
            "root", label=root_disease, color="#e74c3c", size=40,
            title=f"Root Disease: {root_disease}\n"
                  f"Disease module: {tree['root'].get('disease_module_size', '?')} genes",
            shape="diamond",
        )

        drug_nodes_added = set()

        for l2 in tree.get("root", {}).get("layer2", []):
            disease = l2.get("disease", "")
            n_scored = l2.get("n_scored", 0)
            did = f"d_{_slugify(disease)}"
            label = f"{disease[:28]}\n({n_scored} scored)"
            net.add_node(did, label=label, color="#3498db", size=22,
                         title=f"{disease}\nCandidates: {l2.get('n_candidates', 0)}\nScored: {n_scored}",
                         shape="ellipse")
            net.add_edge("root", did, width=2, color="#444")

            for drug in l2.get("drugs", []):
                fs = drug.get("final_score") or 0
                if fs < min_score:
                    continue
                drug_key = drug["drug_name"].lower()
                if drug_key in drug_nodes_added:
                    # Just add edge from this panel disease to existing drug node
                    drug_id = f"drug_{drug_key[:20]}"
                    net.add_edge(did, drug_id, width=max(1, fs * 5), color=score_color(fs))
                    continue
                drug_nodes_added.add(drug_key)
                drug_id = f"drug_{drug_key[:20]}"
                color = score_color(fs)
                sc = drug.get("scores", {})
                caution = " ⚠ CAUTION" if drug.get("caution") else ""
                title_d = (
                    f"{drug['drug_name']}{caution}\n"
                    f"Final Score: {fmt(fs)}\n"
                    f"BBB: {fmt(sc.get('bbb'))}  |  Safety: {fmt(sc.get('safety'))}\n"
                    f"Pathway: {fmt(sc.get('pathway_alignment'))}  |  Network: {fmt(sc.get('network_proximity'))}  |  Phenotypic: {fmt(sc.get('phenotypic_overlap'))}\n"
                    f"Phase: {drug.get('drug_phase', '?')}"
                )
                net.add_node(drug_id, label=f"{drug['drug_name'][:18]}\n{fmt(fs)}",
                             color=color, size=14 + fs * 12, title=title_d, shape="box")
                net.add_edge(did, drug_id, width=max(1, fs * 5), color=color)

        with tempfile.NamedTemporaryFile(delete=False, suffix=".html", mode="w") as f:
            net.save_graph(f.name)
            html = Path(f.name).read_text()
        os.unlink(f.name)
        components.html(html, height=590, scrolling=True)

    except ImportError:
        st.warning("Install pyvis for the interactive tree: `pip install pyvis`")
        st.info("See the Ranked Candidates tab for results in table form.")
    except Exception as e:
        st.error(f"Tree render error: {e}")


def render_evidence_ledger(drug: dict):
    st.markdown(f"### Evidence Ledger — {drug['drug_name']}")
    scores = drug.get("scores", {})
    ev = drug.get("evidence", {})

    cols = st.columns(6)
    labels = [("BBB", "bbb"), ("Safety", "safety"), ("Pathway", "pathway_alignment"),
              ("Network", "network_proximity"), ("Phenotypic", "phenotypic_overlap")]
    for col, (label, key) in zip(cols[:5], labels):
        v = scores.get(key)
        col.metric(label, fmt(v))
    cols[5].metric("**Final**", fmt(drug.get("final_score")))

    if drug.get("caution"):
        st.warning(drug.get("caution_label", "Caution: limited pediatric data"))

    col_a, col_b = st.columns(2)
    with col_a:
        st.markdown("**Shared Pathways**")
        for p in (ev.get("top_shared_pathways") or []):
            st.markdown(f"- {p}")
        if not (ev.get("top_shared_pathways") or []):
            st.caption("None identified")

        st.markdown("**Drug Targets**")
        targets = ev.get("drug_targets") or []
        st.caption(", ".join(targets[:8]) or "—")

    with col_b:
        st.markdown("**Shared Genes with Disease Module**")
        shared = ev.get("shared_genes") or []
        st.caption(", ".join(shared[:8]) or "—")

        dist = ev.get("mean_network_distance")
        st.markdown("**Network Distance**")
        st.caption(f"{dist:.2f} hops from disease module" if dist is not None else "—")

        if ev.get("bbb_experimental_evidence"):
            st.markdown("**CNS Evidence**")
            st.caption(ev["bbb_experimental_evidence"][:200])

    if drug.get("moa"):
        st.caption(f"Mechanism of action: {drug['moa']}")


def render_methodology():
    st.subheader("Methodology")
    st.markdown("""
**Evaluation Pipeline**

```
STAGE 1 — GATING (every candidate)
  Gate 1: BBB Permeability     score >= 0.8
  Gate 2: Pediatric Safety     score >= 0.6

STAGE 2 — SCORING (gate-passing candidates only)
  Criterion 1: Pathway Mechanism Alignment  (0-1)
  Criterion 2: Network Proximity            (0-1)
  Criterion 3: Phenotypic Overlap           (0-1)

  evidence  = [C1, C2, C3]
  reference = [1, 1, 1]          <- ideal drug
  Final Score = clamp((evidence . reference) / (reference . reference), 0, 1)
              = clamp(mean(C1, C2, C3), 0, 1)   <- magnitude-aware cosine
```

**Gate 1 — BBB Score (>= 0.8)**
CNS MPO (Wager et al. 2010, 6 physicochemical descriptors) x Pgp efflux adjustment.
When a drug has experimental CNS evidence with score >= 0.80, that score is used directly
as the authoritative signal (e.g. everolimus=0.90, sirolimus=0.85).
Otherwise: 0.35 x physicochemical + 0.65 x clinical evidence.

**Gate 2 — Pediatric Safety Score (>= 0.6)**
0.40 x FDA pediatric approval status + 0.30 x min approved age + 0.30 x serious AE rate.
Hard disqualifiers: pediatric black box warning or explicit pediatric contraindication.

**Criterion 1 — Pathway Alignment**
0.5 x Jaccard(Reactome pathway sets) + 0.5 x cosine(Enrichr GO-BP score vectors).
Drug targets expanded via STRING 1-hop neighbors (confidence >= 700) when < 5 direct targets.

**Criterion 2 — Network Proximity**
1/(1 + mean shortest-path distance) in STRING PPI subgraph (confidence >= 0.7, 2-hop neighborhood).
Uses direct MOA targets only (not expanded) to preserve mechanistic specificity.

**Criterion 3 — Phenotypic Overlap**
0.4 x HPO Jaccard (gene-based) + 0.3 x gene Jaccard + 0.3 x pathway cosine similarity.

**Data Sources (all free, no API keys)**
- OpenTargets Platform (disease-gene associations, drug-target MOA)
- STRING (protein-protein interactions)
- Enrichr (pathway enrichment: KEGG, Reactome, GO Biological Process)
- Reactome (curated pathway memberships)
- HPO (Human Phenotype Ontology, gene-phenotype annotations)
- PubChem (molecular properties for CNS MPO score)
- openFDA (pediatric safety labels)

**Comparison Panel**
Candidate drugs are sourced from 6 related neurological diseases:
TSC, NF1, SYNGAP1, Angelman, Rett, and Fragile X syndromes.

**Validated on:** everolimus/sirolimus top-2 for TSC (score=0.416);
selumetinib/trametinib top-2 for NF1 (score=0.472).
    """)


# ──────────────────────────────────────────────────────────────────────────────
# Main App
# ──────────────────────────────────────────────────────────────────────────────

def main():
    # Sidebar
    st.sidebar.title("🧬 HackRare 2026")
    st.sidebar.markdown("**Drug Repurposing Navigator**")
    st.sidebar.markdown("---")

    st.sidebar.markdown("#### Enter a rare disease")
    disease_input = st.sidebar.text_input(
        "Disease name",
        placeholder="e.g. Dravet Syndrome, CDKL5 Deficiency...",
        label_visibility="collapsed",
    )

    # Search button
    if st.sidebar.button("🔍 Search", type="primary"):
        if disease_input.strip():
            with st.spinner("Searching OpenTargets..."):
                hits = ot_search_disease(disease_input.strip())
            if hits:
                st.session_state["search_hits"] = hits
                st.session_state["selected_disease"] = None
            else:
                st.sidebar.error("No diseases found. Try a different name.")
                st.session_state["search_hits"] = []
        else:
            st.sidebar.warning("Please enter a disease name.")

    # Show search results
    hits = st.session_state.get("search_hits", [])
    if hits:
        st.sidebar.markdown("**Matches found:**")
        hit_names = [h["name"] for h in hits]
        chosen = st.sidebar.radio("Select:", hit_names, key="disease_radio")
        st.session_state["selected_disease"] = chosen

        if st.sidebar.button("✕ Clear", key="clear_btn"):
            st.session_state["search_hits"] = []
            st.session_state["selected_disease"] = None
            st.session_state.pop("pipeline_tree", None)
            st.rerun()

    st.sidebar.markdown("---")
    min_score_filter = st.sidebar.slider("Min score filter", 0.0, 1.0, 0.0, 0.05)

    selected = st.session_state.get("selected_disease")

    # Run pipeline button
    if selected:
        st.sidebar.markdown("---")
        if st.sidebar.button(f"▶ Run Pipeline for\n{selected}", type="primary"):
            progress = st.progress(0, text="Starting pipeline...")
            def update_progress(frac, msg):
                progress.progress(min(frac, 1.0), text=msg)
            tree = run_pipeline_web(selected, progress_callback=update_progress)
            st.session_state["pipeline_tree"] = tree
            st.session_state["pipeline_disease"] = selected
            progress.empty()
            st.rerun()

    # ── Main content ──────────────────────────────────────────────────────────

    tree = st.session_state.get("pipeline_tree")
    pipeline_disease = st.session_state.get("pipeline_disease", "")

    if tree is None:
        # Landing page
        st.title("🧬 Rare Disease Drug Repurposing Navigator")
        st.markdown("""
        **Find repurposable drugs for any rare disease using mechanistic evidence.**

        This tool scores existing drugs for potential repurposing to rare diseases
        by evaluating three dimensions of mechanistic alignment:

        1. **Pathway Alignment** — Do the drug's targets overlap with the disease's
           biological pathways?
        2. **Network Proximity** — How close are the drug's targets to the disease
           genes in the protein interaction network?
        3. **Phenotypic Overlap** — Do the drug's target genes share phenotypic
           annotations with the disease?

        Drugs must first pass two safety gates: **blood-brain barrier permeability**
        and **pediatric safety profile**.
        """)

        st.markdown("---")

        st.markdown("#### Quick Start")
        st.markdown("""
        1. Type a rare disease name in the sidebar (e.g. *Dravet Syndrome*)
        2. Click **Search** to find matching diseases in OpenTargets
        3. Select the correct match
        4. Click **Run Pipeline** to score candidate drugs
        5. Explore results in the probability tree and ranked table
        """)

        cols = st.columns(3)
        with cols[0]:
            st.markdown("**Try these:**")
            st.markdown("- Dravet Syndrome")
            st.markdown("- CDKL5 Deficiency")
            st.markdown("- Phelan-McDermid Syndrome")
        with cols[1]:
            st.markdown("**Also works for:**")
            st.markdown("- Tuberous Sclerosis")
            st.markdown("- Rett Syndrome")
            st.markdown("- Fragile X Syndrome")
        with cols[2]:
            st.markdown("**Or any disease in:**")
            st.markdown("- [OpenTargets Platform](https://platform.opentargets.org)")
            st.markdown("- 23,000+ diseases indexed")

        st.markdown("---")
        render_methodology()
        return

    # ── Results view ──────────────────────────────────────────────────────────
    st.title(f"🧬 Results — {pipeline_disease}")

    summary = tree.get("summary", {})
    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric("Evaluated", summary.get("n_total_candidates", 0))
    c2.metric("Passed BBB", summary.get("n_gate1_pass", 0))
    c3.metric("Passed Safety", summary.get("n_gate2_pass", 0))
    c4.metric("Scored", summary.get("n_scored", 0))
    c5.metric("Panel Diseases", summary.get("n_layer2_diseases", 0))

    st.markdown("---")

    t1, t2, t3 = st.tabs([
        "🌳 Probability Tree",
        "💊 Ranked Candidates",
        "📖 Methodology",
    ])

    with t1:
        st.subheader("Drug Repurposing Probability Tree")
        st.caption(
            "🟢 score >= 0.7  ·  🟡 0.5–0.7  ·  🟠 0.3–0.5  ·  🔴 < 0.3  ·  "
            "Node size and edge thickness scale with final score"
        )
        render_tree_viz(tree, min_score_filter)

    with t2:
        st.subheader("All Scored Candidates")

        # Collect all scored drugs across panel diseases, deduplicate
        all_drugs = []
        for node in tree.get("root", {}).get("layer2", []):
            for drug in node.get("drugs", []):
                fs = drug.get("final_score") or 0
                if fs < min_score_filter:
                    continue
                drug["_via"] = node.get("disease", "")
                all_drugs.append(drug)

        seen = set()
        unique_drugs = []
        for d in sorted(all_drugs, key=lambda x: x.get("final_score") or 0, reverse=True):
            name = d["drug_name"].lower()
            if name not in seen:
                seen.add(name)
                unique_drugs.append(d)

        if unique_drugs:
            rows = []
            for d in unique_drugs:
                sc = d.get("scores", {})
                rows.append({
                    "Rank": len(rows) + 1,
                    "Drug": d["drug_name"],
                    "Score": d.get("final_score"),
                    "BBB": sc.get("bbb"),
                    "Safety": sc.get("safety"),
                    "Pathway (C1)": sc.get("pathway_alignment"),
                    "Network (C2)": sc.get("network_proximity"),
                    "Phenotypic (C3)": sc.get("phenotypic_overlap"),
                    "Phase": d.get("drug_phase", ""),
                    "Caution": "⚠" if d.get("caution") else "",
                    "Via Disease": d.get("_via", ""),
                })

            df = pd.DataFrame(rows)
            float_cols = [c for c in ["Score", "BBB", "Safety", "Pathway (C1)", "Network (C2)", "Phenotypic (C3)"]
                          if c in df.columns]
            st.dataframe(
                df.style.format({c: "{:.3f}" for c in float_cols}),
                use_container_width=True, height=420,
            )

            st.markdown("---")
            selected_name = st.selectbox("View detailed evidence for:", [d["drug_name"] for d in unique_drugs])
            selected_drug = next((d for d in unique_drugs if d["drug_name"] == selected_name), None)
            if selected_drug:
                render_evidence_ledger(selected_drug)

            st.markdown("---")
            csv = df.to_csv(index=False)
            st.download_button(
                "⬇ Download Results (CSV)", csv,
                f"repurposing_{_slugify(pipeline_disease)}.csv", "text/csv",
            )

            # Also offer JSON download
            json_str = json.dumps(tree, indent=2, default=str)
            st.download_button(
                "⬇ Download Full Tree (JSON)", json_str,
                f"tree_{_slugify(pipeline_disease)}.json", "application/json",
            )
        else:
            st.info("No candidates above the current filter threshold.")

    with t3:
        render_methodology()


if __name__ == "__main__" or True:
    main()
