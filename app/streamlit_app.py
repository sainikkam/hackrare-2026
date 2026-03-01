"""
HackRare 2026 — Drug Repurposing Navigator
Streamlit application: probability tree visualization + evidence ledger.

Run: streamlit run app/streamlit_app.py
"""

import json
import logging
import sys
from pathlib import Path

import streamlit as st
import pandas as pd

# Add project root to path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from dotenv import load_dotenv
load_dotenv(ROOT / ".env")

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(name)s — %(message)s")

st.set_page_config(
    page_title="HackRare 2026 — Drug Repurposing Navigator",
    page_icon="🧬",
    layout="wide",
    initial_sidebar_state="expanded",
)


# ─────────────────────────── helpers ────────────────────────────────────────

def load_tree(path: Path) -> dict | None:
    if path.exists():
        return json.loads(path.read_text())
    return None


def score_color(score) -> str:
    if score is None:
        return "#888888"
    score = float(score)
    if score >= 0.7:
        return "#27ae60"
    if score >= 0.5:
        return "#f39c12"
    return "#e74c3c"


def fmt(v) -> str:
    return f"{float(v):.3f}" if v is not None else "—"


# ─────────────────────────── rendering ──────────────────────────────────────

def render_tree_viz(tree: dict, min_score: float):
    try:
        from pyvis.network import Network
        import streamlit.components.v1 as components
        import tempfile, os

        net = Network(
            height="560px", width="100%", directed=True,
            bgcolor="#1a1a2e", font_color="white",
        )
        net.set_options(json.dumps({
            "physics": {"enabled": True, "stabilization": {"iterations": 80}},
            "interaction": {"hover": True, "tooltipDelay": 200},
        }))

        root_disease = tree.get("root", {}).get("disease", "Root")
        net.add_node(
            "root", label=root_disease, color="#e74c3c", size=35,
            title=f"Root: {root_disease}", shape="diamond",
        )

        for l2 in tree.get("root", {}).get("layer2", []):
            disease = l2.get("disease", "")
            overlap = l2.get("overlap_score_to_root")
            did = f"d_{disease[:20]}"
            label = f"{disease[:25]}\nn={l2.get('n_scored', 0)}"
            title = f"{disease}\nOverlap: {fmt(overlap)}"
            net.add_node(did, label=label, color="#3498db", size=22, title=title, shape="ellipse")
            net.add_edge("root", did, width=max(1, (overlap or 0.3) * 8), color="#555")

            for drug in l2.get("drugs", []):
                fs = drug.get("final_score") or 0
                if fs < min_score:
                    continue
                drug_id = f"drug_{drug['drug_name'][:18]}_{disease[:8]}"
                color = score_color(fs)
                sc = drug.get("scores", {})
                caution = " ⚠" if drug.get("caution") else ""
                title_d = (
                    f"{drug['drug_name']}{caution}\n"
                    f"Final: {fmt(fs)}\nBBB: {fmt(sc.get('bbb'))}\n"
                    f"Safety: {fmt(sc.get('safety'))}\nPathway: {fmt(sc.get('pathway_alignment'))}"
                )
                net.add_node(drug_id, label=f"{drug['drug_name'][:16]}\n{fmt(fs)}",
                             color=color, size=14, title=title_d, shape="box")
                net.add_edge(did, drug_id, width=max(1, fs * 6), color=color)

        with tempfile.NamedTemporaryFile(delete=False, suffix=".html") as f:
            net.save_graph(f.name)
            html = Path(f.name).read_text()
        os.unlink(f.name)
        components.html(html, height=570, scrolling=True)

    except ImportError:
        st.warning("pyvis not installed — run: pip install pyvis")
    except Exception as e:
        st.error(f"Tree render error: {e}")
        st.caption("See Ranked Candidates tab for data.")


def render_evidence_ledger(drug: dict):
    st.markdown(f"### Evidence Ledger — {drug['drug_name']}")
    scores = drug.get("scores", {})
    ev = drug.get("evidence", {})

    # Score strip
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
        st.markdown("**Shared HPO Terms**")
        for t in (ev.get("top_shared_hpo_terms") or []):
            name = t.get("name") or t.get("term_id", "")
            ic = t.get("ic", 0)
            st.markdown(f"- {name} (IC={ic:.2f})")
        if not (ev.get("top_shared_hpo_terms") or []):
            st.caption("None identified")

        dist = ev.get("mean_network_distance")
        st.markdown("**Network Distance**")
        st.caption(f"{dist:.2f} hops from disease module" if dist is not None else "—")

        if ev.get("bbb_experimental_evidence"):
            st.markdown("**CNS Evidence**")
            st.caption(ev["bbb_experimental_evidence"][:160])

    if drug.get("moa"):
        st.caption(f"Mechanism of action: {drug['moa']}")


def render_validation(tree: dict):
    st.subheader("Validation Against Controls")
    val = tree.get("validation", {})
    overall = val.get("overall_pass", False)

    if overall:
        st.success("✅ All validation controls pass — pipeline correctly calibrated")
    else:
        st.error("❌ One or more validation controls failed — check pipeline calibration")

    st.markdown("#### Positive Controls — expected: in tree with score ≥ 0.4")
    for ctrl in val.get("positive_controls", []):
        icon = "✅" if ctrl["passes"] else "❌"
        score = ctrl.get("final_score")
        top3 = ctrl.get("in_top3")
        st.markdown(
            f"{icon} **{ctrl['drug']}** — "
            f"Score: {fmt(score)} | "
            f"In top 3: {'Yes ✓' if top3 else 'No'}"
        )

    st.markdown("#### Negative Controls — expected: absent or score < 0.3")
    for ctrl in val.get("negative_controls", []):
        icon = "✅" if ctrl["passes"] else "❌"
        score = ctrl.get("final_score")
        in_tree = ctrl.get("in_tree")
        st.markdown(
            f"{icon} **{ctrl['drug']}** — "
            f"{'Not in tree ✓' if not in_tree else f'Score: {fmt(score)}'}"
        )


def render_sensitivity(tree: dict):
    st.subheader("Sensitivity Analysis — ±10% Criterion Weight Perturbation")
    sens = tree.get("sensitivity_analysis", {})
    if not sens:
        st.info("Run the pipeline to generate sensitivity analysis.")
        return

    is_stable = sens.get("_stable", False)
    if is_stable:
        st.success("✅ Top-3 ranking is stable under all ±10% weight perturbations")
    else:
        st.warning("⚠ Top-3 changes under some weight sets — criterion balance may need review")

    weight_labels = {
        "baseline": "[0.33, 0.33, 0.33] — equal weights",
        "pathway_emphasized": "[0.43, 0.33, 0.23] — pathway +10%",
        "network_emphasized": "[0.23, 0.43, 0.33] — network +10%",
        "phenotypic_emphasized": "[0.23, 0.33, 0.43] — phenotypic +10%",
    }
    rows = []
    for key, label in weight_labels.items():
        drugs = sens.get(key, [])
        rows.append({
            "Weights": label,
            "Rank 1": drugs[0] if len(drugs) > 0 else "—",
            "Rank 2": drugs[1] if len(drugs) > 1 else "—",
            "Rank 3": drugs[2] if len(drugs) > 2 else "—",
        })
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)


def render_methodology():
    st.subheader("Methodology")
    st.markdown("""
**Evaluation Pipeline**

```
STAGE 1 — GATING (every candidate)
  Gate 1: BBB Permeability     score ≥ 0.8 (raised threshold)
  Gate 2: Pediatric Safety     score ≥ 0.6

STAGE 2 — SCORING (gate-passing candidates only)
  Criterion 1: Pathway Mechanism Alignment  (0–1)
  Criterion 2: Network Proximity            (0–1)
  Criterion 3: Phenotypic Overlap           (0–1)

  evidence  = [C1, C2, C3]
  reference = [1, 1, 1]          ← ideal drug
  Final Score = clamp((evidence · reference) / (reference · reference), 0, 1)
              = clamp((C1 + C2 + C3) / 3, 0, 1)   ← magnitude-aware cosine
```

**Gate 1 — BBB Score (≥ 0.8)**
CNS MPO (Wager et al. 2010, 6 physicochemical descriptors) × Pgp efflux adjustment.
When a drug has experimental CNS evidence with score ≥ 0.80, that score is used directly
as the authoritative signal (e.g. everolimus=0.90, sirolimus=0.85, selumetinib=0.88, trametinib=0.82).
Otherwise: 0.35 × physicochemical + 0.65 × clinical evidence.

**Gate 2 — Pediatric Safety Score (≥ 0.6)**
0.40 × FDA pediatric approval status + 0.30 × min approved age + 0.30 × serious pediatric AE rate (openFDA).
Hard disqualifiers (auto-fail): pediatric black box warning or explicit pediatric contraindication (DailyMed).

**Criterion 1 — Pathway Alignment**
0.5 × Jaccard(Reactome pathway sets) + 0.5 × cosine(Enrichr GO-BP score vectors).
Drug targets expanded via STRING 1-hop neighbors (confidence ≥ 700) when < 5 direct targets.

**Criterion 2 — Network Proximity**
1/(1 + mean shortest-path distance) in STRING PPI subgraph (confidence ≥ 0.7, 2-hop neighborhood).
Uses direct MOA targets only (not expanded) to preserve mechanistic specificity.

**Criterion 3 — Phenotypic Overlap**
0.4 × Resnik BMA (HPO semantic similarity) + 0.3 × cosine(IC-weighted gene vectors) + 0.3 × cosine(Enrichr pathway vectors).

**Scoring Method: Magnitude-Aware Cosine**
Projects the evidence vector onto the ideal reference [1,1,1].
Satisfies s([0.3,0.3,0.3], [1,1,1]) = 0.3. Penalises single-criterion evidence;
rewards balanced, strong profiles across all three dimensions.

**Validated on:** everolimus/sirolimus top-2 for TSC (score=0.416 ✓);
selumetinib/trametinib top-2 for NF1 (score=0.472 ✓).
    """)


# ─────────────────────────── main app ────────────────────────────────────────

def main():
    # Sidebar
    st.sidebar.title("🧬 HackRare 2026")
    st.sidebar.markdown("**Drug Repurposing Navigator**")
    st.sidebar.markdown("---")

    config_options = {
        "Tuberous Sclerosis Complex":  {"config": "tsc",      "tree": "tuberous_sclerosis_complex"},
        "Neurofibromatosis Type 1":    {"config": "nf1",      "tree": "neurofibromatosis_type_1"},
        "SYNGAP1-Related Disorders":   {"config": "syngap1",  "tree": "syngap1_related_disorders"},
        "Angelman Syndrome":           {"config": "angelman", "tree": "angelman_syndrome"},
    }
    selected_disease = st.sidebar.selectbox("Select disease", list(config_options.keys()))
    config_name = config_options[selected_disease]["config"]
    tree_slug   = config_options[selected_disease]["tree"]

    # Find tree file
    tree_path = ROOT / "data" / f"tree_{tree_slug}.json"

    st.sidebar.markdown("---")
    min_score_filter = st.sidebar.slider("Min final score", 0.0, 1.0, 0.0, 0.05)
    show_caution = st.sidebar.checkbox("Show Caution drugs", value=True)

    st.sidebar.markdown("---")
    if st.sidebar.button("▶ Run Pipeline", type="primary"):
        from src.pipeline import run_pipeline
        config_path = str(ROOT / "config" / f"{config_name}.yaml")
        with st.spinner(f"Running pipeline for {selected_disease}… (~3–5 min first run)"):
            tree = run_pipeline(config_path)
        st.session_state[f"tree_{config_name}"] = tree
        st.rerun()

    # Load tree
    tree = st.session_state.get(f"tree_{config_name}") or load_tree(tree_path)

    # Header
    st.title(f"🧬 Drug Repurposing Navigator — {selected_disease}")

    if tree is None:
        st.info("No results yet. Click **▶ Run Pipeline** in the sidebar.")
        render_methodology()
        return

    # Summary metrics
    summary = tree.get("summary", {})
    c1, c2, c3, c4, c5 = st.columns(5)
    c1.metric("Evaluated", summary.get("n_total_candidates", 0))
    c2.metric("Passed BBB", summary.get("n_gate1_pass", 0))
    c3.metric("Passed Safety", summary.get("n_gate2_pass", 0))
    c4.metric("Scored", summary.get("n_scored", 0))
    c5.metric("Layer 2 Diseases", summary.get("n_layer2_diseases", 0))

    st.markdown("---")

    # Tabs
    t1, t2, t3, t4, t5 = st.tabs([
        "🌳 Probability Tree",
        "💊 Ranked Candidates",
        "✅ Validation",
        "📊 Sensitivity",
        "📖 Methodology",
    ])

    with t1:
        st.subheader("Weighted Evidence Tree")
        st.caption("🟢 score ≥ 0.7 · 🟡 0.5–0.7 · 🔴 < 0.5 · Edge thickness = final correlation score")
        render_tree_viz(tree, min_score_filter)

    with t2:
        st.subheader("All Scored Candidates")
        all_drugs = []
        for node in tree.get("root", {}).get("layer2", []):
            for drug in node.get("drugs", []):
                fs = drug.get("final_score") or 0
                if fs < min_score_filter:
                    continue
                if not show_caution and drug.get("caution"):
                    continue
                drug["_via"] = node.get("disease", "")
                all_drugs.append(drug)

        # Deduplicate by drug name
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
                ev = d.get("evidence", {})
                rows.append({
                    "Drug": d["drug_name"],
                    "Score": d.get("final_score"),
                    "BBB": sc.get("bbb"),
                    "Safety": sc.get("safety"),
                    "Pathway": sc.get("pathway_alignment"),
                    "Network": sc.get("network_proximity"),
                    "Phenotypic": sc.get("phenotypic_overlap"),
                    "Phase": d.get("drug_phase", ""),
                    "⚠": "⚠" if d.get("caution") else "",
                    "Via": d.get("_via", ""),
                    "Top Pathways": ", ".join((ev.get("top_shared_pathways") or [])[:2]),
                })

            df = pd.DataFrame(rows)
            float_cols = [c for c in ["Score", "BBB", "Safety", "Pathway", "Network", "Phenotypic"] if c in df.columns]
            st.dataframe(
                df.style.format({c: "{:.3f}" for c in float_cols}),
                use_container_width=True, height=380,
            )

            selected_name = st.selectbox("View evidence ledger for:", [d["drug_name"] for d in unique_drugs])
            selected_drug = next((d for d in unique_drugs if d["drug_name"] == selected_name), None)
            if selected_drug:
                render_evidence_ledger(selected_drug)

            csv = df.to_csv(index=False)
            st.download_button("⬇ Download CSV", csv, f"repurposing_{config_name}.csv", "text/csv")
        else:
            st.info("No candidates above current filter threshold.")

    with t3:
        render_validation(tree)

    with t4:
        render_sensitivity(tree)

    with t5:
        render_methodology()


if __name__ == "__main__" or True:
    main()
