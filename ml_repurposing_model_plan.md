# ML Model for Cross-Disease Drug Repurposing Correlation
## Design & Validation Plan

**Project:** HackRare 2026 — Drug Repurposing Navigator Extension
**Author:** HackRare Engineering
**Date:** 2026-03-01
**Status:** PROPOSED — Awaiting Implementation Approval

---

## 1. Problem Statement

The existing HackRare pipeline uses **fixed equal-weight heuristic scoring** (1/3 pathway + 1/3 network proximity + 1/3 phenotypic overlap) to rank drug candidates against a disease module. This is correct for mechanistic scoring, but has two limitations:

1. **Equal weights are unjustified.** For some diseases, pathway alignment is the dominant signal; for others, phenotypic overlap or network topology matters more. Historical repurposing successes carry the ground truth for which features predict real outcomes.
2. **No cross-compatibility transfer.** If drug D was successfully repurposed for disease A, and A shares mechanistic ground with disease B, the existing pipeline does not amplify D's score for B — it scores D against B from scratch without using the A→D signal.

The ML model proposed here addresses both gaps: it **learns feature weights from historical repurposings** and produces a **cross-compatibility correlation score** that explicitly accounts for mechanistic relatedness between diseases.

---

## 2. Formal Problem Formulation

### 2.1 Primary Task: Drug–Disease Repurposing Link Prediction

Given:
- Drug `d` characterized by molecular profile `x_d`
- Disease `b` characterized by genomic/phenotypic profile `x_b`

Predict:
- `P(repurposed | d, b)` — the probability that drug `d` is therapeutically relevant for disease `b`

This is a **binary link prediction** task on a bipartite drug–disease graph, where positive edges are known clinical repurposings.

### 2.2 Secondary Task: Cross-Compatibility Correlation Score

Given:
- Drug `d` with known repurposed indication `a` (approved or clinical-stage)
- Target disease `b` (rare, related to `a` by biology)

Compute:
- `CrossCompat(d, a→b)` — a correlation score representing how the mechanistic basis of `d`'s efficacy in `a` transfers to `b`

This captures the intuition: *everolimus works in TSC via mTOR inhibition; SYNGAP1 involves downstream mTOR dysregulation; therefore the TSC→SYNGAP1 mechanistic bridge amplifies everolimus's score for SYNGAP1.*

---

## 3. Training Data

### 3.1 Positive Examples — Known Repurposings

| Source | Description | Expected Pairs |
|--------|-------------|---------------|
| **DrugCentral** (drugcentral.org) | Curated approved drug–indication pairs, including off-label repurposings | ~2,500 unique pairs |
| **OpenTargets** (evidence type: `known_drug`) | Phase ≥ 3 drug–disease clinical associations with evidence strings | ~4,000 pairs |
| **TTD** (Therapeutic Target Database) | Approved and clinical trial drug–target–disease triples | ~1,800 pairs |
| **Orphanet** | Rare disease–drug associations (Orphanet Drug entities) | ~600 rare disease pairs |

**Total positive set target:** ~6,000 deduplicated drug–disease positive pairs after mapping to ChEMBL IDs and EFO/OMIM disease identifiers.

### 3.2 Negative Examples

Three strategies, used in proportion 2:1:1:

1. **Random negatives**: drug–disease pairs not in any positive set and with no pathway co-annotation. Ratio 5:1 negative:positive.
2. **Hard negatives (same disease class)**: drug approved for a mechanistically unrelated disease in the same phenotypic category (e.g., imatinib vs. SYNGAP1 — both clinically approved, no molecular overlap). These are the most informative for the model.
3. **Hard negatives (same target class)**: drugs with identical primary targets but no known efficacy in the disease (e.g., a different mTOR inhibitor that failed in Phase 2 for TSC).

### 3.3 Feature Computation for Training Pairs

Each training pair `(d, b)` is featurized as described in Section 4. Features are precomputed and cached in a training matrix `X ∈ ℝ^(N × F)` with label vector `y ∈ {0,1}^N`.

---

## 4. Feature Engineering

### 4.1 Drug Molecular Features

| Feature Group | Dimensionality | Source | Description |
|---|---|---|---|
| **Morgan Fingerprint** | 1024-bit | RDKit / PubChem SMILES | Circular fingerprint radius=2; encodes chemical substructures |
| **MACCS Keys** | 167-bit | RDKit | Structural keys; complements Morgan FP |
| **Target Gene Vector** | ~20,000 (sparse) | ChEMBL / OpenTargets MOA | Binary indicator: which human genes the drug targets (pChEMBL ≥ 7) |
| **Reactome Pathway Vector** | 2,500 (sparse) | Reactome | Pathways containing ≥1 direct drug target |
| **GO-BP Enrichment Vector** | 1,000 (dense) | Enrichr | −log10(p-value) for top GO-BP terms from drug target set; captures functional bias |
| **Target Network Centrality** | 4 scalars | STRING | Mean degree, betweenness, eigenvector centrality, clustering coefficient of drug targets |
| **MOA Category** | 12-dim one-hot | ChEMBL MOA | Inhibitor / agonist / antagonist / activator / etc. |

**Drug feature vector total:** ~24,000 dimensions (mostly sparse; dimensionality reduced via SVD to 256-d dense representation for neural models).

### 4.2 Disease Genomic–Phenotypic Features

| Feature Group | Dimensionality | Source | Description |
|---|---|---|---|
| **Disease Gene Vector** | ~20,000 (sparse) | OpenTargets (score ≥ 0.3) + DisGeNET | Binary indicator over human genes in disease module |
| **Reactome Pathway Vector** | 2,500 (sparse) | Reactome | Pathways annotated to disease gene module |
| **GO-BP Enrichment Vector** | 1,000 (dense) | Enrichr | −log10(p) over GO-BP terms for disease gene set |
| **HPO Phenotype Vector** | 15,000 (sparse) | HPO phenotype_to_genes.txt | Binary indicator over HPO terms annotated to disease genes |
| **OMIM Disease Category** | 8-dim one-hot | OMIM | CNS / metabolic / oncological / cardiac / etc. |
| **Tissue Expression Profile** | 54-dim | GTEx | Mean TPM across 54 tissues for disease module genes |

**Disease feature vector total:** ~38,000 dimensions (reduced to 256-d dense).

### 4.3 Drug–Disease Interaction Features

These are computed pairwise and directly represent cross-compatibility signals:

| Feature | Type | Description |
|---|---|---|
| `pathway_jaccard` | scalar [0,1] | Jaccard of Reactome pathway sets (existing C1-A) |
| `enrichr_cosine` | scalar [0,1] | Cosine of GO-BP enrichment vectors (existing C1-B) |
| `network_proximity` | scalar [0,1] | 1/(1 + mean_shortest_path) in STRING (existing C2) |
| `hpo_jaccard` | scalar [0,1] | Jaccard of HPO term sets (existing C3) |
| `target_gene_overlap` | scalar [0,1] | Jaccard of drug targets ∩ disease module genes |
| `direct_target_in_module` | scalar [0,1] | Fraction of drug targets that are direct disease module genes |
| `string_neighbor_overlap` | scalar [0,1] | Jaccard of 1-hop STRING neighborhoods of targets and module |
| `tissue_expression_overlap` | scalar [0,1] | Cosine similarity of GTEx tissue profiles |
| `moa_pathway_in_disease` | binary | Is the drug's primary MOA pathway a known disease pathway? |
| `clinical_phase` | scalar | Max clinical phase of drug in any indication (proxy for safety/developability) |

**Interaction feature vector:** 10 scalars, concatenated with drug and disease dense representations.

---

## 5. Model Architecture

### 5.1 Overview: Two-Stage Pipeline

```
┌─────────────────────────────────────────────────────────────────────┐
│                         STAGE 1: Embeddings                         │
│                                                                     │
│  Drug Encoder          Disease Encoder        Cross-features        │
│  [Morgan FP +     →    [Gene vector +    →    [Pathway Jaccard,     │
│   Target genes +        HPO vector +           Network proximity,   │
│   GO-BP vector]         GO-BP vector]          HPO Jaccard, ...]    │
│       ↓                     ↓                        ↓              │
│   emb_d ∈ ℝ²⁵⁶         emb_b ∈ ℝ²⁵⁶          feats ∈ ℝ¹⁰          │
└─────────────────────────────────────────────────────────────────────┘
                              ↓
┌─────────────────────────────────────────────────────────────────────┐
│                      STAGE 2: Scoring                               │
│                                                                     │
│  concat([emb_d, emb_b, cosine(emb_d,emb_b), feats]) → XGBoost      │
│                                                                     │
│  Output: P(repurposed | d, b) ∈ [0,1]                               │
└─────────────────────────────────────────────────────────────────────┘
```

### 5.2 Stage 1A: Drug Encoder (MLP)

```python
DrugEncoder(
    input_dim  = 1024 + 167 + 1000 + 4 + 12,   # Morgan + MACCS + GO-BP + centrality + MOA
    hidden_dims = [512, 256],
    output_dim  = 256,
    activation  = ReLU,
    dropout     = 0.3,
    batch_norm  = True,
)
```

Pre-training objective: **drug–drug similarity regression** — given two drugs with known Tanimoto fingerprint similarity T(d1, d2), train encoder to produce embeddings with cosine similarity ≈ T. This grounds the embedding in chemical space before fine-tuning on repurposing labels.

### 5.3 Stage 1B: Disease Encoder (MLP)

```python
DiseaseEncoder(
    input_dim  = 1000 + 256 + 8,  # GO-BP + HPO (SVD-256) + category
    hidden_dims = [512, 256],
    output_dim  = 256,
    activation  = ReLU,
    dropout     = 0.3,
    batch_norm  = True,
)
```

Pre-training objective: **disease–disease phenotypic similarity regression** — using HPO semantic similarity (Lin IC) as supervision signal. Grounds the embedding in phenotypic space.

### 5.4 Stage 1C: Heterogeneous Graph Neural Network (Upgrade Path)

For the full-scale version, replace the MLP encoders with a **Relational Graph Convolutional Network (RGCN)**:

**Graph Schema:**
```
Nodes: Drug | Gene | Disease | Pathway | HPO_Term
Edges:
  drug  --[targets]--> gene              (ChEMBL pChEMBL ≥ 7)
  gene  --[in_pathway]--> pathway        (Reactome)
  gene  --[PPI]--> gene                  (STRING ≥ 0.7)
  disease --[associated_with]--> gene    (OpenTargets ≥ 0.3)
  disease --[has_phenotype]--> hpo_term  (HPO annotations)
  drug  --[repurposed_for]--> disease    (training labels — MASKED during eval)
```

Node initial features: one-hot ID + known attributes (SMILES fingerprint for drugs, sequence embeddings for genes, GO annotations for pathways).

The RGCN learns to propagate information across relation types. The final drug and disease node embeddings capture multi-hop biological context.

**Why RGCN > standalone MLPs:** Message passing allows a drug targeting FKBP1A to "learn" about mTOR pathway nodes two hops away — capturing the STRING expansion logic the current pipeline does manually via heuristic rules.

### 5.5 Stage 2: XGBoost Classifier

**Input features to XGBoost:**
```
[cosine(emb_d, emb_b),          # learned embedding similarity
 emb_d · emb_b,                 # dot product (raw)
 ||emb_d - emb_b||_2,           # L2 distance
 pathway_jaccard,                # from existing C1-A
 enrichr_cosine,                 # from existing C1-B
 network_proximity,              # from existing C2
 hpo_jaccard,                    # from existing C3
 target_gene_overlap,
 direct_target_in_module,
 string_neighbor_overlap,
 tissue_expression_overlap,
 moa_pathway_in_disease,
 clinical_phase]                 # 13 features total
```

**XGBoost hyperparameters (initial):**
```
n_estimators: 500
max_depth: 5
learning_rate: 0.05
subsample: 0.8
colsample_bytree: 0.8
scale_pos_weight: 10   # class imbalance correction
eval_metric: aucpr     # precision-recall AUC (better for imbalanced)
early_stopping_rounds: 50
```

**Why XGBoost as Stage 2:** Handles mixed feature types (continuous + binary), resistant to irrelevant features, outputs calibrated probabilities with Platt scaling, and supports SHAP for per-prediction explanations — critical for scientific credibility.

---

## 6. Cross-Compatibility Correlation Score

### 6.1 Definition

For drug `d` with a known approved/Phase 3 indication `a` (the "source disease"), and a rare target disease `b`:

```
CrossCompat(d, a→b) =
    sigmoid(
        w₁ · P(repurposed | d, b)              # Stage 2 model score
      + w₂ · cosine(embed(b), embed(a))        # Disease similarity in learned space
      + w₃ · pathway_jaccard(a, b)             # Shared biological pathways
      + w₄ · hpo_jaccard(a, b)                 # Shared clinical phenotypes
    )
```

Where weights `w₁...w₄` are learned by logistic regression on a held-out set of known disease-bridging repurposings (e.g., mTOR inhibitors working across TSC, SYNGAP1, NF1, and Rett syndrome).

### 6.2 Multi-Indication Drugs

When drug `d` has multiple known indications `{a₁, a₂, ..., aₖ}`:

```
CrossCompat(d, B) = max over aᵢ [ CrossCompat(d, aᵢ → b) · source_confidence(aᵢ) ]
```

Where `source_confidence(aᵢ)` = clinical phase of drug in indication `aᵢ` / 4 (normalized).

**Rationale:** The best mechanistic bridge wins. If everolimus is approved for TSC (phase=4) and renal cell carcinoma (phase=4), it gets two chances to transfer — the higher-scoring source disease bridge dominates.

### 6.3 Disease Bridge Map

The model produces a **Disease Bridge Map**: for each pair of diseases (a, b), compute:

```
Bridge(a, b) = cosine(embed(a), embed(b))
```

Diseases with `Bridge(a, b) > 0.7` are "mechanistically proximal" and drugs approved for either are scored with a cross-compatibility amplification factor of 1.2×.

This map is precomputed over the Layer 2 disease set and stored in the pipeline config. Example expected results:

| Disease Pair | Expected Bridge Score | Rationale |
|---|---|---|
| TSC ↔ SYNGAP1 | > 0.75 | Shared mTOR, shared synaptic phenotype |
| NF1 ↔ SYNGAP1 | > 0.80 | Both RAS-GAP proteins; RAS/MAPK is shared exact mechanism |
| Costello ↔ NF1 | > 0.72 | Both RAS hyperactivation, Costello via HRAS GOF |
| Rett ↔ Angelman | > 0.65 | Shared synaptic plasticity, overlapping HPO terms |
| TSC ↔ CML | < 0.25 | mTOR vs. BCR-ABL; different system, different phenotype |
| SYNGAP1 ↔ Gaucher | < 0.15 | No shared pathway; lysosomal vs. synaptic signaling |

---

## 7. Integration with Existing HackRare Pipeline

The ML model integrates as a **fourth scoring criterion** (`criterion4_ml`) alongside the existing three:

```
existing pipeline:
  final_score = (c1 + c2 + c3) / 3

upgraded pipeline:
  ml_score     = CrossCompat(d, known_indications → root_disease)
  final_score  = 0.25*c1 + 0.25*c2 + 0.25*c3 + 0.25*ml_score

  (weights replaced by learned XGBoost output in full version)
```

### 7.1 New Module Structure

```
src/
  ml/
    __init__.py
    encoders.py          # DrugEncoder, DiseaseEncoder (MLP)
    graph_builder.py     # Build heterogeneous drug-gene-disease graph
    rgcn_model.py        # RGCN node embedding model
    xgboost_scorer.py    # Stage 2 XGBoost classifier
    cross_compat.py      # CrossCompat() and Disease Bridge Map
    train.py             # Training loop (offline, produces model.pkl)
    inference.py         # score_drug_ml() — called from scorer.py
  training_data/
    drugcentral_pairs.csv
    opentargets_positive.csv
    negative_pairs.csv
    feature_matrix.pkl   # Cached training features
```

### 7.2 New Function Signature

```python
# src/ml/inference.py
def score_drug_ml(
    drug_name: str,
    drug_chembl_id: str,
    drug_direct_targets: list[str],
    drug_expanded_targets: list[str],
    disease_module_genes: list[str],
    disease_pathway_vector: dict[str, float],
    disease_hpo_terms: list[str],
    model_path: str = "models/repurposing_xgb.pkl",
) -> dict:
    """
    Returns:
        score: float [0,1]
        embedding_similarity: float
        cross_compat_sources: list of (indication, bridge_score) tuples
        shap_values: dict of feature -> contribution
    """
```

---

## 8. Validation Framework

### 8.1 Validation Tier 1 — Standard ML Metrics (Training Data)

**Protocol:** Stratified 5-fold cross-validation on training pairs.

| Metric | Target Threshold | Rationale |
|---|---|---|
| AUC-ROC | ≥ 0.85 | Standard classification quality |
| AUC-PR | ≥ 0.40 | Precision-recall; more informative under 10:1 imbalance |
| Hits@10 | ≥ 0.30 | Fraction of true repurposings recovered in top-10 candidates |
| Hits@50 | ≥ 0.60 | Fraction recovered in top-50 |
| Mean Reciprocal Rank | ≥ 0.20 | Average 1/rank of true repurposing |

**Baseline comparisons:**

| Baseline | Description |
|---|---|
| Random ranking | AUC-ROC = 0.50 by definition |
| Drug phase only | Rank by max clinical phase (frequency-based) |
| Heuristic equal-weight | Current pipeline: 1/3 × 3 |
| Pathway Jaccard only | Single-criterion C1-A |
| Network proximity only | Single-criterion C2 |

The ML model must beat all five baselines on AUC-ROC and AUC-PR.

### 8.2 Validation Tier 2 — Known Repurposing Recovery (Held-Out Set)

A held-out set of **historically validated rare disease repurposings** (not in training data):

| Drug | Source Indication | Target Disease | Expected Rank | Validation Signal |
|---|---|---|---|---|
| **everolimus** | Renal cell carcinoma | TSC (tuberous sclerosis) | Top 3 | FDA approved 2012 for TSC |
| **sirolimus** | Organ transplant rejection | TSC | Top 3 | mTOR1 mechanism identical |
| **selumetinib** | Thyroid cancer | NF1 plexiform neurofibromas | Top 3 | FDA approved 2020 for NF1 |
| **trametinib** | Melanoma (BRAF V600E) | NF1 | Top 5 | MEK inhibitor in NF1 trials |
| **rapamycin (sirolimus)** | Transplant | SYNGAP1 (downstream mTOR) | Top 10 | Phase 2 trials ongoing |
| **metformin** | Type 2 diabetes | Various RASopathies | Top 20 | AMPK/mTOR downstream effects |

**Pass criteria:** All drugs in this table must appear in the model's top-N ranked candidates for their respective target disease, where N is the column "Expected Rank."

### 8.3 Validation Tier 3 — Negative Control Exclusion

The model must assign low scores to mechanistically inappropriate drugs:

| Drug | Against Disease | Expected Score | Rationale |
|---|---|---|---|
| **imatinib** (BCR-ABL inhibitor) | SYNGAP1 | < 0.15 | ABL1 is not in SYNGAP1/mTOR/RAS pathway |
| **imatinib** | TSC | < 0.15 | BCR-ABL mechanism orthogonal to mTOR |
| **imiglucerase** (glucocerebrosidase replacement) | SYNGAP1 | < 0.10 | Lysosomal enzyme replacement; no CNS-pathway biology |
| **insulin glargine** | NF1 | < 0.20 | Metabolic hormone; no RAS/MAPK mechanism |
| **adalimumab** (TNF-α inhibitor) | TSC | < 0.15 | Immunological target; no mTOR connection |

**Pass criteria:** All negative controls must score below their "Expected Score" threshold.

### 8.4 Validation Tier 4 — Cross-Compatibility Specific Tests

Tests for the `CrossCompat(d, a→b)` function:

| Test | Drug | Source | Target | Min CrossCompat | Rationale |
|---|---|---|---|---|---|
| CC-1 | sirolimus | TSC | SYNGAP1 | ≥ 0.50 | mTOR identical; strong bridge |
| CC-2 | selumetinib | NF1 | SYNGAP1 | ≥ 0.45 | MEK/RAS shared downstream |
| CC-3 | everolimus | TSC | Costello syndrome | ≥ 0.40 | RAS→mTOR pathway bridge |
| CC-4 | trametinib | NF1 | Costello syndrome | ≥ 0.55 | Both HRAS/NRAS/MEK hyperactivation |
| CC-5 | imatinib | CML | TSC | < 0.20 | Disease bridge TSC↔CML is < 0.25 |
| CC-6 | imiglucerase | Gaucher | SYNGAP1 | < 0.10 | Bridge(Gaucher, SYNGAP1) < 0.15 |

### 8.5 Validation Tier 5 — Sensitivity Analysis

Re-run validation after perturbing individual features:

| Perturbation | Expected Behavior |
|---|---|
| Remove Morgan fingerprint features | AUC-ROC drops ≤ 0.05 (FP = supporting feature, not dominant) |
| Remove network proximity feature | Top-3 ranking stable for TSC/NF1 test cases |
| Remove HPO features | HPO-heavy diseases (Rett, Angelman) may reorder; top picks stable for mTOR diseases |
| Replace GO-BP vector with random noise | AUC-ROC drops > 0.05 (GO-BP should be informative) |
| Swap positive/negative labels | AUC-ROC ≈ 0.50 (sanity check: model learns real signal, not artifacts) |

### 8.6 Validation Tier 6 — SHAP Explainability Check

For each top-ranked positive prediction (e.g., sirolimus → SYNGAP1), SHAP should identify biologically plausible top contributors:

| Drug | Top Disease | Expected Top SHAP Features |
|---|---|---|
| sirolimus | SYNGAP1 | `pathway_jaccard` (mTOR shared), `embedding_cosine`, `network_proximity` |
| selumetinib | SYNGAP1 | `pathway_jaccard` (MAPK shared), `target_gene_overlap` (NF1≈SYNGAP1 module) |
| imatinib | SYNGAP1 | All features low — negative prediction must be coherent |

**Pass criteria:** For positive predictions, the top-3 SHAP features must be biologically interpretable (not noise features like `clinical_phase` alone). A manual review by a domain expert (or the pipeline authors) must confirm the explanation makes biological sense.

---

## 9. Implementation Roadmap

### Phase 1: Data Collection & Feature Matrix (Week 1)

- [ ] Download DrugCentral `drug.target.interaction` table (PostgreSQL export or API)
- [ ] Pull OpenTargets `known_drug` evidence strings via GraphQL (phase ≥ 3)
- [ ] Fetch Orphanet drug–disease XML
- [ ] Map all drugs to ChEMBL IDs; map all diseases to EFO + OMIM IDs
- [ ] Generate negative pairs: random + hard negatives
- [ ] Compute pairwise interaction features using existing `src/criteria/` code
- [ ] Compute drug encoder inputs: SMILES → Morgan FP via RDKit, GO-BP enrichment via Enrichr
- [ ] Compute disease encoder inputs: gene set → HPO via `src/loaders/hpo_loader.py`
- [ ] Serialize feature matrix to `training_data/feature_matrix.pkl`

**Deliverable:** Training matrix `X` (N×F), label vector `y`, train/val/test split indices

### Phase 2: MLP Encoders + XGBoost (Week 2)

- [ ] Implement `DrugEncoder` (MLP, PyTorch)
- [ ] Implement `DiseaseEncoder` (MLP, PyTorch)
- [ ] Pre-train drug encoder on Tanimoto similarity regression
- [ ] Pre-train disease encoder on HPO semantic similarity regression
- [ ] Train XGBoost Stage 2 on concatenated features
- [ ] Run Validation Tier 1 (cross-validation metrics)
- [ ] Run Validation Tier 2 (held-out repurposing recovery)
- [ ] Run Validation Tier 3 (negative controls)

**Deliverable:** `models/drug_encoder.pt`, `models/disease_encoder.pt`, `models/repurposing_xgb.pkl`

### Phase 3: Cross-Compatibility & Integration (Week 3)

- [ ] Implement `CrossCompat(d, a→b)` with Disease Bridge Map
- [ ] Precompute Bridge Map over Layer 2 disease set
- [ ] Implement `score_drug_ml()` inference function
- [ ] Integrate into `src/scoring/scorer.py` as `criterion4_ml`
- [ ] Run Validation Tier 4 (cross-compatibility tests)
- [ ] Add SHAP output to evidence ledger
- [ ] Update Streamlit UI to display ML score and top SHAP contributors

**Deliverable:** Fully integrated ML criterion in pipeline; updated Streamlit app

### Phase 4: RGCN Upgrade (Stretch Goal, Week 4)

- [ ] Build heterogeneous graph using `src/loaders/` data sources
- [ ] Implement RGCN using PyTorch Geometric
- [ ] Fine-tune on repurposing link prediction
- [ ] Compare RGCN embeddings vs. MLP embeddings on Tier 1–4 metrics
- [ ] Deploy best-performing model

**Deliverable:** RGCN model with performance report vs. MLP baseline

---

## 10. Dependencies

### New Python Packages

```
rdkit-pypi>=2023.9         # Morgan fingerprints from SMILES
xgboost>=2.0               # Stage 2 classifier
shap>=0.44                 # SHAP explainability
torch>=2.2                 # MLP encoders
torch-geometric>=2.5       # RGCN (stretch goal)
scikit-learn>=1.4          # Cross-validation, preprocessing
imbalanced-learn>=0.12     # SMOTE for class imbalance handling
```

### New Data Files (not in existing pipeline)

```
training_data/
  drugcentral_pairs.csv        # Drug–indication pairs from DrugCentral
  opentargets_positive.csv     # Phase ≥ 3 clinical associations
  orphanet_drugs.csv           # Rare disease drug annotations
  negative_pairs.csv           # Generated negatives
  hpo_ic_scores.pkl            # HPO information content (for semantic similarity)
  drug_smiles.csv              # ChEMBL SMILES strings
models/
  drug_encoder.pt
  disease_encoder.pt
  repurposing_xgb.pkl
  bridge_map.pkl               # Precomputed disease bridge scores
```

---

## 11. Risk Assessment & Mitigations

| Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|
| Insufficient training data for rare diseases | High | High | Augment with synthetic negatives; transfer from common disease repurposings; down-weight out-of-domain pairs |
| Overfitting to training diseases (e.g. model memorizes TSC) | Medium | High | Leave-disease-out cross-validation; ensure held-out test set contains unseen diseases |
| Morgan FP cannot distinguish stereoisomers | Low | Medium | Add 3D conformer features (RDKit ETKDG); use InChI key-based deduplication |
| API rate limits at feature computation scale | High | Medium | Precompute and cache all features; batch Enrichr/STRING calls; offline computation phase |
| DisGeNET HTML SPA blocking REST access | Known | Low | Already handled in existing pipeline; fallback to OpenTargets + seed genes |
| SHAP values not biologically interpretable | Low | Medium | Include molecular biologist review as explicit validation step; document any non-intuitive SHAP results |
| XGBoost calibration poor at extreme scores | Low | Medium | Apply Platt scaling or isotonic regression post-training; validate calibration plot |

---

## 12. Success Criteria (Go/No-Go Gate)

The model is accepted into the pipeline if and only if **all** of the following hold:

| # | Criterion | Threshold |
|---|---|---|
| S1 | AUC-ROC on held-out test set | ≥ 0.82 |
| S2 | AUC-PR on held-out test set | ≥ 0.38 |
| S3 | everolimus in top-3 for TSC | Required |
| S4 | selumetinib in top-3 for NF1 | Required |
| S5 | imatinib score < 0.15 for SYNGAP1 | Required |
| S6 | imiglucerase score < 0.10 for SYNGAP1 | Required |
| S7 | sirolimus CrossCompat(TSC→SYNGAP1) ≥ 0.50 | Required |
| S8 | SHAP top feature for sirolimus→SYNGAP1 is biologically valid | Manual review required |
| S9 | No regression vs. heuristic pipeline on existing test cases | Scores within ±0.05 |

If any required criterion fails, the model is not merged into the pipeline and a root-cause analysis is filed before retry.

---

## 13. References

1. Menche J, et al. (2015). Uncovering disease-disease relationships through the incomplete interactome. *Science*, 347(6224).
2. Zitnik M, et al. (2018). Modeling polypharmacy side effects with graph convolutional networks. *Bioinformatics*, 34(13).
3. Gysi DM, et al. (2021). Network medicine framework for identifying drug repurposing opportunities for COVID-19. *PNAS*, 118(19).
4. Lotfi Shahreza M, et al. (2018). A review of network-based approaches to drug repositioning. *Briefings in Bioinformatics*, 19(5).
5. Schlichtkrull M, et al. (2018). Modeling relational data with graph convolutional networks. *ESWC 2018*.
6. Lundberg SM & Lee SI (2017). A unified approach to interpreting model predictions. *NeurIPS 2017*.
7. OpenTargets Platform 23.12 documentation — https://platform-docs.opentargets.org
8. DrugCentral (2024) — https://drugcentral.org
9. Orphanet Drug Database — https://www.orpha.net

---

*This document is the validation plan for the ML extension to the HackRare 2026 Drug Repurposing Navigator. Implementation begins upon approval of the go/no-go gates in Section 12.*
