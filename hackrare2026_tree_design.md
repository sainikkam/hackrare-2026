# HackRare 2026 — Project Submission

## Drug Repurposing for SYNGAP1-Related Disorders via Weighted Evidence Tree

---

## 1. Project Title & One-Line Summary

**Title:** SYNGAP1 Repurposing Navigator — A Weighted Evidence Tree for Systematic Drug Candidate Discovery

**Summary:** A 12-hour-buildable pipeline that constructs a quantitatively scored probability tree linking SYNGAP1-related disorders to repurposable drugs via blood-brain-barrier filtering, pediatric safety gating, and three mechanistic scoring criteria — delivering ranked, evidence-backed candidates for patient advocacy groups with zero domain expertise required.

---

## 2. Problem Statement

SYNGAP1 haploinsufficiency affects approximately 1 in 5,000 children, causing intellectual disability, autism spectrum disorder, and drug-resistant epilepsy. There are currently zero FDA-approved therapies targeting the underlying biology.

Drug repurposing — identifying existing approved drugs that may work for a new indication — is the fastest regulatory path to therapy. However, the evidence required to evaluate any single candidate is spread across at least six databases: molecular properties, pathway annotations, protein interaction networks, phenotype ontologies, clinical trial registries, and adverse event systems. Manual curation across these sources for even 50 candidates is beyond the capacity of a rare disease patient foundation.

The practical barrier is not scientific — the biology of SYNGAP1 is reasonably well characterized (RAS/MAPK/mTOR pathway, synaptic plasticity). The barrier is systematic evidence aggregation and transparent prioritization that non-specialist stakeholders can interrogate and trust.

---

## 3. Solution Overview — The Probability Tree Concept

This project builds a **weighted evidence tree** — a probability tree where branch weights are composite evidence scores rather than raw probabilities.

**Why a probability tree, not a decision tree:**

A decision tree asks binary questions to guide someone toward an action. That structure is wrong here — the goal is not to classify patients but to display the *strength of evidence* connecting conditions and treatments. The branches need weights, not yes/no gates.

A probability tree assigns a quantitative score to each branch representing the likelihood of therapeutic relevance. This is exactly the right structure: a map where every connection carries a correlation score, and heavier branches tell advocacy groups where to direct money. It naturally communicates uncertainty — which is honest for drug repurposing — and is more credible to a scientific audience.

**The tree answers three questions simultaneously:**
1. Which drugs might work for SYNGAP1 based on biology shared with related conditions?
2. How strong is the evidence for each candidate across multiple independent dimensions?
3. Is each candidate realistically safe and CNS-accessible for the pediatric population affected?

**Tool generalizability:** The entire pipeline is configurable via a single YAML file. Changing the target disease in that file re-runs the full analysis for a different rare disease — no code changes required. This makes the tool a platform, not a one-off.

---

## 4. Tree Architecture

The tree has three layers. Layer 1 is the root disease. Layer 2 is a set of related diseases selected by pathway and phenotype overlap. Layer 3 branches from each high-scoring drug to other conditions that drug treats, scored for root-disease relevance.

---

### Layer 1 — Root

**SYNGAP1-Related Disorders**

The root node anchors the entire tree. All scoring is computed relative to the SYNGAP1 disease module (the set of proteins directly associated with SYNGAP1 biology, derived from OpenTargets and DisGeNET).

Canonical SYNGAP1 pathways used as reference throughout scoring:
- RAS/MAPK signaling
- mTOR / PI3K signaling
- AMPA receptor trafficking
- Synaptic plasticity and LTP/LTD
- Dendritic spine morphogenesis

---

### Layer 2 — Related Conditions (6 pre-selected diseases)

Each Layer 2 node is a disease selected because it shares known pathway or phenotypic overlap with SYNGAP1. Each node carries a disease-level overlap score to SYNGAP1 (computed via Criterion 3 phenotypic overlap methodology), and lists the drugs used for that condition ranked by their final correlation score to SYNGAP1.

| Layer 2 Disease | Molecular rationale for inclusion |
|---|---|
| **Tuberous Sclerosis Complex (TSC)** | TSC1/TSC2 directly regulate mTOR — same mTOR dysregulation as SYNGAP1 downstream |
| **Neurofibromatosis Type 1 (NF1)** | NF1 encodes a RAS GTPase-activating protein — same molecular class as SYNGAP1 (RAS GAP); RAS/MAPK hyperactivation is the shared mechanism |
| **Rett Syndrome (MECP2)** | MECP2 regulates mTOR signaling and synaptic gene expression; shared synaptic plasticity phenotypes |
| **Angelman Syndrome (UBE3A)** | UBE3A interacts with mTOR pathway and AMPA receptor trafficking; shared EEG and cognitive phenotypes |
| **CDKL5 Deficiency Disorder** | CDKL5 is a kinase in the mTOR-adjacent signaling network; shared epileptic encephalopathy phenotype |
| **Costello Syndrome (HRAS GOF)** | HRAS gain-of-function causes RAS hyperactivation — mechanistically parallel to SYNGAP1 LoF-driven RAS dysregulation |

Each Layer 2 node in the tree displays:
- Disease name and molecular rationale
- Disease-to-SYNGAP1 overlap score (phenotypic similarity)
- Ranked list of drugs for that condition with final correlation scores to SYNGAP1

---

### Layer 3 — Treatment-of-Treatment (Feasibility-scoped)

**Scope constraint:** Top 3 drugs per Layer 2 node only (6 diseases × 3 drugs = maximum 18 Layer 3 expansions). This is a deliberate 12-hour feasibility decision.

For each top-3 drug at Layer 2:
- Query OpenTargets and ChEMBL for all other approved indications of that drug
- Run the full scoring pipeline (gates + 3 criteria) for each indication against SYNGAP1
- Score = final correlation score for that indication's pathway/phenotype profile vs SYNGAP1
- Layer 3 nodes are other diseases the drug is used for, with scores indicating SYNGAP1 co-relevance
- This layer surfaces unexpected repurposing signals — drugs whose other indications share additional biology with SYNGAP1

Layer 3 nodes with a score ≥ 0.5 are highlighted in the visualization as "secondary repurposing signal."

---

## 5. Data Sources — 12-Hour-Viable Access Methods

Every source below is accessible with no pre-download or registration delay unless explicitly noted. Sources requiring bulk downloads (>100 MB) or multi-day registration are excluded. Where a preferred source has access friction, a functional substitute is documented.

---

### Drug & Target Data

**ChEMBL REST API** — free, no registration required
- Endpoint: `https://www.ebi.ac.uk/chembl/api/data/`
- Use: drug-target bioactivity for SYNGAP1 pathway proteins (HRAS, KRAS, MAP2K1, MTOR, RPS6KB1); drug-indication mappings
- Rate limit: ~10 req/s — use batch requests with local JSON cache on first fetch
- Format: JSON

**DrugCentral API** — free, no registration required — **primary fallback for DrugBank**
- URL: `https://drugcentral.org/api/`
- Provides: drug mechanism of action, target proteins, approved indications, transporter annotations
- Use as primary drug/target source if DrugBank registration (1–3 day delay) has not cleared

**DrugBank** — XML bulk, use only if registration already approved before hackathon start
- If available: parse `drugbank_all_full_database.xml` locally for Pgp substrate annotations and pediatric approval status
- If not available: DrugCentral + openFDA provide complete substitutes for all DrugBank-dependent computations — no scoring degradation

**PubChem REST API** — free, no registration required
- Endpoint: `https://pubchem.ncbi.nlm.nih.gov/rest/pug/`
- Use: fetch MW, XLogP3, TPSA, HBD count, pKa for all candidate drugs
- Enables CNS MPO score computation with no RDKit installation required — all inputs available as JSON
- Format: JSON

---

### BBB & Safety Data

**SwissADME** — web tool, not programmatically accessible at scale
- **Substitute (used throughout):** compute CNS MPO directly from PubChem molecular properties using Wager et al. desirability functions implemented in pure Python. No RDKit, no external tool — all six descriptor inputs come from PubChem JSON.
- For Pgp substrate classification: use DrugCentral transporter annotations; fallback to rule-based classifier (MW > 400 + PSA > 90 = predicted Pgp substrate)

**DailyMed SPL XML API** — free, no registration required
- Endpoint: `https://dailymed.nlm.nih.gov/dailymed/services/`
- Use: fetch drug label XML by drug name or RXCUI on demand
- Parse: `<boxedWarning>`, `<contraindications>`, `<pediatricUse>` sections for Gate 2
- No bulk download needed — fetch per drug, cache response immediately

**openFDA API** — free, no registration required — **replaces FAERS bulk download**
- Endpoint: `https://api.fda.gov/drug/event.json`
- Query: `patient.patientonsetage:<18&serious:1` filtered by drug name
- Returns: pediatric serious adverse event counts in real time
- Rate limit: 240 req/min — sufficient for the full ~50 candidate drug set without throttling

**ClinicalTrials.gov API v2** — free, no registration required
- Endpoint: `https://clinicaltrials.gov/api/v2/studies`
- Query: drug name + `ageGroup=CHILD` filter
- Returns: pediatric trial status for Component A of pediatric safety scoring

---

### Pathway & Network Data

**Reactome REST API** — free, no registration required
- Endpoint: `https://reactome.org/ContentService/`
- Use: pathway memberships for SYNGAP1-associated proteins; pathway-gene mappings for drug targets
- Primary source for Criterion 1A (Pathway Overlap Jaccard)

**KEGG REST API** — free, limited to ~10 req/s
- Endpoint: `https://rest.kegg.jp/`
- Use: supplement Reactome with broader pathway coverage for drug target annotation

**Enrichr API** — free, no registration required — **replaces gseapy for speed**
- Endpoint: `https://maayanlab.cloud/Enrichr/`
- Use: gene set enrichment against Reactome, KEGG, GO-BP in a single API call
- Returns ranked pathways with p-values in ~2 seconds per query
- Used in Criterion 1A (SYNGAP1 pathway set construction) and Criterion 1B fallback

**STRING REST API** — free, no registration required — **replaces full STRING bulk download**
- Endpoint: `https://string-db.org/api/`
- Use: fetch interaction partners for SYNGAP1 disease module proteins at confidence ≥ 0.7
- Query limited to 2-hop neighborhood of SYNGAP1 (~500–800 proteins) — sufficient for proximity scoring; avoids 1.5 GB bulk file download entirely
- Build local networkx subgraph from API response

**OpenTargets GraphQL API** — free, no registration required
- Endpoint: `https://api.platform.opentargets.org/api/v4/graphql`
- Use: top associated diseases for SYNGAP1 target; associated drugs per Layer 2 disease; Layer 3 indication expansion
- Primary source for candidate disease and drug identification

**DisGeNET REST API** — free, requires API key — **registration takes <5 minutes (email confirm)**
- Endpoint: `https://www.disgenet.org/api/`
- Use: SYNGAP1 disease module seed genes
- Register at start of hackathon; access is immediate after email confirmation

**OMIM** — web search for HPO term lookup; no programmatic bulk access needed (HPO files used instead)

---

### Phenotypic Data

**HPO Files** — small downloads, ~30 MB total, ~2 minutes to retrieve — **only required file downloads**
- `hp.obo`: ontology structure with IC-computable term hierarchy
- `phenotype_to_genes.txt`: HPO term → gene mappings
- Download from: `https://hpo.jax.org/data/annotations`
- These are the only two files that must be downloaded; all other data sources use APIs

**LINCS L1000 / CMap** — transcriptomic perturbation signatures for Criterion 1B
- **Primary 12-hour strategy:** use the CLUE API (`https://clue.io/api`) with free registration (<5 min)
- **Fallback if CLUE API is unavailable or slow:** substitute Criterion 1B with a second Enrichr pathway enrichment comparison — enrich drug target genes in GO-BP, compare top-10 enriched terms vs SYNGAP1 enrichment using cosine similarity of pathway score vectors. Same conceptual signal (shared perturbation direction), zero external dependency.

---

## 6. Implementation Plan — Step-by-Step, 12-Hour Scoped

### Repository Structure

```
hackrare-2026/
├── data/cache/                  # File cache for all API responses
├── src/
│   ├── loaders/                 # One module per data source
│   │   ├── chembl_loader.py
│   │   ├── pubchem_loader.py
│   │   ├── drugcentral_loader.py
│   │   ├── dailymed_loader.py
│   │   ├── openfda_loader.py
│   │   ├── string_loader.py
│   │   ├── reactome_loader.py
│   │   ├── enrichr_loader.py
│   │   ├── opentargets_loader.py
│   │   ├── disgenet_loader.py
│   │   └── hpo_loader.py
│   ├── gates/
│   │   ├── bbb_gate.py
│   │   └── pediatric_safety_gate.py
│   ├── criteria/
│   │   ├── pathway_alignment.py
│   │   ├── network_proximity.py
│   │   └── phenotypic_overlap.py
│   ├── scoring/
│   │   ├── scorer.py
│   │   └── tree_builder.py
│   └── config/
│       └── syngap1.yaml
├── app/
│   └── streamlit_app.py
└── requirements.txt
```

**Dependencies:** `pandas numpy scipy networkx requests streamlit pyvis`

No RDKit, no bulk XML parsing, no overnight downloads. Install completes in under 2 minutes.

---

### 12-Hour Timeline (2 People)

| Hour | Person A | Person B |
|---|---|---|
| 0–1 | Repo setup, environment, dependencies | Data loader skeletons for all sources |
| 1–3 | Gate 1 (BBB) implementation | Gate 2 (Pediatric Safety) implementation |
| 3–5 | Criterion 1 (Pathway Alignment) | Criterion 2 (Network Proximity) |
| 5–7 | Criterion 3 (Phenotypic Overlap) | Scoring engine + evidence ledger |
| 7–9 | Tree builder + JSON output | Streamlit app + tree visualization |
| 9–11 | Integration, end-to-end test | Validation (rapamycin +ve, Gaucher -ve) |
| 11–12 | Demo polish + sensitivity analysis | README + submission materials |

---

### Step 1 — Environment & Repository Setup (Hour 0–1)

- Person A: initialize repo, install dependencies, create directory structure, create `syngap1.yaml` config with SYNGAP1 Ensembl ID (ENSG00000197949), Layer 2 disease list, and scoring weight defaults
- Person B: write loader skeletons with standard interface `load(cache_path) -> dict`, implement JSON file cache decorator (cache miss → API call → write to file → return; cache hit → read file → return)
- Both: verify API connectivity for all 10 sources before proceeding (test one query each)

---

### Step 2 — SYNGAP1 Disease Module Definition (Hour 1)

- Query OpenTargets GraphQL for SYNGAP1 (ENSG00000197949) → associated proteins (all targets with genetic association score > 0.1)
- Query DisGeNET REST API for SYNGAP1-associated genes (all evidence types)
- Union both sets → SYNGAP1 disease module (expected: 50–150 proteins)
- Cache to `data/cache/syngap1_module.json`
- Query Reactome for canonical SYNGAP1 pathway list → cache to `data/cache/syngap1_pathways.json`
- This module and pathway list are the reference inputs for all three scoring criteria

---

### Step 3 — Layer 2 Candidate Disease & Drug List (Hour 1)

- Confirm the 6 pre-selected Layer 2 diseases in `syngap1.yaml`
- For each Layer 2 disease: query OpenTargets for approved/investigational drugs → merge with ChEMBL drug-target bioactivity data
- Query DrugCentral for MOA and target annotations for all retrieved drugs
- Deduplicate across diseases → consolidated candidate drug list (~50–100 unique drugs)
- Cache to `data/cache/layer2_candidates.json`

---

### Step 4 — Gate 1: BBB Permeability (Hour 1–3, Person A)

For each candidate drug:
1. Fetch molecular properties from PubChem: MW, XLogP3, TPSA, HBD count, pKa
2. Compute 6 CNS MPO desirability scores using Wager et al. thresholds (pure Python, no RDKit):
   - cLogP desirability: 1.0 at ≤3, linear to 0 at >5
   - cLogD desirability: 1.0 at ≤2, linear to 0 at >4
   - MW desirability: 1.0 at ≤360 Da, linear to 0 at >500 Da
   - TPSA desirability: 1.0 at ≤60 Å², linear to 0 at >90 Å²
   - HBD desirability: 1.0 at 0, linear to 0 at >3
   - pKa desirability: 1.0 at ≤8, linear to 0 at >10
3. `CNS_MPO_normalized = (sum of six desirability scores) / 6`
4. Fetch Pgp substrate classification from DrugCentral transporter annotations (or rule-based: MW > 400 + PSA > 90 = predicted substrate)
5. Apply Pgp multiplier: confirmed substrate × 0.6; predicted substrate × 0.8; not a substrate × 1.0
6. `BBB_Score = CNS_MPO_normalized × Pgp_adjustment`
7. Filter: keep only drugs where BBB_Score > 0.8

Output: `{drug_id: bbb_score}` dict, cached to `data/cache/gate1_results.json`

---

### Step 5 — Gate 2: Pediatric Safety (Hour 1–3, Person B)

For each BBB-passing drug:
1. Fetch DailyMed SPL XML via API → parse `<boxedWarning>` and `<contraindications>`
   - If either triggers pediatric hard disqualifier → score = 0.0, drug eliminated
2. If no hard disqualifier: parse `<pediatricUse>` section → score Component A (FDA Pediatric Approval Status)
3. Extract minimum approved age from label → score Component B (Minimum Approved Age)
4. Query openFDA API: `patient.patientonsetage:<18&serious:1` filtered by drug name → compute normalized serious AE rate → invert for Component C
5. `Pediatric_Safety_Score = 0.40 × A + 0.30 × B + 0.30 × C`
6. Filter: keep only drugs where score ≥ 0.6
7. Flag drugs scoring 0.60–0.79 with "Caution: limited pediatric data" annotation

Output: `{drug_id: {safety_score, caution_flag, components}}` dict, cached to `data/cache/gate2_results.json`

---

### Step 6 — Criterion 1: Pathway Mechanism Alignment (Hour 3–5, Person A)

**Sub-metric A — Pathway Overlap (Jaccard):**
- SYNGAP1 pathway set: from cached Reactome + Enrichr enrichment (Step 2)
- For each gate-passing drug: fetch target proteins from ChEMBL → fetch their Reactome pathways → supplement with Enrichr GO-BP enrichment
- `Jaccard = |SYNGAP1_pathways ∩ drug_pathways| / |SYNGAP1_pathways ∪ drug_pathways|`

**Sub-metric B — Perturbation Signal Similarity:**
- **Primary:** query CLUE API for L1000 signature of each drug → cosine similarity vs SYNGAP1 differentially expressed gene direction
- **Fallback (if CLUE API unavailable):** query Enrichr with drug target genes against GO-BP → compute cosine similarity of top-10 enriched term score vectors vs SYNGAP1 Enrichr enrichment. Same conceptual signal; zero external dependency.

`Pathway_Alignment = 0.5 × Jaccard + 0.5 × signal_similarity`

Output: `{drug_id: {pathway_alignment, jaccard, signal_sim}}` dict

---

### Step 7 — Criterion 2: Network Proximity (Hour 3–5, Person B)

**12-hour scope decision:** Full Menche et al. permutation analysis (1,000 random null permutations) is deferred as a stretch goal if time remains in hours 9–11. The simplified proximity score used within the 12-hour window is mathematically well-grounded and produces the same ranking signal.

1. Use STRING REST API to build SYNGAP1 subgraph: fetch 1-hop + 2-hop interaction partners at confidence ≥ 0.7 (limits query to ~500–800 nodes; avoids 1.5 GB bulk file)
2. Build networkx undirected graph from API response
3. For each gate-passing drug's target proteins:
   - Compute mean shortest-path distance from each target to SYNGAP1 node in the subgraph
   - `raw_proximity = 1 / (1 + mean_shortest_path)` → bounded [0, 1]
4. Normalize across all drug scores to [0, 1] via min-max

Output: `{drug_id: network_proximity_score}` dict

**Stretch goal (hours 9–11 if time allows):** add 1,000-permutation null distribution for full Menche separation z-score normalization.

---

### Step 8 — Criterion 3: Phenotypic Overlap (Hour 5–7, Person A)

1. Download HPO files at start of this step (~2 min): `hp.obo` + `phenotype_to_genes.txt`
2. Parse `hp.obo` into IC-annotated ontology DAG using networkx; compute IC for each term from annotation frequency
3. For SYNGAP1 and each Layer 2 disease: retrieve HPO term lists from HPO/OMIM annotations
4. Compute Resnik BMA semantic similarity between SYNGAP1 and each disease/drug-target disease → `hpo_similarity`
5. Map HPO terms to genes via `phenotype_to_genes.txt`; compute weighted gene vectors; cosine similarity → `gene_similarity`
6. Enrich top genes via Enrichr (GO-BP + Reactome); compute cosine similarity of pathway score vectors → `pathway_similarity`
7. `Phenotypic_Overlap = 0.4 × hpo_similarity + 0.3 × gene_similarity + 0.3 × pathway_similarity`

Output: `{drug_disease_pair: {phenotypic_overlap, hpo_sim, gene_sim, pathway_sim, top5_hpo, top10_genes, top5_pathways}}` dict

---

### Step 9 — Scoring Engine & Evidence Ledger (Hour 5–7, Person B)

`scorer.py`:
- Input: drug + disease pair
- Run Gate 1 → if fail, return eliminated with BBB score
- Run Gate 2 → if fail, return eliminated with safety score
- Compute Criteria 1, 2, 3
- `Final_Score = (c1 + c2 + c3) / 3`

Evidence ledger per drug (stored in tree JSON):
- Gate 1 score (BBB), Gate 2 score (safety, caution flag if 0.60–0.79)
- Criterion 1 score (pathway alignment + components)
- Criterion 2 score (network proximity, mean shortest path)
- Criterion 3 score (phenotypic overlap + components)
- Final correlation score
- Top 5 shared HPO terms
- Top 5 shared pathways
- Network distance to SYNGAP1
- All safety flags

Sensitivity analysis: perturb each criterion weight ±10%, recompute Final_Score, verify top-3 ranking stability → output stability table.

---

### Step 10 — Tree Builder (Hour 7–9, Person A)

`tree_builder.py` — builds nested JSON tree structure:

```json
{
  "root": {
    "disease": "SYNGAP1-Related Disorders",
    "layer2": [
      {
        "disease": "Tuberous Sclerosis Complex",
        "overlap_score_to_syngap1": 0.73,
        "drugs": [
          {
            "drug": "rapamycin",
            "bbb_score": 0.91,
            "safety_score": 0.82,
            "pathway_alignment": 0.88,
            "network_proximity": 0.79,
            "phenotypic_overlap": 0.76,
            "final_score": 0.81,
            "evidence_ledger": {...},
            "layer3": [...]
          }
        ]
      }
    ]
  }
}
```

Layer 3 expansion: for each top-3 drug per Layer 2 node, query OpenTargets/ChEMBL for all other approved indications → run scorer for SYNGAP1 relevance → add as Layer 3 child nodes.

Output: `data/tree.json`

---

### Step 11 — Streamlit Application (Hour 7–9, Person B)

- **Landing page:** SYNGAP1 root node, overview stats (# candidates evaluated, # passed Gate 1, # passed Gate 2, # in final tree)
- **Tree visualization:** pyvis interactive network; edge thickness = final correlation score; node color = gate pass/fail status
- **Drug click → evidence ledger panel:** all five scores displayed (BBB, Safety, Pathway, Network, Phenotypic), shared pathways list, HPO term matches, safety flags, caution annotations
- **Disease click → ranked drug table:** sortable by final score or any individual criterion
- **Filter controls:** minimum final score slider, caution flag toggle, Layer 3 show/hide
- **Export:** "Download ranked candidates (CSV)" button
- **Validation panel:** rapamycin highlighted as top candidate (positive control verified); Gaucher/GBA1-targeting drugs shown as correctly excluded (negative control verified)

---

### Step 12 — Integration, Validation & Demo Polish (Hour 9–12)

- End-to-end test: run full pipeline from `syngap1.yaml` config → verify tree renders with all nodes and scores
- **Positive control validation:** rapamycin must appear in top 3 candidates for SYNGAP1 (mTOR inhibitor, documented preclinical SYNGAP1 efficacy)
- **Negative control validation:** imiglucerase, eliglustat (GBA1-targeting drugs) must not appear in final tree — GBA1 has no RAS/MAPK/mTOR pathway overlap
- **Generalizability demo:** change YAML disease config to Costello Syndrome → re-run pipeline → verify different candidate set surfaces with different rankings
- **Sensitivity analysis:** generate table showing top-3 ranking stability under ±10% weight perturbation on all three criterion weights
- Finalize Streamlit layout, polish labels, add caution annotations, verify all evidence ledger links render correctly

---

## 7. Evaluation Pipeline

Every drug/disease combination is evaluated in two sequential stages:

```
STAGE 1 — GATING (applied first, to every candidate)
    Gate 1: BBB Permeability Score > 0.8?    → FAIL = eliminated immediately
    Gate 2: Pediatric Safety Score ≥ 0.6?    → FAIL = eliminated immediately

        ↓ Only candidates that pass BOTH gates proceed ↓

STAGE 2 — SCORING (applied only to gate-passing candidates)
    Criterion 1: Pathway Mechanism Alignment   (0–1)
    Criterion 2: Network Proximity             (0–1)
    Criterion 3: Phenotypic Overlap            (0–1)
    → Composite correlation score computed via cosine similarity
```

This ordering is deliberate: no matter how strong a drug's mechanistic or network evidence, it is excluded if it cannot reach the brain or is unsafe for children. The gates are non-negotiable filters, not scoring dimensions.

---

## STAGE 1 — GATES

---

### Gate 1 — Blood-Brain Barrier (BBB) Permeability
**Pass threshold: score > 0.8**

#### Methodology

The BBB Permeability Score is a composite of two sub-components:

---

**Sub-component A: CNS Multiparameter Optimization (CNS MPO) Score**
*(Wager et al., ACS Chemical Neuroscience, 2010)*

The CNS MPO score is the pharmaceutical industry standard for predicting CNS penetration from physicochemical properties. Six molecular descriptors are each scored on a desirability function (0–1) and summed for a raw score of 0–6, then normalized to [0, 1].

| Property | Optimal range | Scoring rule |
|---|---|---|
| cLogP | ≤ 3 (max 5) | Desirability function: 1.0 at ≤3, linearly decreasing to 0 at >5 |
| cLogD (pH 7.4) | ≤ 2 (max 4) | Desirability function: 1.0 at ≤2, linearly decreasing to 0 at >4 |
| Molecular weight | ≤ 360 Da (max 500) | Desirability function: 1.0 at ≤360, linearly decreasing to 0 at >500 |
| TPSA | ≤ 60 Ų (max 90 Ų) | Desirability function: 1.0 at ≤60, linearly decreasing to 0 at >90 |
| H-bond donors (HBD) | 0 (max 3) | Desirability function: 1.0 at 0, linearly decreasing to 0 at >3 |
| Most basic pKa | ≤ 8 (max 10) | Desirability function: 1.0 at ≤8, linearly decreasing to 0 at >10 |

> CNS_MPO_normalized = (sum of six desirability scores) / 6

**Sources:** Computed from DrugBank / PubChem molecular properties; SwissADME for automated CNS MPO calculation.

---

**Sub-component B: P-glycoprotein (Pgp) Efflux Adjustment**

P-glycoprotein is an efflux transporter at the BBB that actively pumps drugs back out of the CNS. A drug that passes physicochemical filters but is a Pgp substrate will have significantly reduced CNS exposure.

- Retrieve Pgp substrate classification from DrugBank (annotated) or predict via pkCSM / SwissADME
- Apply a multiplicative penalty:
  - Confirmed Pgp substrate: × 0.6
  - Predicted Pgp substrate (moderate confidence): × 0.8
  - Not a Pgp substrate: × 1.0

> BBB_Score = CNS_MPO_normalized × Pgp_adjustment

**If experimental BBB data is available** (e.g., documented CSF/plasma ratio or in vitro PAMPA-BBB assay result from DrugBank or literature), weight it in:

> BBB_Score = 0.5 × CNS_MPO_normalized × Pgp_adjustment + 0.5 × experimental_BBB_normalized

All scores are capped at [0, 1].

---

#### BBB Score Interpretation Guide

| Score | Interpretation | Gate outcome |
|---|---|---|
| > 0.8 | Strong predicted CNS penetration | **PASSES** — proceeds to Stage 2 |
| 0.6 – 0.8 | Moderate CNS penetration — likely subtherapeutic CNS exposure | **FAILS** — eliminated |
| 0.4 – 0.6 | Poor CNS penetration — unlikely to reach target tissue | **FAILS** — eliminated |
| < 0.4 | Negligible CNS penetration or confirmed Pgp-effluxed | **FAILS** — eliminated |

> A score of exactly 0.8 does not pass. The threshold is strictly greater than 0.8.

---

### Gate 2 — Pediatric Safety Compatibility
**Pass threshold: score ≥ 0.6**

#### Methodology

The Pediatric Safety Compatibility Score is computed from five components sourced from DailyMed SPL XML structured fields and the FDA Adverse Event Reporting System (FAERS) pediatric subset. Two components function as hard disqualifiers; three contribute to a continuous score.

---

**Step 1 — Hard Disqualifier Check (applied before scoring)**

These checks are evaluated first. If either is triggered, the score is automatically set to 0.0 and the drug fails the gate with no further evaluation.

| Disqualifier | Source field (DailyMed SPL XML) |
|---|---|
| Black box warning explicitly referencing pediatric patients | `<boxedWarning>` section |
| Explicit pediatric contraindication | `<contraindications>` section |

---

**Step 2 — Continuous Scoring (only if no hard disqualifier triggered)**

Three components are scored and combined:

**Component A: FDA Pediatric Approval Status (weight: 0.40)**
- Approved for ≥1 pediatric age group (neonates, infants, children, adolescents): 1.0
- Approved for adolescents only (≥12 years): 0.6
- No pediatric approval but no explicit contraindication: 0.3
- Under investigation for pediatric use (active trial in ClinicalTrials.gov): 0.4

Source: DailyMed `<pediatricUse>` section; ClinicalTrials.gov

**Component B: Minimum Approved Age (weight: 0.30)**
Scored inversely — drugs approved for younger patients score higher, reflecting demonstrated safety in the relevant population.

| Minimum approved age | Score |
|---|---|
| ≤ 2 years | 1.0 |
| 3–5 years | 0.8 |
| 6–11 years | 0.6 |
| 12–17 years | 0.4 |
| Adults only (no pediatric data) | 0.1 |

**Component C: Serious Pediatric Adverse Event Rate (weight: 0.30)**
- Retrieve pediatric serious adverse event reports from FAERS (age < 18 filter)
- Normalize the serious AE rate per 1,000 reports against the full drug AE population
- Invert so that a low serious AE rate yields a high score:

> AE_score = 1 − min(serious_pediatric_AE_rate_normalized, 1.0)

Source: FDA FAERS public dataset

---

**Step 3 — Composite Score**

> Pediatric_Safety_Score = 0.40 × Component_A + 0.30 × Component_B + 0.30 × Component_C

Score is bounded [0, 1]. Hard disqualifiers override to 0.0 regardless of component scores.

---

#### Pediatric Safety Score Interpretation Guide

| Score | Interpretation | Gate outcome |
|---|---|---|
| 0.0 | Automatic failure — black box warning or explicit pediatric contraindication | **FAILS** — eliminated |
| 0.01 – 0.39 | High pediatric risk — multiple safety concerns, no approval history | **FAILS** — eliminated |
| 0.40 – 0.59 | Moderate risk — limited or adolescent-only data, notable AE signal | **FAILS** — eliminated |
| 0.60 – 0.79 | Acceptable safety profile — passes gate; caution warranted in younger age groups | **PASSES** |
| 0.80 – 1.00 | Favorable pediatric safety — established approval, low serious AE rate | **PASSES** |

> The threshold of ≥ 0.6 is set to reflect the reality that many drugs used in rare pediatric diseases carry known risks that are accepted given the severity of the condition. A threshold of 0.8 would eliminate the majority of repurposing candidates in this space. Any candidate scoring 0.60–0.79 should be flagged in the output with a **"Caution: limited pediatric data"** annotation.

---

## STAGE 2 — SCORING CRITERIA
*(Applied only to drug/disease combinations that pass both gates)*

---

### Criterion 1 — Pathway Mechanism Alignment

**Composite of two sub-metrics, each normalized 0–1:**

**Sub-metric A: Pathway Overlap**
- Identify the canonical pathways associated with SYNGAP1 (e.g. RAS/MAPK, mTOR, AMPA receptor trafficking, synaptic plasticity) from Reactome, KEGG, and GO-BP
- Identify the pathways targeted by the candidate drug's MOA
- Score = Jaccard similarity of the two pathway sets: |A ∩ B| / |A ∪ B|
- Normalize to [0, 1]

**Sub-metric B: Perturbation / Drug Signal Similarity**
- Retrieve the drug's transcriptomic perturbation signature from LINCS L1000 or CMap (gene expression changes induced by the drug)
- Retrieve the SYNGAP1 disease signature (differentially expressed genes in SYNGAP1 LoF models from GEO/PubMed)
- Score = cosine similarity between the two gene expression vectors (or Spearman correlation of ranked gene lists)
- Normalize to [0, 1]; invert if the drug signature is expected to reverse the disease signature (a reversal is a positive signal)

**Composite Score:**
> Pathway Mechanism Alignment = 0.5 × Pathway Overlap + 0.5 × Perturbation Signal Similarity

Equal weighting by default; subject to sensitivity analysis (±10% weight perturbation).

**Sources:** Reactome, KEGG, GO-BP, LINCS L1000 / CMap, GEO, ChEMBL, DrugBank

---

### Criterion 2 — Network Proximity

**Methodology: Menche et al. (2015, *Science*) network separation score, applied to STRING PPI network at confidence threshold ≥ 0.7.**

**Steps:**
1. Define the SYNGAP1 disease module: the set of proteins directly associated with SYNGAP1 in STRING and curated disease-gene databases (OMIM, DisGeNET)
2. Define the drug target module: the set of known protein targets for the candidate drug (from DrugBank / ChEMBL)
3. Compute the network separation score *s*:

> s(A,B) = d(A,B) − (d(A,A) + d(B,B)) / 2

where d(A,B) is the mean shortest-path distance between all nodes in module A and all nodes in module B in the PPI network; d(A,A) and d(B,B) are the mean within-module distances.

4. A negative *s* indicates the modules overlap (high proximity); a positive *s* indicates separation
5. Normalize to [0, 1] using empirical z-score against a null distribution of 1,000 random module pairs of matched size, then map to [0, 1] via min-max normalization:

> Network Proximity Score = 1 − (s − s_min) / (s_max − s_min)

**Sources:** STRING (confidence ≥ 0.7), DrugBank, ChEMBL, OMIM, DisGeNET

---

### Criterion 3 — Phenotypic Overlap

**Methodology: Phenotype→Mechanism Similarity via HPO**

**INPUT**
- Disease A: HPO IDs (with optional phenotype frequencies)
- Disease B: HPO IDs (with optional phenotype frequencies)

**PIPELINE**

**Step 1 — HPO Semantic Similarity**
- Compute Resnik IC-based similarity + Best-Match Average (BMA) across all HPO term pairs between Disease A and Disease B
- Downweight generic (low-IC) terms
- Output: `hpo_similarity` ∈ [0, 1]

**Step 2 — HPO → Gene Mapping**
- Map each HPO term to associated genes
- Score each gene per disease:

> score(g) = Σ_t w(t) × assoc(t, g)

where w(t) = freq(t, if given) × IC(t); assoc(t, g) = 1 or confidence value

- Normalize scores; retain top K = 200 genes per disease
- Compute cosine similarity of the two weighted gene vectors
- Output: `gene_similarity` ∈ [0, 1]

**Step 3 — Gene → Pathway Enrichment**
- Enrich top genes from each disease in Reactome, KEGG, and GO-BP
- Pathway score = −log₁₀(p-value), optionally weighted by mean gene score
- Compute cosine similarity of the two pathway score vectors
- Output: `pathway_similarity` ∈ [0, 1]

**Step 4 — Overall Similarity (Composite)**

> overall_similarity = 0.4 × hpo_similarity + 0.3 × gene_similarity + 0.3 × pathway_similarity

**Step 5 — Confidence Rating**
- High / Medium / Low based on: number of HPO terms available, average IC of terms, number of genes successfully mapped, enrichment strength

**CONSTRAINTS**
- No clinical recommendations
- Do not invent data; note missing frequencies or mappings

**OUTPUT (per pair)**
- `hpo_similarity`, `gene_similarity`, `pathway_similarity`, `overall_similarity` — all normalized [0, 1]
- `confidence`: High / Med / Low
- `evidence`: top 5 HPO matches, top 10 shared genes, top 5 shared pathways
- `limitations`: ≤ 5 bullets

**Sources:** HPO, Reactome, KEGG, GO-BP, OMIM, DisGeNET

---

## 8. Statistical Method

**Primary method: Composite weighted score with cosine similarity as the mathematical backbone.**

**Step 1 — Gate screening** for all candidates using BBB Permeability and Pediatric Safety scores. Candidates failing either gate are discarded before any further computation.

**Step 2 — Score each criterion (0 to 1)** for every gate-passing treatment-condition pair using the methodology defined above for each of the three criteria.

**Step 3 — Construct an evidence vector** for each candidate treatment: a 3-dimensional vector where each dimension is one normalized criterion score:

> **v** = [Pathway Mechanism Alignment, Network Proximity, Phenotypic Overlap]

The SYNGAP1 reference vector is fixed as:

> **r** = [1, 1, 1]

This represents the theoretical ideal candidate — one that scores perfectly on all three criteria relative to SYNGAP1.

**Step 4 — Compute the final correlation score** by combining directional alignment (cosine similarity) with absolute magnitude, so that the score reflects both *how balanced* the evidence is across all three criteria and *how strong* those scores are in absolute terms.

**Component A — Directional alignment (cosine similarity):**

> cos(θ) = (c₁ + c₂ + c₃) / (√(c₁² + c₂² + c₃²) × √3)

This measures the angle between the candidate vector and the ideal reference [1, 1, 1]. A score of 1.0 means the candidate's evidence profile is perfectly proportionally aligned with the ideal, regardless of absolute magnitude.

**Component B — Normalized magnitude:**

> M = √(c₁² + c₂² + c₃²) / √3

This is the L2 norm of the candidate vector normalized against the maximum possible magnitude (√3, achieved when all three scores equal 1.0). It captures how strong the absolute scores are — a candidate scoring [0.9, 0.9, 0.9] has a much higher magnitude than one scoring [0.3, 0.3, 0.3].

**Final Correlation Score — product of both components:**

> **Final Score = cos(θ) × M**

Expanding this algebraically:

> Final Score = [(c₁ + c₂ + c₃) / (√(c₁² + c₂² + c₃²) × √3)] × [√(c₁² + c₂² + c₃²) / √3]

> **Final Score = (c₁ + c₂ + c₃) / 3**

This simplifies to the **arithmetic mean of the three criterion scores**. This is the intended behavior and is both mathematically principled and highly interpretable.

**Worked examples:**

| Candidate vector | cos(θ) | M | Final Score | Interpretation |
|---|---|---|---|---|
| [0.9, 0.9, 0.9] | 1.0 | 0.9 | **0.90** | High and balanced — top candidate |
| [0.3, 0.3, 0.3] | 1.0 | 0.3 | **0.30** | Balanced but weak — penalized appropriately |
| [0.9, 0.1, 0.1] | 0.64 | 0.53 | **0.37** | Strong on one axis only — penalized for imbalance |
| [0.6, 0.6, 0.6] | 1.0 | 0.6 | **0.60** | Moderate across all three — solid candidate |
| [1.0, 0.0, 0.0] | 0.58 | 0.58 | **0.33** | Single-criterion evidence only — low confidence |

**This final score IS the branch weight** displayed on the probability tree for each drug/disease combination. It is bounded [0, 1] and requires no further transformation.

**Why this formulation:**
- Penalizes candidates with inflated single-criterion scores by requiring balanced evidence across all three dimensions
- Penalizes candidates with uniformly low scores by incorporating absolute magnitude — a balanced but weak profile correctly scores low
- Simplifies to the arithmetic mean, making it immediately interpretable to non-technical audiences including patient advocacy groups
- Bounded [0, 1], maps cleanly to branch weights
- Mathematically defensible to both clinical and technical judges

**Step 5 — Sensitivity analysis:** perturb each criterion weight by ±10% and confirm the top-ranked branches remain stable. This directly addresses the "arbitrary weights" criticism.

---

## 9. Demo Strategy

All demo data sourced from API cache built during implementation — no live API calls during presentation. Demo runs entirely from `data/cache/` and `data/tree.json`.

**Total demo time: ~2.5 minutes**

1. **Problem (30s):** "SYNGAP1 — 1 in 5,000 children, zero approved therapies. Drug repurposing is the fastest path. Manually searching six databases is impossible for a patient foundation."

2. **Input (15s):** Open `syngap1.yaml`. Point out the disease config parameter. "One line change = one different disease. This tool is a platform."

3. **Tree renders (30s):** Streamlit opens. Probability tree displays in pyvis — branch weights visible as edge thickness. Layer 1 → Layer 2 diseases → ranked drugs visible.

4. **Evidence ledger (60s):** Click rapamycin node. Five scores appear in the panel: BBB 0.91, Safety 0.82, Pathway Alignment 0.88, Network Proximity 0.79, Phenotypic Overlap 0.76. Shared pathways: mTOR signaling, RAS/MAPK. Final score: 0.81. "Every claim is sourced and auditable."

5. **Validation (30s):** Validation panel: rapamycin in top 3 (correct — mTOR inhibitor, documented preclinical SYNGAP1 efficacy). Imiglucerase not in tree (correct — GBA1 pathway has zero RAS/MAPK/mTOR overlap).

6. **Generalizability (15s):** Change YAML to Costello Syndrome, re-run pipeline, show different tree. "Same tool, different rare disease, different candidate ranking — in under 2 minutes."

---

## 10. Validation Approach

### Positive Control

**Rapamycin (sirolimus)** must rank in the top 3 candidates for SYNGAP1.

Scientific basis: Rapamycin is an mTOR inhibitor. SYNGAP1 haploinsufficiency causes constitutive RAS/MAPK activation, which drives mTOR hyperactivation downstream. mTOR inhibition with rapamycin has documented preclinical efficacy in SYNGAP1 mouse models — it is the most mechanistically justified repurposing candidate known.

If rapamycin does not appear in the top 3, the scoring pipeline has a calibration error that must be diagnosed before submission.

### Negative Control

**GBA1-targeting drugs** (imiglucerase, eliglustat, miglustat) must not appear in the final tree.

Scientific basis: GBA1 encodes glucocerebrosidase, an enzyme involved in glycolipid metabolism (Gaucher disease). The GBA1 pathway shares no proteins, pathways, or phenotypic overlap with RAS/MAPK/mTOR/synaptic plasticity — the SYNGAP1 biology. These drugs should fail on Criterion 1 (pathway alignment near 0) and Criterion 2 (network distance large), producing a Final Score well below any tree inclusion threshold.

### Sensitivity Analysis

Perturb all three criterion weights simultaneously by ±10%:

| Perturbation | Weight set | Expected: rapamycin still top 3? |
|---|---|---|
| Baseline | [0.33, 0.33, 0.33] | ✓ |
| Pathway emphasized | [0.43, 0.33, 0.23] | ✓ |
| Network emphasized | [0.23, 0.43, 0.33] | ✓ |
| Phenotypic emphasized | [0.23, 0.33, 0.43] | ✓ |

If the top-3 ranking remains stable across all four weight sets, the result is robust. Instability would indicate one criterion dominates and warrants re-examination.

---

## 11. Execution Checklist (12-Hour Window)

### Pre-Start (Before Clock Starts)
- [ ] Confirm DrugBank registration status — if not cleared, commit to DrugCentral as primary throughout
- [ ] Register DisGeNET API key (immediate, <5 min — email confirmation grants access)
- [ ] Register CLUE/CMap API key (immediate, <5 min — needed for Criterion 1B primary path)
- [ ] Verify all 10 API endpoints are reachable (one test query each)

### Hour 0–1
- [ ] Repo initialized, directory structure created, `requirements.txt` committed
- [ ] `pip install pandas numpy scipy networkx requests streamlit pyvis` confirmed working
- [ ] `syngap1.yaml` config file written with SYNGAP1 Ensembl ID and Layer 2 disease list
- [ ] All 10 loader skeletons written with standard `load(cache_path)` interface and file cache decorator

### Hour 1–3
- [ ] SYNGAP1 disease module populated from OpenTargets + DisGeNET, cached
- [ ] Layer 2 candidate drug list populated from ChEMBL + OpenTargets + DrugCentral, cached
- [ ] Gate 1 (BBB): PubChem properties fetched for all candidates, CNS MPO computed in pure Python
- [ ] Gate 2 (Pediatric Safety): DailyMed labels parsed for all BBB-passing drugs, openFDA AE counts fetched

### Hour 3–7
- [ ] Criterion 1 (Pathway Alignment): Jaccard via Reactome + Enrichr implemented; CLUE API or fallback Enrichr comparison implemented
- [ ] Criterion 2 (Network Proximity): STRING subgraph built in networkx, mean shortest-path proximity score computed
- [ ] Criterion 3 (Phenotypic Overlap): HPO files downloaded, OBO parsed, Resnik BMA + gene vector + pathway vector all computed
- [ ] Scorer.py integrates all gates and criteria into `{drug_id: full_evidence_ledger}` output

### Hour 7–9
- [ ] `tree_builder.py` outputs valid `tree.json` with all three layers populated
- [ ] Layer 3 expansion complete for top 3 drugs per Layer 2 node
- [ ] Streamlit app renders tree with pyvis, evidence ledger panel, filter controls, and CSV export

### Hour 9–11
- [ ] End-to-end test passes: YAML → pipeline → tree renders
- [ ] Positive control: rapamycin confirmed in top 3
- [ ] Negative control: GBA1 drugs confirmed absent from tree
- [ ] Generalizability demo: Costello Syndrome YAML swap confirmed to produce different tree
- [ ] Sensitivity analysis table generated (top-3 stable under ±10% weight perturbation)

### Hour 11–12
- [ ] All API response caches flushed and pre-loaded for demo (no live calls during presentation)
- [ ] Streamlit layout polished, caution annotations visible, validation panel finalized
- [ ] Sensitivity analysis output visible in Streamlit or printable table
- [ ] README written: problem, solution, setup instructions, YAML config documentation
- [ ] Submission materials assembled

### Fallback Decisions (Pre-committed, No Decision Time Required During Execution)
- [ ] DrugBank unavailable → DrugCentral throughout (no scoring degradation)
- [ ] CLUE API unavailable → Enrichr GO-BP pathway comparison for Criterion 1B (documented substitute)
- [ ] STRING API slow → reduce to 1-hop neighborhood only (still sufficient for ranking signal)
- [ ] Full Menche permutation analysis (1,000 permutations) → deferred to stretch goal in hours 9–11; simplified proximity score used if time does not allow
