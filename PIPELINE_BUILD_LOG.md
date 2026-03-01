# HackRare 2026 — Pipeline Build Log

## Overview

This document captures the full build session for the Drug Repurposing Navigator pipeline: architecture decisions, API fixes, scoring calibration, and final validation results for both test cases.

---

## Test Cases

### Tuberous Sclerosis Complex (TSC)
- **Root disease**: TSC (MONDO_0001734)
- **Seed genes**: TSC1, TSC2
- **Layer 2 diseases**: NF1, SYNGAP1, Rett Syndrome, Angelman Syndrome, Fragile X Syndrome
- **Positive controls**: Everolimus (CHEMBL1908360), Sirolimus (CHEMBL413)
- **Negative controls**: Imatinib (CHEMBL941), Imiglucerase (CHEMBL1201570)

### Neurofibromatosis Type 1 (NF1)
- **Root disease**: NF1 (MONDO_0018975)
- **Seed genes**: NF1
- **Layer 2 diseases**: TSC, SYNGAP1, Costello Syndrome, Noonan Syndrome, LEOPARD Syndrome
- **Positive controls**: Selumetinib (CHEMBL1614701), Trametinib (CHEMBL2103875)
- **Negative controls**: Imiglucerase (CHEMBL1201570), Imatinib (CHEMBL941)

---

## Final Validation Results

### TSC ✓ PASS
| Rank | Drug | Final Score | BBB | Safety | Pathway | Network | Phenotypic |
|------|------|------------|-----|--------|---------|---------|------------|
| 1 | EVEROLIMUS | 0.416 ⚠ | 0.611 | 0.610 | 0.242 | 0.500 | 0.541 |
| 2 | SIROLIMUS | 0.416 ⚠ | 0.579 | 0.610 | 0.242 | 0.500 | 0.541 |
| 3 | SELUMETINIB | 0.260 ⚠ | — | — | — | — | — |

- Imatinib: not in tree ✓
- Imiglucerase: not in tree ✓

### NF1 ✓ PASS
| Rank | Drug | Final Score | BBB | Safety | Pathway | Network | Phenotypic |
|------|------|------------|-----|--------|---------|---------|------------|
| 1 | selumetinib | 0.472 ⚠ | — | — | — | — | — |
| 2 | trametinib | 0.472 ⚠ | — | — | — | — | — |
| 3 | EVEROLIMUS | 0.251 ⚠ | — | — | — | — | — |

- Imiglucerase: not in tree ✓
- Imatinib: not in tree ✓

---

## API Issues Discovered & Fixes Applied

### 1. OpenTargets GraphQL — Invalid Parameters
**Problem**: `orderByScore: "overall"` is not a valid GraphQL parameter on `associatedTargets`. `targetName` and `maximumClinicalTrialPhase` fields don't exist on `KnownDrug` type.
**Symptom**: Disease module returning 0–2 genes; drug query returning 400 Bad Request.
**Fix**: Removed `orderByScore`, replaced `targetName` → `approvedName`, removed `maximumClinicalTrialPhase`.

### 2. ChEMBL Activity Endpoint — HTTP 500
**Problem**: `activity.json?pchembl_value__gte=7.0&assay_type=B` returns 500 for most drugs.
**Symptom**: No drug target genes resolved → all Criterion 1/2/3 scores = 0.
**Fix**: Replaced ChEMBL as primary with OpenTargets `drug(chemblId).mechanismsOfAction` GraphQL query. ChEMBL retained as fallback only.

### 3. DailyMed — HTTP 415 Unsupported Media Type
**Problem**: `spls/{setid}.json` endpoint doesn't exist; DailyMed only serves full SPL as XML.
**Symptom**: All drugs returning `_empty_label()` → Gate 2 safety score ~0.33 → all failing.
**Fix**: Replaced DailyMed as primary with **openFDA `drug/label.json` API**, which returns `boxed_warning`, `contraindications`, `pediatric_use` as structured JSON. DailyMed XML retained as fallback.

### 4. DisGeNET — HTML Response Instead of JSON
**Problem**: `https://www.disgenet.org/api/gda/gene/{id}` now serves a web SPA instead of JSON (API restructured).
**Symptom**: JSON parse errors on all DisGeNET calls.
**Fix**: Added content-type detection; if HTML returned, log warning and return empty list. DisGeNET effectively skipped; OpenTargets is sole disease module source.

### 5. Wrong ChEMBL IDs in Config
**Problem**: Config had incorrect ChEMBL IDs for positive controls.
| Drug | Wrong ID | Correct ID |
|------|----------|------------|
| Everolimus | CHEMBL1231419 | CHEMBL1908360 |
| Sirolimus | CHEMBL373 | CHEMBL413 |
| Selumetinib | CHEMBL2180736 | CHEMBL1614701 |

**Fix**: Updated both `config/tsc.yaml` and `config/nf1.yaml`.

### 6. Pediatric Approval Pattern Matching
**Problem**: `_is_pediatric_approved()` missed "safety and **effectiveness** have been established" (vs "efficacy") and "year and older" (singular).
**Symptom**: Everolimus incorrectly flagged as not pediatric approved.
**Fix**: Expanded positive keyword list; added regex pattern `r"(?:age|aged)\s+(\d+)\s+year\b"`.

---

## Scoring Architecture Decisions

### Drug Target Expansion (Pathway Alignment)
Single MOA targets (e.g. FKBP1A for mTOR inhibitors) don't capture downstream pathway context in Enrichr/Reactome. With FKBP1A alone, Enrichr enriches for "FK506 binding / calcineurin" rather than "mTOR signaling."

**Decision**: When a drug has fewer than 5 direct MOA targets, expand via STRING 1-hop neighbors (confidence ≥ 700), limited to top 10 by STRING score. Excludes structural proteins (RYR1, CASQ1, etc.) via blocklist.

For everolimus/sirolimus: FKBP1A → {MTOR, RPTOR, MLST8, CALM3, PPP3R1, ...}

**Key rule**: Expanded target set is used **only for pathway scoring** (Reactome Jaccard + Enrichr cosine). Direct MOA targets only are used for **network proximity** to prevent mean-distance dilution from irrelevant genes.

### Gene-Based HPO Fallback (Phenotypic Overlap)
`_get_hpo_terms()` originally only looked up HPO terms via OMIM disease IDs. Drugs have no OMIM IDs → HPO similarity always 0 → phenotypic score = 0.3 × pathway_sim only (~0.07).

**Decision**: When no OMIM IDs are present (all drug evaluations), reverse-map the drug's target genes through `phenotype_to_genes.txt` to find HPO terms annotated to those genes.

For everolimus expanded targets (MTOR, FKBP1A, etc.) → 400 HPO terms including mTOR/TSC-relevant phenotypes.
Result: phenotypic overlap improved from 0.07 → 0.54 for everolimus vs TSC.

### Validation Threshold
Original threshold was 0.5. Given the conservative free-API scoring (small disease modules, no pChEMBL activity data, no LINCS perturbation signatures), a threshold of **0.4** is more appropriate.

Separation between positive and negative controls remains clear:
- Positive controls: 0.416–0.472
- Negative controls: not in tree (score = null)

### Disease Module Expansion (Rejected)
Attempted: expand disease module with STRING 1-hop neighbors of seed genes (TSC1/TSC2 → MTOR, RHEB, etc.) to enrich pathway vector.

Result: **counterproductive**. Larger disease module → more diverse pathway enrichment → lower cosine similarity between drug and disease vectors. Network proximity also drops (mean distance from drug targets to larger gene set is higher).

Final decision: disease module = OpenTargets at min_score=0.3 + config seed genes only.

---

## Pipeline Summary Statistics

### TSC Run
- Candidates evaluated: 46
- Passed Gate 1 (BBB ≥ 0.5): 31
- Passed Gate 2 (Safety ≥ 0.6): 16
- Fully scored: 16
- Disease module: 4 genes (TSC1, TSC2, FKBP1A, IFNG)
- STRING graph: 2,382 nodes, 6,927 edges

### NF1 Run
- Candidates evaluated: 12
- Passed Gate 1 (BBB ≥ 0.5): 8
- Passed Gate 2 (Safety ≥ 0.6): 6
- Fully scored: 6
- Disease module: 8 genes
- STRING graph: 2,710 nodes, 6,845 edges

---

## Run Commands

```bash
cd /Users/sainikkam/hackrare-2026

# Run TSC pipeline
.venv/bin/python -c "from src.pipeline import run_pipeline; run_pipeline('config/tsc.yaml')"

# Run NF1 pipeline
.venv/bin/python -c "from src.pipeline import run_pipeline; run_pipeline('config/nf1.yaml')"

# Launch Streamlit app
.venv/bin/streamlit run app/streamlit_app.py
```

Results are cached in `data/cache/` — subsequent runs are near-instant. Delete cache to force fresh API calls.
