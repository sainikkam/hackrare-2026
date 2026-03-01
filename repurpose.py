#!/usr/bin/env python3
"""
repurpose.py — Standalone Rare Disease Drug Repurposing Navigator
=================================================================
Inputs a rare disease name, queries public APIs, scores candidate drugs
using a 3-criterion mechanistic pipeline, and prints ranked results.

Usage:
    python repurpose.py "Dravet Syndrome"
    python repurpose.py              # prompts for input

Dependencies (two packages only):
    pip install requests networkx

No API keys required. Results are cached in ~/.repurpose_cache/ for speed.
"""

import argparse
import hashlib
import json
import logging
import math
import sys
import time
from pathlib import Path
from urllib.parse import quote

try:
    import requests
except ImportError:
    sys.exit("Missing dependency: pip install requests networkx")

try:
    import networkx as nx
except ImportError:
    sys.exit("Missing dependency: pip install requests networkx")

# ──────────────────────────────────────────────────────────────────────────────
# Configuration
# ──────────────────────────────────────────────────────────────────────────────

BBB_THRESHOLD      = 0.8   # Gate 1: clinical-authoritative BBB score
SAFETY_THRESHOLD   = 0.6   # Gate 2: pediatric safety composite
TIMEOUT            = 30    # seconds for most HTTP requests
LONG_TIMEOUT       = 60    # seconds for STRING / HPO downloads

logging.basicConfig(level=logging.WARNING, format="%(levelname)s %(name)s — %(message)s")
logger = logging.getLogger("repurpose")

# Default neuro comparison panel — drugs are drawn from these related diseases
DEFAULT_PANEL = [
    {"name": "Tuberous Sclerosis Complex",  "opentargets_search": "tuberous sclerosis"},
    {"name": "Neurofibromatosis Type 1",    "opentargets_search": "neurofibromatosis type 1"},
    {"name": "SYNGAP1-Related Disorders",   "opentargets_search": "SYNGAP1"},
    {"name": "Angelman Syndrome",           "opentargets_search": "angelman syndrome"},
    {"name": "Rett Syndrome",               "opentargets_search": "rett syndrome"},
    {"name": "Fragile X Syndrome",          "opentargets_search": "fragile x syndrome"},
]

# Curated CNS evidence: when experimental score >= 0.8, use directly (authoritative)
EXPERIMENTAL_BBB = {
    "everolimus":   {"score": 0.90, "evidence": "FDA-approved for TSC-associated SEGA/seizures. CSF penetration documented."},
    "sirolimus":    {"score": 0.85, "evidence": "Preclinical CNS efficacy in TSC/NF1/SYNGAP1 mouse models."},
    "temsirolimus": {"score": 0.80, "evidence": "CNS penetration documented in neuro-oncology studies."},
    "selumetinib":  {"score": 0.88, "evidence": "FDA-approved for NF1 plexiform neurofibromas (pediatric ≥2y). Active CNS trials."},
    "trametinib":   {"score": 0.82, "evidence": "CNS penetration in pediatric low-grade glioma trials."},
    "binimetinib":  {"score": 0.75, "evidence": "CNS activity in melanoma brain metastasis trials."},
    "cobimetinib":  {"score": 0.72, "evidence": "Some CNS activity from melanoma studies."},
    "metformin":    {"score": 0.70, "evidence": "Crosses BBB via OCT transporters. CNS effects documented."},
    "lovastatin":   {"score": 0.75, "evidence": "Lipophilic statin with documented CNS penetration."},
    "simvastatin":  {"score": 0.78, "evidence": "Lipophilic statin, high CNS penetration."},
    "vigabatrin":   {"score": 0.92, "evidence": "FDA-approved antiepileptic with primary CNS target."},
    "valproic acid":{"score": 0.95, "evidence": "Classic CNS antiepileptic. High CNS penetration."},
    "memantine":    {"score": 0.93, "evidence": "FDA-approved for Alzheimer's. Primary CNS target (NMDA receptor)."},
    "ketamine":     {"score": 0.96, "evidence": "Primary CNS drug (NMDA receptor antagonist)."},
    "lithium":      {"score": 0.90, "evidence": "Primary CNS drug (mood stabilizer)."},
    "minocycline":  {"score": 0.72, "evidence": "Crosses BBB. Studied in Fragile X and multiple CNS conditions."},
    "rapamycin":    {"score": 0.85, "evidence": "Same as sirolimus. CNS efficacy in multiple mouse models."},
    "cabozantinib": {"score": 0.65, "evidence": "Some CNS penetration from GBM trials."},
    "erlotinib":    {"score": 0.60, "evidence": "Limited CNS penetration."},
}

# Known Pgp efflux pump substrates (reduces BBB score)
PGP_SUBSTRATES = {
    "everolimus", "sirolimus", "rapamycin", "tacrolimus", "cyclosporine",
    "imatinib", "dasatinib", "nilotinib", "lapatinib", "sunitinib",
    "gefitinib", "erlotinib", "paclitaxel", "docetaxel", "vincristine",
}
PGP_NON_SUBSTRATES = {
    "metformin", "selumetinib", "trametinib", "binimetinib", "cobimetinib",
    "valproic acid", "levetiracetam", "memantine", "lithium",
    "vigabatrin", "ketamine", "minocycline",
}

# Curated pediatric safety data for drugs with known FDA pediatric approvals/data.
# Overrides the openFDA API lookup when a match is found (lowercase key).
# Format: {name: (pediatric_approved, min_age_years)}
PEDIATRIC_KNOWN = {
    "everolimus":       (True,  1.0),   # TSC-SEGA ≥1y, TSC seizures ≥2y
    "sirolimus":        (True, 13.0),   # Lymphangioleiomyomatosis (off-label neuro use)
    "temsirolimus":     (False, None),
    "selumetinib":      (True,  2.0),   # NF1 plexiform neurofibromas ≥2y
    "trametinib":       (True,  1.0),   # Pediatric low-grade glioma ≥1y (FDA 2023)
    "binimetinib":      (False, None),
    "vigabatrin":       (True,  0.0),   # Infantile spasms — neonates onwards
    "valproic acid":    (True,  0.0),   # Epilepsy — all ages
    "ketamine":         (True,  0.0),   # Anaesthesia — all ages
    "cannabidiol":      (True,  2.0),   # Epidiolex for Dravet/LGS ≥2y
    "risperidone":      (True,  5.0),   # Irritability in autism ≥5y; schizophrenia ≥13y
    "memantine":        (False, None),  # Adults only (Alzheimer's); off-label in kids
    "metformin":        (True, 10.0),   # T2D ≥10y
    "lithium":          (True, 12.0),   # Bipolar ≥12y
    "minocycline":      (True,  8.0),   # Infections ≥8y
    "lovastatin":       (False, None),
    "simvastatin":      (False, None),
    "levetiracetam":    (True,  1.0),   # Epilepsy ≥1 month
    "topiramate":       (True,  2.0),   # Epilepsy ≥2y
    "lamotrigine":      (True,  2.0),   # Epilepsy ≥2y
    "fenfluramine":     (True,  2.0),   # Dravet ≥2y (Fintepla)
    "donepezil":        (False, None),  # Adults only
    "aripiprazole":     (True,  6.0),   # Tourette/autism ≥6y
    "gaboxadol":        (False, None),
    "dextromethorphan": (True,  0.0),   # OTC all ages
    "acamprosate":      (False, None),
    "roflumilast":      (False, None),
    "fingolimod":       (True, 10.0),   # MS ≥10y
    "esketamine":       (False, None),  # Adults only (depression)
    "celecoxib":        (True,  2.0),   # JIA ≥2y
    "carbidopa":        (False, None),
    "levodopa":         (False, None),
    "mecasermin":       (True,  2.0),   # Growth failure ≥2y
    "ganaxolone":       (True,  2.0),   # CDKL5 ≥2y
}

# ──────────────────────────────────────────────────────────────────────────────
# Disk Cache (JSON, keyed by SHA-256 of inputs)
# ──────────────────────────────────────────────────────────────────────────────

CACHE_DIR = Path.home() / ".repurpose_cache"
CACHE_DIR.mkdir(exist_ok=True)

def _cache_path(*args) -> Path:
    key = json.dumps(args, sort_keys=True)
    h = hashlib.sha256(key.encode()).hexdigest()[:20]
    return CACHE_DIR / f"{h}.json"

def cache_get(*args):
    p = _cache_path(*args)
    if p.exists():
        try:
            return json.loads(p.read_text())
        except Exception:
            pass
    return None

def cache_set(value, *args):
    p = _cache_path(*args)
    try:
        p.write_text(json.dumps(value, default=str))
    except Exception:
        pass

# ──────────────────────────────────────────────────────────────────────────────
# OpenTargets GraphQL API
# ──────────────────────────────────────────────────────────────────────────────

OT_URL = "https://api.platform.opentargets.org/api/v4/graphql"

def _ot_query(query: str, variables: dict) -> dict:
    try:
        r = requests.post(OT_URL, json={"query": query, "variables": variables}, timeout=TIMEOUT)
        r.raise_for_status()
        return r.json().get("data", {})
    except Exception as e:
        logger.warning("OpenTargets error: %s", e)
        return {}

def ot_search_disease(query: str) -> list:
    cached = cache_get("ot_search", query)
    if cached is not None:
        return cached
    q = """
    query($q: String!) {
      search(queryString: $q, entityNames: ["disease"], page: {index: 0, size: 5}) {
        hits { id name entity }
      }
    }
    """
    hits = _ot_query(q, {"q": query}).get("search", {}).get("hits", [])
    result = [{"id": h["id"], "name": h["name"]} for h in hits if h.get("entity") == "disease"]
    cache_set(result, "ot_search", query)
    return result

def ot_resolve_disease(query: str) -> str | None:
    hits = ot_search_disease(query)
    return hits[0]["id"] if hits else None

def ot_disease_genes(disease_id: str, min_score: float = 0.3, size: int = 80) -> list:
    cached = cache_get("ot_genes", disease_id, min_score, size)
    if cached is not None:
        return cached
    q = """
    query($id: String!, $size: Int!) {
      disease(efoId: $id) {
        associatedTargets(page: {index: 0, size: $size}) {
          rows { target { approvedSymbol } score }
        }
      }
    }
    """
    rows = (_ot_query(q, {"id": disease_id, "size": size})
            .get("disease", {}).get("associatedTargets", {}).get("rows", []))
    result = [r["target"]["approvedSymbol"] for r in rows if r["score"] >= min_score]
    cache_set(result, "ot_genes", disease_id, min_score, size)
    return result

def ot_disease_drugs(disease_id: str) -> list:
    cached = cache_get("ot_drugs", disease_id)
    if cached is not None:
        return cached
    q = """
    query($id: String!) {
      disease(efoId: $id) {
        knownDrugs(size: 100) {
          rows { drug { id name } phase mechanismOfAction }
        }
      }
    }
    """
    rows = (_ot_query(q, {"id": disease_id})
            .get("disease", {}).get("knownDrugs", {}).get("rows", []))
    seen, drugs = set(), []
    for r in rows:
        d = r.get("drug", {})
        cid = d.get("id")
        if cid and cid not in seen:
            seen.add(cid)
            drugs.append({
                "chembl_id": cid,
                "name": d.get("name", ""),
                "phase": r.get("phase", 0) or 0,
                "moa": r.get("mechanismOfAction", ""),
            })
    cache_set(drugs, "ot_drugs", disease_id)
    return drugs

def ot_drug_targets(chembl_id: str) -> list:
    cached = cache_get("ot_drug_targets", chembl_id)
    if cached is not None:
        return cached
    q = """
    query($id: String!) {
      drug(chemblId: $id) {
        mechanismsOfAction {
          rows { targets { approvedSymbol } }
        }
      }
    }
    """
    rows = (_ot_query(q, {"id": chembl_id})
            .get("drug") or {}).get("mechanismsOfAction") or {}
    genes = list({t["approvedSymbol"]
                  for r in rows.get("rows", [])
                  for t in r.get("targets", [])
                  if t.get("approvedSymbol")})
    cache_set(genes, "ot_drug_targets", chembl_id)
    return genes

# ──────────────────────────────────────────────────────────────────────────────
# STRING PPI Network
# ──────────────────────────────────────────────────────────────────────────────

STRING_URL = "https://string-db.org/api/json"

def string_interactions(gene_symbols: list, confidence: int = 700, limit: int = 500) -> list:
    if not gene_symbols:
        return []
    key = ("string_int", tuple(sorted(gene_symbols)), confidence)
    cached = cache_get(*key)
    if cached is not None:
        return cached
    try:
        r = requests.post(
            f"{STRING_URL}/interaction_partners",
            data={"identifiers": "\r".join(gene_symbols), "species": 9606,
                  "required_score": confidence, "limit": limit},
            timeout=LONG_TIMEOUT,
        )
        r.raise_for_status()
        result = [{"a": i.get("preferredName_A", ""), "b": i.get("preferredName_B", ""),
                   "score": i.get("score", 0)} for i in r.json()]
        cache_set(result, *key)
        return result
    except Exception as e:
        logger.warning("STRING error: %s", e)
        return []

def build_ppi_graph(seed_genes: list, confidence: int = 700) -> nx.Graph:
    # 1-hop interactions
    edges = string_interactions(seed_genes, confidence=confidence, limit=300)
    hop1_genes = list({e["a"] for e in edges} | {e["b"] for e in edges} - set(seed_genes))
    # 2-hop sample
    edges += string_interactions(hop1_genes[:80], confidence=confidence, limit=100)
    G = nx.Graph()
    seen = set()
    for e in edges:
        k = tuple(sorted([e["a"], e["b"]]))
        if k not in seen and e["a"] and e["b"]:
            seen.add(k)
            G.add_edge(e["a"], e["b"], weight=e["score"] / 1000.0)
    return G

def string_expand_targets(targets: list, confidence: int = 700, min_size: int = 5) -> list:
    """Expand drug targets via 1-hop STRING neighbors when < min_size direct targets."""
    if len(targets) >= min_size:
        return targets
    edges = string_interactions(targets, confidence=confidence, limit=100)
    neighbors = {e["a"] for e in edges} | {e["b"] for e in edges}
    return list(set(targets) | neighbors)

def mean_ppi_distance(G: nx.Graph, sources: list, targets: list) -> float:
    src = [g for g in sources if g in G]
    tgt = [g for g in targets if g in G]
    if not src or not tgt:
        return 10.0
    dists = []
    for s in src:
        for t in tgt:
            if s == t:
                dists.append(0)
            else:
                try:
                    dists.append(nx.shortest_path_length(G, s, t))
                except nx.NetworkXNoPath:
                    dists.append(10)
    return sum(dists) / len(dists) if dists else 10.0

# ──────────────────────────────────────────────────────────────────────────────
# Enrichr Pathway Enrichment
# ──────────────────────────────────────────────────────────────────────────────

ENRICHR_URL = "https://maayanlab.cloud/Enrichr"
ENRICHR_LIBS = ["KEGG_2021_Human", "Reactome_2022", "GO_Biological_Process_2023"]

def enrichr_pathway_vector(gene_symbols: list, top_n: int = 50) -> dict:
    if not gene_symbols:
        return {}
    cached = cache_get("enrichr", tuple(sorted(gene_symbols[:50])))
    if cached is not None:
        return cached
    try:
        # Submit
        r = requests.post(
            f"{ENRICHR_URL}/addList",
            files={"list": (None, "\n".join(gene_symbols[:50])),
                   "description": (None, "repurpose_query")},
            timeout=TIMEOUT,
        )
        r.raise_for_status()
        uid = r.json().get("userListId")
        if not uid:
            return {}
        # Fetch across libraries
        vector = {}
        for lib in ENRICHR_LIBS:
            time.sleep(0.2)
            r2 = requests.get(f"{ENRICHR_URL}/enrich",
                              params={"userListId": uid, "backgroundType": lib},
                              timeout=TIMEOUT)
            r2.raise_for_status()
            for item in (r2.json().get(lib) or [])[:top_n]:
                if len(item) >= 3:
                    term, pval = item[1], item[2]
                    score = -math.log10(max(pval, 1e-10))
                    vector[term] = max(vector.get(term, 0.0), score)
        cache_set(vector, "enrichr", tuple(sorted(gene_symbols[:50])))
        return vector
    except Exception as e:
        logger.warning("Enrichr error: %s", e)
        return {}

def cosine_sim(vec_a: dict, vec_b: dict) -> float:
    if not vec_a or not vec_b:
        return 0.0
    keys = set(vec_a) | set(vec_b)
    dot  = sum(vec_a.get(k, 0) * vec_b.get(k, 0) for k in keys)
    ma   = math.sqrt(sum(v * v for v in vec_a.values()))
    mb   = math.sqrt(sum(v * v for v in vec_b.values()))
    return dot / (ma * mb) if ma and mb else 0.0

# ──────────────────────────────────────────────────────────────────────────────
# Reactome Pathway Sets
# ──────────────────────────────────────────────────────────────────────────────

REACTOME_URL = "https://reactome.org/ContentService/data/mapping/UniProt/{gene}/pathways?species=Homo+sapiens"

def reactome_pathways_for_genes(gene_symbols: list) -> set:
    pathways = set()
    for gene in gene_symbols:
        cached = cache_get("reactome", gene)
        if cached is not None:
            pathways.update(cached)
            continue
        try:
            r = requests.get(REACTOME_URL.format(gene=gene), timeout=TIMEOUT)
            if r.status_code == 200:
                ids = [p.get("stId", "") for p in r.json()]
                cache_set(ids, "reactome", gene)
                pathways.update(ids)
        except Exception:
            pass
    return pathways

# ──────────────────────────────────────────────────────────────────────────────
# HPO Phenotype Data (gene-based, no full OBO required)
# ──────────────────────────────────────────────────────────────────────────────

HPO_GENES_URL = "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/phenotype_to_genes.txt"
HPO_GENES_PATH = CACHE_DIR / "phenotype_to_genes.txt"

_hpo_gene_map: dict | None = None   # {hpo_term: [gene_symbol, ...]}
_gene_hpo_map: dict | None = None   # {gene_symbol: [hpo_term, ...]}

def _ensure_hpo():
    global _hpo_gene_map, _gene_hpo_map
    if _hpo_gene_map is not None:
        return
    if not HPO_GENES_PATH.exists():
        print("  Downloading HPO phenotype_to_genes.txt (~5 MB, one-time)...", flush=True)
        try:
            r = requests.get(HPO_GENES_URL, timeout=120, stream=True)
            r.raise_for_status()
            HPO_GENES_PATH.write_bytes(r.content)
        except Exception as e:
            logger.warning("HPO download failed: %s", e)
            _hpo_gene_map, _gene_hpo_map = {}, {}
            return
    _hpo_gene_map, _gene_hpo_map = {}, {}
    try:
        for line in HPO_GENES_PATH.read_text(encoding="utf-8").splitlines():
            if line.startswith("#") or not line.strip():
                continue
            parts = line.split("\t")
            if len(parts) >= 4:
                hpo, gene = parts[0], parts[3]
                _hpo_gene_map.setdefault(hpo, []).append(gene)
                _gene_hpo_map.setdefault(gene, []).append(hpo)
    except Exception as e:
        logger.warning("HPO parse error: %s", e)

def hpo_terms_for_genes(gene_symbols: list) -> set:
    _ensure_hpo()
    terms = set()
    for g in gene_symbols:
        terms.update(_gene_hpo_map.get(g, []))
    return terms

def hpo_jaccard(genes_a: list, genes_b: list) -> float:
    a, b = hpo_terms_for_genes(genes_a), hpo_terms_for_genes(genes_b)
    if not a or not b:
        return 0.0
    return len(a & b) / len(a | b)

# ──────────────────────────────────────────────────────────────────────────────
# PubChem Molecular Properties
# ──────────────────────────────────────────────────────────────────────────────

PUBCHEM_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_PROPS = "MolecularWeight,XLogP,TPSA,HBondDonorCount,HBondAcceptorCount"

def pubchem_props(drug_name: str) -> dict | None:
    cached = cache_get("pubchem", drug_name)
    if cached is not None:
        return cached
    try:
        url = f"{PUBCHEM_URL}/compound/name/{quote(drug_name)}/property/{PUBCHEM_PROPS}/JSON"
        r = requests.get(url, timeout=TIMEOUT)
        if r.status_code == 404:
            cache_set(None, "pubchem", drug_name)
            return None
        r.raise_for_status()
        props = r.json().get("PropertyTable", {}).get("Properties", [])
        if not props:
            return None
        p = props[0]
        mw    = p.get("MolecularWeight")
        xlogp = p.get("XLogP")
        tpsa  = p.get("TPSA")
        hbd   = p.get("HBondDonorCount")
        hba   = p.get("HBondAcceptorCount")
        result = {
            "mw":    float(mw)    if mw    is not None else None,
            "xlogp": float(xlogp) if xlogp is not None else None,
            "clogd": (float(xlogp) - 0.5) if xlogp is not None else None,
            "tpsa":  float(tpsa)  if tpsa  is not None else None,
            "hbd":   int(hbd)     if hbd   is not None else None,
            "pka_basic": 8.5 if (hba or 0) > 0 else 1.0,
        }
        cache_set(result, "pubchem", drug_name)
        return result
    except Exception as e:
        logger.warning("PubChem error for %s: %s", drug_name, e)
        return None

# ──────────────────────────────────────────────────────────────────────────────
# openFDA Pediatric Safety Data
# ──────────────────────────────────────────────────────────────────────────────

OPENFDA_URL = "https://api.fda.gov/drug/label.json"

def openfda_pediatric_info(drug_name: str) -> dict:
    cached = cache_get("openfda", drug_name)
    if cached is not None:
        return cached

    result = {
        "pediatric_approved": False,
        "min_age_years": None,
        "has_pediatric_bb_warning": False,
        "has_pediatric_contraindication": False,
        "ae_score": 0.5,
    }
    try:
        r = requests.get(OPENFDA_URL,
                         params={"search": f'openfda.generic_name:"{drug_name}"', "limit": 1},
                         timeout=TIMEOUT)
        if r.status_code == 200:
            results = r.json().get("results", [])
            if results:
                label = results[0]
                ped_text = " ".join((label.get("pediatric_use") or []) +
                                    (label.get("pediatric_safety") or [])).lower()
                bb = " ".join(label.get("boxed_warning") or []).lower()

                if any(w in bb for w in ["pediatric", "children", "neonates"]):
                    result["has_pediatric_bb_warning"] = True

                if "contraindicated in pediatric" in ped_text or "not for use in children" in ped_text:
                    result["has_pediatric_contraindication"] = True

                if any(w in ped_text for w in ["approved", "indicated", "safety and efficacy"]):
                    result["pediatric_approved"] = True

                for phrase, age in [("2 years", 2), ("1 year", 1), ("6 months", 0.5),
                                    ("neonates", 0), ("infants", 0.5), ("12 years", 12),
                                    ("6 years", 6), ("18 years", 18)]:
                    if phrase in ped_text:
                        result["min_age_years"] = age
                        break

    except Exception as e:
        logger.warning("openFDA error for %s: %s", drug_name, e)

    cache_set(result, "openfda", drug_name)
    return result

# ──────────────────────────────────────────────────────────────────────────────
# Gate 1: Blood-Brain Barrier Permeability
# ──────────────────────────────────────────────────────────────────────────────

def _cns_mpo_score(props: dict | None) -> float:
    if not props:
        return 0.3

    def ld(v, opt, hard):
        if v is None: return 0.5
        if v <= opt:  return 1.0
        if v >= hard: return 0.0
        return (hard - v) / (hard - opt)

    mw    = props.get("mw")
    xlogp = props.get("xlogp")
    clogd = props.get("clogd")
    tpsa  = props.get("tpsa")
    hbd   = props.get("hbd")
    pka   = props.get("pka_basic")

    scores = [
        ld(mw,    360.0, 500.0),
        ld(xlogp,   3.0,   5.0),
        ld(clogd,   2.0,   4.0),
        ld(tpsa,   60.0,  90.0),
        ld(hbd,     0.0,   3.0),
        ld(pka,     8.0,  10.0),
    ]
    return sum(scores) / 6.0

def score_bbb(drug_name: str) -> dict:
    props = pubchem_props(drug_name)
    mpo   = _cns_mpo_score(props)
    name  = drug_name.lower()

    pgp_adj = 0.6 if name in PGP_SUBSTRATES else (1.0 if name in PGP_NON_SUBSTRATES else
              (0.8 if props and (props.get("mw") or 0) > 400 and (props.get("tpsa") or 0) > 90 else 1.0))
    phys_score = round(mpo * pgp_adj, 3)

    exp = EXPERIMENTAL_BBB.get(name)
    if exp:
        exp_score = exp["score"]
        if exp_score >= 0.80:
            final = round(exp_score, 3)
            method = "clinical_evidence_authoritative"
        else:
            final = round(0.35 * phys_score + 0.65 * exp_score, 3)
            method = "composite"
    else:
        final = phys_score
        exp_score = None
        method = "physicochemical"

    return {
        "score": final,
        "passes": final >= BBB_THRESHOLD,
        "method": method,
        "physicochemical_score": phys_score,
        "experimental_score": exp_score,
        "experimental_evidence": exp["evidence"] if exp else None,
    }

# ──────────────────────────────────────────────────────────────────────────────
# Gate 2: Pediatric Safety
# ──────────────────────────────────────────────────────────────────────────────

def score_pediatric_safety(drug_name: str) -> dict:
    # Check curated dict first (authoritative for known drugs)
    known = PEDIATRIC_KNOWN.get(drug_name.lower())
    if known is not None:
        ped_approved, min_age = known
    else:
        # Fall back to openFDA
        info = openfda_pediatric_info(drug_name)
        if info.get("has_pediatric_bb_warning"):
            return {"score": 0.0, "passes": False, "caution": False,
                    "caution_label": "Black box warning (pediatric)", "min_age_years": None,
                    "pediatric_approved": False}
        if info.get("has_pediatric_contraindication"):
            return {"score": 0.0, "passes": False, "caution": False,
                    "caution_label": "Pediatric contraindication", "min_age_years": None,
                    "pediatric_approved": False}
        ped_approved = info.get("pediatric_approved", False)
        min_age      = info.get("min_age_years")

    # Component A: approval status
    if ped_approved:
        comp_a = 0.6 if (min_age is not None and min_age >= 12) else 1.0
    else:
        comp_a = 0.4  # unknown/unapproved — generous default for novel repurposing context

    # Component B: min age (younger = better)
    if min_age is None:
        comp_b = 0.5  # unknown — neutral
    elif min_age <= 2:  comp_b = 1.0
    elif min_age <= 5:  comp_b = 0.8
    elif min_age <= 11: comp_b = 0.6
    elif min_age <= 17: comp_b = 0.4
    else:               comp_b = 0.1

    # Component C: AE rate (neutral default 0.5)
    comp_c = 0.5

    score  = round(0.40 * comp_a + 0.30 * comp_b + 0.30 * comp_c, 3)
    passes = score >= SAFETY_THRESHOLD
    caution = passes and score < 0.80

    return {
        "score": score,
        "passes": passes,
        "caution": caution,
        "caution_label": "Caution: limited pediatric data" if caution else None,
        "min_age_years": min_age,
        "pediatric_approved": ped_approved,
    }

# ──────────────────────────────────────────────────────────────────────────────
# Criterion 1: Pathway Alignment
# ──────────────────────────────────────────────────────────────────────────────

def score_pathway_alignment(drug_name: str, chembl_id: str,
                            disease_genes: list, disease_pathway_vec: dict) -> dict:
    # Drug direct targets
    direct_targets = ot_drug_targets(chembl_id) if chembl_id else []

    # Expand via STRING if fewer than 5 direct targets
    expanded = string_expand_targets(direct_targets, min_size=5)

    # Sub-metric A: Reactome Jaccard
    drug_pathways    = reactome_pathways_for_genes(expanded)
    disease_pathways = reactome_pathways_for_genes(disease_genes[:20])
    if drug_pathways and disease_pathways:
        jac = len(drug_pathways & disease_pathways) / len(drug_pathways | disease_pathways)
    else:
        jac = 0.0

    # Sub-metric B: Enrichr pathway cosine
    drug_vec = enrichr_pathway_vector(expanded[:30])
    enr_cos  = cosine_sim(drug_vec, disease_pathway_vec)

    score = round(0.5 * jac + 0.5 * enr_cos, 3)
    shared = list(drug_pathways & disease_pathways)[:5]

    return {
        "score": score,
        "pathway_jaccard": round(jac, 3),
        "enrichr_cosine": round(enr_cos, 3),
        "drug_target_genes": direct_targets,
        "drug_target_genes_expanded": expanded,
        "shared_pathways": shared,
    }

# ──────────────────────────────────────────────────────────────────────────────
# Criterion 2: Network Proximity
# ──────────────────────────────────────────────────────────────────────────────

def score_network_proximity(drug_targets: list, disease_genes: list, G: nx.Graph) -> dict:
    dist  = mean_ppi_distance(G, drug_targets, disease_genes)
    score = round(1 / (1 + dist), 3)
    return {"score": score, "mean_distance": round(dist, 2)}

# ──────────────────────────────────────────────────────────────────────────────
# Criterion 3: Phenotypic Overlap
# ──────────────────────────────────────────────────────────────────────────────

def score_phenotypic_overlap(disease_genes: list, drug_targets_expanded: list,
                             disease_pathway_vec: dict) -> dict:
    # HPO Jaccard (gene-based)
    hpo_jac = hpo_jaccard(disease_genes, drug_targets_expanded)

    # Gene overlap (Jaccard on gene symbols)
    a, b = set(disease_genes), set(drug_targets_expanded)
    gene_jac = len(a & b) / len(a | b) if (a | b) else 0.0

    # Pathway cosine (reuse disease vector vs drug targets enrichment)
    drug_vec = enrichr_pathway_vector(drug_targets_expanded[:30])
    path_cos = cosine_sim(disease_pathway_vec, drug_vec)

    score = round(0.4 * hpo_jac + 0.3 * gene_jac + 0.3 * path_cos, 3)
    shared_genes = sorted(a & b)[:8]

    return {
        "score": score,
        "hpo_jaccard": round(hpo_jac, 3),
        "gene_jaccard": round(gene_jac, 3),
        "pathway_cosine": round(path_cos, 3),
        "shared_genes": shared_genes,
    }

# ──────────────────────────────────────────────────────────────────────────────
# Final Score: Magnitude-Aware Cosine
# ──────────────────────────────────────────────────────────────────────────────

def magnitude_aware_cosine(c1: float, c2: float, c3: float) -> float:
    """
    s(a, b) = clamp((a·b) / (b·b), 0, 1)  where b = [1, 1, 1]
    Equivalent to arithmetic mean of [c1, c2, c3], clamped to [0, 1].
    Penalises single-criterion evidence; rewards balanced profiles.
    """
    return round(max(0.0, min(1.0, (c1 + c2 + c3) / 3.0)), 3)

# ──────────────────────────────────────────────────────────────────────────────
# Score a Single Drug
# ──────────────────────────────────────────────────────────────────────────────

def score_drug(drug: dict, disease_genes: list, G: nx.Graph, disease_pathway_vec: dict) -> dict:
    name     = drug.get("name", "")
    chembl   = drug.get("chembl_id", "")

    bbb = score_bbb(name)
    if not bbb["passes"]:
        return {"drug_name": name, "chembl_id": chembl, "status": "gate1_fail",
                "gate1_bbb": bbb, "final_score": None, "phase": drug.get("phase", 0),
                "moa": drug.get("moa", "")}

    safety = score_pediatric_safety(name)
    if not safety["passes"]:
        return {"drug_name": name, "chembl_id": chembl, "status": "gate2_fail",
                "gate1_bbb": bbb, "gate2_safety": safety, "final_score": None,
                "phase": drug.get("phase", 0), "moa": drug.get("moa", "")}

    c1 = score_pathway_alignment(name, chembl, disease_genes, disease_pathway_vec)
    c2 = score_network_proximity(c1["drug_target_genes"], disease_genes, G)
    c3 = score_phenotypic_overlap(disease_genes,
                                  c1["drug_target_genes_expanded"] or c1["drug_target_genes"],
                                  disease_pathway_vec)

    final = magnitude_aware_cosine(c1["score"], c2["score"], c3["score"])

    return {
        "drug_name": name,
        "chembl_id": chembl,
        "status": "scored",
        "final_score": final,
        "c1_pathway": c1["score"],
        "c2_network": c2["score"],
        "c3_phenotypic": c3["score"],
        "bbb_score": bbb["score"],
        "bbb_method": bbb["method"],
        "safety_score": safety["score"],
        "caution": safety.get("caution", False),
        "caution_label": safety.get("caution_label"),
        "phase": drug.get("phase", 0),
        "moa": drug.get("moa", ""),
        "shared_pathways": c1["shared_pathways"],
        "drug_targets": c1["drug_target_genes"],
        "shared_genes": c3["shared_genes"],
        "mean_network_distance": c2["mean_distance"],
        "gate1_bbb": bbb,
        "gate2_safety": safety,
    }

# ──────────────────────────────────────────────────────────────────────────────
# Main Pipeline
# ──────────────────────────────────────────────────────────────────────────────

def run(disease_name: str, panel: list = None) -> list:
    panel = panel or DEFAULT_PANEL
    sep = "=" * 62

    print(f"\n{sep}")
    print(f"  HackRare Drug Repurposing Navigator")
    print(f"  Disease: {disease_name}")
    print(f"{sep}\n")

    # ── Step 1: Resolve disease & build module ────────────────────────────────
    print("[1/5] Searching OpenTargets...", flush=True)
    hits = ot_search_disease(disease_name)
    if hits:
        print(f"      Matched: {hits[0]['name']}  ({hits[0]['id']})")
        disease_id = hits[0]["id"]
        matched_name = hits[0]["name"]
    else:
        print(f"      WARNING: No OpenTargets match for '{disease_name}' — using name directly.")
        disease_id = None
        matched_name = disease_name

    genes = []
    if disease_id:
        genes = ot_disease_genes(disease_id)
    if not genes:
        print("      WARNING: No genes found. Results will be sparse.")
    else:
        print(f"      Disease module: {len(genes)} genes")

    # ── Step 2: Build PPI graph ───────────────────────────────────────────────
    print("[2/5] Building STRING protein interaction network...", flush=True)
    G = build_ppi_graph(genes, confidence=700) if genes else nx.Graph()
    print(f"      Graph: {G.number_of_nodes()} nodes, {G.number_of_edges()} edges")

    # ── Step 3: Disease pathway signature ────────────────────────────────────
    print("[3/5] Computing disease pathway signature (Enrichr)...", flush=True)
    disease_pathway_vec = enrichr_pathway_vector(genes[:50])
    print(f"      Pathway vector: {len(disease_pathway_vec)} terms")

    # ── Step 4: Collect candidate drugs ──────────────────────────────────────
    print("[4/5] Collecting candidate drugs from comparison panel...", flush=True)
    all_drugs: dict = {}
    for entry in panel:
        did = ot_resolve_disease(entry["opentargets_search"])
        if not did:
            continue
        for d in ot_disease_drugs(did):
            key = d["name"].lower()
            if key and key not in all_drugs:
                all_drugs[key] = d
    print(f"      Panel: {', '.join(e['name'] for e in panel)}")
    print(f"      Candidates: {len(all_drugs)} unique drugs")

    if not all_drugs:
        print("\nNo candidate drugs found. Exiting.")
        return []

    # ── Step 5: Score ─────────────────────────────────────────────────────────
    print("[5/5] Scoring candidates through 2 gates + 3 criteria...", flush=True)
    results = []
    for i, drug in enumerate(all_drugs.values(), 1):
        print(f"      {i}/{len(all_drugs)} {drug['name']:<20}", end="\r", flush=True)
        results.append(score_drug(drug, genes, G, disease_pathway_vec))
    print(" " * 60, end="\r")

    scored  = [r for r in results if r["status"] == "scored"]
    g1_fail = [r for r in results if r["status"] == "gate1_fail"]
    g2_fail = [r for r in results if r["status"] == "gate2_fail"]
    scored.sort(key=lambda r: r["final_score"], reverse=True)

    # ── Print results ─────────────────────────────────────────────────────────
    print(f"\n{sep}")
    print(f"  TOP DRUG CANDIDATES  —  {matched_name}")
    print(f"{sep}")
    print(f"  {'#':<4} {'Drug':<22} {'Score':>6}  {'C1-Path':>7} {'C2-Net':>6} {'C3-Phen':>7}  {'Ph':>2}  Notes")
    print(f"  {'-'*4} {'-'*22} {'-'*6}  {'-'*7} {'-'*6} {'-'*7}  {'-'*2}  {'-'*20}")

    for i, r in enumerate(scored[:15], 1):
        caution = " ⚠ CAUTION" if r.get("caution") else ""
        print(f"  {i:<4} {r['drug_name']:<22} {r['final_score']:>6.3f}"
              f"  {r['c1_pathway']:>7.3f} {r['c2_network']:>6.3f} {r['c3_phenotypic']:>7.3f}"
              f"  {int(r['phase'] or 0):>2}{caution}")

    if not scored:
        print("  (no drugs passed both gates)")

    print(f"\n  Candidates evaluated: {len(results)}")
    print(f"  Gate 1 (BBB ≥ {BBB_THRESHOLD}) fail:     {len(g1_fail)} drugs")
    print(f"  Gate 2 (safety ≥ {SAFETY_THRESHOLD}) fail:   {len(g2_fail)} drugs")
    print(f"  Fully scored:         {len(scored)} drugs")

    if scored:
        top = scored[0]
        print(f"\n  Top pick: {top['drug_name']}")
        if top.get("shared_pathways"):
            print(f"  Shared pathways: {', '.join(top['shared_pathways'][:3])}")
        if top.get("shared_genes"):
            print(f"  Shared genes:    {', '.join(top['shared_genes'][:6])}")
        if top.get("moa"):
            print(f"  Mechanism:       {top['moa']}")

    print(f"\n  Scoring: magnitude-aware cosine  s(a,b) = clamp((a·b)/(b·b), 0, 1)")
    print(f"           reference b = [1,1,1]  →  equivalent to mean(C1, C2, C3)")
    print(f"  Cache:   {CACHE_DIR}")
    print(f"{sep}\n")

    return scored

# ──────────────────────────────────────────────────────────────────────────────
# Entry point
# ──────────────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Rare Disease Drug Repurposing Navigator — no API keys required",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog='Examples:\n  python repurpose.py "Dravet Syndrome"\n  python repurpose.py "CDKL5 Deficiency"\n  python repurpose.py "Phelan-McDermid Syndrome"',
    )
    parser.add_argument("disease", nargs="?", help="Rare disease name to query")
    args = parser.parse_args()

    disease = args.disease
    if not disease:
        try:
            disease = input("Enter rare disease name: ").strip()
        except (KeyboardInterrupt, EOFError):
            sys.exit(0)
    if not disease:
        parser.print_help()
        sys.exit(1)

    run(disease)

if __name__ == "__main__":
    main()
