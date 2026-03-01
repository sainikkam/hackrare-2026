"""
HPO (Human Phenotype Ontology) loader.
Downloads hp.obo and phenotype_to_genes.txt on first run (~30 MB, ~2 min).
Parses ontology into an IC-annotated DAG for Resnik BMA semantic similarity.
"""

import math
import logging
import requests
import networkx as nx
from pathlib import Path
from .cache import file_cache

logger = logging.getLogger(__name__)
HPO_DIR = Path(__file__).resolve().parents[2] / "data" / "hpo"
OBO_PATH = HPO_DIR / "hp.obo"
GENES_PATH = HPO_DIR / "phenotype_to_genes.txt"
ANNOTATIONS_PATH = HPO_DIR / "phenotype.hpoa"

OBO_URL = "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/hp.obo"
GENES_URL = "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/phenotype_to_genes.txt"
ANNOTATIONS_URL = "https://github.com/obophenotype/human-phenotype-ontology/releases/latest/download/phenotype.hpoa"


def ensure_hpo_files():
    """Download HPO files if not already present."""
    HPO_DIR.mkdir(parents=True, exist_ok=True)

    for path, url in [(OBO_PATH, OBO_URL), (GENES_PATH, GENES_URL)]:
        if not path.exists():
            logger.info("Downloading %s ...", path.name)
            resp = requests.get(url, timeout=120, stream=True)
            resp.raise_for_status()
            path.write_bytes(resp.content)
            logger.info("Downloaded %s (%.1f MB)", path.name, path.stat().st_size / 1e6)


def parse_obo(path: Path) -> dict[str, dict]:
    """
    Parse hp.obo into dict: {term_id: {id, name, parents: [ids], obsolete: bool}}
    """
    terms = {}
    current = None

    with open(path, encoding="utf-8") as f:
        for line in f:
            line = line.rstrip()
            if line == "[Term]":
                current = {"parents": [], "obsolete": False}
            elif line == "[Typedef]" or line == "":
                if current and "id" in current:
                    terms[current["id"]] = current
                current = None
            elif current is not None:
                if line.startswith("id: "):
                    current["id"] = line[4:]
                elif line.startswith("name: "):
                    current["name"] = line[6:]
                elif line.startswith("is_a: "):
                    parent_id = line[6:].split(" ! ")[0].strip()
                    current["parents"].append(parent_id)
                elif line == "is_obsolete: true":
                    current["obsolete"] = True

    if current and "id" in current:
        terms[current["id"]] = current

    # Remove obsolete terms
    return {k: v for k, v in terms.items() if not v.get("obsolete")}


def build_dag(terms: dict[str, dict]) -> nx.DiGraph:
    """Build a directed acyclic graph (child -> parent) from HPO terms."""
    G = nx.DiGraph()
    for term_id, term in terms.items():
        G.add_node(term_id, name=term.get("name", ""))
        for parent_id in term.get("parents", []):
            if parent_id in terms:
                G.add_edge(term_id, parent_id)  # child -> parent
    return G


def parse_phenotype_to_genes(path: Path) -> dict[str, list[str]]:
    """
    Parse phenotype_to_genes.txt.
    Returns: {hpo_term_id: [gene_symbol, ...]}
    """
    mapping = {}
    with open(path, encoding="utf-8") as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) >= 4:
                hpo_id = parts[0]
                gene_symbol = parts[3]
                mapping.setdefault(hpo_id, []).append(gene_symbol)
    return mapping


@file_cache("hpo")
def get_disease_hpo_terms(omim_id: str) -> list[str]:
    """
    Get HPO terms annotated to a disease by OMIM ID from phenotype.hpoa.
    Falls back to empty list if file not available.
    """
    ann_path = ANNOTATIONS_PATH
    if not ann_path.exists():
        # Try to download
        try:
            resp = requests.get(ANNOTATIONS_URL, timeout=60, stream=True)
            resp.raise_for_status()
            ann_path.write_bytes(resp.content)
        except Exception as e:
            logger.error("Could not download phenotype.hpoa: %s", e)
            return []

    terms = []
    omim_prefix = f"OMIM:{omim_id}"
    try:
        with open(ann_path, encoding="utf-8") as f:
            for line in f:
                if line.startswith("#"):
                    continue
                parts = line.strip().split("\t")
                if len(parts) >= 4 and parts[0] == omim_prefix:
                    hpo_term = parts[3]
                    if hpo_term.startswith("HP:"):
                        terms.append(hpo_term)
    except Exception as e:
        logger.error("HPO annotation parse error for OMIM:%s: %s", omim_id, e)

    return list(set(terms))


class HPOSimilarity:
    """Compute Resnik IC-based semantic similarity between HPO term sets."""

    def __init__(self):
        ensure_hpo_files()
        logger.info("Loading HPO ontology...")
        self.terms = parse_obo(OBO_PATH)
        self.dag = build_dag(self.terms)
        self.gene_map = parse_phenotype_to_genes(GENES_PATH)
        self.ic = self._compute_ic()
        logger.info("HPO loaded: %d terms, %d term-gene mappings", len(self.terms), len(self.gene_map))

    def _compute_ic(self) -> dict[str, float]:
        """Compute information content for each term via annotation frequency propagation."""
        # Count direct annotations per term from gene_map
        freq = {term_id: 0 for term_id in self.terms}
        for term_id in self.gene_map:
            if term_id in freq:
                freq[term_id] += len(self.gene_map[term_id])

        # Propagate counts to ancestors
        propagated = dict(freq)
        for term_id in nx.topological_sort(self.dag):
            for parent in self.dag.successors(term_id):  # parent in child->parent graph
                propagated[parent] = propagated.get(parent, 0) + propagated.get(term_id, 0)

        total = max(propagated.values()) if propagated else 1
        ic = {}
        for term_id, count in propagated.items():
            if count > 0:
                ic[term_id] = -math.log2(count / total)
            else:
                ic[term_id] = 0.0
        return ic

    def _get_ancestors(self, term_id: str) -> set[str]:
        """Get all ancestors of a term (including itself)."""
        if term_id not in self.dag:
            return {term_id}
        return set(nx.ancestors(self.dag, term_id)) | {term_id}

    def resnik_sim(self, term_a: str, term_b: str) -> float:
        """Resnik similarity = IC of Most Informative Common Ancestor (MICA)."""
        anc_a = self._get_ancestors(term_a)
        anc_b = self._get_ancestors(term_b)
        common = anc_a & anc_b
        if not common:
            return 0.0
        return max(self.ic.get(t, 0.0) for t in common)

    def bma_similarity(self, terms_a: list[str], terms_b: list[str]) -> float:
        """Best-Match Average semantic similarity between two HPO term sets. Returns [0,1]."""
        terms_a = [t for t in terms_a if t in self.dag or t in self.terms]
        terms_b = [t for t in terms_b if t in self.dag or t in self.terms]
        if not terms_a or not terms_b:
            return 0.0

        scores_ab = [max((self.resnik_sim(a, b) for b in terms_b), default=0.0) for a in terms_a]
        scores_ba = [max((self.resnik_sim(b, a) for a in terms_a), default=0.0) for b in terms_b]

        bma = (sum(scores_ab) / len(scores_ab) + sum(scores_ba) / len(scores_ba)) / 2
        max_ic = max(self.ic.values()) if self.ic else 1.0
        return min(bma / max_ic, 1.0) if max_ic > 0 else 0.0

    def get_genes_for_terms(self, hpo_terms: list[str]) -> dict[str, float]:
        """
        Map HPO terms to genes with IC-weighted scores.
        Returns: {gene_symbol: weighted_score}
        """
        gene_scores = {}
        for term_id in hpo_terms:
            term_ic = self.ic.get(term_id, 0.0)
            genes = self.gene_map.get(term_id, [])
            for gene in genes:
                gene_scores[gene] = gene_scores.get(gene, 0.0) + term_ic
        return gene_scores

    def top_shared_hpo_terms(self, terms_a: list[str], terms_b: list[str], top_k: int = 5) -> list[dict]:
        """Return top-k HPO terms shared between two diseases, ranked by IC."""
        shared = set(terms_a) & set(terms_b)
        ranked = sorted(shared, key=lambda t: self.ic.get(t, 0.0), reverse=True)[:top_k]
        return [{"term_id": t, "name": self.terms.get(t, {}).get("name", t), "ic": self.ic.get(t, 0.0)} for t in ranked]
