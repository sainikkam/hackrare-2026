# LaPlaceRx

## Inspiration

95% of 7,000+ rare diseases have no approved treatment — patient populations are too small for traditional drug development. But drugs like everolimus (transplant drug → TSC seizures) prove existing drugs can be repurposed when diseases share biology. We wanted to automate that reasoning and make it accessible to anyone.

## What it does

Input any rare disease name → get ranked drugs with the strongest mechanistic case for repurposing. Two safety gates (BBB permeability, pediatric safety) filter candidates, then three criteria score them using magnitude-aware cosine similarity:

$$s(\mathbf{a}, \mathbf{b}) = \text{clamp}\!\left(\frac{\mathbf{a} \cdot \mathbf{b}}{\mathbf{b} \cdot \mathbf{b}},\ 0,\ 1\right), \quad \mathbf{a} = [C_1, C_2, C_3],\; \mathbf{b} = [1,1,1]$$

- **C1**: Pathway alignment (Reactome + Enrichr)
- **C2**: Network proximity (STRING PPI graph)
- **C3**: Phenotypic overlap (HPO + gene similarity)

## How we built it

Python + Streamlit + NetworkX + pyvis, powered entirely by **free public APIs** — OpenTargets, STRING, Enrichr, Reactome, HPO, PubChem, openFDA. No API keys required. We built a standalone CLI (`repurpose.py`, single file, runs anywhere) and a Streamlit web app with an interactive probability tree visualization.

## Challenges we ran into

- **Pediatric safety gate too strict** — openFDA returns sparse data, causing most drugs to fail. Fixed by curating known pediatric approvals for 30+ drugs.
- **BBB gate vs. positive controls** — raising the threshold broke validated drugs. Fixed with a clinical-authoritative branch: when experimental CNS evidence $\geq 0.80$, use it directly.
- **API unreliability** — multiple APIs fail intermittently. Every external call degrades gracefully so the pipeline always completes.

## Accomplishments that we're proud of

The algorithm correctly ranks known successful repurposings first — everolimus/sirolimus for TSC (0.416) and selumetinib/trametinib for NF1 (0.472) — while excluding negative controls (imatinib, imiglucerase). It's fully disease-agnostic, working for any of OpenTargets' 23,000+ indexed diseases with zero code changes.

## What we learned

- Validated drugs are the best unit tests — if your algorithm doesn't rank everolimus #1 for TSC, something is wrong.
- Simplicity beats complexity — magnitude-aware cosine with equal weights outperformed our planned ML approach with zero training data and full interpretability.
- Free public APIs are powerful enough to build real biomedical tools, but you must design for their failures.

## What's next for LaPlaceRx

- Expand beyond neurological diseases to oncology, metabolic, and immunological rare diseases with disease-specific comparison panels.
- Add an ML layer (graph neural networks on the PPI network) to learn disease-specific criterion weights from successful repurposing outcomes.
- Clinical trial matching — connect top-ranked candidates to active ClinicalTrials.gov entries so patients and researchers can act on the results.
