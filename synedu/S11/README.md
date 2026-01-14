# S10 · Automated Rule Library Construction via **SynKit**

<div class="alert alert-block alert-info">
<b>Welcome to SynEdu.</b><br>
This talktorial concludes the <b>SynEdu</b> series by demonstrating how the concepts developed in
<b>S01–S09</b> come together in a realistic, end-to-end workflow: constructing a <b>reaction rule library</b>
from experimental data. We use the <b>USPTO-50K</b> benchmark and rely primarily on <b>SynKit</b>.
</div>

<div class="alert alert-block alert-success">
<b>What you will gain.</b><br>
You will learn an end-to-end, reproducible pipeline:
<ul>
  <li>Atom-map raw reactions (or load cached mappings).</li>
  <li>Canonicalize mapped reactions to remove map-ID / ordering artifacts.</li>
  <li>Construct <b>ITS graphs</b> and extract <b>reaction centers</b>.</li>
  <li>Cluster reaction centers by <b>graph isomorphism</b> (WL prefilter + exact check).</li>
  <li>Export a compact rule library with representatives and members.</li>
</ul>
</div>

<div class="alert alert-block alert-warning">
<b>Prerequisites.</b><br>
You should be comfortable with typed molecular graphs and morphisms (S01),
subgraph matching and MCS (S02), and atom mapping / ITS ideas (S03–S04).
</div>



## Aim of this talktorial

Automated construction of reaction-rule libraries from data relies on two principles:

1. **Canonicalization** — represent equivalent reactions by a unique, deterministic form
   (stable atom identity, stable ordering, stable graph labels).
2. **Equivalence** — identify when two records express the same transformation using
   **graph invariants** and **typed graph isomorphism**, rather than string similarity.

We start from a dataset of reaction records

$$
\mathcal{D} = \{(\mathrm{id}_i,\; c_i,\; \rho_i)\}_{i=1}^{N},
$$

where

$$
\rho_i
$$

is a reaction string (reaction SMILES) and

$$
c_i
$$

is a class label.

We then build a rule library by transforming each record through the pipeline

$$
\rho \;\longmapsto\; \rho^{\mathrm{map}} \;\longmapsto\; \rho^{\mathrm{canon}}
\;\longmapsto\; T \;\longmapsto\; \mathrm{RC}(T),
$$

with the following objects:

$$
\rho^{\mathrm{map}} \;:\; \text{atom-mapped reaction (explicit atom identity)},
$$

$$
\rho^{\mathrm{canon}} \;:\; \text{canonical mapped reaction (stable under map-ID relabeling)},
$$

$$
T \;:\; \text{ITS graph (encodes pre/post bond labels as one typed graph)},
$$

$$
\mathrm{RC}(T) \;:\; \text{reaction-center graph (the core transformation)}.
$$

Rule **equivalence** is defined by typed graph isomorphism

$$
\mathrm{RC}(T_1) \cong \mathrm{RC}(T_2),
$$

after a fast WL-hash prefilter. Each equivalence class becomes one rule (reaction-center pattern),
storing a representative together with member metadata.

---

## Learning outcomes

After completing this talktorial, you will be able to:

- Load and sanity-check a real reaction dataset (USPTO-50K).
- Produce atom-mapped reactions (or load them from cache) and explain why mapping is needed.
- Canonicalize mapped reactions to stabilize identity across map-ID relabelings.
- Construct ITS graphs, extract reaction centers, and compute graph signatures.
- Cluster reaction centers by WL prefiltering and exact typed graph isomorphism.
- Export a compact rule library suitable for downstream rule application.



## Outline

0. **Setup & data**
1. **Load USPTO-50K reactions**
2. **Atom mapping (RXNMapper)**
3. **Canonicalize mapped reactions**
4. **ITS construction and reaction-center extraction**
5. **Reaction-center clustering (WL prefilter + isomorphism)**
6. **Export a compact rule library**
7. **Discussion**
8. **References**



## 0. Setup & data


## 1. Load local USPTO-50K reactions

The provided CSV contains:

- `id`: document/patent id
- `class`: reaction class label (1..10 in many USPTO-50K variants)
- `reactions`: reaction SMILES, typically `reactants>>products` (unmapped)

We normalize to a table with columns:

- `rxn` : unmapped reaction SMILES
- `id`, `class`



**Basic dataset sanity check**

We ensure:
- reactions contain an arrow (`>>` or `>...>`),
- RDKit can parse at least the product side for a small sample.

(We do not attempt to “fix” dataset issues here; rule libraries inherit dataset quality.)



### Exercise Q1 — Sanity-check reaction strings (splitting + basic counts)

**Task.** Add three columns:

- `n_reactants`: number of dot-separated molecules on the left side,
- `n_products`: number of dot-separated molecules on the right side,
- `ok_split`: whether the reaction can be split into two non-empty sides.

Then report the fraction of valid reactions and show the distribution of `n_reactants` / `n_products`.

<details>
<summary><b>Solution</b></summary>

```python
def _count_side(side: str) -> int:
    side = side.strip()
    if not side:
        return 0
    return len([x for x in side.split(".") if x.strip()])

def _split_ok(rxn: str) -> bool:
    try:
        r, _, p = split_rxn(rxn)
        # print(r)
        return bool(r.strip()) and bool(p.strip())
    except Exception:
        return False


df["ok_split"] = df['rxn'].astype(str).apply(_split_ok)
df["n_reactants"] = df['rxn'].astype(str).apply(lambda s: _count_side(split_rxn(s)[0]) if _split_ok(s) else 0)
df["n_products"]  = df['rxn'].astype(str).apply(lambda s: _count_side(split_rxn(s)[2]) if _split_ok(s) else 0)

print("valid fraction:", df["ok_split"].mean())

display(df[["n_reactants","n_products"]].describe())
```
</details>

---


## 2. Atom mapping with RXNMapper

We start from an unmapped reaction SMILES.

$$
\mathrm{rxn} = R > G > P
$$

RXNMapper converts it into an atom-mapped reaction.

$$
\mathrm{rxn}^{\#} = R^{\#} > G > P^{\#}
$$

Mapped reactions are cached to avoid recomputation and ensure determinism.

> **Tip.** Start with  
> $$
> N = 2000
> $$  
> then scale to  
> $$
> N = \lvert \mathrm{USPTO\text{-}50K} \rvert
> $$



### Exercise Q2 — Cache atom-mapped reactions

**Goal.**  
Avoid recomputing atom mappings by caching results to disk.

**Task.**  
Implement a function `load_or_map_reactions(df, cache_path)` that:

1. Checks whether `cache_path` exists.
2. If it exists, loads atom-mapped reactions from disk.
3. Otherwise:
   - runs **RXNMapper** on `df["rxn"]`,
   - saves the mapped reactions to `cache_path`,
   - returns the mapped list.

Assume:
- reactions are stored in column `rxn`,
- caching uses `json.gz`,
- mapping is done in batches.

---

<details>
<summary><b>Solution</b></summary>

```python
import json
import gzip
from pathlib import Path
from typing import List
from tqdm.auto import tqdm

def load_or_map_reactions(
    df,
    cache_path: Path,
    *,
    rxn_col: str = "rxn",
    batch_size: int = 1,
) -> List[str]:
    """
    Load atom-mapped reactions from cache if available;
    otherwise compute them using RXNMapper and cache the result.
    """
    rxns = df[rxn_col].astype(str).tolist()

    # ---- Load from cache ----
    if cache_path.exists() and cache_path.stat().st_size > 0:
        with gzip.open(cache_path, "rt", encoding="utf-8") as f:
            payload = json.load(f)
        return payload["rxn_mapped"]

    # ---- Compute mapping ----
    from rxnmapper import RXNMapper
    mapper = RXNMapper()

    mapped: List[str] = []
    for i in tqdm(range(0, len(rxns), batch_size), desc="RXNMapper"):
        batch = rxns[i : i + batch_size]
        out = mapper.get_attention_guided_atom_maps(batch)
        mapped.extend([rec.get("mapped_rxn", "") for rec in out])

    # ---- Save cache ----
    payload = {"rxn_mapped": mapped}
    with gzip.open(cache_path, "wt", encoding="utf-8") as f:
        json.dump(payload, f)

    return mapped


# Usage
cache_path = OUTDIR / "uspto50k_mapped.json.gz"
df["rxn_mapped"] = load_or_map_reactions(df.iloc[:,:], cache_path)
df.head()
```


## 3. Canonicalize mapped reactions

Atom-mapped reactions are not unique due to variable component order
and atom-map numbering.
We therefore enforce a canonical representation.

$$
\mathrm{rxn}^{\#} \;\longmapsto\; \mathrm{canon}(\mathrm{rxn}^{\#})
$$

Canonicalization ensures that chemically equivalent reactions admit
an identical atom-mapped form and can be compared structurally.

We use SynKit’s canonicalization based on
Weisfeiler–Lehman refinement, which yields a fast, deterministic but
approximate canonical form; for exact canonicalization, SynKit also
supports nauty-like algorithms with full isomorphism guarantees.


### Exercise Q3 — Canonicalization and deduplication rate

**Goal.**  
Quantify how canonicalization reduces representational redundancy in
atom-mapped reactions.

**Task.**

1. Compute how many mapped reactions collapse after canonicalization.
2. Deduplicate the dataset using the canonical representation.

Report:

- `n_mapped`: total number of mapped reactions,
- `n_unique_canon`: number of unique canonical mapped reactions,
- `compression = n_unique_canon / n_mapped`.

---

<details>
<summary><b>Solution</b></summary>

```python
# Work on a copy to avoid side effects
tmp = df.copy()

# Deduplication statistics
n_mapped = len(tmp)
n_unique_canon = tmp["rxn_mapped_canon"].nunique()
compression = n_unique_canon / max(n_mapped, 1)

print("n_mapped:", n_mapped)
print("n_unique_canon:", n_unique_canon)
print("compression:", compression)

# Deduplicate by canonical atom-mapped reaction
df = tmp.drop_duplicates(subset="rxn_mapped_canon")
print("Shape after deduplication:", df.shape)
```



## 4. ITS construction and reaction-center extraction

Each canonical atom-mapped reaction is represented as an
**Imaginary Transition State (ITS)** graph, which encodes bond changes
directly without explicitly materializing reactants and products.

$$
\mathrm{canon}(\mathrm{rxn}^{\#}) \;\longmapsto\; \mathrm{ITS}
$$

An atom or bond belongs to the **reaction center** if its label differs
between the reactant and product states.

$$
\mathrm{RC} \;\subseteq\; \mathrm{ITS}
$$

The reaction center therefore captures precisely *what changed* during
the reaction.

We prefer SynKit’s ITS construction, which preserves
typed atom and bond annotations and supports direct extraction of the
reaction center.

---


### Exercise Q4 — Inspect an ITS and its reaction center for one example

**Task.** Pick one canonical mapped reaction, construct its ITS graph, and print:

- number of ITS nodes / edges,
- number of reaction-center edges (core),
- a small list of edges whose bond labels changed.

<details>
<summary><b>Solution</b></summary>

```python
from synkit.IO import rsmi_to_its

aam = df.loc[0, "rxn_mapped_canon"]
T = rsmi_to_its(aam, core=False)
RC = rsmi_to_its(aam, core=True)

print("ITS nodes:", T.number_of_nodes(), "edges:", T.number_of_edges())
print("RC  nodes:", RC.number_of_nodes(), "edges:", RC.number_of_edges())

changed = []
for u, v, ed in T.edges(data=True):
    # convention: ITS edge attributes include pre/post labels; adjust keys if needed
    pre = ed.get("bond_pre") or ed.get("bR") or ed.get("pre")
    post = ed.get("bond_post") or ed.get("bP") or ed.get("post")
    if pre != post:
        changed.append((u, v, pre, post))

print("changed edges (u,v,pre,post) sample:", changed[:10])
```
</details>

---


## 5. Reaction-center clustering (WL prefilter + exact isomorphism)

Reaction centers are grouped into **reaction rules** by testing
**typed graph isomorphism**.
Two reaction centers belong to the same rule if they are structurally
identical up to relabeling.

To scale to thousands of reactions, we use a two-stage strategy.

---

### Clustering strategy

1. **WL hash prefilter**  
   Reaction centers are first grouped by a Weisfeiler–Lehman hash,
   which is invariant under graph isomorphism.

2. **Exact typed isomorphism**  
   Within each WL bucket, exact graph isomorphism is tested using
   typed node and edge matches.

3. **Export clusters**  
   Each cluster corresponds to one reaction rule, represented by
   a canonical reaction-center graph and its member indices.

Formally, clustering identifies equivalence classes under isomorphism.

$$
\mathrm{RC}_i \;\sim\; \mathrm{RC}_j
\quad\Longleftrightarrow\quad
\mathrm{RC}_i \cong \mathrm{RC}_j
$$

---


### Exercise Q5 — Implement a minimal reaction-center clustering (WL prefilter + isomorphism)

**Task.** For a small subset of reaction centers, build clusters where two RC graphs are in the same cluster
iff they are isomorphic under your typing predicates.

Return a list of clusters (each cluster is a list of row indices).

<details>
<summary><b>Solution</b></summary>

```python
from synkit.Graph.Feature.wl_hash import WLHash
from networkx.algorithms import isomorphism as iso

def rc_signature(G):
    return WLHash(graph=G).lehman_graph_hash(graph)

def cluster_rc_graphs(rc_graphs):
    # Step 1: bucket by fast signature
    buckets = {}
    for i, G in enumerate(rc_graphs):
        buckets.setdefault(rc_signature(G), []).append(i)

    # Step 2: within each bucket, do exact isomorphism clustering
    clusters = []
    for _, idxs in buckets.items():
        reps = []  # representative index for each cluster
        for i in idxs:
            placed = False
            for rep in reps:
                GM = iso.GraphMatcher(
                    rc_graphs[i], rc_graphs[rep],
                    node_match=node_match,
                    edge_match=edge_match,
                )
                if GM.is_isomorphic():
                    # append to existing cluster
                    for cl in clusters:
                        if rep in cl:
                            cl.append(i)
                            break
                    placed = True
                    break
            if not placed:
                reps.append(i)
                clusters.append([i])
    return clusters

# Example on a small sample
sample = df["RC"].iloc[:200].tolist()
clusters = cluster_rc_graphs(sample)
print("clusters:", len(clusters))
print("cluster sizes (top 10):", sorted([len(c) for c in clusters], reverse=True)[:10])
```
</details>


## 6. Export a compact rule library

We export:

- `rule_library.pkl.gz`: one record per cluster
  - representative core graph (nodes + edges + labels)
  - member indices
  - optional metadata (class distribution, example ids)


### Exercise Q6 — Export and validate a compact rule library

**Task.**  
From your deduplicated dataframe (after canonicalization and clustering), export a compact
rule library and validate basic invariants:

1. Sample at most n rules per class (stratified).
2. Export to a CSV file with columns: `id`, `class`, `rxn_mapped_canon`, `rule_id`.
3. Verify that `rule_id` is unique in the exported library.

<details>
<summary><b>Solution</b></summary>

```python
# Assumes you already have a dataframe like `df_rules` with one row per rule
# and a stable identifier `rule_id` (e.g., cluster id) plus the canonical mapped reaction.
# If you used `data` for clustered RC graphs, adapt accordingly.

MAX_PER_CLASS = 200

df_rules = df.copy()

# Example: use `cluster_id` if available, otherwise fall back to index
if "cluster_id" not in df_rules.columns:
    df_rules["cluster_id"] = range(len(df_rules))
df_rules["rule_id"] = df_rules["cluster_id"]

# Stratified sampling per reaction class
df_lib = stratified_random_sample(
    df_rules,
    by="class",
    n=MAX_PER_CLASS,
    random_state=0,
)

# Keep a minimal schema
keep = ["id", "class", "rxn_mapped_canon", "rule_id"]
df_lib = df_lib[keep].copy()

# Validate uniqueness
assert df_lib["rule_id"].is_unique, "rule_id must be unique in the exported library"

OUT_PATH = "rule_library_compact.csv"
df_lib.to_csv(OUT_PATH, index=False)

print("exported:", len(df_lib), "rules ->", OUT_PATH)
print("classes:", df_lib["class"].nunique())
```
</details>



## 7. Discussion

- Determinism vs. speed: WL signatures speed up clustering but must be validated by exact isomorphism.
- Atom mapping quality dominates: mapping errors propagate into ITS and rule identity.
- Canonicalization is essential: without it, the same transformation can appear under many map-ID permutations.
- Practical tip: cache mappings and intermediate artifacts (canonical AAM, ITS) to make experiments reproducible.



### Exercise Q7 — Matching semantics: what changes when you change typing?

**Task.**  
In this notebook, rule equivalence is defined by isomorphism of typed graphs (ITS / RC graphs).
Consider the following two modifications and predict their effect on clustering:

1. Remove `formal_charge` from the node labels used in matching.
2. Ignore bond order in the edge labels (treat single/double as identical).

For each case, explain whether you expect the number of clusters to increase, decrease,
or remain similar, and why.

<details>
<summary><b>Solution (expected reasoning)</b></summary>

- Removing `formal_charge` makes the typing **coarser**.
  More nodes become compatible under `node_match`, so more graphs become isomorphic.
  **Expected effect:** the number of clusters tends to **decrease** (more merges).

- Ignoring bond order also makes the typing **coarser**.
  Many transformations that differ only by order changes become indistinguishable under matching.
  **Expected effect:** the number of clusters tends to **decrease**, possibly sharply,
  because distinct reaction centers may collapse to the same unlabeled topology.

In both cases, the clustering becomes less chemically specific: it may improve recall
for noisy data but risks merging chemically distinct rules.
</details>



## 8. References and further reading

The following resources provide additional context and deeper technical details
related to graph-based reaction modeling, rule extraction, and reaction templates.

- **SynKit**  
  A graph-theoretic toolkit for modeling chemical reactions using typed graphs,
  atom mappings, and rule-based transformations.  
  *JCIM (2025).*  
  https://pubs.acs.org/doi/10.1021/acs.jcim.5c02123

- **SynTemp**  
  A framework for extracting, analyzing, and applying reaction templates with
  an emphasis on structural changes and reaction centers.  
  *JCIM (2024).*  
  https://pubs.acs.org/doi/10.1021/acs.jcim.4c01795

**Suggested background reading**

- Ehrig, H. *et al.* **Fundamentals of Algebraic Graph Transformation** —  
  the foundational reference for DPO rewriting.
- Willett, P. **Chemical similarity searching** —  
  background on substructure and MCS concepts in chemistry.
- Schneider, N. *et al.* **RXNMapper** —  
  attention-based atom mapping for chemical reactions.

These references collectively connect the formal graph-rewriting perspective
used in this notebook with practical cheminformatics workflows.
