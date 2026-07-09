# S07: From Atom-Mapped Reactions to DPO Rules

This talktorial builds the bridge from atom-mapped reaction records to a reusable DPO reaction-rule library. It combines standardization, ensemble mapping, ITS construction, WL prefiltering, graph isomorphism, and optional MØD export [\[1\]](#6.-References), [\[2\]](#6.-References), [\[3\]](#6.-References).



## Aim of this talktorial

1. Standardize and canonicalize mapped reaction strings so equivalent maps can be compared reliably.
2. Construct high-confidence reaction centers from ensemble atom maps and ITS graphs.
3. Cluster, deduplicate, and export reaction centers as DPO rules.

---

## Learning outcomes

After completing this talktorial, you will be able to:

- standardize and canonicalize mapped reaction SMILES,
- construct an **ensemble atom map** by comparing multiple mapper outputs,
- build an **ITS graph** from a mapped reaction,
- extract a **reaction center** from an ITS graph,
- cluster reaction centers using **WL hashing** and exact **graph isomorphism**,
- convert a reaction center into a **DPO rule** $L \leftarrow K \rightarrow R$, and
- explain why consensus mapping and deduplication improve the reliability of a rule library.

---

## Outline

- [0. Setup & data](#0.-Setup-&-data)
- [1. Reaction preprocessing](#1.-Reaction-preprocessing)
- [2. Ensemble atom-to-atom mapping agreement](#2.-Ensemble-atom-to-atom-mapping-agreement)
- [3. Rule library construction](#3.-Rule-library-construction)
- [4. Discussion](#4.-Discussion)
- [5. Quiz](#5.-Quiz)
- [6. References](#6.-References)


<p>
  The complete S07 workflow can be summarized as a linear pipeline.
  Each step progressively reduces mapping ambiguity and improves the reproducibility
  of graph-based reaction rule extraction.
</p>

<figure style="text-align: left;">
  <img src="../../docs/_static/S07/flow_chart.svg"
       alt="S07 workflow for consensus-driven reaction rule extraction"
       style="width: 100%; max-width: 900px;">
  <figcaption>
    <b>Figure 1.</b> Overview of the S07 workflow for extracting graph-based reaction
    rules from atom-to-atom mapped reactions. Reactions are first standardized and
    canonicalized, filtered by ensemble agreement, converted into ITS graphs, and
    reduced to reaction centers. Reaction centers are clustered to identify unique
    representatives, completed with hydrogen information, and finally exported as
    GML reaction rules. Discard branches indicate cases rejected because of
    insufficient mapping consensus or ambiguous hydrogen completion.
  </figcaption>
</figure>


## 0. Setup & data



This notebook uses **SynKit** [\[1\]](#6.-References) to accelerate development: fast access to molecular graphs, rule application (DPO/ITS), pattern matching, and utilities for building reaction-centric workflows.

we can easily install synkit via pypi:
```bash
pip install synkit
```


The dataset incorporates three distinct atom mapping methodologies: `rxn_mapper` [\[4\]](#6.-References), `graphormer` mapper [\[5\]](#6.-References), and `local_mapper` [\[6\]](#6.-References).


## 1. Reaction preprocessing

1. Validate & standardize SMILES using `Standardize` [\[7\]](#6.-References) (drop invalid/unparsable).

2. Canonicalize atom-map numbers with `CanonRSMI` (from **S06**) for identical mappings.


## 2. Ensemble atom-to-atom mapping agreement

Compare the canonicalized mapper outputs by exact string equality; rows where all canonical maps match are high-confidence AAM assignments [\[8\]](#6.-References).


**Definition (Canonical Reaction SMILES).**  
Two atom-mapped reaction SMILES $s_1$ and $s_2$ are *equivalent* if their corresponding ITS graphs $\Gamma_1$ and $\Gamma_2$ are isomorphic as labeled graphs: $\Gamma_1 \cong \Gamma_2$. The *canonical reaction SMILES* $\mathrm{canon}(s)$ is the lexicographically smallest string in the equivalence class, obtained by applying a canonical atom-map permutation $\pi: \mathbb{N} \to \mathbb{N}$ computed via WL refinement (approximation) or IR (exact).

**Definition (Ensemble Agreement).**  
Given $n$ atom-mapping functions $\varphi_1, \ldots, \varphi_n$ applied to the same reaction $\varrho$, the ensemble *agrees* on $\varrho$ if all canonical strings are identical:

$$
\mathrm{canon}(\varphi_1(\varrho)) = \mathrm{canon}(\varphi_2(\varrho)) = \cdots = \mathrm{canon}(\varphi_n(\varrho))
$$

When the ensemble disagrees, *refinement* (exact isomorphism check) is used to resolve the discrepancy.

**Definition (WL Reaction Hash).**  
The *WL hash* of a reaction center ITS $\Gamma_{\mathrm{RC}}$ is the string representation of the stable WL partition histogram:

$$
h_{\mathrm{WL}}(\Gamma_{\mathrm{RC}}) = \mathrm{encode}\!\left(\mathcal{P}^*(\Gamma_{\mathrm{RC}})\right)
$$

Reactions with identical WL hashes are *candidate isomorphs*; exact isomorphism (graph matching) is then used to confirm.



Weisfeiler-Lehman [\[2\]](#6.-References) based canonicalization is fast but not exact. Use WL/canonical strings as a cheap filter, then run exact ITS isomorphism on rows where canonical strings disagree or ambiguity remains.

In **SynKit**, we expose `AAMValidator` to compare directly atom-to-atom map by converting to ITS and checking isomorphism


### Ensemble mapper agreement distribution

Three mappers (RXNMapper, Graphormer, Local Mapper) each produce an atom map for each reaction. The ensemble step classifies reactions by how many mappers agree. High-confidence reactions (all three agree) form the primary training signal; disagreements are resolved by WL-based isomorphism checking or discarded.


**Q1 — Ensemble mapping**

Given a DataFrame with several mapper columns, `ensemble_maps(df, map_cols)`:
- canonicalizes each mapper column with `std_canon(...)`,
- accepts rows where all canonical values are identical,
- runs `refinement(...)` on the leftovers,
- returns `(solved_df, unsolved_df)` where `solved_df` keeps only metadata + one `aam` column.
- combine all high-confident maps and drop duplicate if possible


---

<details> <summary><b>Solution:</b></summary>

```python
from typing import List, Tuple
import pandas as pd
from tqdm.auto import tqdm

tqdm.pandas() 

def ensemble_maps(df: pd.DataFrame, map_cols: List[str]) -> Tuple[pd.DataFrame, pd.DataFrame]:
    missing = [c for c in map_cols if c not in df.columns]
    if missing:
        raise KeyError(f"missing columns: {missing}")

    canon_cols = []
    for col in map_cols:
        tgt = f"{col}_canon"
        canon_cols.append(tgt)
        tqdm.pandas(desc=f"Canonicalization {col}") 
        df[tgt] = df[col].progress_apply(std_canon)

    match, unmatch = extract_aam(df, keys=tuple(canon_cols))
    solve, unsolve = refinement(unmatch, aam_cols=tuple(canon_cols))

    out = pd.concat([match, solve], axis=0).reset_index(drop=True)
    out = out.drop_duplicates(subset="aam", keep="first").reset_index(drop=True)
    return out, unsolve

    
```
</details>


## 3. Rule library construction

With high-confidence `aam` maps in hand, convert each mapped reaction into an ITS graph, extract the reaction center, and turn that center into a reusable rule entry for the library.


### 3.1. Graph construction


We now use **synedu.Utils** to convert each atom-mapped reaction to an ITS graph and extract the reaction center.



### 3.2. Core refinement


#### Graph Clustering

Now we have 8,962 reaction centers, but some of them can be isomorphic [\[9\]](#6.-References), we can develop a simple pairwise grouping algorithm

With a set $S=\{\Gamma_1,\Gamma_2,\dots,\Gamma_n\},\qquad n\ge 2$, we want a partition $T$ of $S$ into isomorphism classes.

**Pseudocode**

$$
\begin{array}{l}
T \leftarrow \varnothing\\[4pt]
S \leftarrow \{\Gamma_1,\Gamma_2,\dots,\Gamma_n\}\\[6pt]
\text{while } S \neq \varnothing :\\
\quad q^{*} \leftarrow \text{arbitrary element of } S\\
\quad S \leftarrow S \setminus \{q^{*}\}\\
\quad C \leftarrow \{q^{*}\}\\
\quad \text{for each } q \in \mathrm{copy}(S) :\\
\quad\quad \text{if }\operatorname{VF2\_iso}(q,q^{*}) \text{ then}\\
\quad\quad\quad C \leftarrow C \cup \{q\}\\
\quad\quad\quad S \leftarrow S \setminus \{q\}\\
\quad\quad \text{end if}\\
\quad \text{end for}\\
\quad T \leftarrow T \cup \{C\}\\
\text{end while}\\[4pt]
\text{return } T
\end{array}
$$



In SynKit, this algorithm is implemented in `GraphCluster`


**Q2 — Unique rules**

You have a list of dictionaries where each item represents a reaction-center record and contains a `class` label (the isomorphism class). Implement `get_unique_rc(...)` to extract one representative reaction center per class.

---

<details> <summary><b>Solution:</b></summary>

```python
from typing import List, Dict, Any

def get_unique_rc(
    items: List[Dict[str, Any]],
    class_key: str = "class",
    rc_key: str = "RC",
    keep: str = "first",
) -> List[Any]:
    if keep not in ("first", "last"):
        raise ValueError("keep must be 'first' or 'last'")

    reps: dict = {}  # class -> (idx, rc)
    for i, it in enumerate(items):
        cls = it.get(class_key)
        if cls is None:
            continue
        rc = it.get(rc_key)
        if keep == "first":
            if cls not in reps:
                reps[cls] = (i, rc)
        else:  # keep == "last"
            reps[cls] = (i, rc)

    # sort by stored index and return rc values
    ordered = sorted(reps.values(), key=lambda x: x[0])
    return [rc for _, rc in ordered]


    
```
</details>


Alternatively, we can leverage Synkit’s utility functions


#### Hydrogen completion in reaction centers

Implicit hydrogens are stored as a `hcount = (h_r, h_p)` tuple on each ITS node.
When `h_r ≠ h_p` a hydrogen atom is *transferred* between molecules.
Graph Modelling Language format requires **explicit H nodes** so the rule engine can track every atom.

**`h_expand_its` — four-step pipeline**

| Step | Action |
|------|--------|
| ① | Decompose ITS → separate reactant / product graphs (convert tuple attrs → scalars) |
| ② | Identify *donor* atoms (`h_r > h_p`) and *acceptor* atoms (`h_r < h_p`) |
| ③ | Pair donors ↔ acceptors (`n_shared = min(|donors|, |acceptors|)`); add one **shared H node** per pair |
| ④ | Rebuild ITS — the shared H node gets **two edges** |

| H edge | `order` | GML block | Colour |
|--------|---------|-----------|--------|
| H – donor    | `(1, 0)` | `left`  | 🔴 broken |
| H – acceptor | `(0, 1)` | `right` | 🟢 formed |

> **Only balanced transfers are expanded.**  
> If donors and acceptors do not have the same count (e.g. H is gained from an unmapped
> solvent molecule), the surplus is left implicit.  Adding single-ended H nodes (broken
> *or* formed only) to a DPO rule would introduce atoms with no counterpart on the other
> side — chemically meaningless and rejected by MØD.

> **Ambiguous cases** — when multiple donors and acceptors exist, different pairings can
> produce **non-isomorphic** reaction centers (see §3.2).


#### Ambiguous hydrogen expansion

When an RC contains **multiple donor atoms and multiple acceptor atoms**,
more than one valid donor-acceptor pairing exists.  Different pairings can produce
**non-isomorphic** reaction centers, each representing a different mechanistic picture.

The reaction below involves **trichloroacetonitrile reacting with an *N*-hydroxy compound**.
Its RC has two donors and two acceptors:

| Atom | Role | hcount change |
|------|------|---------------|
| **O6** (N–O–H oxygen) | donor  | (1→0) loses H |
| **C12** (aromatic C–H) | donor  | (1→0) loses H |
| **N5** (amine nitrogen) | acceptor | (0→1) gains H |
| **N18** (nitrile nitrogen) | acceptor | (0→1) gains H |

This gives **two distinct pairings**:

| Pairing | O6 donates to | C12 donates to | Topology |
|---------|---------------|----------------|----------|
| A (natural) | **N5** | **N18** | H stays near its original bond partner |
| B (crossed) | **N18** | **N5**  | H crosses — different bridge topology |

`h_expand_all` enumerates every permutation and returns only the **structurally distinct** outcomes.



### 3.3. Rule selection


#### Rule library size vs dataset size

A key question: does the rule vocabulary *saturate* (good for generalization)
or grow indefinitely (poor coverage at test time)?
We subsample the dataset at increasing sizes and plot unique rule count.


Applying the same reduction logic to isomorphic ITS graphs is feasible, though computationally more expensive than focusing solely on the reaction center.
To optimize the computationally intensive ITS reduction, we can leverage Weisfeiler-Lehman [\[2\]](#6.-References) pre-filtering (per **S06**). This ensures that expensive isomorphism checks are only performed on high-probability candidates


#### Rule frequency and cluster size distribution

The WL hash partitions reactions by their reaction-centre topology. The **left panel** shows the distribution of cluster sizes on a log scale — a long tail is typical (a few very common rules, many rare ones). The **right panel** shows the 20 most frequent rules. Singleton clusters represent unique reaction patterns seen only once in the dataset.


#### Top-10 rule classes — reaction center gallery

The most frequent WL hash classes represent the most common reaction mechanisms in the dataset.  
We visualize the reaction center ITS graph for the representative of each top-10 class.



While `wl_hash` effectively yields the same equivalence classes as graph isomorphism in this instance, it is important to note that the WL hash is an approximation. For absolute precision, it should be supplemented or replaced by a formal isomorphism check.


**Q3 — ITS uniqueness**

Now integrate a function named `cluster_with_wl_prefilter(df, graph_key, wl_key)` to prefilter the bucket before

---

<details> <summary><b>Solution:</b></summary>

```python
import pandas as pd
from synedu.Utils import cluster_its_graphs

def cluster_with_wl_prefilter(
    df: pd.DataFrame,
    graph_key: str = "RC",
    wl_key: str = "wl_hash",
) -> pd.DataFrame:
    df2 = df.copy()
    if wl_key not in df2.columns:
        raise KeyError(f"{wl_key} not found in DataFrame")

    results = []
    next_class = 0
    for wl_val, idxs in df2.groupby(wl_key).groups.items():
        bucket = df2.loc[idxs]
        graphs = bucket[graph_key].tolist()
        labels = cluster_its_graphs(graphs)
        bucket = bucket.copy()
        bucket["class"] = [lbl + next_class for lbl in labels]
        next_class = int(bucket["class"].max()) + 1
        results.append(bucket)

    return pd.concat(results, ignore_index=True) if results else pd.DataFrame()
```
</details>



#### DPO rules

Each representative reaction center in `rc_h` (the H-expanded list built in §3.1)
is now converted to a **MØD GML rule** — the formal DPO representation $L \leftarrow K \rightarrow R$.

| GML block | Meaning | ITS edges included |
|-----------|---------|-------------------|
| `left`    | bonds **broken** in the reaction | `order = (b_r, 0)` — exist only in reactant |
| `right`   | bonds **formed** in the reaction | `order = (0, b_p)` — exist only in product |
| `context` | bonds **unchanged** + all atoms | `order = (b, b)` with `b > 0` |

Because we use `rc_h`, every transferred hydrogen is represented as an **explicit H node**
with a broken bond to its donor (`left`) and a formed bond to its acceptor (`right`).
This satisfies the MØD requirement that all atoms — including hydrogens — are tracked explicitly.



**Building `its_to_gml` — step by step**

The conversion has three sub-tasks:

| Sub-task | Input | Output |
|----------|-------|--------|
| ① Classify each edge | `order = (b_r, b_p)` tuple | `left` / `right` / `context` |
| ② Label each node | `element`, `formal_charge` attrs | atom label string (e.g. `"N+"`) |
| ③ Label each bond | scalar bond order `b` | GML bond string (`"1"`, `"2"`, `"#"`, `"3"`) |

Consider the small RC below.
Each edge carries an `order = (b_r, b_p)` tuple:

| Edge | `order` | GML block | Meaning |
|------|---------|-----------|---------|
| C – O | `(1, 0)` | `left` | bond **broken** |
| C – N | `(0, 1)` | `right` | bond **formed** |
| C – C | `(1, 1)` | `context` | bond **preserved** |
| N – H | `(0, 1)` | `right` | H **transferred** (formed on acceptor) |
| O – H | `(1, 0)` | `left` | H **transferred** (broken on donor) |



##### Step ① — helper functions

Two small helpers underpin `its_to_gml`:

| Helper | Input | Output | Example |
|--------|-------|--------|---------|
| `_bond_label(b)` | scalar float bond order | GML label string | `2.0 → "2"`, `1.5 → "#"` |
| `_node_label(attrs)` | node attribute dict | atom label string | `{"element":"N","formal_charge":1} → "N+"` |



##### Step ② — assembling `its_to_gml`

With the helpers in place, the full function is a single pass over edges:

```
for each edge (u, v):
    br, bp = order
    if br > 0 and bp == 0  →  left  (broken)
    if br == 0 and bp > 0  →  right (formed)
    if br > 0 and bp > 0
        br ≠ bp             →  left + right (order changes)
        br == bp            →  context (preserved, optional)

all nodes always go into context
```



**Q5 — Implement `its_to_gml`**

Using the edge classification logic above, implement `its_to_gml(its, rule_id, include_context_edges)` that converts an ITS/RC graph to a MØD GML string.

Requirements:
- `_bond_label(order)` maps a scalar float bond order to its GML string: `1.0→"1"`, `1.5→"#"`, `2.0→"2"`, `3.0→"3"`.
- `_node_label(attrs)` returns the element symbol with `+` or `-` suffix if `formal_charge ≠ 0`; handle both scalar and `(fc_r, fc_p)` tuple forms.
- For each edge `(b_r, b_p)`: if `b_r ≠ b_p`, put the reactant side (`b_r > 0`) in `left` and the product side (`b_p > 0`) in `right`; if `b_r == b_p > 0` and `include_context_edges` is True, put it in `context`.
- Always emit every node in the `context` block.

---

<details> <summary><b>Solution:</b></summary>

```python
import networkx as nx
from typing import Union


def _bond_label(order: float) -> str:
    o = round(float(order), 2)
    return {1.0: "1", 1.5: "#", 2.0: "2", 3.0: "3"}.get(o, str(o))


def _node_label(attrs: dict) -> str:
    elem = attrs.get("element", "*")
    chg_raw = attrs.get("formal_charge", 0)
    chg = int((chg_raw[0] if isinstance(chg_raw, tuple) else chg_raw) or 0)
    if chg > 0: return f"{elem}+"
    if chg < 0: return f"{elem}-"
    return str(elem)


def its_to_gml(
    its: nx.Graph,
    rule_id: Union[int, str] = 0,
    include_context_edges: bool = True,
) -> str:
    left_edges, right_edges, ctx_edges = [], [], []
    for u, v, d in its.edges(data=True):
        br, bp = d.get("order", (0.0, 0.0))
        br, bp = float(br), float(bp)
        if br != bp:
            if br > 0:
                left_edges.append((u, v, _bond_label(br)))
            if bp > 0:
                right_edges.append((u, v, _bond_label(bp)))
        elif include_context_edges and br > 0:
            ctx_edges.append((u, v, _bond_label(br)))

    lines = [f'rule [', f'  ruleID "{rule_id}"', '  left [']
    for u, v, lbl in left_edges:
        lines.append(f'    edge [ source {u} target {v} label "{lbl}" ]')
    lines.append('  ]')

    lines.append('  context [')
    for n in sorted(its.nodes()):
        lines.append(f'    node [ id {n} label "{_node_label(its.nodes[n])}" ]')
    for u, v, lbl in ctx_edges:
        lines.append(f'    edge [ source {u} target {v} label "{lbl}" ]')
    lines.append('  ]')

    lines.append('  right [')
    for u, v, lbl in right_edges:
        lines.append(f'    edge [ source {u} target {v} label "{lbl}" ]')
    lines.append('  ]')
    lines.append(']')
    return '\n'.join(lines)
```

Test on the demo RC:
```python
print(its_to_gml(_demo, rule_id="demo", include_context_edges=True))
```
</details>



**Example GML output** (one rule, L/K/R encoded as sub-blocks):

```text
rule [
  ruleID "0"
  left [
    edge [ source 5 target 19 label "1" ]
  ]
  context [
    node [ id 5 label "N" ]
    node [ id 19 label "H" ]
    node [ id 23 label "N" ]
  ]
  right [
    edge [ source 23 target 19 label "1" ]
  ]
]
```

**Interpretation**

| Block | DPO component | Meaning |
|-------|---------------|---------|
| `context` (K) | preserved subgraph | atoms that exist in both reactant and product; all nodes always appear here |
| `left` (L∖K) | deleted bonds | bonds present only in the reactant — broken during the reaction |
| `right` (R∖K) | created bonds | bonds present only in the product — formed during the reaction |

In the example above, the H node (id 19) bridges from `N` (atom 5) in the reactant to `N` (atom 23) in the product — a direct N→N proton transfer encoded as an explicit H node with one broken and one formed bond.

Because `rc_h` contains H-expanded reaction centers, every transferred hydrogen already appears as an explicit node. The resulting GML rules are directly compatible with **MØD** [\[3\]](#6.-References) without further post-processing.



**Q4 — Data pipeline**

Implement a single high-level function that ingests a pandas DataFrame containing multiple atom-mapping candidates per reaction and produces a deduplicated library of reaction centers / DPO rules (optionally exported to GML), plus the rows that remain ambiguous.

---

<details> <summary><b>Solution:</b></summary>

```python
from typing import List
from synedu.Utils import rsmi_to_its, cluster_its_graphs, stratified_sample, wl_hash

def rule_lib_pipeline(
        data: pd.DataFrame,
        map_cols: List[str],
        export_gml: bool = False):

    # 1. get high-confidence aam
    data, _ = ensemble_maps(data, map_cols)

    # 2. convert to RC
    tqdm.pandas(desc="ITS construction")
    data["ITS"] = data["aam"].progress_apply(lambda s: rsmi_to_its(s, core=False))

    tqdm.pandas(desc="Reaction center construction")
    data["RC"] = data["aam"].progress_apply(lambda s: rsmi_to_its(s, core=True))

    # 3. compute wl_hash for prefilter
    tqdm.pandas(desc="WL hash")
    _wl = lambda g: wl_hash(g, node_attrs=["element", "formal_charge"], edge_attrs=["order"])
    data["wl_hash"] = data["RC"].progress_apply(_wl)

    # 4. clustering
    rc_graphs = data["RC"].tolist()
    data["class"] = cluster_its_graphs(rc_graphs,
                                       node_attrs=["element", "formal_charge"],
                                       edge_attrs=["order"])
    result = data.to_dict("records")

    # 5. get unique reaction centers and expand H
    rc_records = stratified_sample(result, key="class")
    rc = [r["RC"] for r in rc_records]
    rc_h = [h_expand_its(g) for g in rc]

    # 6. convert to GML
    if export_gml:
        return [its_to_gml(g, rule_id=i, include_context_edges=False)
                for i, g in enumerate(rc_h)]
    return rc_h
```
</details>



<a id="4-discussion"></a>

## 4. Discussion

- Atom-mapping quality can be significantly improved through ensemble techniques, specifically by performing isomorphism checks on generated ITS graphs.
- While WL-based canonicalization is effective for pre-filtering and deduplication, it remains an approximation and cannot fully replace exact isomorphism checks for definitive verification.
- Clustering reaction centers in large-scale datasets is computationally intensive; however, implementing a WL hashing pre-filter significantly accelerates the process by pruning the search space.
- Rules can be persisted as NetworkX graph objects for deep computational tasks or exported in GML format for a lightweight, human-readable, and explainable representation.


<a id="5-quiz"></a>

## 5. Quiz

1. Why is reaction standardization needed before comparing atom-mapped reactions from different mappers?
2. How does comparing several mapper outputs help identify high-confidence atom mappings?
3. Why is WL hashing useful for clustering reaction centers, and why is exact isomorphism still needed?
4. What information must a DPO rule preserve so it can be exported and reused as a reaction-rule library entry?



## 6. References

1. Phan, T.-L.; González Laffitte, M. E.; Weinbauer, K.; Merkle, D.; Andersen, J. L.; Fagerberg, R.; Gatter, T.; Stadler, P. F. *SynKit: A Graph-Based Python Framework for Rule-Based Reaction Modeling and Analysis.* *Journal of Chemical Information and Modeling* (2025). https://doi.org/10.1021/acs.jcim.5c02123
2. Shervashidze, N.; Schweitzer, P.; van Leeuwen, E. J.; Mehlhorn, K.; Borgwardt, K. M. *Weisfeiler-Lehman graph kernels.* *Journal of Machine Learning Research* **12**, 2539-2561 (2011). https://www.jmlr.org/papers/v12/shervashidze11a.html
3. Andersen, J. L.; Flamm, C.; Merkle, D.; Stadler, P. F. *A software package for chemically inspired graph transformation.* In *International Conference on Graph Transformation*, 73-88 (Springer, 2016). https://doi.org/10.1007/978-3-319-40530-8_5
4. Schwaller, P. *et al.* *Extraction of organic chemistry grammar from unsupervised learning of chemical reactions (RXNMapper).* *Science Advances* **7**, eabe4166 (2021). https://doi.org/10.1126/sciadv.abe4166
5. Nugmanov, R.; Dyubankova, N.; Gedich, A.; Wegner, J. K. *Bidirectional graphormer for reactivity understanding: neural network trained to reaction atom-to-atom mapping task.* *Journal of Chemical Information and Modeling* **62**(14), 3307-3315 (2022). https://doi.org/10.1021/acs.jcim.2c00344
6. Chen, S.; An, S.; Babazade, R.; Jung, Y. *Precise atom-to-atom mapping for organic reactions via human-in-the-loop machine learning.* *Nature Communications* **15**, 2250 (2024).
7. Schneider, N.; Stiefl, N.; Landrum, G. A. What's What: The (Nearly) Definitive Guide to Reaction Role Assignment. *Journal of Chemical Information and Modeling* **56**, 2336-2346 (2016). https://doi.org/10.1021/acs.jcim.6b00564
8. Phan, T.-L.; Weinbauer, K.; González Laffitte, M. E.; Pan, Y.; Merkle, D.; Andersen, J. L.; Fagerberg, R.; Flamm, C.; Stadler, P. F. *SynTemp: Efficient Extraction of Graph-Based Reaction Rules from Large-Scale Reaction Databases.* *Journal of Chemical Information and Modeling* (2025). https://doi.org/10.1021/acs.jcim.4c01795
9. Cordella, L. P.; Foggia, P.; Sansone, C.; Vento, M. A (Sub)Graph Isomorphism Algorithm for Matching Large Graphs. *IEEE Transactions on Pattern Analysis and Machine Intelligence* **26**, 1367-1372 (2004). https://doi.org/10.1109/TPAMI.2004.75
