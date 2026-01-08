# S02 · Subgraph Isomorphism and MCS: Symmetry-aware Matching

<div class="alert alert-block alert-info">
<b>Welcome back to SynEdu.</b><br>
This talktorial continues from <b>S01</b>. We build on typed molecular graphs to study
<b>pattern matching</b> and <b>maximum common substructure (MCS)</b> with careful treatment of symmetry.
</div>

<div class="alert alert-block alert-success">
<b>What you will gain.</b><br>
You will learn how to compute <b>subgraph isomorphisms</b> (pattern → host),
how to <b>deduplicate</b> symmetry-inflated matches using automorphism orbits,
and how to interpret <b>MCS</b> results in both RDKit and NetworkX.
</div>

<div class="alert alert-block alert-warning">
<b>Prerequisites.</b><br>
You should be comfortable with the definitions and code utilities from <b>S01</b>
(typed graphs, morphisms, isomorphism, automorphisms).
</div>


## Aim of this talktorial

In **S01**, we built the foundations: **typed molecular graphs** and **typed morphisms**
(isomorphisms, automorphisms) as the matching backbone for graph transformation.

This talktorial (**S02**) moves to the *workhorse alignments* that appear throughout reaction modeling
and rule-based systems (e.g., DPO-style rewriting):

1. **Subgraph matching**: finding when a *pattern* graph occurs inside a *host* graph
   via **typed subgraph monomorphisms** (a.k.a. subgraph isomorphisms in practice).
2. **MCS alignment**: using **maximum common substructure (MCS)** as a pragmatic alignment primitive
   when an exact pattern is unknown or when molecules differ by edits.

We keep implementations minimal and transparent, using:

- **NetworkX** for explicit, attribute-based subgraph morphisms and symmetry-aware deduplication,
- **RDKit** for chemistry-aware substructure/MCS behavior (sanitization, aromaticity, valence, MCS).

**Data example:** `data/molecules.csv`

---

## Learning outcomes

After completing this talktorial, you will be able to:

- Formulate **typed subgraph matching** as an **injective typed graph morphism**
  $$
  \varphi : V(P) \hookrightarrow V(H),
  $$
- Compute **subgraph isomorphisms** (pattern → host) in **NetworkX**, and interpret multiple matches.
- Apply **symmetry-aware deduplication** of matches using automorphism or orbit ideas.
- Use **RDKit substructure matching** and understand when it differs from a pure graph matcher.
- Compute **MCS with RDKit** and turn it into an **alignment map** between two molecules.
- Compare **RDKit vs NetworkX** matching behavior and attribute choices, and decide which is appropriate for
  later steps (rule extraction, reaction center localization, DPO rule application).

---

## Outline

0. **Setup & data**
1. **NetworkX subgraph isomorphism (pattern → host)**
2. **Symmetry and deduplicating equivalent matches**
3. **RDKit substructure matching: chemistry-aware behavior**
4. **MCS with RDKit: maximum common substructure as alignment**
5. **From MCS to atom maps: building a practical correspondence**
6. **RDKit vs NetworkX morphisms: comparison and pitfalls**
7. **Discussion, quiz, and references**



## 0. Setup & Data


## 1. Subgraph isomorphism (pattern → host)

In rule-based reaction modeling, we repeatedly solve the **pattern-in-host** query:

$$
\text{Does a pattern graph } P \text{ occur inside a host graph } G?
\quad \text{If yes, what are the embeddings?}
$$

Formally, a **subgraph isomorphism** is an **injective, label-preserving graph morphism**

$$
f : V(P) \hookrightarrow V(G)
$$

such that:

- **Atom (node) labels are preserved**
  $$a_G(f(v)) = a_P(v)\quad \forall v \in V(P)$$

- **Bond existence and bond types are preserved**
  $$uv \in E(P)\ \Rightarrow\ f(u)f(v) \in E(G), \quad b_G(f(u)f(v)) = b_P(uv)$$

Intuitively, \(f\) is a **typed monomorphism**: it embeds the pattern into the host
without collisions (injective), while respecting chemical identity (types).

### Practical note: many matches due to symmetry
Even when the chemical occurrence is “the same”, symmetric graphs can admit many
equally valid embeddings:

- **Host symmetry** (automorphisms of \(G\)) produces multiple placements.
- **Pattern symmetry** (automorphisms of \(P\)) produces multiple equivalent mappings.
- If both are symmetric, matches can multiply combinatorially.

NetworkX exposes this via:

- `GraphMatcher.subgraph_isomorphisms_iter()` — enumerates all injective embeddings
  that satisfy `node_match` and `edge_match`.

For downstream tasks (reaction center extraction, rule application, deduplication),
we often need to **post-process** these matches to remove symmetry-equivalent
embeddings, typically by orbit-based canonicalization or choosing a canonical
representative embedding.


## 2. Deduplication of subgraph embeddings

Raw subgraph matches often contain many symmetry-equivalent embeddings.  
We give a minimal, renderer-friendly description of the **host-image**
deduplication strategy.

Let
$$
m : V(P) \longrightarrow V(G)
$$
be a mapping from pattern nodes to host nodes.
The **image** of \( m \) is the set of host nodes occupied by the embedding:

$$
\mathrm{img}(m) = \{\, m(v) \mid v \in V(P) \,\} \subseteq V(G).
$$

We deduplicate embeddings by grouping all mappings that share the same image.
Operationally, we use the sorted tuple of \( \mathrm{img}(m) \) as a canonical key:

$$
\mathrm{key}(m) = \mathrm{tuple}\!\left(\mathrm{sorted}\big(\mathrm{img}(m)\big)\right).
$$

This collapses symmetry variants that differ only by permutations of pattern
nodes (i.e., different bijections with the same image) and yields the distinct
*placements* of the pattern inside the host.



### Q1 — Deduplicate modulo pattern automorphisms

**Goal.**  
Implement a function `dedup_by_pattern_image_with_orbits(matches, pattern_autos, ...)` that groups/filters pattern→host mappings **modulo pattern automorphisms**.

Intuitively, two mappings \(m_1,m_2: P \to H\) are equivalent if there exists a pattern automorphism \(\varphi \in \mathrm{Aut}(P)\) such that
$$
m_2 \;=\; m_1 \circ \varphi .
$$
Equivalently, \(m_1\) and \(m_2\) differ only by a permutation of the pattern nodes.


<details> <summary><b>Solution:</b></summary>

```python
from collections import defaultdict
from typing import Dict, Iterable, List, Tuple, Optional, Union

def dedup_by_pattern_image_with_orbits(
    matches: Iterable[Dict[int,int]],
    *,
    pattern_autos: Iterable[Dict[int,int]],
    pattern_node_order: Optional[Iterable[int]] = None,
    return_groups: bool = False,
) -> Union[List[Dict[int,int]], Dict[Tuple[int,...], List[Dict[int,int]]]]:
    matches = list(matches)
    if not matches:
        return {} if return_groups else []

    # deterministic pattern node order
    if pattern_node_order is None:
        pattern_node_order = tuple(sorted(matches[0].keys()))
    else:
        pattern_node_order = tuple(pattern_node_order)

    # pre-list autos for repeated use
    autos = list(pattern_autos)

    # helper: canonical tuple for mapping m under all pattern automorphisms
    def canonical_tuple_for_mapping(m: Dict[int,int]) -> Tuple[int, ...]:
        best = None
        for phi in autos:
            tup = tuple(m[phi[p]] for p in pattern_node_order)
            if best is None or tup < best:
                best = tup
        return best

    # bucket by canonical tuple
    buckets: Dict[Tuple[int,...], List[Dict[int,int]]] = defaultdict(list)
    for m in matches:
        key = canonical_tuple_for_mapping(m)
        buckets[key].append(m)

    # deterministic ordering of keys
    ordered_keys = sorted(buckets.keys())

    if return_groups:
        return {k: buckets[k] for k in ordered_keys}

    # select deterministic representative per bucket:
    # pick mapping with smallest tuple according to pattern_node_order
    def ordering_tuple(m: Dict[int,int]) -> Tuple[int, ...]:
        return tuple(m[p] for p in pattern_node_order)

    representatives: List[Dict[int,int]] = []
    for k in ordered_keys:
        group = buckets[k]
        rep = min(group, key=ordering_tuple)
        representatives.append(rep)

    return representatives

```


### Q2 — Detect ethanol as a subgraph in a molecule table

**Goal.**  
Given a table of molecules represented by SMILES strings, determine which
molecules contain **ethanol** as a subgraph, and store the result in a new
boolean column `is_subgraph_ethanol`.

This exercise connects **subgraph isomorphism** with a practical,
dataframe-based workflow.

---

**Step-by-step tasks**

1. **Convert SMILES to RDKit molecules**  
   Create a column `mol` by parsing the SMILES strings with RDKit.

2. **Convert molecules to typed graphs**  
   Create a column `graph` by converting each RDKit molecule to a
   NetworkX graph using `mol_to_graph`.

3. **Select ethanol as the pattern graph**  
   Use the graph of ethanol (row 0) as the subgraph pattern.

4. **Test subgraph containment**  
   For each molecule, check whether the ethanol graph is a subgraph of
   the molecule graph using `nx_subgraph_matches`.

5. **Store the result**  
   Add a boolean column `is_subgraph_ethanol` indicating whether at least
   one subgraph embedding exists.


<details> <summary><b>Solution:</b></summary>

```python
from rdkit import Chem

# 1. Convert SMILES to RDKit molecules
df["mol"] = df["smiles"].apply(Chem.MolFromSmiles)

# 2. Convert molecules to typed NetworkX graphs
df["graph"] = df["mol"].apply(mol_to_graph)

# 3. Select ethanol as the pattern graph
ethanol_G = df.loc[0, "graph"]

# 4. Subgraph test: ethanol ⊆ host
def has_ethanol_subgraph(host_G):
    matches = nx_subgraph_matches(
        host_G=host_G,
        pattern_G=ethanol_G,
        invert=True,  # pattern_idx -> host_idx
    )
    return len(matches) > 0

# 5. Add boolean result column
df["is_subgraph_ethanol"] = df["graph"].apply(has_ethanol_subgraph)

```


## 3. RDKit subgraph search


RDKit performs pattern-in-host queries via **substructure matching**:

```python
host_mol.GetSubstructMatches(pattern_mol)
```
Each returned match is a mapping **pattern atom index → host atom index**.

**Practical limitations:**
- Atom indices are fixed by RDKit’s internal molecule ordering and cannot be
  freely relabeled or manipulated.
- Node/edge matching criteria are not user-definable; they are hard-coded into
  RDKit’s chemistry model.
- Automorphisms and orbit structure are not exposed, making symmetry handling implicit.
- Fine-grained control over which atom or bond attributes participate in matching
  is limited compared to explicit graph-based approaches.

As a result, RDKit substructure matching is efficient and chemically robust,
but less flexible for algorithmic experimentation and symmetry-aware workflows.



### Q3 — Detect ethanol subgraphs using RDKit and compare with NetworkX

**Goal.**  
Repeat **Q2** using **RDKit substructure matching** instead of NetworkX,
and compare the resulting boolean column with the NetworkX-based result.

This exercise highlights practical differences between:
- explicit graph-based subgraph isomorphism (NetworkX), and
- chemistry-native substructure matching (RDKit).

---

**Step-by-step tasks**

0. **Reuse the RDKit molecule column**  
   Use the existing `mol` column created from SMILES.

1. **Select ethanol as the pattern molecule**  
   Use the RDKit molecule corresponding to ethanol (row 0).

2. **Test RDKit substructure containment**  
   For each molecule, check whether RDKit finds at least one substructure
   match of ethanol using `GetSubstructMatches`.

3. **Store the result**  
   Add a boolean column `is_subgraph_ethanol_rdkit`.

4. **Compare with the NetworkX result**  
   Identify molecules where RDKit and NetworkX disagree.



<details> <summary><b>Solution:</b></summary>


```python
from rdkit import Chem

# 1. Select ethanol as RDKit pattern molecule
ethanol_mol = df.loc[0, "mol"]

# 2. RDKit-based subgraph test
def has_ethanol_subgraph_rdkit(host_mol):
    matches = rdkit_subgraph_matches(
        host_mol=host_mol,
        pattern_mol=ethanol_mol,
        uniquify=True,
    )
    return len(matches) > 0

# 3. Add RDKit result column
df["is_subgraph_ethanol_rdkit"] = df["mol"].apply(has_ethanol_subgraph_rdkit)

# 4. Comparison
comparison = df[
    df["is_subgraph_ethanol"] != df["is_subgraph_ethanol_rdkit"]
][["smiles", "name", "is_subgraph_ethanol", "is_subgraph_ethanol_rdkit"]]

comparison
```


## 4. MCS (Maximum Common Substructure) — RDKit + NetworkX views

The **Maximum Common Substructure (MCS)** problem asks for the largest subgraph that two molecular graphs share.
It is a core *alignment primitive* used in similarity, substructure transfer, and as a common starting point for **atom mapping**.

---

### 4.1 Formal definition (typed molecular graphs)

Let
$$
G = (V_G, E_G, a_G, b_G)
\quad \text{and} \quad
H = (V_H, E_H, a_H, b_H)
$$
be **typed molecular graphs**, where
\( a(\cdot) \) assigns atom labels (element, charge, aromaticity) and
\( b(\cdot) \) assigns bond labels (bond order, aromaticity).

A graph \( S \) is a **common subgraph** of \( G \) and \( H \) if there exist
**injective typed morphisms** (subgraph embeddings)


$$
f: S \hookrightarrow G, \qquad g: S \hookrightarrow H
$$

such that labels and bonds are preserved (same predicates as subgraph isomorphism).

An **MCS** is any common subgraph \(S^\*\) that maximizes a size objective:

$$
S^* \in \arg\max_{S}
\left(
|V(S)| \;\text{or}\; w_V |V(S)| + w_E |E(S)|
\right)
\quad
\text{subject to } S \hookrightarrow G \text{ and } S \hookrightarrow H.
$$


**Notes.**
- The MCS need not be unique (multiple maximum solutions may exist).
- Constraints (ring-only matching, atom types, bond types, chirality) change the feasible set and thus the MCS.

---

### 4.2 RDKit view: MCS as a SMARTS pattern + match lists

RDKit provides a practical MCS solver:

- `rdFMCS.FindMCS([mol1, mol2], ...)` returns an object whose `smartsString`
  encodes a **query substructure** \(Q\) (SMARTS) intended to represent a largest shared substructure under constraints.

We then compute embeddings of the SMARTS query into each molecule:

$$
\mathrm{Match}_G(Q) = \{\, f_i : V(Q)\hookrightarrow V(G)\,\}, \qquad
\mathrm{Match}_H(Q) = \{\, g_j : V(Q)\hookrightarrow V(H)\,\}.
$$

In practice:
- `mol.GetSubstructMatches(query)` returns many embeddings (symmetry variants).
- When multiple matchings exist, choosing a canonical representative (or a best-scoring one) is a separate step.

**Interpretation.**  
RDKit’s MCS output gives you:
1) a *candidate* maximum common substructure (as a query), and  
2) potentially many embeddings in each molecule.

---

### 4.3 NetworkX view: MCS as a maximum common *subgraph isomorphism*

When we convert molecules to NetworkX graphs, MCS corresponds to finding a largest typed graph \(S\)
that is simultaneously subgraph-isomorphic to both graphs.

Conceptually:

$$
S \subseteq G,\; S \subseteq H
\quad\Longleftrightarrow\quad
\exists\, f: S\hookrightarrow G,\; \exists\, g: S\hookrightarrow H.
$$

In the NX world, this is typically attacked by:
- searching over candidate node/bond subsets, or
- using MCS heuristics (e.g., expand from seeds; branch-and-bound; constraint propagation),
because exact MCS is NP-hard.

**Why still use NX here?**
- You can enforce *your* exact chemical typing predicates (`node_match`, `edge_match`).
- You can integrate symmetry handling (automorphism orbits) and custom constraints.
- You can expose intermediate states for teaching (what gets pruned, what expands).

---

### 4.4 Symmetry and non-uniqueness (important for mapping)

Both RDKit and NetworkX share the same caveats:

- **Multiple maximum solutions:** the same maximum subgraph size can be achieved by different subgraphs.
- **Many embeddings:** even one fixed maximum subgraph can have many matches due to molecular symmetry.
- **Implication:** MCS is a useful starting point for atom mapping, but not sufficient on its own — a tie-breaking or scoring rule is still required.

---

### 4.5 Practical pipeline (RDKit ↔ NX)

A common didactic workflow is:

1. **RDKit MCS**
   - compute a SMARTS query \(Q\) using `rdFMCS.FindMCS`.
2. **RDKit embeddings**
   - enumerate `GetSubstructMatches(Q)` for each molecule.
3. **NX analysis**
   - convert selected embeddings to NX node correspondences,
   - optionally deduplicate symmetry-equivalent matches using orbits,
   - use the resulting correspondence as a seed for atom mapping or reaction-center extraction.

This makes MCS an excellent bridge between *chemistry-native toolkits* (RDKit) and *graph-theoretic control* (NetworkX).


### Q4 — Compute the MCS between ethanol and another molecule

**Goal.**  
Use **RDKit MCS** to compute the maximum common substructure between **ethanol**
and each molecule in the table, and record the **size of the MCS**.

This exercise shows how RDKit abstracts chemical similarity beyond strict
subgraph containment.

---

**Step-by-step tasks**

1. **Select ethanol as the reference molecule**  
   Use the RDKit molecule corresponding to ethanol (row 0).

2. **Compute RDKit MCS**  
   For each molecule, compute the MCS SMARTS between ethanol and the molecule
   using `rdkit_mcs_smarts`.

3. **Convert SMARTS to a molecule**  
   Turn the SMARTS into an RDKit query molecule.

4. **Measure MCS size**  
   Record the number of atoms in the MCS (or `0` if no MCS is found).

5. **Store the result**  
   Add a column `mcs_size_ethanol_rdkit`.



<details> <summary><b>Solution:</b></summary>

```python
ethanol_mol = df.loc[0, "mol"]

def rdkit_mcs_size_with_ethanol(mol):
    smarts = rdkit_mcs_smarts([ethanol_mol, mol])
    if not smarts:
        return 0
    mcs_mol = Chem.MolFromSmarts(smarts)
    if mcs_mol is None:
        return 0
    return mcs_mol.GetNumAtoms()

df["mcs_size_ethanol_rdkit"] = df["mol"].apply(rdkit_mcs_size_with_ethanol)
df
```


### Q5 — Compute the MCS between ethanol and another molecule (NetworkX view)

**Goal.**  
Compute the **Maximum Common Substructure (MCS)** between **ethanol** and each
molecule in the table using the **NetworkX-based MCS implementation**, and
record the size of the resulting common subgraph.

This exercise emphasizes a **strict graph-theoretic view** of MCS, in contrast
to the chemistry-aware RDKit approach in Q4.

---

**Step-by-step tasks**

1. **Select ethanol as the reference graph**  
   Use the NetworkX graph corresponding to ethanol (row 0).

2. **Compute NetworkX MCS**  
   For each molecule graph, compute an MCS using `nx_mcs`.

3. **Measure MCS size**  
   Record the number of nodes (atoms) in the resulting MCS graph.

4. **Store the result**  
   Add a new column `mcs_size_ethanol_nx` to the dataframe.



<details> <summary><b>Solution:</b></summary>

```python
# Reference graph: ethanol
ethanol_G = df.loc[0, "graph"]

def nx_mcs_size_with_ethanol(host_G):
    S = nx_mcs(
        ethanol_G,
        host_G,
        node_match=node_match,
        edge_match=edge_match,
    )
    if S is None:
        return 0
    return S.number_of_nodes()

# Add NetworkX MCS size column
df["mcs_size_ethanol_nx"] = df["graph"].apply(nx_mcs_size_with_ethanol)
```


## 5. Discussion

- **Subgraph isomorphism** is the core operation for rule application later in SynEdu.
- **MCS** is a chemistry-aware alignment primitive, but it is heuristic and sometimes non-unique; always log settings and timeouts.
- **RDKit vs NetworkX**:
  - RDKit: chemistry-aware SMARTS matching, built-in `uniquify`.
  - NetworkX: full control over attributes and morphism semantics; you manage deduplication and interpretation.


## 6. Quiz

Answer the following questions using **both chemical intuition and formal graph language**.

---

### 1. Deduplicating subgraph matches  
Subgraph matching often returns many equivalent matches.

- Explain how using **sets of host atom indices** can be used to
  **deduplicate** equivalent matches.
- Why does this work even when atom ordering differs?

---

### 2. RDKit vs NetworkX (practice vs theory)  
Give **one practical advantage** of each approach:

- **RDKit** substructure matching  
- **NetworkX** graph matching  

In which situations would you prefer one over the othe


## 7. References and further reading

- RDKit documentation: https://www.rdkit.org/docs/  
- RDKit Book: https://www.rdkit.org/docs/Book.html  
- NetworkX documentation: https://networkx.org/documentation/stable/  
- NetworkX isomorphism: https://networkx.org/documentation/stable/reference/algorithms/isomorphism.html  
- RDKit MCS (rdFMCS): https://www.rdkit.org/docs/source/rdkit.Chem.rdFMCS.html
