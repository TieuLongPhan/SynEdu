# S01 · From Molecules to Typed Graphs: Fundamentals

<div class="alert alert-block alert-info">
<b>Welcome to SynEdu.</b><br>
This talktorial is part of <b>SynEdu</b>, a lightweight, research-oriented teaching series built around the
<b>Syn</b> ecosystem and <b>RDKit</b> for practical, reproducible chemical and reaction modeling.
</div>

<div class="alert alert-block alert-success">
<b>What you will gain.</b><br>
You will learn how to represent molecules as <b>typed graphs</b> and how to reason about
structure-preserving matches that underpin graph-based reaction modeling: <b>graph morphisms</b>,
<b>isomorphism</b>, and <b>automorphisms</b>.
</div>

<div class="alert alert-block alert-warning">
<b>Reproducibility note.</b><br>
This notebook is self-contained and meant to be run top-to-bottom. If you are following the full SynEdu series,
continue with <b>S02</b> for <i>subgraph isomorphism</i> and <i>MCS</i>.
</div>


## Aim of this talktorial

Reaction modeling via **graph transformation** (e.g. DPO rules) relies on two pillars:

1. **Representation** — converting molecules into graphs with chemically meaningful labels (typed graphs).
2. **Matching** — reliably identifying when two graphs (or parts of graphs) correspond, via **typed graph morphisms**:
   - isomorphisms,
   - automorphisms,
   - subgraph isomorphisms,
   - MCS-based alignments.

This talktorial establishes both pillars with minimal, transparent implementations using:

- **RDKit** as the chemical ground truth (sanitization, aromaticity, valence, MCS),
- **NetworkX** as the generic graph engine for matching, automorphism analysis, and later rewriting.

**Data example:** `data/molecules.csv`  

---

## Learning outcomes

After completing this talktorial, you will be able to:

- Convert SMILES into a **typed NetworkX molecular graph** (atoms → nodes, bonds → edges).
- Perform a **round-trip**: RDKit → NetworkX → RDKit and check what is preserved.
- State and use the formal definition of a **typed graph morphism** and its special cases:
  - homomorphism, monomorphism, isomorphism, automorphism.
- Distinguish and apply:
  - **graph isomorphism** (same structure under relabeling),
  - **automorphisms** (symmetries of one graph),
  - **subgraph isomorphism** (pattern inside host),
  - **MCS** (maximum common substructure) as a practical alignment primitive.
- Compare **RDKit morphisms** (SMARTS-based, chemistry-aware) to **NetworkX morphisms** (pure structure+attributes).
- Understand why **symmetry** and **attribute choices** strongly influence rule extraction and rule application later in SynEdu.

---

## Outline

0. **Setup & data**
1. **Theory: typed graphs and morphisms**
2. **RDKit → NetworkX typed graphs**
3. **Round-trip: NetworkX → RDKit and sanity checks**
4. **Isomorphism & automorphisms**
5. **Subgraph isomorphism (pattern → host)**
6. **MCS alignment with RDKit**
7. **RDKit vs NetworkX morphisms: comparison**
8. **Discussion, quiz, and references**




## 0. Setup & Data


## 1. Theory: Typed Graphs and Morphisms

### 1.1 Typed molecular graphs

In computational reaction modeling, we represent molecules as **typed graphs** so that “matching” respects
chemical identity (elements, charges, bond orders), not just connectivity.

A **typed graph** is a quadruple

$$
G = (V, E, \tau_V, \tau_E),
$$

where:

- **Vertices** $V$ represent **atoms**.
- **Edges** $E \subseteq \{\{u,v\}\mid u,v\in V,\ u\neq v\}$ represent **bonds** (finite, undirected, simple: no loops, no parallel edges).
- $\tau_V: V \to \mathcal{A}_V$ assigns **atom attributes** (chemical labels).
- $\tau_E: E \to \mathcal{A}_E$ assigns **bond attributes** (chemical labels).

We often write $V(G)$ and $E(G)$ for the vertex and edge sets of $G$. For a vertex $v\in V(G)$:

- neighbourhood:
  $$
  N_G(v)=\{w\in V(G)\mid vw\in E(G)\},
  $$
- degree:
  $$
  \deg_G(v)=|N_G(v)|.
  $$

#### Labels (typed graphs)

“Types” are encoded as labelling maps

$$
\ell_V: V(G)\to L_V,\qquad \ell_E: E(G)\to L_E,
$$

where $L_V$ and $L_E$ are finite, non-empty label sets.
For molecular graphs we use the chemistry-specific notation:

$$
a_G: V(G)\to L_V \quad\text{(atom labels)},\qquad
b_G: E(G)\to L_E \quad\text{(bond labels)}.
$$

Let $\mathcal{G}$ denote the class of all labelled molecular graphs equipped with $(a_G,b_G)$.
In chemistry, $a_G(v)$ encodes *what atom this is* (element, charge, aromaticity, …), 
while $b_G(uv)$ encodes *what bond this is* (order, aromaticity, ring status, …).

---

### 1.2 Graph morphisms

Let $G,H\in\mathcal{G}$ be labelled molecular graphs with atom- and bond-labelling functions
$(a_G,b_G)$ and $(a_H,b_H)$.

A **(labelled) graph morphism** from $G$ to $H$ is a map

$$
\varphi:V(G)\to V(H)
$$

that preserves **chemical identity at atoms**, and preserves **bonds and their types**.
Formally, $\varphi$ must satisfy:

#### (M1) Atom-label preservation
For every atom $v\in V(G)$,

$$
a_H(\varphi(v)) = a_G(v).
$$

> **Chemistry meaning**  
> $\varphi$ never maps a carbon to a nitrogen, or a neutral atom to a charged atom, if those are encoded in $a_\cdot$.

#### (M2) Adjacency (bond existence) preservation
For every bond $uv\in E(G)$,

$$
\varphi(u)\varphi(v)\in E(H).
$$

> **Chemistry meaning**  
> Bonded atoms in $G$ must map to bonded atoms in $H$ (no bond can “disappear” under the map).

#### (M3) Bond-label preservation
For every bond $uv\in E(G)$,

$$
b_H\!\big(\varphi(u)\varphi(v)\big) = b_G(uv).
$$

> **Chemistry meaning**  
> If $uv$ is a double bond in $G$, its image must be a double bond in $H$ (and likewise for aromaticity if included).

> **Summary**  
> A morphism $\varphi$ is a label-preserving embedding of the local chemical graph of $G$ into $H$.  
> It preserves “what atoms are” and “how they are connected”.

---

#### Compatibility predicates (practical generalization)

In practice, chemoinformatics representations may differ (e.g. aromatic vs Kekulé, resonance conventions).
We therefore sometimes relax strict equality into **compatibility**:

- atom compatibility:
  $$
  a_H(\varphi(v)) \sim_V a_G(v),
  $$
- bond compatibility:
  $$
  b_H(\varphi(u)\varphi(v)) \sim_E b_G(uv),
  $$

where $\sim_V$ and $\sim_E$ encode allowed correspondences (e.g. “aromatic bond” compatible with alternating single/double under a chosen model).

> **For chemists**  
> This is where you decide whether two representations should be considered “the same chemistry”.  
> For example, strict matching distinguishes aromatic vs Kekulé; compatibility matching can treat them as equivalent.

---

#### Standard special cases (matching tasks)

- **Monomorphism (injective morphism)**  
  $\varphi$ is injective (distinct atoms in $G$ map to distinct atoms in $H$).  
  $\rightarrow$ This is the formal object behind **substructure search / subgraph isomorphism**.

- **Isomorphism**  
  $\varphi$ is bijective and $\varphi^{-1}$ is also a morphism.  
  $\rightarrow$ We write $G\simeq H$ (same molecule up to relabeling).

- **Automorphism**  
  An isomorphism $\varphi:G\to G$.  
  $\rightarrow$ Encodes **molecular symmetry**, which can multiply equivalent matches and motivates deduplication.

> **Interpretation for reaction rules**  
> Applying a graph-rewrite rule begins by finding a **monomorphism** from the rule LHS into the host molecule.  
> Automorphisms (symmetry) can generate many equivalent embeddings; later SynEdu notebooks introduce symmetry-aware deduplication.


## 2. RDKit ⇄ NetworkX: Typed Molecular Graphs

In SynEdu, **RDKit** and **NetworkX** play complementary roles:

- **RDKit** is the chemical authority: sanitization, valence rules, aromaticity perception, and canonicalization.
- **NetworkX** provides an explicit, inspectable graph representation used for matching, symmetry analysis,
  and later graph rewriting.

To ensure that graph-based operations remain chemically meaningful, we require a **reversible interface**
between the two representations.


### Exercise: Round-trip accuracy (RDKit ⇄ NetworkX)

The goal of this exercise is to verify that converting

RDKit → NetworkX → RDKit

preserves the **chemical information we care about**.

You should treat the two functions provided above as a black box.

---

#### Q1 — Heavy-atom SMILES invariance

Write a function `roundtrip_smiles_equal(smiles)` that:

1. parses a SMILES string into an RDKit molecule,
2. converts it to a typed graph using `mol_to_graph`,
3. reconstructs a molecule using `graph_to_mol`,
4. compares the **canonical heavy-atom SMILES** of the original and reconstructed molecules.

The function should return `True` if the two SMILES are identical, and `False` otherwise.

---

#### Q2 — Count invariants

Extend your check in **Q1** to also verify that:

- the number of **heavy atoms** is preserved,
- the number of **heavy-atom bonds** is preserved.

Return `True` only if *all* invariants are satisfied.

> Hint: use `Chem.RemoveHs(mol)` before counting atoms or bonds.

---



<details>
<summary><b>Solution:</b></summary>

### Q1–Q2: Round-trip checker (heavy SMILES + count invariants)

```python
from rdkit import Chem
from rdkit.Chem import rdmolops  # optional: useful for extra invariants

def canonical_heavy_smiles(m: Chem.Mol) -> str:
    """Return canonical SMILES after removing H (heavy-atom skeleton)."""
    return Chem.MolToSmiles(Chem.RemoveHs(m), canonical=True)

def heavy_counts(m: Chem.Mol) -> tuple[int, int]:
    """Return (n_heavy_atoms, n_heavy_bonds) after removing H."""
    mh = Chem.RemoveHs(m)
    return mh.GetNumAtoms(), mh.GetNumBonds()

def roundtrip_ok(smiles: str, verbose: bool = True) -> bool:
    """
    RDKit → NetworkX → RDKit round-trip check.

    Criteria:
    1) canonical heavy-atom SMILES preserved
    2) heavy atom count preserved
    3) heavy bond count preserved
    """
    m1 = Chem.MolFromSmiles(smiles)
    if m1 is None:
        if verbose:
            print("Parse failed:", smiles)
        return False

    G = mol_to_graph(m1, include_implicit_h=True)
    m2 = graph_to_mol(G, make_explicit_h=False)

    s1, s2 = canonical_heavy_smiles(m1), canonical_heavy_smiles(m2)
    c1, c2 = heavy_counts(m1), heavy_counts(m2)

    ok = (s1 == s2) and (c1 == c2)

    if verbose and not ok:
        print("FAIL:", smiles)
        print(" heavy SMILES:", s1, "vs", s2)
        print(" counts:", c1, "vs", c2)

    return ok
```

### Quick test (run on a small subset)

```python
n_test = min(50, len(df))
fails = []

for s in df["smiles"].head(n_test):
    if not roundtrip_ok(s, verbose=False):
        fails.append(s)

print("Checked:", n_test)
print("Failures:", len(fails))
if fails:
    print("Example failures:", fails[:5])

# Optional: inspect one failure in detail
if fails:
    _ = roundtrip_ok(fails[0], verbose=True)
```



## 3. Isomorphism

We connect **formal graph-morphism definitions** to concrete
`networkx` matchers used in practice.

### Attribute predicates

We denote the **vertex** and **edge** attribute predicates as:

$$
\Phi_V : V(G) \times V(H) \rightarrow \{\text{true}, \text{false}\}
$$

$$
\Phi_E : E(G) \times E(H) \rightarrow \{\text{true}, \text{false}\}
$$

In code, these predicates are implemented as Python functions:

- `node_match(n1, n2)` — compares atom (vertex) attributes
- `edge_match(e1, e2)` — compares bond (edge) attributes

### Compatibility used in S01

For **S01**, we adopt a *strict-but-minimal* chemical compatibility model:

- same atom `symbol` (element),
- same `formal_charge`,
- same `aromatic` flag,
- same bond `order`.



### Exercise: Isomorphism
**Q3 — Fix the matcher**

Implement `node_match` that requires matching `symbol` **and** either `total_h` or `formal_charge` (or both). Replace the existing `node_match` with your function and re-run the demo so that:

- `benzene` still matches, and  
- `aniline` (`c1ccccc1N`) **does not** match `anilinium` (`c1ccccc1[NH3+]`).

> Hint: `mol_to_graph(..., include_implicit_h=True)` stores H as `total_h`. Use `n.get("total_h",0)` or `n.get("formal_charge",0)`.




<details> <summary><b>Solution:</b></summary>

```python
# Solution: enhanced matcher that checks symbol + (total_h OR formal_charge)
def enhanced_node_match(n1, n2):
    return (
        n1.get("symbol") == n2.get("symbol")
        and (
            int(n1.get("total_h", 0)) == int(n2.get("total_h", 0))
            or int(n1.get("formal_charge", 0)) == int(n2.get("formal_charge", 0))
        )
    )

# run the demo with the enhanced matcher (uses existing `pairs`, `graphs`, `edge_match`, `iso_and_count`)
print("=== enhanced matcher (symbol + total_h/charge) ===")
for name in pairs:
    G1 = graphs[f"{name}_a"]; G2 = graphs[f"{name}_b"]
    iso_flag, n_maps = iso_and_count(G1, G2, enhanced_node_match, edge_match)
    print(f"{name:8} | isomorphic: {int(iso_flag):1d} | mappings: {n_maps}")

# quick instructor checks
assert iso_and_count(graphs["benzene_a"], graphs["benzene_b"], enhanced_node_match, edge_match)[0]
assert not iso_and_count(graphs["aniline_a"], graphs["aniline_b"], enhanced_node_match, edge_match)[0]
```
<details>


## 4. Automorphisms & orbits

**Observation.** In the benzene example you enumerated **12 mappings** — these are the automorphisms of the benzene heavy-atom graph (the dihedral group \(D_6\), where \(|D_6| = 12\)).

**Definition.** An automorphism is a graph isomorphism from the graph to itself:

$$
f : G \longrightarrow G.
$$

The automorphism group is

$$
\mathrm{Aut}(G).
$$

The **orbit** of a vertex \(v\) is the set of images of \(v\) under all automorphisms:

$$
\mathrm{Orbit}(v)=\{\psi(v)\;|\;\psi\in\mathrm{Aut}(G)\}.
$$

**Facts.** For benzene:

$$
|\mathrm{Aut}(G)| = |D_6| = 12,
$$

and all six carbon atoms lie in a single orbit.

**Why it matters.** Symmetric hosts produce many equivalent embeddings → duplicate matches and wasted work.

**Simple remedies.**
- Deduplicate by host-atom set: use `frozenset(mapping.values())`.  
- Use orbit representatives (e.g. choose the $\min$ index per orbit).  
- Accept only a canonical mapping (WL/lexicographic tie-break).

**Practical tips.**
- Include chemical attributes (`total_h`, `formal_charge`, stereochemistry) in matchers to reduce false symmetry.  
- Pre-filter with cheap signatures (degree, label counts, WL hashes) before enumerating automorphisms.




### Exercise: Automorphisms of a Molecular Graph

### Q4 — Develop a function `enumerate_automorphisms` to enumerate all automorphisms of a graph

**Hint:** An automorphism of a graph \(G\) is an isomorphism from \(G\) to itself.

Equivalently, the automorphism group satisfies

$$
\mathrm{Aut}(G) \subseteq \mathrm{Iso}(G, G).
$$




<details> <summary><b>Solution:</b></summary>

```python
def enumerate_automorphisms(G: nx.Graph):
    GM_self = iso.GraphMatcher(G, G, node_match=node_match, edge_match=edge_match)
    return list(GM_self.isomorphisms_iter())
```


We can now analyse the symmetry of a molecular graph by computing the
**orbits induced by its automorphism group**.

Under the natural action of the automorphism group on the vertex set,
two vertices belong to the same orbit if there exists an automorphism
mapping one to the other.

$$
\text{For } u, v \in V(G), \quad
u \sim v
\;\Longleftrightarrow\;
\exists\, \varphi \in \mathrm{Aut}(G)
\text{ such that }
\varphi(u) = v .
$$

Each orbit therefore represents a set of **symmetry-equivalent atoms**.



## 5. Discussion
- A **typed graph morphism** formalizes structure- and attribute-preserving maps between graphs.
- Our **typed molecular graphs** use a minimal attribute schema (`symbol`, `formal_charge`, `aromatic`, `order`) to define what “same” means.
- **Round-trip conversion** (RDKit → NetworkX → RDKit) is valuable for debugging and peer review; we preserve heavy-atom topology, but exact RDKit internal state may differ.
- **Automorphisms** describe symmetries; they inflate match enumeration. Deduplicate (e.g. by host-atom set) to prevent combinatorial explosion.


## 6. Quiz · Graph Matching Fundamentals

Answer the following questions using **both chemical intuition and formal graph language**.

---

### 1. Typed graph morphism

In one or two sentences, define a **typed graph morphism**.

$$
f : G \rightarrow H
$$

- What objects does f map?
- Which **atom** and **bond** properties must be preserved?
- Give one example of a mapping that would be **invalid** in chemistry.


---

### 2. Isomorphism  
What **additional requirement** must a graph morphism satisfy to become an  
**isomorphism**?

- How does this relate to the idea of *two molecules having the same structure*?

---

### 3. Automorphism and symmetry  
What is an **automorphism** of a molecular graph?

- Why do symmetric molecules (e.g. benzene) have **many automorphisms**?
- Why do automorphisms cause **duplicate subgraph matches** during matching?

---

### 4. Deduplicating subgraph matches  
Subgraph matching often returns many equivalent matches.

- Explain how using **sets of host atom indices** can be used to
  **deduplicate** equivalent matches.
- Why does this work even when atom ordering differs?

---

### 5. RDKit vs NetworkX (practice vs theory)  
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
