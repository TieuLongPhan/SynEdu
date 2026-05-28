# S03: Maximum Common Substructure in Reaction Informatics

This talktorial introduces maximum common substructure (MCS) as a practical molecular-alignment tool for reaction informatics. We compare chemistry-aware RDKit MCS with explicit graph-level reasoning, then use the same idea for reaction rebalancing [\[1\]](#6.-References), [\[2\]](#6.-References), [\[3\]](#6.-References).



## Aim of this talktorial

1. Use **MCS-based molecular alignment** to compare molecules that differ by edits, substitutions, or missing components.
2. Interpret where RDKit MCS and graph-based MCS agree or differ.
3. Apply MCS reasoning to simple **reaction rebalancing** examples and connect it to SynRBL.

---

## Learning outcomes

After completing this talktorial, you will be able to:

- define **maximum common substructure** as a shared labeled subgraph,
- compute and interpret MCS results with **RDKit**,
- convert an MCS into an atom-level alignment between molecules,
- explain why MCS results may be non-unique or sensitive to matching constraints,
- use MCS-based alignment to reason about partial or noisy reactions, and
- apply MCS ideas to simple **reaction rebalancing** examples.

---

## Outline

- [0. Setup & Data](#0.-Setup-&-Data)
- [1. Maximum Common Substructure](#1.-Maximum-Common-Substructure)
- [2. Rule-based reaction rebalancing](#2.-Rule-based-reaction-rebalancing)
- [3. MCS-based reaction rebalancing](#3.-MCS-based-reaction-rebalancing)
- [4. Discussion](#4.-Discussion)
- [5. Quiz](#5.-Quiz)
- [6. References](#6.-References)


## 0. Setup & Data


## 1. Maximum Common Substructure

The **Maximum Common Substructure (MCS)** [\[2\]](#6.-References), [\[4\]](#6.-References) problem seeks the largest subgraph
shared by two molecular graphs. It serves as a core **alignment primitive** for
molecular similarity, substructure transfer, and as a common starting point for
**atom mapping**.

---

### 1.1 Theory

Let
$$
G = \bigl( V_G, E_G, a_G, b_G \bigr), \qquad
H = \bigl( V_H, E_H, a_H, b_H \bigr)
$$
be **labeled molecular graphs**, where
$$
a(\cdot) \text{ denotes atom labels}, \qquad
b(\cdot) \text{ denotes bond labels}.
$$

A graph
$$
S = (V_S, E_S, a_S, b_S)
$$
is a **common subgraph** of \(G\) and \(H\) if there exist **injective graph
morphisms**
$$
f : S \hookrightarrow G, \qquad
g : S \hookrightarrow H,
$$
that preserve adjacency and all labels.

An **MCS** is a common subgraph
$$
S^\star
$$
that maximizes the number of vertices, i.e.
$$
|V(S^\star)| \;\ge\; |V(S)| \qquad
\forall\, S \text{ such that } S \subseteq G \text{ and } S \subseteq H .
$$

**Practical notes.**
- MCSs are generally **not unique**.
- Imposing chemical constraints (atom and bond types, aromaticity, ring
  membership, chirality) restricts the admissible morphisms and can change the
  resulting MCS.


**Definition (Maximum Common Substructure / MCS).**  
Given two labeled graphs $G_1$ and $G_2$, a *common subgraph* is a graph $C$ together with label-preserving injective morphisms $\iota_1: C \hookrightarrow G_1$ and $\iota_2: C \hookrightarrow G_2$.  
The *maximum common subgraph* (MCS) is a common subgraph $C^*$ of maximum cardinality $|V(C^*)|$ (maximum common substructure by atoms) or maximum $|E(C^*)|$ (by bonds).

**Remark.** MCS computation is NP-hard in general [\[4\]](#6.-References). RDKit uses an FMCS heuristic that may return a *locally* maximal but not globally maximal common subgraph. NetworkX-based MCS uses exact backtracking, which is exponential in the worst case but exact.

**Definition (Normalized MCS size).**  
Given an MCS $C^*$ of two molecules, the *normalized MCS size* (Tanimoto-like) is:

$$
\tau_{\text{MCS}}(G_1, G_2) = \frac{|V(C^*)|}{|V(G_1)| + |V(G_2)| - |V(C^*)|}
$$

This is 1 if the molecules are isomorphic and 0 if they share no atom.


### 1.2. RDkit MCS

In **RDKit**, the **maximum common substructure (MCS)** can be computed using
`rdkit.Chem.rdFMCS` [\[5\]](#6.-References). 

RDKit MCS is commonly used for molecular alignment, reaction analysis, and as a
basis for atom mapping.


**Stored in `synedu.Utils`** — the MCS helpers are deposited for later reaction-alignment tasks:
```python
from synedu.Utils.mcs import rdkit_mcs, nx_mcs, mcs_size, mcs_class, map_edges
```
Later notebooks can reuse these functions without carrying the teaching implementation forward.


**Q1 — Compute the MCS between ethanol and another molecule**

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
   using `rdkit_mcs`.

3. **Convert SMARTS to a molecule**  
   Turn the SMARTS into an RDKit query molecule.

4. **Measure MCS size**  
   Record the number of atoms in the MCS (or `0` if no MCS is found).

5. **Store the result**  
   Add a column `mcs_size_ethanol_rdkit`.

<details class="synedu-solution">
<summary><strong>Solution</strong></summary>

```python
ethanol_mol = df.loc[0, "mol"]

def rdkit_mcs_size_with_ethanol(mol):
    smarts, _ = rdkit_mcs([ethanol_mol, mol])
    if not smarts:
        return 0
    mcs_mol = Chem.MolFromSmarts(smarts)
    if mcs_mol is None:
        return 0
    return mcs_mol.GetNumAtoms()

df["mcs_size_ethanol_rdkit"] = df["mol"].apply(rdkit_mcs_size_with_ethanol)
df
```

</details>



### 1.3. NetworkX MCS


### MCS alignment gallery

Each row shows a molecule pair. The **highlighted (orange) atoms and bonds** form the maximum common substructure. Larger MCS → more shared scaffold; smaller MCS → structurally divergent molecules despite similar element composition.


**Q2 — Compute the MCS between ethanol and another molecule**

Compute the **MCS** between **ethanol** and each
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

<details class="synedu-solution">
<summary><strong>Solution</strong></summary>

```python
# Reference graph: ethanol
df['graph'] = df['mol'].apply(mol_to_graph)
ethanol_G = df.loc[0, "graph"]

def nx_mcs_size_with_ethanol(host_G):
    S, _ = nx_mcs(
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

</details>



Now we inspect the difference between rdkit and networkx

```python
df[df['mcs_size_ethanol_rdkit'] != df['mcs_size_ethanol_nx']]
```


Most discrepancies originate from **ring systems**. RDKit MCS is *chemistry-first*, enforcing ring–chain consistency, whereas NetworkX MCS is *graph-first* and performs purely topological matching. With `RingMatchesRingOnly=True`, RDKit allows atom matches only when both the pattern and host atoms belong to rings, preventing acyclic–ring correspondences. Consequently, cyclic hydrocarbons yield zero MCS with ethanol under RDKit but non-zero matches under unconstrained graph matching.


### MCS size vs Tanimoto similarity

Tanimoto similarity (fingerprint-based) and MCS size (graph-based) both measure molecular overlap, but they are **not equivalent** [\[4\]](#6.-References). The scatter plot below samples molecules from the dataset and compares both metrics against ethanol. Ring-containing molecules can have high graph MCS under RDKit's chemistry-aware matcher yet low Tanimoto — or vice versa.


## 2. Rule-based reaction rebalancing

We utilize idea from *[Reaction rebalancing: a novel approach to curating reaction databases](https://link.springer.com/article/10.1186/s13321-024-00875-4)* [\[6\]](#6.-References). 

---

### 2.1 Reaction decomposition

Given a reaction SMILES
$$
r \;=\; R_1 + R_2 + \cdots \;\longrightarrow\; P_1 + P_2 + \cdots ,
$$
define the reactant and product multisets
$$
\mathcal{R}=\{R_i\}, \qquad \mathcal{P}=\{P_j\}.
$$

---

### 2.2 Formula representation

Let $\mathcal{E}$ be the set of chemical elements.
Each molecule $M$ is represented by its element-count vector
$$
\phi(M) \in \mathbb{N}^{|\mathcal{E}|}.
$$
Total formulas are
$$
\Phi_{\mathcal{R}}=\sum_{R_i\in\mathcal{R}}\phi(R_i), \qquad
\Phi_{\mathcal{P}}=\sum_{P_j\in\mathcal{P}}\phi(P_j).
$$

---

### 2.3 Imbalance detection

Define the imbalance vector
$$
\Delta \;=\; \Phi_{\mathcal{R}}-\Phi_{\mathcal{P}} \;\in\; \mathbb{Z}^{|\mathcal{E}|}.
$$
The reaction is balanced iff $\Delta=\mathbf{0}$.

---

### 2.4 Auxiliary-species imputation

Let $\mathcal{L}=\{A_k\}$ be a library of auxiliary species
(e.g. $\mathrm{H_2O}$) with formulas $\phi(A_k)$.
We seek coefficients $c_k\in\mathbb{N}$ such that
$$
\Delta = \sum_k c_k\,\phi(A_k).
$$

---

### 2.5 Reaction update

If a solution exists, rebalance the reaction by adding
$c_k A_k$ to the deficient side:
$$
r' =
\begin{cases}
\mathcal{R}\cup\{c_kA_k\} \rightarrow \mathcal{P}, & \Delta>0,\\
\mathcal{R}\rightarrow \mathcal{P}\cup\{c_kA_k\}, & \Delta<0.
\end{cases}
$$

In SMILES, this corresponds to appending the auxiliary species
(e.g. `O` for $\mathrm{H_2O}$).


### Rebalancing workflow

The five-stage pipeline below converts an imbalanced reaction SMILES into a balanced one. Library-based imputation handles small inorganic species; MCS-based imputation handles missing carbon fragments.


Now consider the following simple reaction. 
By comparing the total molecular formulas of the reactant and product sides,
we can immediately detect a **missing water molecule on the product side**.


To automate this imputation, we introduce a small set of utility functions
for parsing reactions, computing formula imbalances, and resolving them
using predefined auxiliary species.


**Stored in `synedu.Utils`** — reaction-balance helpers are now packaged for later curation and rule-extraction steps:
```python
from synedu.Utils.reaction import parse_reaction_smiles, reaction_imbalance
```
This keeps the formula-balance convention consistent across later tasks.


**Q3 — Auxiliary imputation feasibility**

You are given an element-wise imbalance vector  
$$
\Delta \in \mathbb{Z}^{|\mathcal E|}
$$
and an auxiliary species $A$ with molecular formula  
$$
\phi(A) \in \mathbb{N}^{|\mathcal E|}.
$$

The goal is to decide whether the imbalance can be resolved by adding
$c$ copies of $A$, i.e.
$$
\Delta = c \cdot \phi(A), \qquad c \in \mathbb{N}.
$$

---

**Hint:**

Think element-wise.  
If such a $c$ exists, then for every element $e$,
the ratio $\Delta_e / \phi_e(A)$ must be:
1. an integer,
2. the same for all elements,
3. and strictly positive.

Also ensure that no extra elements remain unexplained.



---

<details class="synedu-solution">
<summary><strong>Solution</strong></summary>

```python
from collections import Counter
from typing import Optional


def can_impute_with_aux(
    delta: Counter[str],
    aux_formula: Counter[str],
) -> Optional[int]:
    """
    Check whether Δ can be written as c · φ(A).
    """
    if not delta:
        return 0

    ratios = []

    for elem, count in aux_formula.items():
        if elem not in delta:
            return None
        if delta[elem] % count != 0:
            return None
        ratios.append(delta[elem] // count)

    if not ratios or len(set(ratios)) != 1:
        return None

    c = ratios[0]
    if c <= 0:
        return None

    remainder = delta - Counter({k: c * v for k, v in aux_formula.items()})
    if remainder:
        return None

    return c
```

</details>



**Stored in `synedu.Utils`** — auxiliary-species imputation is available as a reusable utility:
```python
from synedu.Utils.reaction import AUX_LIBRARY, can_impute_with_aux, impute_rsmi_rule
```
Later notebooks can focus on when to impute, rather than redefining the same bookkeeping.


Now try with the above example


### Before / after rebalancing

Imbalanced reactions are missing one or more species. The gallery below shows each reaction before and after library-based imputation — atoms and bonds are unchanged; only the missing auxiliary species (water, HCl) are appended to the correct side.


**Q4 — Rebalance some reactions**

Extend the auxiliary-species library and use it to rebalance these reactions:

```text
CC(C)O>>CC(C)Cl.O
CC(C)Cl>>CC(C)N.Cl
CC(=O)O>>CC(=O)OC.O
```

<details class="synedu-solution">
<summary><strong>Solution idea</strong></summary>

add common auxiliary species such as water, hydrochloric acid, ammonia, hydroxide, and methanol to the library. For each reaction, compute the imbalance with `reaction_imbalance`, call `impute_rsmi`, and keep the balanced reaction when imputation succeeds.

</details>



## 3. MCS-based reaction rebalancing

For **non-carbon compounds**, missing species can often be imputed from a
finite library. In contrast, for **carbon-containing compounds**, the space of
possible structures grows exponentially, making library-based enumeration
impractical at scale.

Instead, we use **reaction alignment** to identify missing **structural motifs**
via MCS and apply **rule-based reasoning** to complete and rebalance reactions.

In this section, we directly use functionality from **SynKit** [\[6\]](#6.-References)
(`pip install synkit`) and explore an established tool for reaction rebalancing,
**SynRBL** [\[7\]](#6.-References) (`pip install synrbl`).

The focus is on applying **MCS-based alignment** to detect missing species,
infer stoichiometry, and rebalance reactions in a chemically consistent manner.


The missing compound has the formula **C₈H₆O₂**, which corresponds to many possible structures, making full library enumeration infeasible.  
Instead, by extracting the **maximum common subgraph** between reactants and products, the remaining unmatched portion directly reveals the **missing fragment**.


This mapping corresponds to a `pattern_to_host` alignment, i.e. from a
smaller molecule to a larger one.  
In this example, the mapping is from the **product** (pattern) to the
**reactant** (host).

The unmatched region in the host graph therefore identifies the
**missing subgraph**, which must be imputed to complete and rebalance
the reaction.


**Stored in `synedu.Utils`** — induced-subgraph extraction is available for later rule-context construction:
```python
from synedu.Utils.graph import extract_induced_subgraph
```


**Q5 — Detect valence violations in a molecular subgraph**

Given a **molecular graph** represented as a **NetworkX graph**, detect atoms
whose **total valence violates basic chemical constraints**, using only
**local graph information** (node and edge attributes).

This exercise emphasizes a **purely graph-theoretic definition of valence**,
independent of RDKit sanitization, atom typing, or implicit chemistry rules.

---

**Step-by-step tasks**

1. **Define allowed valences**  
   Specify an extensible mapping
   $$
   \mathcal{V} : \text{Element} \;\to\; \{\text{allowed valences}\}
   $$
   for common elements (C, N, O, S, P).

2. **Compute total valence per atom**  
   For each atom $i$, compute the total valence
   $$
   v_i
   \;=\;
   h_i
   \;+\;
   \sum_{(i,j)\in E} b_{ij},
   $$
   where $h_i$ is the hydrogen count and $b_ij$ is the bond order of the
   bond between atoms $i$ and $j$.

3. **Check against valence constraints**  
   Verify whether
   $$
   v_i \in \mathcal{V}\bigl(a(i)\bigr),
   $$
   where $a(i)$ denotes the element type of atom $i$.

4. **Collect violations**  
   Record all atoms $i$ for which
   $$
   v_i \notin \mathcal{V}\bigl(a(i)\bigr).
   $$

5. **Report diagnostics**  
   For each violating atom, report:
   $$
   (\text{atom index},\; \text{element},\; v_i,\; \mathcal{V}(a(i))).
   $$

---

<details class="synedu-solution">
<summary><strong>Solution</strong></summary>

```python
_ALLOWED_VALENCE = {
    "C": {4},
    "N": {3, 4},  
    "O": {2},
    "S": {2, 4, 6},
    "P": {3, 5},
}

import networkx as nx
from typing import Dict


def find_valence_violations(G: nx.Graph) -> Dict[int, Dict]:
    violations = {}

    for node, data in G.nodes(data=True):
        element = data.get("element")
        hcount = data.get("hcount", 0)

        if element not in _ALLOWED_VALENCE:
            continue  # skip unknown elements

        bond_valence = 0.0
        for nbr in G.neighbors(node):
            bond_valence += G.edges[node, nbr].get("order", 1.0)

        total_valence = hcount + bond_valence
        allowed = _ALLOWED_VALENCE[element]

        if total_valence not in allowed:
            violations[node] = {
                "element": element,
                "hcount": hcount,
                "bond_valence": bond_valence,
                "total_valence": total_valence,
                "allowed": sorted(allowed),
            }

    return violations


violations = find_valence_violations(sub)

for idx, info in violations.items():
    print(
        f"Atom {idx} ({info['element']}): "
        f"total={info['total_valence']} "
        f"(h={info['hcount']} + bonds={info['bond_valence']}), "
        f"allowed={info['allowed']}"
    )
```

</details>



Then we need to mark this atom need to be merged with some subgraphs. We introduce wildcard here


**Stored in `synedu.Utils`** — wildcard repair and graph gluing helpers are packaged for later MCS-based completion tasks:
```python
from synedu.Utils.graph import add_wildcard, find_wildcards, merge_graphs_on_wildcard
from synedu.Utils.reaction import impute_smi_mcs
```


**Stored in synedu.Utils** — The I/O helpers used in this section (`smiles_to_graph`, `graph_to_smi`, `rsmi_to_graph`, `graph_to_rsmi`) are now packaged in `synedu.Utils`:
```python
from synedu.Utils.conversion import smiles_to_graph, graph_to_smi
from synedu.Utils.reaction  import rsmi_to_graph, graph_to_rsmi
```
Later notebooks import from these modules directly instead of relying on SynKit.


We now merge the wildcard subgraph with the hydroxyl subgraph,
using the wildcard nodes as anchors for attachment.


**Q6 — Reaction completion via MCS-based imputation**

You are given a reaction SMILES `rsmi` with a detected imbalance and a
missing species identified via MCS (`missing_smi`).

---

**Tasks**

1. **Impute missing structure**  
   Use `impute_smi_mcs` to add `missing_smi` to the appropriate side
   of the reaction.

2. **Re-check balance**  
   Recompute the imbalance using `reaction_imbalance`.

3. **Decide outcome**  
   - Balanced → accept reaction  
   - Unbalanced → apply further imputation


---

<details class="synedu-solution">
<summary><strong>Solution</strong></summary>

```python
def impute_smi_mcs(
    rsmi: str,
    status: str,
    missing_smi: str | None = None,
) -> str | None:

    # --- trivial cases ---
    if status == "balanced" or missing_smi is None:
        return rsmi

    if status == "missing_in_both_sides":
        return None

    # --- split reaction ---
    lhs, rhs = rsmi.split(">>")

    # --- impute on the correct side ---
    if status == "missing_in_products":
        rhs = f"{rhs}.{missing_smi}" if rhs else missing_smi
    elif status == "missing_in_reactants":
        lhs = f"{lhs}.{missing_smi}" if lhs else missing_smi
    else:
        return None

    return f"{lhs}>>{rhs}"



# imput missing smiles
mcs_rsmi  = impute_smi_mcs(rsmi, status, missing_smi=missing_smi)

# re-checking
delta, status = reaction_imbalance(mcs_rsmi)
print(delta, status)

# rule-based imputation
rsmi_balanced = impute_rsmi_rule(mcs_rsmi, delta, status, library=AUX_LIBRARY)
print(rsmi_balanced)
```

</details>



Here we demonstrate the same rebalancing task using **SynRBL** [\[6\]](#6.-References), a dedicated
rule-based reaction rebalancer. Starting from an imbalanced reaction SMILES,
SynRBL automatically detects the missing fragment via MCS-based alignment and
completes the reaction by imputing the appropriate substructure. The output
indicates that the reaction is successfully solved, reports the applied rules,
and assigns a high confidence score, illustrating how MCS-based reasoning can
robustly recover missing species without manual intervention.


## 4. Discussion

- **MCS** is a useful alignment primitive but is heuristic and may yield
  non-unique solutions.
- Many reactions in databases are **incomplete**, most often missing water
  or other small molecules.
- Library-based imputation is feasible for **non-carbon species** but becomes
  intractable for carbon-containing compounds.
- **MCS-based methods** help localize missing carbon fragments but still require
  chemical rules and domain knowledge for valid completion.
- **SynRBL** combines MCS alignment with rule-based reasoning to perform
  practical reaction rebalancing.


<a id="5-quiz"></a>

## 5. Quiz

1. Why is MCS useful for molecular alignment, and why can MCS be non-unique or heuristic in practical chemistry toolkits?
2. Which missing species are easiest to impute in incomplete reaction records, and what chemical information makes them easier?
3. Why does library-based imputation scale poorly for carbon-containing missing fragments?
4. What information can MCS provide when rebalancing reactions with missing carbon fragments, and where can it still fail?



## 6. References

1. RDKit documentation. https://www.rdkit.org/docs/
2. RDKit MCS (`rdFMCS`) documentation. https://www.rdkit.org/docs/source/rdkit.Chem.rdFMCS.html
3. Phan, T.-L. *et al.* Reaction rebalancing: a novel approach to curating reaction databases. *Journal of Cheminformatics* **16**, 82 (2024). https://doi.org/10.1186/s13321-024-00875-4
4. Raymond, J. W.; Willett, P. Maximum common subgraph isomorphism algorithms for the matching of chemical structures. *Journal of Computer-Aided Molecular Design* **16**, 521-533 (2002). https://doi.org/10.1023/A:1021271615909
5. NetworkX documentation. https://networkx.org/documentation/stable/
6. Bonchev, D.; Rouvray, D. H., eds. *Chemical Graph Theory: Introduction and Fundamentals*. Abacus Press (1991).
7. Diestel, R. *Graph Theory*, 5th ed. Springer (2017). https://doi.org/10.1007/978-3-662-53622-3
8. Phan, T.-L. *et al.* SynKit: A graph-based framework for rule-based reaction modeling. *Journal of Chemical Information and Modeling* (2025). https://doi.org/10.1021/acs.jcim.5c02123
