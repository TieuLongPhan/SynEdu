# S04 · From Atom-Mapped Reactions to DPO Rules: ITS, Reaction Centers, and Fast Clustering

<div class="alert alert-block alert-info">
<b>Welcome to SynEdu.</b><br>
This talktorial is part of <b>SynEdu</b>, a lightweight, research-oriented teaching series built around the
<b>Syn</b> ecosystem and <b>RDKit</b> for practical, reproducible reaction modeling.
</div>

<div class="alert alert-block alert-success">
<b>What you will gain.</b><br>
You will learn a complete, <b>mapping → ITS → reaction center → clustering → rule</b> pipeline:
starting from <b>1000 pre-mapped reactions (from a local JSON/JSON.GZ dataset)</b>, we assume atom maps are already present (generated e.g. with RXNMapper upstream), construct an
<b>ITS (Imaginary Transition State)</b> graph, extract a <b>reaction center</b>, cluster centers by
<b>typed isomorphism</b>, accelerate clustering with a <b>Weisfeiler–Lehman (WL) hash prefilter</b>,
and finally convert each center into a <b>DPO (double-pushout) graph transformation rule</b>.
</div>

<div class="alert alert-block alert-warning">
<b>Reproducibility note.</b><br>
Atom-mapping quality and rule extraction depend on tool versions, sanitization, charge handling,
and the precise definition of “reaction center” (bonds only vs bonds + atoms + radius context).
Always record package versions and settings, and persist intermediate artifacts:
mapped reaction SMILES, ITS graphs, centers, hashes, and extracted rules.
</div>

## Aim of this talktorial

Build a chemically interpretable and mathematically clean pipeline from *mapped reactions* to *graph rewriting rules*.

## Learning outcomes

After completing S04, you can:

1. Map a batch of reaction SMILES using <b>RXNMapper</b>.
2. Construct an <b>ITS graph</b> from a mapped reaction (atoms as nodes, (reactant,product) bond labels as edges).
3. Extract a <b>reaction center</b> and (optionally) a <b>radius-r context</b>.
4. Cluster reaction centers using <b>typed graph isomorphism</b>.
5. Speed up clustering via <b>WL hashing</b> as an isomorphism-invariant *prefilter*.
6. Convert a reaction center into a <b>DPO rule</b> \(L \leftarrow K \rightarrow R\) (delete/add bonds; keep context).

## Roadmap

- **0. Setup & data** (1000 mapped reactions from `./data/smart.json.gz`)
- **1. Theory** (mapped reactions, ITS, center, isomorphism, WL, DPO)
- **2. (Optional) Atom-mapping** (RXNMapper, if you start from unmapped reactions)
- **3. ITS construction**
- **4. Reaction center extraction**
- **5. Center clustering via isomorphism**
- **6. WL hash prefilter acceleration**
- **7. Convert centers to DPO rules + export**
- **8. Exercises & quiz**



## Authors and contributions

- Tieu-Long Phan (SynEdu / Syn ecosystem)
- (Add contributors here)

> If you reuse code snippets or figures from this notebook, please cite SynEdu appropriately and pin versions in your supplement.



## 0. Setup & data



### 0.0 Note on atom mapping (RXNMapper)

In **S01/S04** we often start from *unmapped* reaction SMILES and run **RXNMapper** to obtain atom-map IDs.
In this S04 run, we **skip mapping** because our dataset file `./data/smart.json.gz` already contains **atom-mapped** reactions.

> If you ever start from unmapped reactions, you can still plug RXNMapper into the same pipeline:
> produce a column `mapped_rxn` and continue with ITS construction below.



### 0.1 Reaction data (1000 mapped reactions)

We load a local dataset from:

- `./data/smart.json.gz`

The file may be either:

- **plain JSON** (often misnamed as `.gz`), or
- **gzip-compressed JSON**.

Each entry is expected to contain a **mapped reaction string** (reaction SMILES/SMARTS) with atom-map IDs like `:12`.
In our dataset, the mapped reaction is stored under the key **`"smart"`** (and optionally an identifier like `"R-id"`).

We will load **1000 mapped reactions** and store them in a DataFrame `df_map` with at least:

- `mapped_rxn` — mapped reaction string (`reactants >> products`)
- `rxn_id` — optional id (if available)



## 1. Theory (formal, but chemist-friendly)

### 1.1 Atom-mapped reactions as typed structures

An unmapped reaction SMILES is a string

$$
\texttt{rxn} = \texttt{R} \;\gg\; \texttt{P},
$$

where $\texttt{R}$ (reactants and reagents) and $\texttt{P}$ (products) are dot-separated
lists of molecules.

An **atom-mapped reaction** augments each atom with an integer **map identifier**.
Conceptually, atom mapping defines a *partial correspondence* between atoms on the
reactant and product sides.

Let

$$
A_R \quad\text{and}\quad A_P
$$

be the sets of atoms on the reactant and product sides, respectively, and let

$$
M \subset \mathbb{N}
$$

be the set of map identifiers appearing in the reaction.
The atom-mapped reaction induces two functions

$$
\mu_R : A_R \to M,
\qquad
\mu_P : A_P \to M,
$$

where atoms sharing the same map identifier are interpreted as representing
“the same atom through the transformation”.

In practice, atom mapping may be imperfect:
missing map identifiers, duplicated identifiers, or unmapped leaving groups may occur.
All downstream constructions are therefore designed to be **robust**:
invalid or inconsistent cases are detected, warned about, and skipped.

---

### 1.2 ITS (Imaginary Transition State) graph

The **Imaginary Transition State (ITS)** represents the reaction as a single graph
by explicitly encoding bond changes between reactants and products.

Let $B$ be the set of allowed bond labels
(e.g. $\{1,2,3,1.5\}$ for single, double, triple, and aromatic bonds),
and introduce a special symbol $\bot$ denoting **bond absence**.

The ITS of an atom-mapped reaction is defined as a **typed graph**

$$
G_{\mathrm{ITS}} = (V, E, a, b_R, b_P),
$$

where:

- $V \subseteq M$ is the set of atom map identifiers used as node identifiers;
- $a : V \to \Sigma$ assigns atom types (e.g. element symbol, optionally charge);
- $E$ is the set of unordered edges between atom identifiers;
- $b_R : E \to B \cup \{\bot\}$ assigns the **reactant-side bond label**;
- $b_P : E \to B \cup \{\bot\}$ assigns the **product-side bond label**.

Chemical interpretation:

- $b_R(e) = \bot$ and $b_P(e) \neq \bot$ corresponds to **bond formation**;
- $b_R(e) \neq \bot$ and $b_P(e) = \bot$ corresponds to **bond breaking**;
- $b_R(e) \neq b_P(e)$ with both present corresponds to a **bond order or aromaticity change**.

Thus, the ITS compactly encodes the full local transformation in a single object.

---

### 1.3 Reaction center and radius-$r$ context

The **reaction center** is defined as the part of the ITS where changes occur.

First, define the set of changed edges

$$
\Delta E = \{\, e \in E \mid b_R(e) \neq b_P(e) \,\}.
$$

The **core atom set** is then

$$
C_0
=
\{\, v \in V \mid \exists e \in \Delta E \text{ incident to } v \,\},
$$

i.e. the atoms directly involved in bond changes.
Optionally, atoms with changed intrinsic states
(e.g. formal charge changes) may also be included in $C_0$.

To incorporate local chemical context, we define the **radius-$r$ expansion**

$$
C_r
=
\{\, u \in V \mid \mathrm{dist}_{G_{\mathrm{ITS}}}(u, C_0) \le r \,\},
$$

where $\mathrm{dist}_{G_{\mathrm{ITS}}}$ denotes the graph distance in the ITS.
The induced subgraph

$$
G_{\mathrm{ITS}}[C_r]
$$

is referred to as the **reaction-center graph with radius $r$**.

---

### 1.4 Clustering reaction centers by typed isomorphism

Two reaction-center graphs $G$ and $H$ are **typed isomorphic** if there exists a bijection

$$
\varphi : V(G) \to V(H)
$$

such that:

$$
a_G(v) = a_H(\varphi(v))
\quad\text{for all } v \in V(G),
$$

and

$$
\bigl(b_R, b_P\bigr)_G(\{u,v\})
=
\bigl(b_R, b_P\bigr)_H(\{\varphi(u), \varphi(v)\})
\quad\text{for all } \{u,v\} \in E(G).
$$

Chemically, this means that the two reactions share the **same local transformation
pattern**, independent of atom numbering or molecular size.

Reaction-center clustering groups reactions by this equivalence relation.

---

### 1.5 Weisfeiler--Lehman hash as a fast prefilter

Exact graph isomorphism testing becomes expensive for large collections of
reaction-center graphs.
To accelerate clustering, we use the **Weisfeiler--Lehman (WL) graph hash**.

The WL procedure iteratively refines node labels by hashing:

$$
\text{current node label}
\;+\;
\text{multiset of neighboring node and edge labels}.
$$

The final WL hash is an **isomorphism invariant**:
isomorphic graphs always share the same hash.
However, the invariant is **not complete**:
non-isomorphic graphs may occasionally collide.

Therefore, WL hashing is used only as a **prefilter**:

1. group reaction-center graphs by identical WL hash;
2. perform exact typed isomorphism checks only within each hash bucket.

This strategy preserves correctness while dramatically reducing runtime.

---

### 1.6 From reaction centers to DPO rules

A **double-pushout (DPO) rule** in graph rewriting is defined as a span

$$
p
=
\left(
L \xleftarrow{\;\ell\;} K \xrightarrow{\;r\;} R
\right),
$$

where:

- $L$ is the **left-hand side** pattern to be matched;
- $R$ is the **right-hand side** pattern produced by the rewrite;
- $K$ is the **context graph** that is preserved.

In the setting of atom-mapped chemical reactions:

- nodes of $K$ correspond to atoms preserved through the reaction;
- edges of $K$ correspond to bonds whose type does not change;
- edges present in $L$ but not in $R$ are **deleted bonds**;
- edges present in $R$ but not in $L$ are **created bonds**.

Atom-state changes (e.g. formal charge) are naturally represented as
node attributes that may differ between $L$ and $R$,
while $K$ retains only invariant identity information
(e.g. element symbol).

In practice, we extract rules as follows:

- $L$: reactant-side bonds induced on $C_r$;
- $R$: product-side bonds induced on $C_r$;
- $K$: intersection of bonds on $C_r$ whose labels are unchanged.

This construction makes the common intuition of
“reaction centers as rewriting rules”
explicit in precise graph-theoretic terms.



## 2. Atom mapping (already provided)

The downstream ITS / reaction-center pipeline assumes each atom has a **map id** that is consistent
between reactant and product sides.

In this notebook we **do not run RXNMapper**.
Instead, we **validate** that the loaded reactions look mapped:
- the reaction string contains `:n` map annotations,
- RDKit parses the molecules,
- at least one mapped atom appears on each side.



### 2.1 Inspect a mapped reaction

Our dataset stores a reaction SMILES where atoms carry map ids like `:12`.  
We will parse these map ids and build our ITS.



## 3. Construct the ITS graph



### 3.1 Implementation idea (from mapped SMILES → ITS)

We parse the mapped reaction into RDKit molecules on each side.

For each atom (with map id \(m\)) we record (at least):

- element symbol
- formal charge on reactant side and product side

For each bond between map ids \((u,v)\) we record:

- \(b_R(u,v)\): reactant bond order (or \(\bot\))
- \(b_P(u,v)\): product bond order (or \(\bot\))

We then take the **union** of reactant and product bonds as edges of the ITS.

> Chemistry convention used here:
> aromatic bonds are treated as order \(1.5\).



### 3.2 Visualize an ITS as a typed graph

For readability, we draw each edge label as `(r_order → p_order)` where `None` means \(\bot\) (absent bond).
Changed bonds are the mechanistic “action”.



## 4. Extract the reaction center



### 4.1 Center definition used in this notebook

We define the *core* center atoms \(C_0\) as atoms incident to **changed bonds**:

\[
C_0 = \{v \mid \exists e \text{ incident to } v \text{ with } b_R(e)\neq b_P(e)\}.
\]

Optionally, we also include atoms with **formal charge changes** (reactant vs product),
because for some chemistries the bond pattern alone does not capture the transformation.

Finally, we optionally expand to radius \(r\) context and return the induced subgraph.



### 4.2 Visualize the reaction center on the molecules

We highlight center atoms (by map id) on both sides of the mapped reaction.



## 5. Reaction center clustering via typed isomorphism



### 5.1 Typed matching predicates

We treat each center as a typed graph with:

- node labels: \((\text{symbol}, r\_charge, p\_charge)\)
- edge labels: \((r\_order, p\_order)\)

Two centers are equivalent if there exists an isomorphism preserving these labels.

> You can relax or strengthen these predicates depending on your use-case:
> - ignore charge,
> - treat aromatic \(1.5\) specially,
> - include radius-1 context, etc.



### 5.2 Inspect representative centers

A cluster representative is the first center graph assigned to that cluster.



## 6. Accelerate clustering with WL graph hashing (prefilter)



### 6.1 WL hash prefilter

We compute a WL hash for each center graph using label strings:

- node label: `symbol | r_charge → p_charge`
- edge label: `r_order → p_order` (with `⊥` for absent)

Then we group by hash and run expensive isomorphism only within each group.

> **Reminder:** WL hashing is an isomorphism-invariant heuristic, not a proof.
> We still run exact isomorphism within each hash bucket to obtain correct clusters.



### 6.2 Sanity check: WL-prefilter + exact isomorphism equals exact isomorphism alone

Because WL is only a prefilter, the final clustering should be identical to direct isomorphism clustering
(up to cluster id renaming).



## 7. Convert reaction centers into DPO rules



### 7.1 Practical rule extraction used here

Given an ITS graph \(G_{\mathrm{ITS}}\) and a chosen context node set \(C_r\), we build:

- \(L\): **reactant pattern** on \(C_r\) (keep only bonds with \(b_R\neq\bot\))
- \(R\): **product pattern** on \(C_r\) (keep only bonds with \(b_P\neq\bot\))
- \(K\): **preserved context** on \(C_r\) (keep only bonds with \(b_R=b_P\neq\bot\))

Nodes:
- \(L\) nodes carry `symbol` and `charge=r_charge`
- \(R\) nodes carry `symbol` and `charge=p_charge`
- \(K\) nodes carry `symbol` only (identity type)

This yields a DPO-style span \(L \leftarrow K \rightarrow R\) that is directly interpretable as:
“delete bonds in \(L\setminus K\), keep bonds in \(K\), create bonds in \(R\setminus K\)”.

> For production systems, you may also:
> - include stereochemistry / aromatic flags explicitly,
> - attach atom environments,
> - require connected centers,
> - handle reagents vs true reactants,
> - normalize charges, etc.



### 7.2 Extract rules for all reactions and export

We extract one DPO rule per reaction (for a chosen center radius), then optionally export:

- all rules
- cluster representatives (one per isomorphism class of centers)

The export format is a simple JSON with nodes/edges for L,K,R.



## 8. Exercises & quiz

<div class="alert alert-block alert-info">
<b>Cross-referencing.</b> Concepts used here build on:
- <b>S01</b>: typed graphs, isomorphism, automorphisms
- <b>S03</b>: (if you have it) reaction-side graph construction and label design choices
</div>

### Exercise 1 — Center radius as a hyperparameter

Implement `radius=2` and compare:
- cluster counts,
- representative centers,
- rule sizes (nodes/edges).

**Chemistry intuition:** larger radii incorporate more “environment”, making rules more specific.

### Exercise 2 — Relax the typing

Modify the matcher to ignore charge:

- node label uses only `symbol`
- edge label still uses `(r_order, p_order)`

Does the number of clusters decrease? Which chemistries merge?

### Exercise 3 — WL iterations vs collisions

Change `wl_iters` in `cluster_centers_wl_prefilter` and measure:
- runtime,
- number of WL buckets,
- collision rate (bucket size distribution).

### Quiz (conceptual)

1. What would go wrong if we define \(K\) as empty for every rule?
2. How does a bond order change (single→double) appear in the \(L,K,R\) decomposition?
3. In practice, why is WL hashing only a prefilter and not a proof of isomorphism?
4. If RXNMapper assigns a wrong map id, which step in the pipeline is most sensitive, and why?



# References and further reading

- RXNMapper: Schwaller et al., transformer-based atom mapping (see the RXNMapper repository / publication).
- ITS graphs in reaction modeling: classic “imaginary transition state” representations in reaction informatics.
- Double pushout (DPO) graph rewriting: Ehrig et al., algebraic graph transformation (intro texts); chemical applications in rule-based reaction systems.
- Weisfeiler–Lehman (WL) refinement and hashing: WL test / color refinement; widely used in cheminformatics fingerprints and GNN theory.
- RDKit documentation and book: https://www.rdkit.org/docs/
- NetworkX documentation: https://networkx.org/documentation/stable/

> For a publication-quality method section, replace these bullets with full BibTeX entries and pin exact tool versions.
