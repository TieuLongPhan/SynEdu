# S06 : Canonicalizing Atom-Mapped Reactions and Rules.


## Aim of this talktorial

This talktorial (**S06**) explains why atom-mapped reactions and extracted rules need **canonicalization** before they can be compared, counted, or stored in a reproducible library.

Atom-map numbers are identifiers, not chemistry. The same reaction can be mapped with different but equivalent numbers because of atom ordering, molecular symmetry, or mapper-specific conventions. Without canonicalization, duplicate reactions and duplicate rules can look artificially different.

Concretely, we focus on:

1. **Partition and refinement**  
   Using graph partitions and Weisfeiler-Lehman-style refinement to distinguish atoms by local structural context.

2. **Individualization and exact methods**  
   Resolving ambiguous symmetric atoms when refinement alone is insufficient.

3. **Atom-map canonicalization**  
   Reindexing mapped reactions deterministically while preserving the underlying chemistry.

4. **Rule-level deduplication**  
   Measuring how canonicalization changes the number of unique reaction centers and rules.

---

## Learning outcomes

After completing this talktorial, you will be able to:

- explain why atom maps are non-unique and why map numbers should be canonicalized,
- define graph partitions, equitable refinement, and WL-style color refinement,
- implement a deterministic atom-map reindexing $\pi: \mathbb{N} \to \mathbb{N}$ based on map-invariant structural ranks,
- canonicalize mapped reaction SMILES by molecule order and map IDs without changing chemistry,
- identify cases where refinement leaves unresolved symmetry, and
- quantify how canonicalization affects the number of unique centers or rules using hashes and isomorphism checks.

---

## Outline

<ul class="synedu-outline">
  <li><a href="#0-setup--data">0. Setup &amp; Data</a></li>
  <li><a href="#1-partition-refinement-and-approximation">1. Partition, Refinement, and Approximation</a></li>
  <li><a href="#2-individualization-refinement-and-exact-methods">2. Individualization, Refinement, and Exact Methods</a></li>
  <li><a href="#3-atom-mapped-canonicalization">3. Atom-Mapped Canonicalization</a></li>
  <li><a href="#4-quiz">4. Quiz</a></li>
  <li><a href="#5-discussion">5. Discussion</a></li>
  <li><a href="#6-references">6. References</a></li>
</ul>


## 0. Setup & data



## 1. Partition, Refinement and Approximation

### 1.1. Partition

- Let $G=(V,E)$ be a (vertex-)colored graph. A **partition** is $\Pi=\{C_1,\dots,C_k\}$ with $\bigcup_i C_i=V$ and $C_i\cap C_j=\varnothing$ for $i\neq j$.
- $\Pi$ is **equitable** (WL-stable) iff for all cells $C_i,C_j$ there exists a constant $c_{ij}$ such that
  $$
  \forall v\in C_i:\quad |N(v)\cap C_j| = c_{ij}.
  $$


**Definition (Vertex Partition).**  
A *partition* of $V$ is a collection $\mathcal{P} = \{C_0, C_1, \ldots, C_k\}$ of non-empty disjoint sets (*cells*) with $\bigcup_i C_i = V$. A partition is *discrete* if each cell has exactly one element (i.e., $k = |V| - 1$).

**Definition (Equitable Partition).**  
A partition $\mathcal{P}$ of $V(G)$ is *equitable* (or *stable*) if for every pair of cells $C_i, C_j$ and every vertex $v \in C_i$, the number of neighbors of $v$ in $C_j$ depends only on the cell $C_i$ — not on the specific vertex $v$:

$$
|N(v) \cap C_j| = d_{ij} \quad \text{for all } v \in C_i
$$

**Definition (WL Color Refinement).**  
Starting from the initial partition $\mathcal{P}^{(0)}$ induced by node labels $\mathbf{a}$, *Weisfeiler–Lehman (WL) refinement* iterates:

$$
\mathcal{P}^{(t+1)} \;=\; \text{refine}\!\left(\mathcal{P}^{(t)},\, G\right)
$$

by splitting each cell $C \in \mathcal{P}^{(t)}$ whenever two vertices in $C$ have different *neighborhood signatures* — i.e., different multisets of $(\text{cell label, edge label})$ pairs among their neighbors. Refinement terminates at the *stable partition* $\mathcal{P}^* = \mathcal{P}^{(t^*)}$ where no further splitting occurs.

**Theorem (WL Soundness).** If $G_1 \cong G_2$ (as labeled graphs), then WL refinement produces the same stable partition histogram. The converse fails in general: two non-isomorphic graphs can yield identical WL histograms.



- **`C0: [1]` — methyl carbon (unique)**  
  Node `1` is a *non-aromatic carbon* with degree `1` (the methyl substituent). Its tuple `(element='C', aromatic=False, degree=1)` is unique in the graph, so it forms its own color class. 

- **`C2: [2]` — ipso aromatic carbon (unique)**  
  Node `2` is the *ipso* aromatic carbon (the ring carbon bonded to the methyl group) with degree `3`. The triple `(C, aromatic=True, degree=3)` differs from the ring CH positions, so it also forms its own class.

- **`C1: [3,4,5,6,7]` — remaining ring carbons (ambiguous)**  
  Nodes `3–7` are aromatic carbons with degree `2`, so the chosen attributes make them locally indistinguishable. WL groups them into one big class. These are the ring CH positions that, from the local viewpoint (element, aromatic flag, neighbor count), look the same while the global symmetry (reflections/rotations) determines which are actually interchangeable.



### 1.2 Refinement

The refinement operator $R(\Pi)$ returns the *coarsest equitable partition* that refines $\Pi$. A partition $\Pi=\{C_1,\dots,C_k\}$ is equitable iff for every pair of cells $C_i,C_j$ and every $u,v\in C_i$,
$$
|N(u)\cap C_j| \;=\; |N(v)\cap C_j|.
$$
Operationally we compute, for each vertex $v$, the signature
$$
\sigma(v) \;=\; \big(|N(v)\cap C_1|,\;|N(v)\cap C_2|,\;\dots,\;|N(v)\cap C_k|\big),
$$
and split cells by identical signatures until nothing splits.


**Chemical intuition**

- **Signature = local chemical context.**  
  The signature $\sigma(v)$ is a compact count-vector that records *how many* neighbors $v$ has in each *color class* (each class corresponds to a type of chemical environment defined so far). Chemically this is like saying: *“how many of my neighbors are methyl-like, ipso-like, aromatic-CH-like, etc.”* — not which exact atoms they are.

- **Splitting = discovering positional roles.**  
  When a cell $C$ splits into subcells, WL has discovered that atoms that looked the same under the current descriptors actually occupy *different positions* in the molecular graph (e.g., ortho vs. meta vs. para). Splitting is not an abstract algebraic trick — it corresponds to chemically meaningful distinctions (bonding pattern relative to substituents).

- **Equitable partition = indistinguishability under chosen observables.**  
  If a partition is equitable, then *within each cell* all atoms have identical counts of neighbors in every environment-class. That means *no additional local topological information* (based on current coloring) can tell those atoms apart; they are locally indistinguishable under the chosen invariants.


**Toluene example**

Start with an initial coloring by $(\text{element},\ \text{aromatic},\ \deg)$:

```text
C0 = [1] # methyl carbon (non-aromatic, deg=1)
C1 = [3,4,5,6,7] # ring CH positions (aromatic, deg=2)
C2 = [2] # ipso carbon (aromatic, deg=3)
```


Order cells as $(C_0,C_1,C_2)$. Compute signatures $\sigma(v)=(|N(v)\cap C_0|,|N(v)\cap C_1|,|N(v)\cap C_2|)$:

- $\sigma(3)=(0,1,1)$ because node 3 neighbors are $\{2,4\}$: one in $C_1$ (node 4) and one in $C_2$ (node 2).  
- $\sigma(4)=(0,2,0)$ because node 4 neighbors $\{3,5\}$ both lie in $C_1$.  
- $\sigma(5)=(0,2,0)$ (neighbors $\{4,6\}$).  
- $\sigma(6)=(0,2,0)$ (neighbors $\{5,7\}$).  
- $\sigma(7)=(0,1,1)$ (neighbors $\{6,2\}$).

So $C_1$ splits into two signature groups:
- $C_{1a}=[3,7]$ with $\sigma=(0,1,1)$ (these are *ortho* relative to the substituent), and  
- $C_{1b}=[4,5,6]$ with $\sigma=(0,2,0)$ (positions not directly adjacent to ipso).

Refine again using $[C_0,C_{1a},C_{1b},C_2]$:
- Within $C_{1b}$ the new signatures differentiate the para position (node 5) from its neighbors: $C_{1b}\to[4,6]$ and $[5]$.  

Final $R(\Pi)$ for this example:
```text
[[1], [3,7], [4,6], [5], [2]]
```



### 1.3 WL refinement for chemical graphs

At each iteration replace each node’s color by  
`(old_color, sorted_multiset((neighbor_color, edge_label)))`, re-rank tuples lexicographically, stop when stable.

---

**Formal algorithm**

```text
Inputs:
  G = (V,E) with node attrs a_V(v) and edge attrs a_E(u,v)
  max_iter = 50

Initialize:
  color0(v) := tuple(a_V(v))                # base node signature
  map distinct color0 → integer ids (lexicographic)

For t = 0 .. max_iter-1:
  for each v in V:
    neigh_sig := sorted([ ( color_t(u), a_E(v,u) ) for u in N(v) ])
    sig_v := ( color_t(v), tuple(neigh_sig) )
  rank distinct sig_v lexicographically → color_{t+1}(v)
  if color_{t+1} == color_t for all v: break

Output:
  partition := classes of nodes with equal final color
  colors := mapping node → final color id
```

**Chemical intuition**

- Base color = atom identity (element, aromatic flag, H-count, charge, …).

- Refinement = who you are connected to and how (neighbor color + bond label).

- Result = an *equitable* partition: atoms in same class are locally indistinguishable under chosen features.

- Use as: orbit approximation, node features for ML, seed for exact automorphism.


**Toluene**

```text
              3
            /   \
          4       2 --1
          |       |
          5       7
            \   /
              6
```

**t = 0 — initial signatures (by node attrs)**

- `s0(1) = (C, False, h=3)`  → **color A**  — methyl (unique)  
- `s0(2) = (C, True,  h=0)`  → **color B**  — ipso (unique)  
- `s0(3,4,5,6,7) = (C, True, h=1)` → **color C**  — ring carbons (initially same)

---

**t = 1 — incorporate neighbor signatures**

- `sig(1) = (A, [ (B,1.0) ])`                      → `{1}`  
- `sig(2) = (B, [ (A,1.0),(C,1.5),(C,1.5) ])`      → `{2}`  
- `sig(3) = (C, [ (B,1.5),(C,1.5) ])`              → `{3,7}`  
- `sig(4) = (C, [ (C,1.5),(C,1.5) ])`              → `{4,5,6}`  

**Partition after t = 1:** `{1}, {2}, {3,7}, {4,5,6}`

---

**t = 2 — refine within ring classes**

Using t=1 colors, neighbors of `{4,5,6}` differ → split:  
**New partition:** `{1}, {2}, {3,7}, {4,6}, {5}`

---

**t = 3 — stable (no further splits)**

**Final partition:**

- `{1}` — methyl (CH₃)  
- `{2}` — ipso (attached to methyl)  
- `{3,7}` — ortho carbons (equivalent)  
- `{4,6}` — meta carbons (equivalent)  
- `{5}` — para carbon (unique)



#### WL refinement iteration strip

Each column shows the partition class (colour) assigned to each node at one WL iteration. Nodes with the same colour are in the same equivalence class at that step. As refinement proceeds, symmetric atoms split into finer classes until the partition is **stable** (no cell changes between two consecutive columns). Atoms that remain in the same class throughout are **WL-indistinguishable**.


**Figure — WL coloring progression on toluene.**  
Each panel shows the graph after one refinement step. Nodes with the same colour share the same partition class (WL-indistinguishable at that iteration). Iteration 0 is the initial partition by element type; the partition stabilises once no cell splits further.


#### WL cell count convergence across molecules

More symmetric molecules need fewer iterations to reach a stable partition
(benzene stabilizes early because all carbons are equivalent).
Less symmetric molecules keep splitting cells for more iterations.




### 1.4. Approximation

Now we have the **stable WL partition** which is obtained after Weisfeiler–Lehman (WL) color refinement converges.
Each inner list represents a color class, i.e. a set of nodes that WL cannot
distinguish further based on their attributes and local neighborhoods.

In the toluene example, this partition already has a clear chemical meaning:
the methyl carbon, the ipso carbon, the para position, and two symmetric pairs
(ortho and meta positions).

From the partition above, a WL-consistent block ordering may be written as

```text
[1, (3,7), 5, (4,6), 2]
```

Only the two parenthesized blocks are ambiguous. Enumerating their internal
permutations yields \(2! \times 2! = 4\) candidate linearizations. Evaluating
all of them and selecting the lexicographically smallest graph signature
produces a canonical index that is strictly stronger than plain WL
canonicalization, while remaining inexpensive.

This approach is **not a full automorphism solver**, but a WL-seeded,
locally exact approximation that is well sui


## 2. Individualization, Refinement, and Exact Methods

- **Individualization**: for $v\in C$ write $I_v(\Pi)$ for the partition obtained by replacing $C$ with $\{v\}$ and $C\setminus\{v\}$.
- The **IR search tree** has nodes (partitions) and children of a node $\Pi$ given by $R(I_v(\Pi))$ for chosen $v$ in a non-singleton cell.
- During IR Nauty discovers automorphisms $\sigma\in\operatorname{Aut}(G)$ and stores a generating set $\mathcal{G}$. If a discovered $\sigma$ maps a search node to an already visited node, Nauty **prunes** that branch.


**Individualization-Refinement (IR) search tree**

When the stable WL partition $\mathcal{P}^*$ is not fully discrete, we use
*individualization*: pick one ambiguous vertex, fix its color, then refine again.
This creates a binary tree of partition states; the canonical ordering
is the leaf with the lexicographically smallest certificate.


**Individualization depth across molecules**

After WL stabilises, the number of **non-singleton cells** tells us
how many atoms remain ambiguous — and therefore how many individualisation
steps the IR tree needs.  
Highly symmetric molecules (benzene) need the deepest tree;
fully asymmetric ones (alanine) need none at all.



## 3. Atom-Mapped Canonicalization

Now we combine all of them


### Rule deduplication before vs after canonicalization

Two reactions that differ only in atom-map numbering are the **same rule**.
Canonical atom-map numbering collapses this redundancy and reduces the
apparent rule vocabulary.


**Figure — two differently-mapped versions of the same reaction, before and after canonicalization.**  
The raw SMILES differ in atom-map numbers; after `canon_aam` both collapse to the same string, enabling exact deduplication.


### Canonicalization comparison

Different atom-mapping methods (RXNMapper, Graphormer, Local Mapper) may assign different atom-map numbers to the same reaction. After canonicalization, equivalent maps collapse to the **same canonical SMILES**. The table below checks whether all three methods produce identical canonical forms for one example reaction. Green = agreement; red = mismatch.


## 4. Quiz

1. **WL sensitivity to labels**  
   In `wl1_partition_nx`, change `node_attrs` / `edge_attrs` and observe how the stable partition changes.
   - What happens if you drop `edge_attrs=("order",)`?
   - What happens if you add `"degree"` into `node_attrs`?

2. **A maximally symmetric molecule**  
   Build a mapped benzene graph and run labeled WL:
   - Why does WL stop before a discrete partition?
   - Which step in IR is needed to fully canonicalize?

3. **Reaction-level canonicalization**  
   Use `canon_aam` on:
   - the toy example (already provided), and
   - a real mapped reaction from `data` (if available).  
   Confirm that **permuting atom-map ids** yields the **same canonical reaction**.

4. **Challenge**  
   Extend `canon_aam` to canonicalize/sort **agents** too (if they are mapped).



## 5. Discussion


### What you should take away

- **Atom maps are not unique**: even for the *same* molecule/reaction, symmetric atoms can be renumbered without changing chemistry.
- **WL refinement is fast but approximate**: it produces an equitable partition; if it is not fully discrete, symmetry remains.
- **IR (individualization + refinement) gives a true canonical labeling** (at higher computational cost): it resolves remaining symmetries by branching and uses refinement to keep the search manageable.
- For **reactions**, canonicalization must be **reaction-level**: map ids are *global* across molecules, so you need a global object (e.g., an ITS-like graph) to define a stable renumbering.

### Practical notes

- For most chemically labeled graphs, **labeled WL** almost-discretizes the graph, so the IR search is shallow.
- Truly symmetric cases (e.g., benzene, cubane, or highly symmetric ions) are where IR matters most.
- In production pipelines, you typically combine:
  - strong initial labels (chemistry-aware),
  - WL refinement,
  - IR with pruning and orbit-based symmetry breaking (Nauty/Bliss-style ideas).



## 6. References

1. <a id="ref-1"></a> Weisfeiler, B.; Lehman, A.  
   *A reduction of a graph to a canonical form and an algebra arising during this reduction.*  
   (1968).  *(The classic WL refinement report.)*

2. <a id="ref-2"></a> McKay, B. D.; Piperno, A.  
   *Practical Graph Isomorphism, II.*  
   *Journal of Symbolic Computation* **60**, 94–112 (2014).  *(Nauty/Traces family; IR + pruning.)*

3. <a id="ref-3"></a> Morgan, H. L.  
   *The generation of a unique machine description for chemical structures—A technique developed at Chemical Abstracts Service.*  
   *Journal of Chemical Documentation* **5**(2), 107–113 (1965).  *(Early canonicalization ideas in chem-informatics.)*

4. <a id="ref-4"></a> Phan, T.-L.; González Laffitte, M. E.; Weinbauer, K.; Merkle, D.; Andersen, J. L.; Fagerberg, R.; Gatter, T.; Stadler, P. F.  
   *SynKit: A Graph-Based Python Framework for Rule-Based Reaction Modeling and Analysis.*  
   *Journal of Chemical Information and Modeling* (2025).
