# S03 · Reaction Rules as Graph Rewriting (DPO) — Chemistry-first, ITS-native

<div class="alert alert-block alert-info">
<b>SynEdu (S03).</b><br>
This notebook continues S01 (typed molecular graphs & matching) and upgrades from “static graphs” to
<b>transformations</b>: reaction rules, graph rewriting, and chemistry-flavoured examples (Diels–Alder).
</div>

<div class="alert alert-block alert-success">
<b>Aim.</b><br>
Represent a reaction as a <b>rule</b> and apply it to molecules in a way that is:
<ul>
  <li><b>mathematically precise</b> (typed morphisms, DPO rewriting)</li>
  <li><b>chemistry interpretable</b> (bond breaking/forming, conserved atoms)</li>
  <li><b>implementation-ready</b> (RDKit ⇄ NetworkX, atom-mapped reaction SMARTS)</li>
</ul>
</div>

<div class="alert alert-block alert-warning">
<b>Key design decision in S03.</b><br>
For export and comparison we construct an <b>ITS-like change graph</b> <i>directly</i>, then decompose it into
reactant/product graphs. This makes atom correspondence <b>identity-by-node-id</b> (no extra mapping step).
</div>



## Learning outcomes

By the end of this talktorial you can:

1. Write a reaction rule as a **span** \(L \leftarrow K \rightarrow R\) and explain each part chemically.
2. Formally define a **match** \(m: L \hookrightarrow G\) and enumerate all applications in a host molecule.
3. Construct an **ITS graph** \(T\) capturing bond changes (pre/post labels) *without first materializing the product*.
4. Decompose \(T\) back into reactant/product graphs (projection), and export **atom-mapped reaction SMARTS**.
5. Implement a worked Diels–Alder example and inspect **formed/broken/changed** bonds.

## Roadmap

- **0. Setup** (imports, versions, helpers)
- **1. Formalism** (typed graphs, morphisms, DPO rewriting, correctness conditions)
- **2. RDKit ⇄ NetworkX** (typed molecular graphs)
- **3. DPO kernel** (matching + rewrite) and why chemistry needs *wildcards*
- **4. ITS-first rewriting** (construct ITS directly; decompose into reactant/product)
- **5. Export** (atom-mapped reaction SMARTS from ITS)
- **6. Example**: Diels–Alder (butadiene + ethene → cyclohexene)
- **7. Exercises**



## Example molecules used in this notebook

We will use a classic *pericyclic* transformation:

- **Diene**: butadiene (SMILES: `C=CC=C`)
- **Dienophile**: ethene (SMILES: `C=C`)
- Combined as a disconnected reactant mixture: `C=CC=C.C=C`

> Chemistry caveat: this is a toy model (no stereochemistry, no substituents, no endo/exo control).
The point is to demonstrate **bond reorganisation** under a rule-based graph rewriting view.



## 1. Formalism: chemistry as typed graph rewriting

### 1.1 Typed molecular graphs (recap from S01)

A (heavy-atom) molecule is a **typed graph**

\[
G = (V_G, E_G, a_G, b_G),
\]

where:

- \(V_G\) = atoms (vertices),
- \(E_G \subseteq \{\{u,v\}\mid u\neq v\}\) = bonds (undirected edges),
- \(a_G: V_G \to \Sigma_V\) = **atom typing** (element, charge, aromaticity, chirality, …),
- \(b_G: E_G \to \Sigma_E\) = **bond typing** (order, aromaticity, ring flag, …).

Chemistry interpretation:
- vertex labels enforce that carbon maps to carbon, etc.
- edge labels enforce that a double bond is not confused with a single bond unless you *allow* it.

### 1.2 Label-preserving morphisms

A **typed graph morphism** \(\varphi: G \to H\) is a function on vertices
\(\varphi_V: V_G \to V_H\) such that:

1. **Atom label preservation**  
\[
a_H(\varphi_V(v)) = a_G(v)\quad \forall v\in V_G.
\]

2. **Adjacency and bond label preservation**  
For every \(\{u,v\}\in E_G\),
\[
\{\varphi_V(u),\varphi_V(v)\}\in E_H
\quad\text{and}\quad
b_H(\{\varphi_V(u),\varphi_V(v)\}) = b_G(\{u,v\}).
\]

In practice we often restrict to **injective** morphisms (embeddings), written
\[
m: G \hookrightarrow H,
\]
which is what NetworkX returns in subgraph matching.

### 1.3 Reaction rules as spans \(L \leftarrow K \rightarrow R\)

A reaction rule is encoded by three typed graphs:

- \(L\) (**left**) — the pattern to be found in a reactant
- \(R\) (**right**) — what replaces it in the product
- \(K\) (**context/interface**) — what is preserved (the “atom identity backbone”)

Formally a DPO rule is a span of morphisms

\[
L \xleftarrow{\ell} K \xrightarrow{r} R,
\]

where \(\ell\) and \(r\) embed \(K\) into \(L\) and \(R\).

Chemistry interpretation:
- \(V_K\) are the atoms that keep their identities through the reaction (the “mapped atoms”).
- \(E_L \setminus E_K\): bonds that must be **broken** (or bond types that must disappear).
- \(E_R \setminus E_K\): bonds that must be **formed**.

A bond **order change** is naturally represented by
- old bond in \(L\) but not in \(K\),
- new bond in \(R\) but not in \(K\),
while the endpoints remain in \(K\).

### 1.4 Where can a rule apply? matches \(m: L \hookrightarrow G\)

Given a reactant molecule graph \(G\), a match is an injective morphism

\[
m: L \hookrightarrow G.
\]

Chemically: “find a substructure in the reactant that looks like the reacting pattern”.

NetworkX’s `GraphMatcher(G, L).subgraph_isomorphisms_iter()` enumerates all such matches.

### 1.5 Correctness conditions you should care about (chemistry view)

DPO rewriting is “correct” when two constraints hold (informally):

1. **Dangling condition** (no half-bonds):  
if you delete an atom (node), you must also delete all incident bonds.  
Chemistry: you cannot leave a bond to a missing atom.

2. **Identification condition** (don’t accidentally merge distinct atoms):  
the match must not identify two different pattern atoms onto the same host atom.  
Chemistry: two atoms in the reacting pattern cannot map to one atom in the molecule.

In this notebook we focus on the common chemistry case: **atom-conserving rules**
(\(V_L = V_K = V_R\)) and only **bond changes**.



### 1.6 The categorical DPO diagram (what “double pushout” means)

In category language, rewriting is defined by two pushouts:

\[
\require{AMScd}
\begin{CD}
K @>{r}>> R \\
@V{\ell}VV @VVV \\
L @>>> D
\end{CD}
\qquad
\begin{CD}
L @>{m}>> G \\
@VVV @VVV \\
D @>>> H
\end{CD}
\]

- The first square constructs the **pushout complement** \(D\): it “removes” \(L\setminus K\) from \(G\) (subject to the dangling condition).
- The second square glues in \(R\setminus K\) to produce the product \(H\).

**Chemistry translation**
- \(K\) are the mapped atoms.
- \(L\setminus K\) are deleted atoms/bonds (leaving groups, bond cleavage).
- \(R\setminus K\) are created atoms/bonds (new substituents, bond formation).

In this notebook we implement the same effect algorithmically (delete-then-add), and additionally
construct an ITS graph that records the pre/post state in one object.



### 1.7 Practical correctness checks (what to verify in code)

When you implement rules in practice, you usually want to validate:

1. **Interface consistency:** \(K\) must embed into both \(L\) and \(R\).  
   (Here we enforce this by using shared node ids for \(K\) inside \(L\) and \(R\).)

2. **Dangling condition (node deletion):** if a host atom is deleted, it must not have bonds to atoms outside the matched subgraph.

3. **Typing consistency:** labels used in `node_match`/`edge_match` must be consistent across RDKit graphs and rule graphs.

We will implement (2) as an optional filter for matches (useful for leaving group rules).



## 2. RDKit ⇄ NetworkX: typed molecular graphs

We reuse the S01 “typed molecular graph” representation.

**Node attributes** (atoms):
- `symbol`, `formal_charge`, `aromatic`, `chiral_tag`
- optionally `total_h`

**Edge attributes** (bonds):
- `order` (1/2/3), `aromatic`, `in_ring`

This is enough for most pedagogical examples.



## 3. DPO rewriting kernel (chemistry-friendly)

### 3.1 Why chemistry rules need wildcards

A “reaction rule” almost never wants to pin down *every* atom attribute.
Examples:
- a Diels–Alder rule wants “carbon, non-aromatic”, but usually not fixed `chiral_tag`.
- generic rules often allow any `formal_charge` on spectator atoms.
- aromaticity handling depends on the representation (Kekulé vs aromatic bonds).

So we use **pattern-side wildcards**:
- if an attribute is missing in the pattern, we ignore it;
- if present but `None`, we also ignore it;
- otherwise we require equality.

This corresponds to defining a compatibility predicate
\(\Phi_V(d_G, d_L)\) and \(\Phi_E(d_G, d_L)\).



### 3.2 Matches in NetworkX: beware mapping direction

If you build:

```python
GM = GraphMatcher(G, L)
```

then `subgraph_isomorphisms_iter()` yields dictionaries mapping

\[
V(G)\to V(L)
\]

i.e. **host node → pattern node**.

For rewriting we want **pattern node → host node**, so we invert each mapping.



## 4. ITS-first rewriting (DPO-like)

### 4.1 What is an ITS graph here?

We define an ITS-like graph \(T\) that stores **both** pre- and post- bond types
on the *same atom identities*.

- Nodes represent atom identities (plus possibly created/deleted atoms).
- Each edge carries a pair:
  \[
  (b_{\text{pre}}, b_{\text{post}}),
  \]
  where \(b=0\) means “no bond”.

Chemistry reading:
- formed bond: \(0 \to 1\) (or \(0\to 2\))
- broken bond: \(1 \to 0\) (or \(2\to 0\))
- order change: \(1 \to 2\) etc.

### 4.2 Why construct ITS directly?

If you construct a product graph first, you must still maintain an atom mapping.
In DPO, the “identity backbone” is precisely \(K\), so the natural carrier object is:

> **the change graph** on conserved identities, plus any created/deleted identities.

By building ITS directly, the mapping is trivial:
- reactant atom \(v\) maps to product atom \(v\) **by node id**.

This makes export to atom-mapped SMARTS clean and deterministic.



## 5. Export: atom-mapped reaction SMARTS (ITS-native)

### 5.1 DPO-like “identity = node id”

Because we build ITS \(T\) on conserved identities, a node \(v\) refers to:

- the same atom in the reactant graph \(\pi_{\text{pre}}(T)\),
- the same atom in the product graph \(\pi_{\text{post}}(T)\),

whenever \(v\) is present in both (`present_pre=True` and `present_post=True`).

So atom mapping becomes **identity-by-node-id**:
\[
\text{map}(v) := v+1.
\]

New atoms (if created) naturally get fresh ids \(> \max(V_G)\), therefore also fresh map numbers.

### 5.2 Implementation

We:
1. project ITS → reactant graph \(G\) and product graph \(H\),
2. convert each graph to RDKit Mol,
3. set atom maps by graph node id,
4. export `MolToSmarts(reactant)>>MolToSmarts(product)`.

This avoids any post-hoc “mapping inference”.



### 5.3 Practical detail: RDKit atom reordering vs. mapping stability

RDKit may reorder atoms internally during sanitization or SMILES/SMARTS generation.
That is **not a problem** as long as:

- each atom carries its `AtomMapNum`,
- the same identity uses the same map number on both sides.

Because we set `map = node_id + 1` after converting from graphs, the mapping remains stable even if RDKit
chooses a different atom ordering for the final SMARTS string.



## 6. Worked example: Diels–Alder as a DPO rule

### 6.1 Chemistry view

Diels–Alder (4+2 cycloaddition) is a **concerted** pericyclic reaction:

- two \(\pi\) bonds disappear,
- two new \(\sigma\) bonds appear,
- one \(\pi\) bond remains in the cyclohexene product.

We encode this as a bond-rewriting rule that conserves all 6 carbon atoms:

- \(V_L = V_K = V_R\) (atom-conserving)
- only edges/bond orders change

This is the “cleanest” DPO case, and ideal for educational purposes.



### 6.2 Apply rule → ITS → (reactant, product) → mapped SMARTS

We:
- convert RDKit reactant mixture to a typed graph \(G\),
- apply the rule to get ITS graphs \(T\),
- decompose the first ITS to reactant/product,
- export mapped reaction SMARTS.



**Sanity check:** Diels–Alder is atom-conserving (no node deletion), so the dangling condition is trivially satisfied.
For leaving-group rules, `strict_dangling=True` helps prevent chemically nonsensical deletions that cut bonds outside the matched substructure.


### 6.3 Inspect the ITS: formed / broken / changed bonds

For each ITS edge \(uv\), we look at \((\text{pre\_order}, \text{post\_order})\).

- formed: \(0 \to >0\)
- broken: \(>0 \to 0\)
- order change: \(>0 \to >0\) but different



### 6.3b A compact ITS “transition fingerprint”

A very lightweight mapping-free reaction descriptor is the multiset of bond transitions:

\[
\{(b_{\text{pre}}, b_{\text{post}})\}_{uv\in E(T)}.
\]

For example:
- many \(2\to 1\) transitions indicate \(\pi\)-bonds turning into \(\sigma\)-bonds,
- \(0\to 1\) counts are bond formations,
- \(1\to 0\) counts are bond breakages.

This is not as expressive as full ITS isomorphism, but it is cheap and often useful for clustering/filtering.



Chemistry check (Diels–Alder expectation):

- **broken**: the dienophile double bond \(e=f\) becomes single (often represented as an order change),
  and one diene double bond changes.
- **formed**: two new C–C \(\sigma\) bonds close the ring.

Because our minimal rule encodes “remove old double bonds, add new single bonds”, you should see:
- two **formed** edges,
- and **order_change** edges accounting for \(\pi\to\sigma\) shifts.



### 6.4 Verify projections: ITS → reactant graph and product graph

We can reconstruct RDKit molecules from the projected graphs and compare SMILES.



## 7. Technical notes (what you will want in real chemistry pipelines)

### 7.1 Aromaticity and Kekulé form
RDKit may represent aromatic systems with aromatic bonds (`BondType.AROMATIC`) or Kekulé double/single.
A rule library must decide on one representation:
- either match on `aromatic=True` bonds/nodes,
- or kekulize all molecules before rewriting.

### 7.2 Stereochemistry (endo/exo, R/S)
Pericyclic reactions have stereochemical outcomes. Encoding this requires:
- stereo atom attributes (chiral tags),
- double bond stereo (E/Z),
- and rule-side constraints (not just plain graph rewriting).

### 7.3 Atom addition/deletion (leaving groups, proton transfers)
Many reactions are not atom-conserving in the heavy-atom graph, especially if you omit spectators.
The ITS design here supports this via:
- `present_pre`, `present_post` on nodes,
- `pre_order=0` or `post_order=0` on edges.

### 7.4 Multiple matches and automorphisms
Symmetric patterns lead to multiple subgraph isomorphisms.
For chemistry, you often want **unique products**:
- we deduplicate ITS graphs by isomorphism on `(pre,post)` edge labels.



## 8. Exercises

Try these without changing the core kernels.

### Exercise 1 — Diels–Alder variants
Change the Diels–Alder rule so that the remaining double bond is **a–b** instead of **b–c**.
- How do the ITS “order_change” edges differ?
- Does RDKit sanitize the product?

### Exercise 2 — Add a wildcard
Modify matching keys so that the rule matches carbon atoms regardless of `formal_charge`.
(Hint: omit `formal_charge` from `node_keys`.)

### Exercise 3 — Build a tiny rule for “halogen substitution”
Encode a rule that replaces `C–Cl` with `C–O` (as in your earlier example) and export mapped SMARTS via ITS.
- Verify the mapping is identity-by-node-id for conserved atoms.



<details>
<summary><b>Solutions (sketch)</b></summary>

**Ex1.** In `diels_alder_rule_minimal`, change which edge in `R` is double.
You must also ensure the corresponding edges in `L` that should disappear are present (and not in `K`).

**Ex2.** Call:

```python
Ts = apply_rule_to_its(G, rule_da, node_keys=("symbol","aromatic"), edge_keys=("order","aromatic"))
```

**Ex3.** Use \(K\) with the carbon node only; in ITS, carbon keeps its node id so map numbers are consistent.

</details>
