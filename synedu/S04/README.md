# S03 · Atom Mapping as Graph Morphism: MCS, RXNMapper, and ITS Equivalence

<div class="alert alert-block alert-info">
<b>Welcome to S03.</b><br>
This talktorial extends <b>S01–S02</b> from single-molecule matching to the
<b>reaction setting</b>. We move from pattern alignment to
<b>atom mapping</b>, treating mappings as explicit graph morphisms between
reaction sides.
</div>

<div class="alert alert-block alert-success">
<b>What you will gain.</b><br>
You will learn how to construct atom maps from MCS alignments, understand
their failure modes, apply a modern ML-based mapper (RXNMapper), and compare
atom maps rigorously using <b>ITS graph isomorphism</b>.
</div>

<div class="alert alert-block alert-warning">
<b>Prerequisites.</b><br>
You should be comfortable with:
<ul>
  <li>Typed molecular graphs and morphisms (S01)</li>
  <li>Subgraph isomorphism, symmetry, and MCS (S02)</li>
  <li>Basic reaction SMILES notation</li>
</ul>
</div>



## Aim of this talktorial

In **S01**, we formalized molecules as typed graphs and defined
label-preserving morphisms.
In **S02**, we studied subgraph matching and MCS as alignment primitives.

This talktorial (**S03**) brings these ideas into the **reaction domain**.
We view a chemical reaction as two multisets of molecular graphs
(reactants → products) and treat **atom mapping** as a
*label-preserving partial isomorphism* between them.

Concretely, we focus on:

1. **MCS-based atom mapping**  
   Building an atom correspondence by aligning maximum common substructures
   between reactants and products.

2. **Failure modes of naive MCS mapping**  
   Understanding why MCS alone is insufficient:
   symmetry, multi-component ambiguity, and unbalanced reactions.

3. **Attention-guided atom mapping (RXNMapper)**  
   Running RXNMapper to obtain mapped reaction SMILES and extracting atom maps.

4. **ITS graph construction**  
   Encoding a mapped reaction as an **Imaginary Transition State (ITS)** graph.

5. **Map comparison via graph isomorphism**  
   Comparing different atom maps by testing **ITS isomorphism**, making the
   comparison invariant to atom-map IDs.

---

## Learning outcomes

After completing this talktorial, you will be able to:

- Interpret an **atom map** as a label-preserving (partial) graph morphism.
- Construct a simple **MCS-based atom map** in code.
- Identify when and why MCS-based mapping fails.
- Run **RXNMapper** and parse mapped reaction SMILES.
- Build an **ITS graph** from a mapped reaction.
- Decide whether two atom maps are equivalent by checking
  **ITS graph isomorphism**.

---

## Outline

0. **Setup & reaction data**
1. **Atom mapping as a graph morphism**
2. **MCS-based atom mapping (baseline)**
3. **Failure cases: symmetry, components, imbalance**
4. **RXNMapper: attention-based atom mapping**
5. **ITS construction from mapped reactions**
6. **Comparing atom maps via ITS isomorphism**
7. **Discussion, quiz, and references**


## 0. Setup & Data

We use:
- **RDKit** for molecules + MCS (`rdFMCS`)
- **NetworkX** for ITS graphs + isomorphism tests
- **rxnmapper** (optional) for a strong baseline mapper

If `rxnmapper` is not installed, the notebook still runs and will skip those parts.


## 1. Formal description

### 1.1 Molecular graphs (S01 recap)

In **S01**, a molecule is modeled as a **typed labeled graph**

$$
G = (V, E, a, b),
$$

with the following data:

$$
V \;\text{is a finite set of atoms},
\qquad
E \subseteq \{\{u,v\} \mid u,v \in V,\; u \neq v\},
$$

$$
a : V \longrightarrow \mathcal{A},
\qquad
\mathcal{A} = \text{atom label space},
$$

$$
b : E \longrightarrow \mathcal{B},
\qquad
\mathcal{B} = \text{bond label space}.
$$

A **label-preserving graph morphism**

$$
\varphi : V(G) \longrightarrow V(H)
$$

satisfies:

$$
\text{(M1)} \quad
a_H(\varphi(v)) = a_G(v)
\quad \forall v \in V(G),
$$

$$
\text{(M2)} \quad
\{u,v\} \in E(G)
\;\Rightarrow\;
\{\varphi(u), \varphi(v)\} \in E(H),
$$

$$
\text{(M3)} \quad
b_H(\{\varphi(u), \varphi(v)\})
=
b_G(\{u,v\})
\quad \text{(optional)}.
$$

An **isomorphism** is a bijective morphism whose inverse is also a morphism.

---

### 1.2 Reactions and atom mapping as partial isomorphisms

A reaction SMILES is written as

$$
R \;\gg\; P,
$$

where \(R\) and \(P\) are multisets of molecules.
Each side is represented as a **disjoint union of molecular graphs**:

$$
G_R = \biguplus_i G_{R_i},
\qquad
G_P = \biguplus_j G_{P_j}.
$$

In general,

$$
G_R \not\cong G_P,
$$

because reactions may be unbalanced and some atoms may not participate.

An **atom map** is therefore modeled as a **partial isomorphism**:
find subgraphs

$$
G_R' \subseteq G_R,
\qquad
G_P' \subseteq G_P,
$$

and an isomorphism

$$
\varphi : V(G_R') \longrightarrow V(G_P')
$$

that maximizes

$$
|V(G_R')|.
$$

---

### 1.3 Maximum Common Substructure (MCS) and its role (S02)

As formalized in **S02**, a **Maximum Common Substructure (MCS)** between
\(G_R\) and \(G_P\) is a graph

$$
K = (V_K, E_K, a_K, b_K)
$$

together with injective morphisms

$$
\iota_R : V_K \hookrightarrow V(G_R),
\qquad
\iota_P : V_K \hookrightarrow V(G_P),
$$

such that

$$
K \in
\arg\max_{S}
\left(
|V(S)|
\;\middle|\;
S \hookrightarrow G_R,\;
S \hookrightarrow G_P
\right).
$$

The embeddings \(\iota_R\) and \(\iota_P\) induce an atom map

$$
\varphi
=
\iota_P \circ \iota_R^{-1}
:
V(G_R') \longrightarrow V(G_P'),
$$

where \(G_R' = \iota_R(K)\) and \(G_P' = \iota_P(K)\).

---

### 1.4 From MCS to ITS and map equivalence

Given an atom map \(\varphi\), the **Imaginary Transition State (ITS)** graph
encodes bond changes via paired labels

$$
b_{\mathrm{ITS}} : E_{\mathrm{ITS}} \longrightarrow \mathcal{B}_R \times \mathcal{B}_P.
$$

Two atom maps are considered **equivalent** if their corresponding ITS graphs
are **isomorphic** under label-preserving graph isomorphism, as defined in **S01**.



## 2. Utilities: parsing and RDKit helpers

Conventions:
- mapped atoms are those with `atom.GetAtomMapNum() > 0`
- our teaching MCS mapper allows **partial mapping** (unbalanced OK)


## 3. Atom mapping via MCS (teaching version)

### 3.1 Problem with the naive idea

A naive MCS mapper assumes **one reactant vs one product**.
It fails on:
- multi-component reactions (assignment problem),
- symmetry (many equally valid embeddings),
- unbalanced reactions (extra/missing atoms).

### 3.2 Patch for teaching

We extend the mapper to:
- choose a component assignment maximizing total MCS size,
- allow partial mapping (unmatched atoms remain unmapped),
- keep deterministic choices (first embedding) for reproducibility.


### 3.3 Demo set (balanced, multi-component, symmetry, unbalanced)


### 3.4 Symmetry deep dive: counting MCS embeddings

In symmetric systems, MCS has many matches (S01: automorphisms).
We count how many embeddings RDKit finds for benzene → chlorobenzene ring conservation.


## 4. Attention-guided atom mapping (RXNMapper)

RXNMapper is an *unsupervised* atom-mapping method that extracts atom correspondences from
attention patterns learned by a Transformer model pretrained on reaction SMILES with a
masked-language objective (ALBERT/BERT-style), i.e. without requiring gold atom maps.

Let a reaction be written as a token sequence
$$
x = (x_1,\dots,x_T),
\qquad x = R \; \gg \; P,
$$
where a subset of tokens correspond to atoms on the reactant side and product side.
During pretraining, a random subset of tokens is masked and the model is trained to recover them:
$$
\min_\theta \; \mathbb{E}\big[ -\log p_\theta(x_{\mathcal{M}} \mid x_{\setminus\mathcal{M}})\big].
$$

At inference time, we extract an attention-based alignment matrix between reactant-atom tokens
and product-atom tokens:
$$
A \in \mathbb{R}_{\ge 0}^{m \times n},
\qquad
A_{ij} := \text{Attn}(r_i \rightarrow p_j),
$$
where $r_i$ is the $i$-th reactant atom token and $p_j$ is the $j$-th product atom token
(typically using a selected layer/head and a scalar sharpening factor).

Atom mapping is then posed as a maximum-weight assignment:
$$
\mu = \arg\max_{\mu \in \Pi} \sum_{i=1}^{m} A_{i,\mu(i)},
$$
where $\Pi$ is the set of valid (partial) matchings subject to basic chemical constraints
(e.g., element-type consistency). Because $\mu$ is inferred from attention, bond-order
preservation is **not** imposed, so bond formation/cleavage and order changes are allowed.

RXNMapper returns:
- `mapped_rxn`: reaction SMILES annotated with atom-map indices induced by $\mu$,
- `confidence`: a scalar summary of alignment consistency (higher typically indicates more reliable maps).


We follow the Molecular Transformer framework for reaction modeling[^mt] and use RXNMapper for attention-guided atom mapping[^rxn].

[^mt]: Schwaller *et al.* (2019), *Molecular Transformer: A Model for Uncertainty-Calibrated Chemical Reaction Prediction*, ACS Cent. Sci. DOI: 10.1021/acscentsci.9b00576. https://pubs.acs.org/doi/10.1021/acscentsci.9b00576
[^rxn]: Schwaller *et al.* (2021), *Extraction of organic chemistry grammar from unsupervised learning of chemical reactions* (RXNMapper), Sci. Adv. DOI: 10.1126/sciadv.abe4166. https://www.science.org/doi/10.1126/sciadv.abe4166



## 5. ITS graph construction

ITS nodes are map numbers. ITS edges carry labels `(br, bp)`:
- `br`: bond order in reactants (0 if absent)
- `bp`: bond order in products (0 if absent)

Reaction center edges are those with `br != bp`.


## 6. Comparing atom maps via ITS isomorphism (map-number invariant)

We use `networkx.GraphMatcher` (S01) with:
- `node_match`: compare node labels (`R_*`, `P_*`)
- `edge_match`: compare `(br,bp)`


## 7. End-to-end table: MCS mapper vs RXNMapper

We compute:
- `mcs_mapped`
- `rxnmapper_mapped` + `conf` (if available)
- `ITS_iso(MCS,RXNMapper)` when both exist


## 8. Bonus: reaction center from ITS

Reaction center edges:
\[
RC = \{ (u,v)\in E(G_{ITS}) : b_R(u,v)\neq b_P(u,v)\}.
\]

This is a lightweight extractor that doesn't need templates.


## 9. Exercises

### Exercise 1 (S01 link): make the isomorphism “too permissive”
Modify `node_match` in `its_isomorphic` to ignore charge and compare results on `sn2_substitution`.

### Exercise 2: symmetry
Force a different MCS embedding for benzene chlorination (swap which match you pick) and show ITS-isomorphism stays `True`.

### Exercise 3: scaling
Why does brute-force component assignment explode? Estimate complexity for `n` components.

### Exercise 4: aromatic handling
Encode aromatic bonds with a special integer (e.g., 15) and see if comparisons change.


## 10. Takeaways

- Atom mapping can be formalized as a **partial isomorphism** between reaction-side graphs (S01 morphisms).
- MCS is a principled baseline but needs help for multi-component, symmetry, and unbalanced cases.
- RXNMapper provides a strong practical mapping + confidence.
- ITS graphs give a compact transformation representation.
- ITS isomorphism is a clean map-number-invariant way to compare atom maps.
