# S01: From Molecules to Labeled Graphs

This talktorial introduces molecular representations for SynEdu: RDKit molecules, SMILES strings, and explicit NetworkX labeled graphs. The focus is on what chemical information is preserved, what can be lost, and why representation choices matter for later graph-matching and reaction-rule tasks [\[1\]](#6.-References), [\[2\]](#6.-References), [\[3\]](#6.-References), [\[4\]](#6.-References).



## Aim of this talktorial

1. Build a compact **RDKit** foundation for parsing, inspecting, and canonicalizing molecules.
2. Understand **SMILES** as a graph encoding with choices around aromaticity, hydrogens, and stereochemistry.
3. Convert molecules into **labeled molecular graphs** and test round-trip behavior between RDKit and NetworkX.

---

## Learning outcomes

After completing this talktorial you will be able to:

- Parse SMILES with **RDKit** and understand what sanitization does.  
- Produce canonical SMILES.
- Convert SMILES into **labeled molecular graphs** (atoms → nodes, bonds → edges).  
- Perform a **round-trip conversion** between RDKit and NetworkX and identify which chemical details are preserved or lost.  
- Explain why **symmetry** and **label design** are critical for reaction rule discovery and application.

---

## Outline

<ul class="synedu-outline">
  <li><a href="#0-setup--data">0. Setup &amp; data</a></li>
  <li><a href="#s01-rdkit">1. Introduction to RDKit</a></li>
  <li><a href="#s01-smiles">2. SMILES</a></li>
  <li><a href="#s01-graph">3. Molecular Graph representation</a></li>
  <li><a href="#s01-discussion">4. Discussion</a></li>
  <li><a href="#5-quiz">5. Quiz</a></li>
  <li><a href="#6.-References">6. References</a></li>
</ul>


<a id="0-setup--data"></a>

## 0. Setup & data


<a id="s01-rdkit"></a>

## 1. Introduction to RDKit


Use `IPythonConsole` to render RDKit atom indices on molecular structures.


RDKit’s `Draw.MolsToGridImage` allows you to display a list of molecules
as a grid for rapid visual inspection and comparison.


We now explore basic information encoded in the RDKit representation of osimertinib.


Using `rdkit.Chem.Descriptors`, we can compute common physicochemical
properties directly from the RDKit molecular representation of osimertinib.


**Q1 — LogP**
 
Compute the octanol–water partition coefficient (logP) of osimertinib using RDKit.

<details> <summary><b>Solution</b></summary>

```python
from rdkit.Chem.Crippen import MolLogP
logp = MolLogP(osimertinib)
```
</details>


In medicinal chemistry, it is important to prioritize the most promising molecules in order to reduce experimental and computational costs. One of the simplest and most widely used heuristics for early-stage compound selection is [Lipinski’s Rule of Five](https://en.wikipedia.org/wiki/Lipinski%27s_rule_of_five), which describes four physicochemical criteria associated with favorable oral bioavailability.

## Lipinski’s Rule of Five

<div style="border:1px solid #d0d7de; padding:1rem; border-radius:6px;">

**Lipinski’s Rule of Five** defines a set of empirical physicochemical criteria commonly used to assess whether a small molecule is likely to exhibit acceptable oral drug-like properties.

</div>

<br>

<table>
  <thead>
    <tr>
      <th align="left">Parameter</th>
      <th align="left">Recommended threshold</th>
      <th align="left">Physicochemical rationale</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <td><strong>Hydrogen bond donors</strong></td>
      <td>≤ 5</td>
      <td>Excessive donor capacity may reduce passive membrane permeability.</td>
    </tr>
    <tr>
      <td><strong>Hydrogen bond acceptors</strong></td>
      <td>≤ 10</td>
      <td>High acceptor count is often associated with increased polarity and reduced intestinal absorption.</td>
    </tr>
    <tr>
      <td><strong>Molecular weight</strong></td>
      <td>&lt; 500 Da</td>
      <td>Larger molecules generally show reduced passive diffusion across biological membranes.</td>
    </tr>
    <tr>
      <td><strong>logP</strong></td>
      <td>&lt; 5</td>
      <td>Excessive lipophilicity may impair aqueous solubility and pharmacokinetic behavior.</td>
    </tr>
  </tbody>
</table>

<br>

<div style="border-left:3px solid #666; padding:0.75rem 1rem;">

<strong>Interpretation:</strong> Molecules satisfying these criteria are more likely to possess physicochemical properties compatible with oral absorption. However, the rule should be considered a prioritization heuristic rather than an absolute filter. Deviations may occur for specific chemical classes, including natural products, macrocycles, peptides, and compounds relying on active transport mechanisms.

</div>

<br>

<div class="hover-figure caption-hover" style="max-width:800px; width:100%; margin:0 auto;">
  <img
    src="../../docs/_static/images/talks/RO5.svg"
    alt="Visualization of Lipinski’s Rule of Five"
    style="width:100%; height:auto;"
  >
  <div class="figure-caption">
    Figure: Visualization of Lipinski’s Rule of Five, summarizing four physicochemical thresholds associated with oral drug-likeness.
  </div>
</div>


**Q2**
 
Compute the number of hydrogen bond donors and acceptors for each molecule in
`mols`.

<details> <summary><b>Solution</b></summary>

```python
from rdkit.Chem.Descriptors import NumHAcceptors, NumHDonors

nh_acc_mols = [NumHAcceptors(mol) for mol in mols]
nh_do_mols  = [NumHDonors(mol) for mol in mols]
```
</details>


We can now use the `Lipinski` function to filter molecules based on
Lipinski’s Rule of Five.


**Q3 — Lipinski filtering**
 
Apply the Lipinski filter to the dataset and annotate the DataFrame with a
pass/fail flag.

<details>
<summary><b>Solution</b></summary>

```python
from rdkit import Chem

df["mol"] = df["smiles"].apply(Chem.MolFromSmiles)
df["Lipinski"] = df["mol"].apply(passes_lipinski)

```
</details>


### Dataset-level property distributions

The four Lipinski properties computed across the full 1,000-molecule dataset.
Dashed black lines mark the Rule of Five thresholds; bars are coloured by pass/fail status.


<a id="s01-smiles"></a>

## 2. SMILES

<div style="padding: 1rem; border-left: 6px solid #4D96FF; background: #F3F8FF; border-radius: 10px;">

**SMILES** stands for **Simplified Molecular Input Line Entry System** [\[2\]](#6.-References).

It is a compact, human-readable text notation that encodes a molecular **graph**, where:

| Graph concept | Molecular meaning | SMILES representation |
|---|---|---|
| **Nodes** | Atoms | `C`, `N`, `O`, `Cl`, `[NH4+]` |
| **Edges** | Bonds | implicit single bonds, `=`, `#` |
| **Branches** | Side chains | `(...)` |
| **Cycles** | Rings | matching digits such as `1...1` |
| **Geometry** | Stereochemistry | `@`, `@@`, `/`, `\` |

</div>


---




### 2.1. Overview

SMILES encodes molecular structure using a small set of compact symbols.

<div style="display: grid; grid-template-columns: repeat(2, minmax(0, 1fr)); gap: 1rem; margin-top: 1rem;">

<div style="padding: 1rem; border: 1px solid #DDD; border-radius: 10px; background: #FAFAFA;">
<b>Atoms</b>

<code>C</code>, <code>N</code>, <code>O</code>, <code>Cl</code>

<p>Atoms are written using element symbols.</p>
</div>

<div style="padding: 1rem; border: 1px solid #DDD; border-radius: 10px; background: #FAFAFA;">
<b>Bonds</b>

<code>CC</code>, <code>C=C</code>, <code>C#N</code>

<p>Single bonds are usually implicit.</p>
</div>

<div style="padding: 1rem; border: 1px solid #DDD; border-radius: 10px; background: #FAFAFA;">
<b>Branches</b>

<code>CC(O)C</code>

<p>Parentheses create side chains.</p>
</div>

<div style="padding: 1rem; border: 1px solid #DDD; border-radius: 10px; background: #FAFAFA;">
<b>Rings</b>

<code>C1CCCCC1</code>

<p>Matching digits close rings.</p>
</div>

</div>

<br>

#### Token legend

| Token type | Example | Meaning |
|---|---|---|
| <span style="color:#4D4D4D;"><b>Aliphatic atom</b></span> | `C`, `N`, `O`, `Cl` | Standard non-aromatic atoms |
| <span style="color:#1F77B4;"><b>Aromatic atom</b></span> | `c`, `n`, `o` | Aromatic atoms, usually lowercase |
| <span style="color:#D62728;"><b>Double bond</b></span> | `=` | Explicit double bond |
| <span style="color:#8C564B;"><b>Triple bond</b></span> | `#` | Explicit triple bond |
| <span style="color:#FF7F0E;"><b>Branch</b></span> | `(...)` | Side chain |
| <span style="color:#E377C2;"><b>Ring closure</b></span> | `1`, `2`, `%10` | Connects two atoms to form a ring |
| <span style="color:#9467BD;"><b>Stereochemistry</b></span> | `@`, `@@`, `/`, `\` | Chirality or double-bond geometry |
| <span style="color:#2CA02C;"><b>Bracketed atom</b></span> | `[NH4+]`, `[13C]`, `[O-]` | Explicit hydrogens, isotopes, charges, or uncommon valence |

<div style="padding: 1rem; border-left: 5px solid #999; background: #F7F7F7; border-radius: 8px;">

**Key rule:**  
SMILES is not just a string. It is a compact way to encode a molecular graph.

</div>


---




### 2.2. Atoms & bonds

Atoms are represented by element symbols such as `C`, `N`, `O`, and `Cl`.

Single bonds are usually **implicit**, while double and triple bonds are written explicitly.

| SMILES | Visual interpretation | Meaning |
|---|---|---|
| `C` | C | methane carbon with implicit hydrogens |
| `CCO` | C–C–O | ethanol-like fragment |
| `C=O` | C=O | carbonyl double bond |
| `C#N` | C≡N | triple bond |
| `[Na+]` | Na⁺ | explicit sodium cation |

```text
C      # methane carbon with implicit hydrogens
CCO    # ethanol fragment, C–C–O
C=O    # formaldehyde-like carbonyl
C#N    # hydrogen cyanide or nitrile-like triple bond
[Na+]  # explicit sodium cation
```

<div class="hover-figure caption-hover">
  <img src="../../docs/_static/images/S01/carbonyl.svg" alt="Carbonyl group — C=O double bond">
  <div class="figure-caption">Figure: Carbonyl group — C=O double bond</div>
</div>



### 2.3. Branches

Parentheses create branches or side chains.

```text
CCO       # linear chain
CC(O)C    # branch: hydroxyl on the middle carbon
C(C)(C)C  # tert-butyl-like branching
```
<div class="hover-figure caption-hover">
  <img src="../../docs/_static/images/S01/ethanol.svg" alt="Ethanol — branch notation CC(O)C">
  <div class="figure-caption">Figure: Ethanol — branch notation CC(O)C</div>
</div>



### 2.4. Rings

Ring closures use matching digits to connect two positions in the SMILES string.

```text
C1CCCCC1      # cyclohexane
c1ccccc1      # benzene, aromatic SMILES
C1=CC=CC=C1   # benzene-like Kekulé notation
```
For ring labels above `9`, use `%`.

```text
C%10CCCCCCCCC%10
```

<div class="hover-figure caption-hover">
  <img src="../../docs/_static/images/S01/cyclohexane.svg" alt="Cyclohexane — ring closure C1CCCCC1">
  <div class="figure-caption">Figure: Cyclohexane — ring closure C1CCCCC1</div>
</div>



### 2.5. Aromaticity

Aromatic atoms are commonly written using lowercase letters.

| SMILES | Meaning |
|---|---|
| `c1ccccc1` | benzene |
| `n1ccccc1` | pyridine-like aromatic ring |

```text
c1ccccc1   # benzene, aromatic
n1ccccc1   # pyridine-like aromatic ring
```
<div class="hover-figure caption-hover">
  <img src="../../docs/_static/images/S01/benzene.svg" alt="Benzene — aromatic SMILES c1ccccc1">
  <div class="figure-caption">Figure: Benzene — aromatic SMILES c1ccccc1</div>
</div>



### 2.6. Stereochemistry

SMILES can encode stereochemical information.

| Notation | Used for | Example |
|---|---|---|
| `@` | tetrahedral chirality | <code>C&#91;C@H&#93;(O)Cl</code> |
| `@@` | opposite tetrahedral configuration | <code>C&#91;C@@H&#93;(O)Cl</code> |
| `/` and `\` | double-bond geometry | <code>C/C=C\C</code> |

```text
C&#91;C@H&#93;(O)Cl      # chiral centre with explicit stereochemistry
C/C=C\C          # defined double-bond geometry
```
RDKit preserves stereochemical flags when parsing and writing SMILES with:

```python
Chem.MolToSmiles(mol, isomericSmiles=True)
```

<div class="hover-figure caption-hover">
  <img src="../../docs/_static/images/S01/chiral_example.svg" alt="Chiral centre — tetrahedral stereo C&#91;C@H&#93;(O)Cl">
  <div class="figure-caption">Figure: Chiral centre — tetrahedral stereo C&#91;C@H&#93;(O)Cl</div>
</div>



### 2.7. Bracketed atoms & special cases

Brackets are used when the default SMILES rules are not enough.

Use brackets for:

| Case | Examples |
|---|---|
| explicit hydrogens | `[NH2]`, `[nH]` |
| charges | `[O-]`, `[NH4+]` |
| isotopes | `[13C]` |
| unusual valence states | `[Fe+2]`, `[Se]` |
| explicit atom specification | `[Cl]`, `[Na+]` |

```text
c1ccccc1[NH2]     # aniline with explicit NH2
c1cc[nH]c1        # aromatic nitrogen with explicit H
CC(=O)[O-]        # acetate anion
[NH4+]            # ammonium cation
[NH4+].[Cl-]      # ammonium chloride, ionic pair
[13CH3]C(=O)O     # acetic acid with 13C-labelled methyl carbon
C[Se]C            # dimethyl selenide
```

<div class="hover-figure caption-hover">
  <img src="../../docs/_static/images/S01/acetate.svg" alt="Acetate anion — CC(=O)[O-]">
  <div class="figure-caption">Figure: Acetate anion — CC(=O)[O-]</div>
</div>



### 2.8. Common beginner mistakes

<details>
<summary><b>Click to expand common SMILES mistakes</b></summary>

<br>

| Mistake | Example | Why it matters |
|---|---|---|
| Confusing uppercase and lowercase atoms | `C` vs `c` | `C` is aliphatic, `c` is aromatic |
| Forgetting ring digit pairs | `C1CCCC` | ring closure is incomplete |
| Misplacing branches | `CC(O)C` vs `CCC(O)` | branches attach to the atom before `(` |
| Ignoring stereochemistry | <code>C&#91;C@H&#93;(O)Cl</code> vs <code>CC(O)Cl</code> | stereochemistry may change molecular identity |
| Overusing brackets | `[C][C][O]` | most common atoms do not need brackets |
| Missing charge notation | `[O-]` vs `O` | charged and neutral atoms are different |

</details>


---




### 2.9. Summary

<div style="padding: 1rem; border-left: 6px solid #4D96FF; background: #F3F8FF; border-radius: 10px;">

SMILES is compact because it relies on a small set of rules:

| Feature | Main notation |
|---|---|
| Atoms | `C`, `N`, `O`, `Cl` |
| Aromatic atoms | lowercase, e.g. `c`, `n` |
| Single bonds | usually implicit |
| Double bonds | `=` |
| Triple bonds | `#` |
| Branches | `(...)` |
| Rings | matching digits, e.g. `1...1` |
| Charges/isotopes/explicit H | brackets, e.g. `[NH4+]`, `[13C]` |
| Stereochemistry | `@`, `@@`, `/`, `\` |

</div>


Now experiment with some examples

The cell below color-codes each SMILES token by its syntactic role — useful for parsing unfamiliar strings.
Each color corresponds to a token class (aromatic atom, branch, ring closure, stereo, etc.).


SMILES strings can have many valid variants for the same molecule [\[5\]](#6.-References).



The
simplest and most reliable way to compare or normalise SMILES is to convert
them to a **canonical form** in RDKit [\[1\]](#6.-References), [\[5\]](#6.-References), [\[6\]](#6.-References).



Graph canonicalization (canonical SMILES, toolkit differences) will be covered in later talktorials (**S06**).



<a id="s01-graph"></a>

## 3. Molecular Graph representation

In computational reaction modeling, we represent molecules as **labeled graphs** so that any notion of
“matching” respects **chemical identity**, such as element type, charge, and bond order, rather than bare
connectivity alone [\[4\]](#6.-References).

A **labeled molecular graph** is a quadruple

$$
G = (V, E, a, b),
$$

where:

- **Vertices** $V$ represent **atoms**.
- **Edges** $E \subseteq \{\{u,v\}\mid u,v\in V,\ u\neq v\}$ represent **bonds**
  (finite, undirected, simple: no loops, no parallel edges).
- $a$: atom labelling function.
- $b$: bond labelling function.

We often write $V(G)$ and $E(G)$ for the vertex and edge sets of $G$. For a vertex $v\in V(G)$:

- neighbourhood:
  $$
  N_G(v)=\{w\in V(G)\mid vw\in E(G)\},
  $$
- degree:
  $$
  \deg_G(v)=|N_G(v)|.
  $$

These graph-theoretic notions correspond chemically to an atom’s bonded neighbours and coordination number,
abstracting away geometry while retaining connectivity.

---

### 3.1. Labels and chemical types

Attributes are encoded via explicit labelling maps

$$
\ell_V: V(G)\to L_V,\qquad \ell_E: E(G)\to L_E,
$$

where $L_V$ and $L_E$ are finite, non-empty label sets.
For molecular graphs, we use the chemistry-specific notation:

$$
a_G: V(G)\to L_V \quad\text{(atom labels)},\qquad
b_G: E(G)\to L_E \quad\text{(bond labels)}.
$$

Let $\mathcal{G}$ denote the class of all labelled molecular graphs equipped with $(a_G,b_G)$.
In chemistry, $a_G(v)$ encodes *what atom this is* (element, charge, aromaticity, hydrogen count, …), 
while $b_G(uv)$ encodes *what bond this is* (order, aromaticity, ring status, …).

All subsequent notions of equivalence, symmetry, and matching in this talktorial are defined **relative to these labels**.

---

### 3.2. Graph representations in practice

In SynEdu, **RDKit** and **NetworkX** play complementary roles:

- **RDKit** is the chemical authority: sanitization, valence rules, aromaticity perception, and canonicalization.
- **NetworkX** provides an explicit, inspectable graph representation used for matching, symmetry analysis,
  and later graph rewriting.

To ensure that graph-based operations remain chemically meaningful, we require a **reversible interface**
between the two representations:
a molecule converted from RDKit to a labeled NetworkX graph must be convertible back without loss of
the chemical information encoded in $(a_G,b_G)$.

This reversible interface forms the foundation for all later notions—graph isomorphism, automorphisms,
and eventually reaction rules—introduced in subsequent SynEdu notebooks.



**Definition (Labeled Molecular Graph).**  
A *labeled molecular graph* is a 4-tuple

$$
G = (V,\, E,\, \mathbf{a},\, \mathbf{b})
$$

where $V$ is the set of **heavy atoms** (nodes), $E \subseteq V \times V$ is the set of **covalent bonds** (edges), $\mathbf{a}: V \to \mathcal{A}$ assigns each atom a tuple of **node attributes** (element symbol, formal charge, aromaticity, hydrogen count), and $\mathbf{b}: E \to \mathcal{B}$ assigns each bond a tuple of **edge attributes** (bond order, aromaticity).

A *labeled graph morphism* $\varphi: G_1 \to G_2$ is a pair of maps $\varphi_V: V_1 \to V_2$, $\varphi_E: E_1 \to E_2$ that (i) preserves adjacency: $\varphi_E(\{u,v\}) = \{\varphi_V(u), \varphi_V(v)\}$, and (ii) preserves labels: $\mathbf{a}_2(\varphi_V(v)) = \mathbf{a}_1(v)$ and $\mathbf{b}_2(\varphi_E(e)) = \mathbf{b}_1(e)$ for all $v \in V_1$, $e \in E_1$.

**Remark.** The attribute schema used throughout SynEdu is  
`a(v) = (element, formal_charge, aromatic, hcount)` and `b(e) = (order, aromatic)`.  
Hydrogen atoms are stored implicitly in `hcount` to keep graphs small.



#### Graph attribute inspection

The DataFrames below show the concrete realization of the label functions $a_G$ (node attributes) and $b_G$ (edge attributes) for phenol.
This makes explicit what information is stored — and therefore what is available for downstream matching and rewriting.


> **Stored in `synedu.Utils`** — both `mol_to_graph` and `graph_to_mol` (defined above) are packaged in `synedu.Utils.conversion` so every later talktorial can import them without redefining them:
> ```python
> from synedu.Utils.conversion import mol_to_graph, graph_to_mol
> ```
> `smiles_to_graph` and `graph_to_smi` (thin wrappers around these) are available there too.


#### Kekulé vs aromatic — graph label differences

The same molecule (benzene) written in two SMILES notations produces **different bond-order labels** in the graph:
- **Kekulé** (`C1=CC=CC=C1`): alternating single (1.0) and double (2.0) bonds.
- **Aromatic** (`c1ccccc1`): uniform aromatic bonds (order 1.5), with an `aromatic=True` flag on all nodes and edges.

This distinction matters for graph isomorphism and reaction rule matching in later notebooks.


### 3.3. Matrix representations

A labeled molecular graph can be encoded as several **matrices**, each capturing a different structural aspect [\[4\]](#6.-References), [\[7\]](#6.-References), [\[8\]](#6.-References):

| Matrix | Shape | Entry meaning |
|---|---|---|
| **Adjacency** $A$ | $n \times n$ | $A_{ij}$ = bond order between atoms $i$ and $j$ |
| **Distance** $D$ | $n \times n$ | $D_{ij}$ = shortest-path length (in bonds) between $i$ and $j$ |
| **Incidence** $B$ | $n \times m$ | $B_{ij} = 1$ if atom $i$ participates in bond $j$ |
| **Bond-Electron** $BE$ | $n \times n$ | off-diag = bond order; **diagonal = free (non-bonding) electrons** |

The BE matrix generalises the adjacency matrix by placing lone-pair electron counts on the diagonal.
For a reaction, the **difference** $\Delta BE = BE_\text{products} - BE_\text{reactants}$ encodes the full electron reorganisation — this is developed in **S04**.


#### Adjacency

The **adjacency matrix** $A \in \mathbb{R}^{n \times n}$ encodes pairwise atom connectivity.
Three variants are informative:

| Panel | Entry $A_{ij}$ | When useful |
|---|---|---|
| **Binary** | 1 if bond exists, else 0 | Connectivity-only algorithms (BFS, diameter) |
| **Weighted — aromatic** | Bond order from aromatic SMILES (1.5 for ring bonds) | Aromatic-aware descriptors |
| **Weighted — Kekulé** | Bond order from Kekulé SMILES (1.0 or 2.0) | Reaction-rule matching, BE matrix |

Key properties shared by all variants:
- **Symmetric**: $A_{ij} = A_{ji}$ (undirected bonds)
- **Zero diagonal**: no self-loops ($A_{ii} = 0$)
- Diagonal entries of $A + D$ give the degree (valence) of each atom, where $D$ is the degree matrix


#### Distance matrix

The **topological distance matrix** $D \in \mathbb{N}^{n \times n}$ records the
shortest-path length (measured in bonds) between every pair of atoms.
$D_{ij}$ is the minimum number of bonds to traverse to reach atom $j$ from atom $i$,
computed via Floyd–Warshall or breadth-first search on the unweighted graph.

Key properties:
- **Zero diagonal**: $D_{ii} = 0$
- **Symmetric**: $D_{ij} = D_{ji}$
- The **graph diameter** $\max_{i,j} D_{ij}$ is the longest shortest path — a
  compact measure of molecular "stretch"
- Distance-based **Wiener index** $W = \tfrac{1}{2}\sum_{i,j} D_{ij}$ correlates
  with boiling points for alkanes [\[7\]](#6.-References)


#### Incidence matrix

The **node-edge incidence matrix** $B \in \{0,1\}^{n \times m}$ maps atoms to bonds:
$B_{ij} = 1$ if atom $i$ participates in bond $j$, and 0 otherwise.

Key properties:
- Each **column** has exactly two 1s (every bond connects exactly two atoms)
- Each **row sum** equals the heavy-atom degree (valence) of that atom
- $B B^\top = \Delta + A$, where $\Delta$ is the diagonal degree matrix and $A$ is
  the binary adjacency matrix — a fundamental identity in algebraic graph theory
- Incidence matrices appear in spectral graph theory and in the cycle-space
  formulation of Kirchhoff's current laws


#### Bond-electron matrix

The **Bond-Electron (BE) matrix** $M \in \mathbb{R}^{n \times n}$ extends the
weighted adjacency matrix by encoding electron counts on the diagonal [\[8\]](#6.-References):

$$
M_{ij} = \begin{cases}
  b_{ij} & i \neq j \quad (\text{bond order between atoms } i \text{ and } j) \\
  v_i - q_i - \displaystyle\sum_j b_{ij} - h_i & i = j \quad (\text{free / non-bonding electrons on atom } i)
\end{cases}
$$

where $v_i$ is the valence electron count, $q_i$ the formal charge, and $h_i$ the
implicit hydrogen count.

**Aromatic vs Kekulé** — two panels below highlight a key practical concern:
- In the **aromatic** form, ring bond orders are 1.5, giving fractional diagonal values
  that have no physical meaning.
- In the **Kekulé** form, all bond orders are integers, and the diagonal is a
  well-defined electron count.

For this reason, reaction-informatics tools always work with the **Kekulé BE matrix**.
The reaction-level version — the **ΔBE matrix** — is developed in **S04**.


> **Stored in `synedu.Utils`** — `build_be_matrix` is also packaged for later graph-representation tasks:
> ```python
> from synedu.Utils.graph import build_be_matrix
> ```
> Later notebooks can reuse the same bond-electron matrix convention instead of redefining it.



**Q4 — Conversion**

**Goal**
Convert the SMILES below to a graph and back to SMILES.  
Explain **why the round-trip fails**.

```python
smiles = "c1cc[nH]c1"
```

<br>

<details>
<summary><b>Solution</b></summary>

```python
from rdkit import Chem

def naive_roundtrip(smiles: str) -> str:
    mol = Chem.MolFromSmiles(smiles)
    G = mol_to_graph(mol)
    mol2 = graph_to_mol(G)
    return Chem.MolToSmiles(mol2, canonical=True)

naive_roundtrip("c1cc[nH]c1")
>> 'c1ccnc1'
```
Why it fails

- `[nH]` is an aromatic nitrogen with an **explicit hydrogen**
- **Heavy-atom graphs discard hydrogens**
- RDKit cannot recover this on rebuild: `[nH] → [n]`
- This is **information loss by representation**, not a bug

> **Conclusion:** Heavy-atom graphs are **not information-complete** for SMILES round-trip.


</details> 


Since the hydrogen count is already stored in the graph (`hcount`), it can be
propagated back during molecule reconstruction:

```python
for node, data in G.nodes(data=True):
    n_h = int(data.get("hcount", 0))
    heavy_atom.SetNoImplicit(True)
    heavy_atom.SetNumExplicitHs(n_h)
```
This logic is exposed via:
```python
graph_to_mol(G, use_hcount=True)
```
With `use_hcount=True`, explicit hydrogens (e.g. `[nH]`) are preserved, making the
SMILES round-trip reversible.



**Q5 — Round-trip conversion**

**Goal.** Verify that SMILES in a DataFrame survive a  
SMILES → RDKit → graph → RDKit round-trip **without losing heavy-atom information**.

**Input.** A DataFrame `df` with:
- `smiles`: SMILES string  
- `name`: molecule name  

**Task.**
1. Implement `roundtrip_smiles_equal(smiles)` that:
   - parses a SMILES with RDKit,
   - converts it to a labeled graph (`mol_to_graph`),
   - reconstructs a molecule (`graph_to_mol`),
   - compares **canonical heavy-atom SMILES** (hydrogens removed).
2. Apply it to `df["smiles"]`.
3. Report molecules that fail the round-trip (if any).

---

<details>
<summary><b>Solution</b></summary>

```python
def standardize(smiles):
    return Chem.CanonSmiles(smiles, useChiral=False)


def roundtrip_smiles_equal(smiles: str) -> bool:
    new_smiles = roundtrip(smiles)
    return (
        standardize(new_smiles) == standardize(smiles)) 

df["ok"] = df["smiles"].apply(roundtrip_smiles_equal)
df.loc[~df["ok"], ["name", "smiles"]]
```
</details> 


<a id="s01-discussion"></a>

## 4. Discussion

- A **labeled graph morphism** provides a precise abstraction for
  structure- and attribute-preserving mappings, making the notion of
  “sameness” explicit and inspectable.

- Our **labeled molecular graph** adopts a deliberately minimal attribute
  schema (`symbol`, `formal_charge`, `aromatic`, `order`). This strikes a
  balance between chemical faithfulness and algorithmic tractability:
  too few labels induce spurious symmetries, while too many hinder matching
  and reuse.

- **Round-trip conversion** (RDKit → NetworkX → RDKit) serves as a practical
  validation tool. While exact RDKit internal states need not be preserved,
  maintaining **heavy-atom topology and labels** ensures semantic
  equivalence and reproducibility.


<a id="5-quiz"></a>

## 5. Quiz

Answer briefly using what you learned about **RDKit**, **SMILES**, and **molecular graphs**.

1. What does RDKit molecule sanitization check, and why is it useful before graph conversion?
2. Why can one molecule have multiple valid SMILES strings, and what problem does canonical SMILES solve?
3. In a molecular graph, what do nodes and edges represent? Name two atom labels and one bond label that matter for matching.
4. When converting RDKit molecules to NetworkX graphs, what information is preserved, what may be simplified, and why is that acceptable for the later SynEdu tasks?



## 6. References

1. RDKit documentation. https://www.rdkit.org/docs/
2. Weininger, D. SMILES, a chemical language and information system. 1. Introduction to methodology and encoding rules. *Journal of Chemical Information and Computer Sciences* **28**, 31-36 (1988). https://doi.org/10.1021/ci00057a005
3. NetworkX documentation. https://networkx.org/documentation/stable/
4. Bonchev, D.; Rouvray, D. H., eds. *Chemical Graph Theory: Introduction and Fundamentals*. Abacus Press (1991).
5. Weininger, D.; Weininger, A.; Weininger, J. L. SMILES. 2. Algorithm for generation of unique SMILES notation. *Journal of Chemical Information and Computer Sciences* **29**, 97-101 (1989). https://doi.org/10.1021/ci00062a008
6. Morgan, H. L. The generation of a unique machine description for chemical structures: a technique developed at Chemical Abstracts Service. *Journal of Chemical Documentation* **5**, 107-113 (1965). https://doi.org/10.1021/c160017a018
7. Wiener, H. Structural determination of paraffin boiling points. *Journal of the American Chemical Society* **69**, 17-20 (1947). https://doi.org/10.1021/ja01193a005
8. Dugundji, J.; Ugi, I. An algebraic model of constitutional chemistry as a basis for chemical computer programs. *Topics in Current Chemistry* **39**, 19-64 (1973).
