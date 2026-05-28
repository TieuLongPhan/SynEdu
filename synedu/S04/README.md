# S04: Atom Mapping as Graph Morphism

This talktorial connects molecular alignment to atom-to-atom mapping. We build transparent MCS-based maps, compare them with RXNMapper, and use Imaginary Transition State (ITS) graphs as a map-number-invariant representation of reaction change [\[1\]](#6.-References), [\[2\]](#6.-References).



## Aim of this talktorial

1. Construct **MCS-based atom maps** between reactants and products.
2. Run **RXNMapper** and inspect attention information for student-facing interpretation.
3. Build and compare **ITS graphs** so atom-map equivalence can be checked by graph isomorphism.

---

## Learning outcomes

After completing this talktorial, you will be able to:

- construct a simple **MCS-based atom map** in code,
- run **RXNMapper** and parse mapped reaction SMILES,
- extract and visualize mapper attention for student-facing interpretation,
- build an **ITS graph** from a mapped reaction,
- separate the reactant and product views of an ITS graph for visualization, and
- decide whether two atom maps are equivalent by checking **ITS graph isomorphism**.

---

## Outline

<ul class="synedu-outline">
  <li><a href="#0-setup--data">0. Setup &amp; Data</a></li>
  <li><a href="#1-atom-mapping">1. Atom mapping</a></li>
  <li><a href="#2-imaginary-transition-state">2. Imaginary Transition State</a></li>
  <li><a href="#3-its-equivalence">3. ITS equivalence</a></li>
  <li><a href="#4-discussion">4. Discussion</a></li>
  <li><a href="#5-quiz">5. Quiz</a></li>
  <li><a href="#6.-References">6. References</a></li>
</ul>


## 0. Setup & Data


## 1. Atom mapping
### 1.1. Alignment method


Now we convert them into reactant and product graphs using `rsmi_to_graph`.


We now develop `mcs_networkx` to identify the Maximum Common Substructure (MCS), using the graph-matching perspective introduced in S03 [\[3\]](#6.-References).
Note that in a reaction context, **bonds may be formed or broken** between
reactants and products. Therefore, bond attributes should **not** be used
as matching constraints (`edge_attrs`), and the MCS is computed based on
atom-level corr



Now we can just add `atom_maps` to node attribute


**Q1 — Assign atom maps from reactant–product alignment**

You are given a node mapping from a **reactant graph** `r` to a
**product graph** `p`, obtained from an MCS-based alignment.
Your task is to assign **atom-map numbers** to the corresponding atoms
on both sides of the reaction.

The convention is that each matched atom pair shares the **same
`atom_map` label**, which is taken from the reactant node index.

---

<details> <summary><b>Solution:</b></summary>

```python
def assign_atom_map_r_to_p(r, p, mapping):
    """
    Assign atom_map attributes on both reactant and product graphs
    using an r -> p node mapping.

    Parameters
    ----------
    r : nx.Graph
        Reactant graph
    p : nx.Graph
        Product graph
    mapping : dict[int, int]
        Mapping r_index -> p_index
    """
    for r_idx, p_idx in mapping.items():
        if r_idx not in r:
            raise KeyError(f"Node {r_idx} not found in reactant graph")
        if p_idx not in p:
            raise KeyError(f"Node {p_idx} not found in product graph")

        # assign the same atom-map label to matched atoms
        r.nodes[r_idx]["atom_map"] = r_idx
        p.nodes[p_idx]["atom_map"] = r_idx

```

</details>


#### Atom mapping colour coding

After MCS alignment, matched atom pairs share the same colour across the reactant and product graphs. Atoms that could not be matched (spectators or unmapped atoms) appear grey. Comparing the two panels reveals which atoms participate in the reaction center.


Now we combine into 1 function


**Q2 — Atom mapping via MCS: full vs partial**

Consider the following reaction SMILES:

```text
CCCl.[OH-]>>CCO.[Cl-]
c1ccccc1.ClCl>>c1ccccc1Cl.Cl
CCO>>C=C
```
For each reaction, can we recover a full atom map using MCS-based alignment, or only a partial atom map?
Explain why.

---

<details> <summary><b>Solution:</b></summary>

```python
aams = []
rsmis = ['CCCl.[OH-]>>CCO.[Cl-]', 'c1ccccc1.ClCl>>c1ccccc1Cl.Cl', 'CCO>>C=C']
for rsmi in rsmis:
    aam = mcs_aam(rsmi)
    print(aam)
    aams.append(aam)
    
```

Output

```text
['C([CH3:1])[Cl:3].[OH-:4]>>C([CH3:1])[OH:4].[Cl-:3]']
['Cl[Cl:7].[cH:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1>>Cl[c:6]1[cH:1][cH:2][cH:3][cH:4][cH:5]1.[ClH:7]', '[Cl:7][Cl:8].[cH:1]1[cH:2][cH:3][cH:4][cH:5][cH:6]1>>Cl[c:6]1[cH:1][cH:2][cH:3][cH:4][cH:5]1.[ClH:8]']
['O[CH2:2][CH3:1]>>[CH2:1]=[CH2:2]']
```

All three reactions yield **partial atom maps only**.

This is because:
- **Connectivity changes** (bond breaking/forming) prevent full graph isomorphism.
- **Stoichiometric imbalance** (missing species) leaves some atoms unmappable.

</details>


### 1.2. Attention-guided atom mapping

**RXNMapper** extracts atom correspondences from **Transformer attention** learned on reaction SMILES with a masked-language objective (BERT family). It is **unsupervised** and does not require gold atom maps.

Write a reaction as a token sequence $x=(x_1,\dots,x_T)$ with reactant and product atom tokens. During pretraining the model minimizes the masked-token objective:

$$
\min_\theta \; \mathbb{E}\!\left[-\log p_\theta\big(x_{\mathcal{M}}\mid x_{\setminus\mathcal{M}}\big)\right].
$$

At inference we extract an attention alignment matrix between reactant and product atom tokens:

$$
A \in \mathbb{R}_{\ge 0}^{m\times n},\qquad
A_{ij} := \mathrm{Attn}(r_i \to p_j),
$$

(typically from a chosen head/layer and optionally sharpened). Atom mapping is then posed as a maximum-weight assignment over admissible matchings $\Pi$:

$$
\mu = \arg\max_{\mu\in\Pi}\sum_{i=1}^{m} A_{i,\mu(i)},
$$

where $\Pi$ enforces lightweight chemical constraints (element consistency, optional charge checks). Because the assignment is driven by attention, **bond-order preservation is not enforced**, so RXNMapper tolerates bond formation/cleavage and order changes.

Typical RXNMapper outputs:

- `mapped_rxn`: reaction SMILES annotated with atom-map indices induced by \(\mu\),  
- `confidence`: a scalar summary of alignment consistency (higher → more reliable).

We adopt this attention-guided strategy following the Molecular Transformer framework [\[4\]](#6.-References) and RXNMapper [\[2\]](#6.-References).



#### RXNMapper attention alignment heatmap

The atom map is selected from an attention matrix whose rows are product atom tokens and whose columns are reactant atom tokens. Bright cells indicate product atoms that strongly attend to a reactant atom. The final atom map is the high-weight assignment through this matrix, so visualising it helps connect the Transformer output to the atom-map numbers we use later.



## 2. Imaginary Transition State
The **Imaginary Transition State (ITS)** [\[2\]](#6.-References) (or Condensed Graph of the Reaction [\[5\]](#6.-References)) is a compact, chemistry-oriented way to represent *what changes* in a reaction by **superimposing reactants and products via an atom-atom map**.

Think of the ITS as a single graph whose **nodes are atom-map labels** (one node per mapped atom) and whose **edges record the bond before and after the reaction**. Reading the ITS tells you, at a glance, which bonds are preserved, broken or formed.

---


Given a balanced reaction and an atom map $\alpha$ (a bijection between reactant and product atoms):

1. Let the ITS node set be the atom-map labels
   $$
   V=\{a_1,\dots,a_M\}
   $$

2. For each unordered pair of nodes \((i,j)\) read the bond order in reactants and in products:
   $$
   br_{ij}=\text{bond order of atoms mapped to }a_i,a_j\ \text{in reactants (or }0\text{ if absent)},
   $$
   $$
   bp_{ij}=\text{bond order of atoms mapped to }a_i,a_j\ \text{in products (or }0\text{ if absent)}.
   $$

3. Add an ITS edge between \(i\) and \(j\) labelled with the pair
   $$
   (br_{ij},\,bp_{ij}).
   $$

4. The **reaction center** is simply the set of edges with a change:
   $$
   E_{\mathrm{rc}}=\{(i,j):\; br_{ij}\ne bp_{ij}\}.
   $$

---

How to read the ITS (chemical rules of thumb)

- An edge labelled (1,1) — bond preserved (single → single).  
- (1,0) — bond **broken** (present in reactants, absent in products).  
- (0,1) — bond **formed** (new bond in products).  
- (2,1) or (1.5,1.5) — bond order **changed** (double→single, or aromatic preserved).  
- Nodes keep atom identity via labels (element, charge, hybridization if desired) so you can tell *which atom* changed connectivity.

---



**Definition (ITS Graph).**  
Given a reaction with atom mapping $\mu: V_R \to V_P$ (a bijection between reactant and product atoms), the *Imaginary Transition State graph* (ITS) is the labeled graph

$$
\Gamma = (V,\, E,\, \mathbf{a},\, \mathbf{b}_{\text{ITS}})
$$

where $V = V_R = V_P$ (shared atom set via $\mu$), $E = E_R \cup E_P$ (union of reactant and product bonds), and the ITS edge attribute is:

$$
\mathbf{b}_{\text{ITS}}(e) = (b_r(e),\, b_p(e))
\quad b_r, b_p \in \{0, 1, 1.5, 2, 3\}
$$

where $b_r(e)$ is the bond order of edge $e$ in the reactants (0 if absent) and $b_p(e)$ is the bond order in the products (0 if absent).

**Definition (Reaction Center).**  
The *reaction center* $\mathrm{RC}(\Gamma)$ is the subgraph of the ITS induced by edges where $b_r \neq b_p$:

$$
E_{\mathrm{RC}} = \{e \in E(\Gamma) \mid b_r(e) \neq b_p(e)\}
$$

These are the bonds that are **broken** ($b_r > 0, b_p = 0$), **formed** ($b_r = 0, b_p > 0$), or **changed** ($b_r \neq b_p$, both $> 0$).

**Definition (ΔBE Entry).**  
For each atom pair $(i, j)$ with $i \neq j$: $\Delta\mathrm{BE}[i,j] = b_p(i,j) - b_r(i,j)$.  
Positive values indicate bond formation; negative values indicate bond breaking.  
The diagonal entries represent changes in free-electron count (valence electrons not in bonds).



### 2.1. ITS graph — edge coloring by bond-change type


For drawing, the same ITS can be placed on different molecular coordinate systems.  
`coordinate="reactant"` projects the ITS onto the reactant-side bonds, while `coordinate="product"` projects it onto the product-side bonds. The edge colors and labels still show the full `(b_r, b_p)` change, but the picture is easier to read because the atom positions come from a chemically valid molecule instead of a generic graph layout.
`its_to_side_graph(...)` and `its_to_side_mol(...)` expose the reactant/product projections for later tasks. `its_coordinate_layout(...)` returns these coordinates directly, which is useful when you want to draw a reaction-center subgraph on the same frame as the full ITS.

Each edge in the ITS carries a pair $(b_r, b_p)$ encoding the bond order in reactants and products.  
We color edges by change type:

| $(b_r, b_p)$ | Interpretation | Color |
|---|---|---|
| $b_r = b_p > 0$ | **Preserved** bond | gray |
| $b_r > 0,\; b_p = 0$ | **Broken** bond | red |
| $b_r = 0,\; b_p > 0$ | **Formed** bond | green |
| $b_r \neq b_p$, both $> 0$ | **Changed** (e.g. single→double) | orange |



**Q3 — Reaction center**

Now develop function `get_reaction_center` to extract reacton center

---

<details> <summary><b>Solution:</b></summary>

```python
import networkx as nx


def get_reaction_center(
    its: nx.Graph,
    *,
    node_view: bool = False,
):
    rc_edges = []
    rc_nodes = set()

    for u, v, d in its.edges(data=True):
        br, bp = d.get("order", (0.0, 0.0))
        if br != bp:
            rc_edges.append((u, v))
            rc_nodes.update([u, v])

    if node_view:
        return rc_nodes

    rc = nx.Graph()
    rc.add_nodes_from((n, its.nodes[n]) for n in rc_nodes)
    rc.add_edges_from(
        (u, v, its.edges[u, v]) for u, v in rc_edges
    )

    return rc

```
</details>


Now we combine to 1 function `rsmi_to_its`


### 2.2. ITS reaction center and ΔBE matrix

The left panel shows the full ITS: **red edges** are broken bonds, **green edges** are formed bonds, **black edges** are preserved. The right panel isolates the **reaction center** — the subgraph where $b_r \ne b_p$.

The ΔBE matrix encodes the same information numerically: $\Delta BE_{ij} = b^P_{ij} - b^R_{ij}$. Red cells indicate bond cleavage; blue cells indicate bond formation. The pattern of signed changes uniquely fingerprints the reaction type.


### 2.3. ΔBE matrix heatmap

The Bond–Electron (BE) matrix introduced in **S01** can be computed for both the reactant
and product graphs. Their difference $\Delta\mathrm{BE} = \mathrm{BE}_P - \mathrm{BE}_R$
captures exactly the same information as the ITS edge labels:

- **Positive** off-diagonal: bond formed
- **Negative** off-diagonal: bond broken
- **Zero**: atom-pair unchanged


## 3. ITS equivalence

Now try with 3 reactions below


Even for experienced chemists, it is non-trivial to determine whether three atom mappings are equivalent.


While atom-map equivalence is hard to judge directly, the corresponding ITS graphs are equivalent if and only if they are isomorphic, a concept we already encountered in **S02** [\[6\]](#6.-References).


**MCS vs RXNMapper — ITS isomorphism comparison**

Two atom mappings are **equivalent** if and only if their ITS graphs are isomorphic. The table below runs both methods on three test reactions and uses ITS isomorphism to check whether they agree. Disagreement flags a case where one method makes a chemically different assignment.


**Q4 — Atom-map equivalence**

Consider the dataset `data/reaction_maps.csv`
Each reaction is associated with three atom-mapped SMILES, generated by different atom-mapping methods (`map1`, `map2`, `map3`).


Tasks:
1. Write a function that determines whether three atom-mapped reactions are equivalent, by converting each to an ITS graph and checking graph isomorphism.
2. Apply this function to all reactions in the dataset.
3. Record, for each reaction, whether the three mappings are equivalent.

---

<details> <summary><b>Solution:</b></summary>

```python
import pandas as pd
def is_maps_eq(
    maps,
    *,
    node_attrs=("element", "aromatic", "charge", "hcount"),
    edge_attrs=("order",),
):
    """
    Check whether multiple mapped reactions are equivalent
    by comparing their ITS graphs up to isomorphism.
    """
    its_list = [rsmi_to_its(m) for m in maps]

    ref = its_list[0]
    for its in its_list[1:]:
        if not its_isomorphic(
            ref,
            its,
            node_attrs=node_attrs,
            edge_attrs=edge_attrs,
        ):
            return False
    return True

data = pd.read_csv('./data/reaction_maps.csv').to_dict('records')

for entry in data:
    maps = [entry["map1"], entry["map2"], entry["map3"]]
    entry["is_eq"] = is_maps_eq(maps)
data = pd.DataFrame(data)
data.head()
```

</details>



<a id="4-discussion"></a>

## 4. Discussion

- **MCS** is a principled alignment baseline, but struggles with multi-component reactions, molecular symmetry, and unbalanced equations.
- **RXNMapper** offers a strong, practical solution by leveraging learned attention patterns and providing confidence estimates.
- **ITS graphs** compactly encode bond changes, abstracting away atom-map labeling details.
- **ITS isomorphism** provides a clean, map-number-invariant criterion for comparing atom mappings across methods.


<a id="5-quiz"></a>

## 5. Quiz

1. Why can maximum common substructure mapping fail for reactions with multiple components, symmetry, or large structural rearrangements?
2. What information does RXNMapper use that is not explicitly enforced by a graph-only MCS alignment?
3. In one or two sentences, explain what an ITS graph represents chemically.
4. Why can two atom-mapped reaction SMILES look different but still describe equivalent chemistry, and how does ITS isomorphism help compare them?



## 6. References

1. Phan, T.-L. *et al.* SynKit: A graph-based framework for rule-based reaction modeling. *Journal of Chemical Information and Modeling* (2025). https://doi.org/10.1021/acs.jcim.5c02123
2. Schwaller, P. *et al.* Extraction of organic chemistry grammar from unsupervised learning of chemical reactions. *Science Advances* **7**, eabe4166 (2021). https://doi.org/10.1126/sciadv.abe4166
3. Dugundji, J.; Ugi, I. An algebraic model of constitutional chemistry as a basis for chemical computer programs. *Topics in Current Chemistry* **39**, 19-64 (1973).
4. Schwaller, P. *et al.* Molecular Transformer: A Model for Uncertainty-Calibrated Chemical Reaction Prediction. *ACS Central Science* (2019). https://doi.org/10.1021/acscentsci.9b00576
5. Phan, T.-L. *et al.* SynTemp: Efficient Extraction of Graph-Based Reaction Rules from Large-Scale Reaction Databases. *Journal of Chemical Information and Modeling* (2025). https://doi.org/10.1021/acs.jcim.4c01795
6. Bonchev, D.; Rouvray, D. H., eds. *Chemical Graph Theory: Introduction and Fundamentals*. Abacus Press (1991).
7. RDKit documentation. https://www.rdkit.org/docs/
8. NetworkX documentation. https://networkx.org/documentation/stable/
9. Nugmanov, R. I. *et al.* CGRtools: Python Library for Molecule, Reaction, and Condensed Graph of Reaction Processing. *Journal of Chemical Information and Modeling* **59**, 2516-2521 (2019). https://doi.org/10.1021/acs.jcim.9b00102
