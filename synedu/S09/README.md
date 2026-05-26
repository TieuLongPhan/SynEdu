# S09: Context Graph Expansion

This talktorial introduces **hierarchical context expansion** as a principled way to
control the generality–specificity trade-off in reaction-center templates.
By expanding the reaction center with progressively larger chemical neighbourhoods
(radii *r* = 0 … 5), the template library becomes either more general (small *r*,
high coverage) or more specific (large *r*, high precision).
We first build the core BFS expansion by hand, then scale it with SynKit's
`HierContext`, and evaluate multi-radius libraries on a real dataset.


## Aim of this talktorial

This talktorial (**S09**) studies **context graph expansion**: how much chemical neighbourhood should be included around a reaction center when building reaction templates.

A minimal reaction center is highly general and can match many substrates, but it may also be too permissive. A larger context is more specific and often more precise, but it may fail to generalize. We explore this trade-off by expanding the reaction-center context with increasing radius and evaluating how template behavior changes.

Concretely, we focus on:

1. **The r-hop context around a reaction center**  
   Defining $N_r(C)$ on an ITS graph and interpreting the resulting context subgraph.

2. **Manual BFS expansion**  
   Implementing context expansion from scratch so the algorithm is transparent.

3. **SynKit hierarchical context extraction**  
   Scaling the same idea with `HierContext` for multi-radius template libraries.

4. **Generality-specificity evaluation**  
   Measuring how recall, precision, and F1 change as the radius increases.

---

## Learning outcomes

After completing this talktorial, you will be able to:

- define the *r*-hop context neighbourhood $N_r(C)$ on an ITS graph,
- explain how context controls the generality-specificity trade-off of a DPO rule,
- implement BFS-based context expansion and verify it against SynKit's `HierContext`,
- extract and visualize multi-radius template libraries for a rule family,
- evaluate forward-prediction recall, precision, and F1 at each radius, and
- choose a practical context radius for a dataset using evidence rather than guesswork.

---

## Outline

<ul class="synedu-outline">
  <li><a href="#0-setup--data">0. Setup &amp; data</a></li>
  <li><a href="#1-what-is-context">1. What is context?</a></li>
  <ul>
    <li><a href="#11-the-reaction-center--recap">1.1 The reaction center — recap</a></li>
    <li><a href="#12-defining-the-r-hop-context">1.2 Defining the r-hop context</a></li>
    <li><a href="#13-implement-expand_context-from-scratch">1.3 Implement expand_context from scratch</a></li>
    <li><a href="#14-visualise-r--0-1-2-on-one-its">1.4 Visualise r = 0, 1, 2 on one ITS</a></li>
    <li><a href="#15-the-context-subgraph-as-the-dpo-k-graph">1.5 The context subgraph as the DPO K-graph</a></li>
  </ul>
  <li><a href="#2-rule-family-deep-dive--radii-05">2. Rule family deep-dive — radii 0-5</a></li>
  <ul>
    <li><a href="#21-pick-the-working-family">2.1 Pick the working family</a></li>
    <li><a href="#22-extract-multi-radius-templates-with-hiercontext">2.2 Extract multi-radius templates with HierContext</a></li>
    <li><a href="#23-template-count-vs-radius">2.3 Template count vs. radius</a></li>
    <li><a href="#24-visual-tour--k-graphs-at-r--0-1-2">2.4 Visual tour — K-graphs at r = 0, 1, 2</a></li>
    <li><a href="#25-multi-radius-evaluation-loop">2.5 Multi-radius evaluation loop</a></li>
    <li><a href="#26-trade-off-curve">2.6 Trade-off curve</a></li>
  </ul>
  <li><a href="#3-generalising--family-dependent-optima">3. Generalising — family-dependent optima</a></li>
  <li><a href="#4-choosing-a-radius-in-practice">4. Choosing a radius in practice</a></li>
  <li><a href="#5-quiz">5. Quiz</a></li>
  <li><a href="#6-discussion">6. Discussion</a></li>
</ul>


## 0. Setup & data


**Split the dataset**

We use an 80 / 20 train / test split with the same seed as S08 so that results
are comparable. The training split is used to build template libraries at every
radius; the test split is reserved for evaluation.


**Build the rule library — ITS, WL hash, and clustering**

We reuse the exact pipeline from S07/S08:

1. Convert each training reaction to its full ITS graph and to its
   reaction-center (RC) graph.
2. Compute a WL hash of each RC for fast grouping.
3. Cluster with `GraphCluster` to identify rule families.

This is shown quickly here; refer to S07 for a step-by-step explanation.


**Standardise the test set**

We remove atom mapping from test reactions so that the predictor never sees the
ground-truth correspondence. Standardisation also normalises ring perception and
aromaticity, which is required for fair chemical comparison.


## 1. What is context?

The *reaction center* captures the minimum subgraph that changes in a reaction.
**Context expansion** augments that minimal representation by pulling in the
chemical neighbourhood around it.  Before using SynKit's ready-made tool, we
build the expansion by hand so the mechanism is fully transparent.


### 1.1 The reaction center — recap

From S04: an edge (u, v) in the ITS graph belongs to the **reaction center**
E_rc if its bond orders on the reactant and product sides differ:

$$E_{\text{rc}} = \{(u,v) \in E \mid b_r(u,v) \neq b_p(u,v)\}$$

The **reaction-center nodes** $V_{rc}$ are all endpoints of edges in $E_{rc}$.
At radius *r* = 0, $V_{rc}$ is the starting seed for context expansion.


### 1.2 Defining the r-hop context

**Definition.** Let C ⊆ V be the reaction-center atoms of an ITS graph G.
The *r-hop context* is:

$$N_r(C) = \{ v \in V \mid \text{dist}(v,\, C) \leq r \}$$

where dist(v, C) = min_{c ∈ C} d_G(v, c) is the shortest-path distance in G.

The **context subgraph** K_r is the subgraph of G induced by N_r(C).

| r | What is included |
|---|---|
| 0 | Only the atoms that directly break or form bonds |
| 1 | Those atoms + all their immediate ITS neighbours |
| 2 | r=1 context + one more shell of neighbours |
| k | All atoms within k hops of any reaction-center atom |

At r = 0, K_r is the smallest possible rule (maximum generality).
Each additional hop adds chemical context, increasing specificity.


### 1.3 Implement `expand_context` from scratch

**Q1.** Before reading the solution below, try to implement `expand_context(its, r)`
yourself: given an ITS graph and a radius *r*, return the induced subgraph of all nodes
within *r* hops of the reaction center.

*Hint*: start with a frontier set equal to the reaction-center nodes.
Each iteration, extend the frontier by one shell of neighbours and merge into
the accumulated set.


### 1.4 Visualise r = 0, 1, 2 on one ITS


### 1.5 The context subgraph as the DPO K-graph

In S05 we defined a **DPO span** L ← K → R where K is the *interface*
subgraph — the part of the rule that is not rewritten.

The context subgraph K_r built above is exactly this interface graph:

- **L** = the reactant-side view of K_r (with b_r bond orders)
- **R** = the product-side view of K_r (with b_p bond orders)
- **K** = the shared atoms and unchanged bonds

A larger K_r carries more chemical context → the rule pattern-matches fewer
(but more precisely chosen) substrates → higher precision, lower coverage.


## 2. Rule family deep-dive — radii 0–5

We now focus on a **single rule family** — the largest cluster — and study
how the template library and prediction quality change as we vary *r* from 0 to 5.


### 2.1 Pick the working family


Let us look at a few representative reactions from this family to confirm
they share the same reaction-center topology.


### 2.2 Extract multi-radius templates with `HierContext`

`HierContext` does exactly what our hand-built `expand_context` does,
but applied to every training reaction simultaneously and with automatic
deduplication at each radius level.

```
HierContext(max_radius=5).fit(family_train, its_key="ITS")
    → (demo,  templates)
       demo        : updated training records with hierarchical cluster indices
       templates   : list of length max_radius+1
                     templates[r] = list of unique K-graphs at radius r
                     each entry: {"K": nx.Graph, "class": int, ...}
```


### 2.3 Template count vs. radius


### 2.4 Visual tour — K-graphs at r = 0, 1, 2

We draw one representative K-graph per radius using `visualize_its`
(which handles ITS-style edge attributes automatically).


### 2.5 Multi-radius evaluation loop

We now apply each radius's template library to the test set and measure:

| Metric | What it captures |
|---|---|
| **Coverage** (recall) | Fraction of test reactions where the ground truth appears in the candidate set |
| **Recognition** (precision) | Fraction of *unique* generated candidates that match the ground truth |
| **F1** | Harmonic mean of Coverage and Recognition |

We use the same `SynReactor` and `_compute_metrics` tools introduced in S08.


### 2.6 Trade-off curve

The two-panel figure below is the central diagnostic of S09:

- **Left panel** — Coverage vs. Recognition scatter, one point per radius.
  Moving right and up is better; the best operating point is the point
  closest to the top-right corner.
- **Right panel** — Template count and wall-clock time vs. radius.
  These capture the *cost* of higher specificity.


## 3. Generalising — family-dependent optima

The optimal radius found for our working family may not apply universally.
Different rule families have different symmetry profiles, ring environments,
and neighbour diversity, so the coverage–precision trade-off curve will differ.

The code below repeats the evaluation across the **top 3 families** at radii
r = 0, 1, 2 and assembles the results into a **heatmap** (family × radius → F1).
Running it confirms that the best radius is family-dependent.


## 4. Choosing a radius in practice

The trade-off curves and heatmap together inform a practical decision workflow.
The schematic below summarises the main decision points.


**Practical guidelines (summary)**

| Situation | Action |
|---|---|
| Need to identify core transformation types | Start at r = 0 |
| Symmetry-induced ambiguous matches | Increase to r = 1 |
| Retrosynthesis planning (high precision) | r = 2 or 3 |
| Coverage drops sharply above some r | Normalise charges / tautomers upstream |
| Latency is a hard constraint | r ≤ 2 |
| F1 plateaus | You are in the saturation regime — larger r gives no gain |

**Never tune the radius on the test set.** Use a validation split for radius
selection and reserve the test set for final reporting only.


**Q2.** At r = 0, two reactions that differ only in a remote substituent
(e.g. –F vs. –Cl on the meta position of an aromatic ring) are treated as
**identical** templates because the reaction center is the same.

Give a concrete example where this causes an incorrect forward prediction.
What is the *minimum* r needed to distinguish the two reactions, assuming the
substituent is located exactly 3 bonds from the nearest reaction-center atom?


## 5. Quiz

### Q1 — Template count saturation

If you increase r from 2 to 3 for a given family, do you always get *more*
unique templates? Explain why or why not, referring to the deduplication step
performed by `HierContext`.

*Hint*: think about what happens when every training reaction in the family
shares exactly the same 3-hop neighbourhood.

---

### Q2 — When does r = 0 fail? (see Section 4)

### Q3 — Ring-aware context (implementation)

Implement `ring_aware_context(its, r)`: identical to `expand_context` but,
once any atom of a ring is included in the context, **all** atoms of that ring
are also included (even if they are farther than r hops away).

```python
def ring_aware_context(its: nx.Graph, r: int) -> nx.Graph:
    ...
```

*Hint*: use `nx.cycle_basis(its)` to enumerate rings.
After building the r-hop context with ordinary BFS, iterate over all rings:
if any ring member is already in the context, add the whole ring.


## 6. Discussion

### Key takeaways

**The generality–specificity trade-off is fundamental.**
Any pattern-matching system over structured data faces it. The graph-rewriting
formalism makes it explicit, controllable, and measurable via standard
information-retrieval metrics.

**The optimal radius is family-dependent.**
Families whose reactions share a highly symmetric neighbourhood saturate
at r = 1.  Families involving ring-fused or stereocentre-bearing contexts
may need r = 3 or higher.

**`HierContext` is just principled BFS at scale.**
The hand-built `expand_context` in §1 and the SynKit `HierContext` are
algorithmically equivalent; the latter adds deduplication and multi-level
indexing across an entire training corpus.

**What this series does not cover** — and where to go next:

| Extension | Notes |
|---|---|
| Stereocentre-aware templates | Requires encoding CIP descriptors in node attributes |
| Tautomer normalisation | Pre-process with RDKit's MolStandardize before extracting ITS |
| Learning-based candidate ranking | Train a graph neural network to rank candidates from large *r* |
| End-to-end retrosynthesis | Use templates as moves in an MCTS search tree |
| Reaction-centre definition | Replace bond-order change with charge, radical, or isotope change |

### The full pipeline in one sentence

SMILES → ITS → reaction center → WL hash → cluster by family →
`HierContext` per family → multi-radius template library → `SynReactor` →
`_compute_metrics` → choose *r*.

This is the conceptual endpoint of the SynEdu series.
