# S09: Context Radius Expansion

This talktorial asks how much local neighbourhood should be kept around a reaction center. We expand reaction-center templates by radius, apply the resulting rule libraries, and use recall, enrichment, and F1 to choose a practical context size [\[1\]](#6.-References), [\[2\]](#6.-References), [\[3\]](#6.-References).



## Aim of this talktorial

1. Define the **r-hop context** around a reaction center on an ITS graph.
2. Implement transparent BFS-based context expansion and build deduplicated template libraries for several radii.
3. Evaluate the recall-enrichment trade-off and choose a practical radius from validation evidence.

---

## Learning outcomes

After completing this talktorial, you will be able to:

- define the *r*-hop context neighbourhood $N_r(C)$ on an ITS graph,
- explain how context controls the generality-specificity trade-off of a DPO rule,
- implement BFS-based context expansion,
- build deduplicated template libraries at several radii,
- evaluate forward prediction with recall@K and enrichment@K, and
- choose a practical context radius using evidence rather than guesswork.

---

## Outline

<ul class="synedu-outline">
  <li><a href="#0-setup--data">0. Setup &amp; data</a></li>
  <li><a href="#1-what-is-context">1. What is context?</a></li>
  <li><a href="#2-template-radius-deep-dive--radii-05">2. Template radius deep-dive — radii 0-5</a></li>
  <li><a href="#3-choosing-a-radius-in-practice">3. Choosing a radius in practice</a></li>
  <li><a href="#4-discussion">4. Discussion</a></li>
  <li><a href="#5-quiz">5. Quiz</a></li>
  <li><a href="#6.-References">6. References</a></li>
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


**Standardize the test set**

We remove atom mapping from test reactions so that the predictor never sees the
ground-truth correspondence. Standardization also normalizes ring perception and
aromaticity, which is required for fair chemical comparison.


## 1. What is context?

The *reaction center* captures the minimum subgraph that changes in a reaction.
**Context expansion** augments that minimal representation by pulling in the
chemical neighbourhood around it.  Here, we
build the expansion by hand so the mechanism is fully transparent.


### 1.1 Recap reaction center

From **S04**: an edge (u, v) in the ITS graph belongs to the **reaction center**
$E_{rc}$ if its bond orders on the reactant and product sides differ:

$$E_{\text{rc}} = \{(u,v) \in E \mid b_r(u,v) \neq b_p(u,v)\}$$

The **reaction-center nodes** $V_{rc}$ are all endpoints of edges in $E_{rc}$ [\[1\]](#6.-References), [\[3\]](#6.-References).
At radius *r* = 0, $V_{rc}$ is the starting seed for context expansion.


### 1.2 Defining the r-hop context

**Definition.** Let  $C \subseteq V$  be the reaction-center atoms of an ITS graph G.
The *r-hop context* is:

$$
N_r(C) = \left\{ v \in V \mid \operatorname{dist}(v, C) \leq r \right\},
$$
where

$$
\operatorname{dist}(v, C) = \min_{c \in C} d_G(v, c)
$$

is the shortest-path distance in G.

The **context subgraph** $K_r$ is the subgraph of G induced by $N_r(C)$.

| r | What is included |
|---|---|
| 0 | Only the atoms that directly break or form bonds |
| 1 | Those atoms + all their immediate ITS neighbours |
| 2 | r=1 context + one more shell of neighbours |
| k | All atoms within k hops of any reaction-center atom |

At r = 0, $K_r$ is the smallest possible rule (maximum generality).
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

The context subgraph $K_r$ built above is exactly this interface graph:

- **L** = the reactant-side view of $K_r$ (with $b_r$ bond orders)
- **R** = the product-side view of $K_r$ (with $b_p$ bond orders)
- **K** = the shared atoms and unchanged bonds

A larger $K_r$ carries more chemical context. The rule therefore pattern-matches fewer
substrates and usually generates a more focused candidate list. In evaluation, this often
appears as a trade-off: recall@K may decrease, while enrichment can improve because fewer
irrelevant candidates need to be inspected.



## 2. Template radius deep-dive — radii 0–5

We now use the context expansion function from Section 1 to build template libraries at
several radii. To keep the lesson simple, we study one representative rule family rather
than comparing many families at once.



### 2.1 Pick one rule family


Let us look at a few representative reactions from this family to confirm
they share the same reaction-center topology.


### 2.2 Build templates at each radius

For each training reaction in the selected family, we extract the r-hop context subgraph [\[4\]](#6.-References)
with `expand_context(its, r)`. We then deduplicate the resulting K-graphs with a WL hash.

The result is a plain list:

```python
_templates[r] = [
    {"K": context_graph, "hash": wl_hash, "source_idx": training_index},
    ...
]
```

This is intentionally simple: no extra hierarchy, just one template set per radius.



### 2.3 Template count vs. radius


We draw one representative K-graph per radius using `visualize_its`
(which handles ITS-style edge attributes automatically).


### 2.4 Evaluate each radius

We now apply each radius's template library to test reactions from the **same reaction-center family**. This keeps the comparison local: each radius is trying to solve the same class of chemistry, rather than being averaged across unrelated rule classes.

| Metric | What it captures |
|---|---|
| **Recall@K** | Fraction of same-family test reactions where the standardized ground truth appears in the first K unique candidates |
| **Enrichment@K** | Mean hit density inside the first K unique candidates; with one ground truth this is `hit / n_candidates` |
| **F1@K** | Harmonic mean of recall@K and enrichment@K, used here as the radius trade-off score |

We use `SynReactor` for generation, then standardize and deduplicate predictions before computing recall@K, enrichment@K, and F1@K. In this section, K is fixed to 50 so the radius comparison stays easy to read.



### 2.5 Trade-off curve

The two-panel figure below is the central diagnostic of S09:

- **Left panel** — Recall@K, enrichment@K, and their F1 score as radius increases.
  This shows whether larger context is helping both recovery and focus.
- **Right panel** — Recall@K versus enrichment@K, one point per radius.
  Moving right means the answer appears more often; moving up means the candidate list is more focused when the answer appears.

The highlighted radius maximizes `tradeoff_score`, defined here as the harmonic mean of recall@K and enrichment@K. This is not a universal objective; it is a compact visual guide for comparing radii inside one reaction-center family.



## 3. Choosing a radius in practice

The single-family curve above is useful for learning, but radius selection should be checked on more than one reaction-center class. Here we run a small validation-style exercise:

1. choose three reaction-center classes that have both training and test reactions;
2. build class-specific template libraries at **r = 0, 1, 2**;
3. apply each class's templates to all test reactions from the same class;
4. compare recall, enrichment, and F1 using **all generated candidates** within each class, then average across classes.

This keeps the lesson simple while avoiding a common mistake: choosing a radius from one large class and assuming it behaves the same everywhere.

For this setup experiment we set `k=None`: every unique generated candidate is inspected. This gives a clearer view of whether a radius loses the true product entirely, without mixing in a rank cutoff.



### 3.1 Trade-off illustration across three classes

The next figure separates two questions. Because this check uses all generated candidates, the enrichment value directly reflects how diluted the answer is in the full generated set:

- the left panel shows whether each radius is stable across classes;
- the right panel shows the macro-average operating point for each radius.

The highlighted macro point is the best practical default for this small check. It is a default, not a law: a production workflow would repeat the same logic on a larger validation split.



**Practical guidelines**

| Situation | Action |
|---|---|
| r = 0 has high all-candidate recall but low enrichment | Add local context and test r = 1 |
| r = 1 improves all-candidate F1 across most classes | Use r = 1 as the first practical default |
| r = 2 improves only one class but hurts the macro mean | Keep r = 1 globally and override per class only when justified |
| Recall drops sharply from r = 1 to r = 2 | The larger template is becoming too specific |
| Runtime is too high | Cap the maximum radius or pre-filter candidate templates |

In practice, choose a radius from a validation set, report the all-candidate macro-average trade-off, and inspect per-class curves before treating the radius as a default.



<a id="4-discussion"></a>

## 4. Discussion

### Key takeaways

**Context radius controls generality.**
Small radii make broad templates. Larger radii add local chemical information and usually make templates more specific.

**Deduplication matters.**
Many training reactions can produce the same context graph at a given radius. Deduplicating those graphs keeps the template library compact.

**Recall@K and enrichment@K should be read together.**
Recall@K asks whether the correct reaction appears. Enrichment@K asks how focused the candidate list is. In this talktorial, `tradeoff_score` is the F1 score between them.

**Evaluate within a comparable class.**
The radius decision is easier to interpret when all test reactions share the same reaction-center family.

**The simple pipeline**

SMILES -> ITS -> reaction center -> r-hop context -> deduplicated templates -> `SynReactor` -> standardized recall@K, enrichment@K, and F1@K -> choose *r*.

This is the conceptual endpoint of the SynEdu series.



<a id="5-quiz"></a>

## 5. Quiz

1. If you increase radius *r* from 2 to 3 for one reaction family, do you always get more unique templates? Explain the role of deduplication.
2. At *r* = 0, why can two reactions with the same reaction center but different remote substituents collapse to the same template?
3. When two radii have similar recall but different enrichment, which radius is more useful for a focused search workflow and why?
4. Why should the radius trade-off be evaluated within comparable reaction-center families rather than across unrelated classes?



## 6. References

1. Phan, T.-L. *et al.* SynKit: A Graph-Based Python Framework for Rule-Based Reaction Modeling and Analysis. *Journal of Chemical Information and Modeling* (2025). https://doi.org/10.1021/acs.jcim.5c02123
2. Shervashidze, N.; Schweitzer, P.; van Leeuwen, E. J.; Mehlhorn, K.; Borgwardt, K. M. Weisfeiler-Lehman Graph Kernels. *Journal of Machine Learning Research* **12**, 2539-2561 (2011).
3. Ehrig, H.; Ehrig, K.; Prange, U.; Taentzer, G. *Fundamentals of Algebraic Graph Transformation*. Springer (2006). https://doi.org/10.1007/3-540-31188-2
4. Coley, C. W.; Green, W. H.; Jensen, K. F. RDChiral: An RDKit Wrapper for Handling Stereochemistry in Retrosynthetic Template Extraction and Application. *Journal of Chemical Information and Modeling* **59**, 2529-2537 (2019). https://doi.org/10.1021/acs.jcim.9b00286
5. RDKit documentation. https://www.rdkit.org/docs/
6. Daylight Theory Manual. *Reaction SMILES and SMARTS*.
7. Schneider, N.; Stiefl, N.; Landrum, G. A. What's What: The (Nearly) Definitive Guide to Reaction Role Assignment. *Journal of Chemical Information and Modeling* **56**, 2336-2346 (2016). https://doi.org/10.1021/acs.jcim.6b00564
