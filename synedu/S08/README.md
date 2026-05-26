# S08: One-Step Reaction Prediction

This talktorial introduces a **graph-based, template-driven workflow** for
**single-step reaction prediction** and **single-step retrosynthesis**.


## Aim of this talktorial

In this talktorial, we construct a compact and interpretable **one-step reaction
prediction pipeline** from **mapped reaction data**. Starting from a curated reaction
set, we extract **reaction-center templates**, apply them to unseen molecules,
and evaluate the generated candidates under two complementary notions of correctness:
the **chemical level** and the **mechanistic level**. We then reverse the same rule
library for **backward prediction**, showing how a forward reaction model can also be
used to propose precursor hypotheses.

This notebook highlights an important principle in
reaction informatics: a reaction can be viewed not only as a string or record in a
database, but as a **local graph transformation**. The template library therefore acts as
a compact representation of reaction knowledge. Instead of memorizing full reactions, the
system stores **generalized local edit patterns** and reapplies them to new substrates.

---

## Learning outcomes

After completing this talktorial, you will be able to:

- derive a compact template library from mapped reactions using **reaction centers**,
- perform **forward one-step prediction** by applying all templates to a substrate,
- evaluate predictions using **recall** and **enrichment** under two equivalence notions,
- distinguish **chemical correctness** from **mechanistic correctness**, and
- reuse the same template library for **backward prediction** while recognizing why
  retrosynthetic generation is usually more ambiguous.

---

## Outline

<ul class="synedu-outline">
  <li><a href="#0-setup--data">0. Setup &amp; data</a></li>
  <li><a href="#1-forward-prediction">1. Forward prediction</a></li>
  <li><a href="#2-evaluation">2. Evaluation</a></li>
  <li><a href="#3-backward-prediction">3. Backward prediction</a></li>
  <li><a href="#4-discussion">4. Discussion</a></li>
  <li><a href="#5-quiz">5. Quiz</a></li>
  <li><a href="#6-references">6. References</a></li>
</ul>


## 0. Setup & data

We work with a subset of **mapped reaction SMILES** and divide the workflow into two
roles:

- the **training split** is used to extract **reaction rules**, and
- the **test split** is converted into an **unmapped, standardized input set** for
  prediction and evaluation.

At inference time, the predictor does **not** replay full training reactions; instead, it reuses the **generalized local
transformation patterns** distilled from the training subset. In that sense, template
extraction functions as a form of **knowledge compression**.

> **Conceptual note**
>
> A mapped reaction contains two kinds of information at once:
> 1. the **overall chemical transformation**, and
> 2. the **atom correspondence** between reactants and products.
>
> The second piece is what allows us to localize the changed region and extract a
> reaction-center rule.


**Split the dataset**

We first divide the reactions into a **training split** and a **test split**. The
training portion is used to build the **template library**, while the test portion is
reserved for evaluating whether those templates can reproduce unseen reactions.

This is a minimal but important experimental design. Without a split, the notebook would
only show that rules can be extracted and replayed on the same examples. With a held-out
test subset, we instead ask a more meaningful question:

> *Can local transformations learned from one set of reactions be transferred to new
> molecular contexts?*

Although this is still a simplified benchmark, it already captures two central issues in
reaction prediction: **template coverage** and **generalization to unseen substrates**.


**Extract a compact template library**

Each training reaction is converted into a **reaction-center graph**, that is, a local
representation of the atoms and bonds involved in the transformation. We then compute a
**Weisfeiler–Lehman (WL) hash** to group equivalent or near-equivalent local patterns and
keep one representative template per class.

This serves two purposes:

1. it removes unnecessary redundancy from the library, and
2. it makes downstream prediction faster and easier to interpret.

Scientifically, this step emphasizes that many recorded reactions differ only in their
global molecular context, while sharing the **same local transformation logic**.
Grouping by reaction-center identity therefore moves the representation from
**instance-level data** toward **rule-level knowledge**.

> **Why this matters**
>
> A raw reaction dataset may contain many examples of the same transformation class.
> Keeping all of them would enlarge the search space without increasing the true diversity
> of available chemistry.


**Standardize the evaluation set**

For prediction, we remove atom mapping from the test reactions and keep only standardized
reactant and product strings. This makes the inference task more realistic: the predictor
must generate the correct chemistry **without** being handed the ground-truth mapping.

Standardization is also essential for fair evaluation. The same chemistry may admit
multiple syntactically different reaction strings because of atom ordering, component
ordering, or mapping variation. Canonicalization reduces these superficial differences so
that comparisons reflect **chemical content** rather than formatting artifacts.

At this point, the notebook has two clearly separated objects:

- a **template library** extracted from mapped training reactions, and
- an **evaluation set** expressed in a standardized but prediction-oriented form.


## 1. Forward prediction

In **forward prediction**, we start from the **reactant side** and ask which products can
be generated by applying a library of reaction-center templates. This is the most direct
graph-rewriting view of a chemical reaction: a substrate graph is matched against a local
rule, and the corresponding edit is applied to produce candidate product graphs.

Pedagogically, it is useful to begin with a **single example** before scaling to the full
dataset. This lets us inspect what the reactor is doing, identify duplicate or
mapping-equivalent predictions, and build intuition for why a seemingly simple one-step
task can still generate a large candidate set.


### 1.1. Inspect one substrate-template application

Before scaling up, it is helpful to examine a **single substrate** in detail. We
visualize one test example, apply a template, and inspect the generated reactions to see
whether the ground-truth transformation appears among the candidates.

This micro-level inspection answers several educational questions:

- What exactly is being matched on the substrate?
- Why can one template produce multiple mapped outcomes?
- When are two generated reactions genuinely different, and when are they merely
  different representations of the same chemistry?

Working through one example first makes the later dataset-level evaluation much easier to
interpret.


**Q1 — Inspect reaction uniqueness**

The reactions generated by the reactor currently include **atom-to-atom mapping (AAM)**.
At this level, several predictions may appear distinct because their atom indices differ,
even when they encode the **same overall chemical transformation**.

Your task is to:

1. **remove atom mapping** from each generated reaction,
2. **standardize** the resulting reaction representation, and
3. check whether there are still **5 unique reactions**.

If the number of unique reactions decreases, explain **what type of redundancy** has been
collapsed.

---

<details>
<summary><b>Why this question matters</b></summary>

This exercise teaches an important lesson in reaction informatics: **string diversity is
not the same as chemical diversity**. A candidate list can look large simply because the
same transformation is represented with different atom maps or ordering conventions.

</details>

---

<details>
<summary><b>Hint</b></summary>

Use the reaction standardization utility:

```python
from collections import defaultdict
from synkit.Chem.Reaction.standardize import Standardize

std = Standardize()

grouped = defaultdict(list)

for smarts in reactor.smarts_list:
    canonical = std.fit(smarts)
    if canonical is not None:
        grouped[canonical].append(smarts)

grouped = dict(grouped)
grouped
```

Think about whether two mapped reactions can differ only in their **atom correspondence**
while still representing the same transformation after mapping is removed.

</details>

---

<details>
<summary><b>Solution</b></summary>

Removing atom maps and standardizing reduces the set to **4 unique reactions**. The drop
occurs because part of the apparent diversity comes from **mapping-level redundancy**:
two mapped variants differ in atom correspondence, but after removing AAM they collapse
to the **same standardized chemical transformation**.

This is exactly why candidate post-processing is necessary before evaluation: otherwise,
the search space would appear more diverse than it really is.
</details>


### Template match count distribution

For each test reaction, we count how many template applications fired
(the **branching factor**). Reactions with zero matches are **coverage failures**;
reactions with many matches face **disambiguation** pressure.


**Q2 — Scale up the forward reaction generation**

So far, we have applied templates to **one example substrate**. The next step is to turn
this into a **dataset-level prediction routine**.

Your task is to:

1. write a function that applies **all templates** to a single test entry,
2. collect all generated reaction SMARTS into `value["fw"]`,
3. run this function for the **entire test set**, and
4. speed up the process using **parallel processing**.

Aim for a solution that is clean, reusable, and suitable for larger benchmark runs.

---

<details>
<summary><b>Why this question matters</b></summary>

A one-step predictor is not just a chemistry object; it is also a **search procedure**.
Once the template library grows, throughput and redundancy become practical concerns.
This exercise therefore connects **chemical reasoning** with **computational scaling**.

</details>

---

<details>
<summary><b>Hint</b></summary>

Encapsulate the logic for one entry first. Once the per-entry workflow is stable, use
`joblib.Parallel` or a similar tool to distribute the workload across CPU cores.

</details>

---

<details>
<summary><b>Solution</b></summary>

```python
from copy import deepcopy
from joblib import Parallel, delayed
import logging


def silence_logging() -> None:
    logging.disable(logging.CRITICAL)


def generate_forward_reactions(entry, templates):
    silence_logging()  # apply inside each worker
    value = deepcopy(entry)
    r = value["r"]
    result = []

    for t in templates:
        reactor = SynReactor(
            substrate=r,
            template=t["RC"],
            explicit_h=True,
            implicit_temp=False,
            automorphism=True,
        )
        result.extend(reactor.smarts_list)

    value["fw"] = result
    return value


silence_logging()

test_processed = Parallel(n_jobs=4, verbose=1)(
    delayed(generate_forward_reactions)(entry, templates) for entry in test
)
```

This pattern separates **domain logic** (template application) from **execution logic**
(parallelization), which makes the notebook easier to scale and easier to maintain.
</details>


## 2. Evaluation

Raw rule application usually produces a **large and redundant candidate set**. Two main
sources of redundancy are common:

1. **exact duplicates**, where the same reaction SMARTS is generated more than once, and
2. **mapping variants**, where different atom maps encode the same chemistry.

To interpret model behavior correctly, we therefore evaluate the predictions under **two
equivalence notions**:

- **with atom mapping**: strict, mechanism-aware comparison;
- **without atom mapping**: transformation-level chemical comparison.

These two views answer different scientific questions. The mapped setting asks whether the
system recovers not only the right product, but also a consistent **atom flow**. The
standardized setting asks whether the system recovers the correct **net chemical change**,
even if the exact mapping representation differs.


**Definition (Template-Based One-Step Prediction).**  
Given a template library $\mathcal{T} = \{p_1, \ldots, p_m\}$ of DPO rules and a query molecule $q$, the *one-step forward prediction* generates the candidate set:

$$
\hat{\mathcal{P}}(q) = \bigcup_{p_i \in \mathcal{T}} \{H \mid q \Rightarrow_{p_i} H \text{ via some valid match}\}
$$

**Definition (Recall@K).**  
Let $q^*$ be the ground-truth product of query $q$. *Recall at $K$* is:

$$
\mathrm{Recall@}K(q) = \mathbf{1}\!\left[q^* \in \hat{\mathcal{P}}_K(q)\right]
$$

where $\hat{\mathcal{P}}_K(q)$ is the top-$K$ candidates ranked by, e.g., template frequency. The *dataset Recall@K* is the mean over all test queries.

**Definition (Mechanistic vs Chemical Evaluation).**  
- *Mechanistic (mapped) evaluation*: two reactions are equal iff their canonical atom-mapped SMARTS strings are identical — i.e., $\Gamma_1 \cong \Gamma_2$.
- *Chemical (unmapped) evaluation*: two reactions are equal iff their product canonical SMILES are identical (atom-map numbers stripped).

Chemical recall is always $\geq$ mechanistic recall since it ignores the atom correspondence.

**Definition (Branching Factor).**  
The *branching factor* of a query $q$ is $|\hat{\mathcal{P}}(q)|$ — the total number of distinct candidates generated. A high branching factor indicates low template specificity; a low branching factor (possibly 0) indicates poor coverage.



### Recall@K curves — forward and backward prediction

We sweep $K$ from 1 to 50 and measure how often the ground-truth reaction appears
in the top-$K$ candidates. A steep initial rise means the correct answer is
consistently ranked near the top.


**Why two evaluation settings?**

Reaction prediction can be judged at different levels of strictness. A model may recover
the right **chemical transformation** even when the **atom mapping** is not identical to
the reference. Conversely, a mapped agreement is stronger evidence that the model has
captured a mechanistically consistent transformation pattern.

Both views are informative, so we keep them separate throughout the analysis.

> **Interpretation guide**
>
> - A hit **without atom mapping** indicates recovery of the correct **reaction outcome**.
> - A hit **with atom mapping** indicates recovery of both the outcome and a compatible
>   **atom correspondence**.
>
> The gap between the two can reveal where the system already knows the chemistry, but
> still struggles with representation-level or mechanism-level precision.


### Forward prediction gallery — ground truth vs candidates

For one example substrate, the template library generates a set of candidate reactions. The ground-truth reaction is highlighted with a green check mark. Comparing candidates to the ground truth directly shows whether the correct chemistry appears in the top predictions and which template was responsible.


A first cleanup step is to remove **duplicate predictions** before computing metrics.
This prevents repeated candidates from artificially inflating the denominator in
enrichment calculations.

For one test example, the reactor may generate **hundreds of candidate reactions**.
However, many of them are not genuinely different transformations. Deduplication lets us
separate **search breadth** from **true reaction diversity**.

From an educational standpoint, this is a good reminder that evaluation metrics are only
as meaningful as the representation they operate on. Poor post-processing can make a
predictor look worse—or occasionally better—than it actually is.


### 2.1. With atom map

At the mapped level, we keep the **atom-to-atom correspondence** and canonicalize each
reaction SMARTS before comparison. This gives a **strict equivalence test**:
a prediction is counted as correct only if it matches both the **overall transformation**
and the **underlying atom flow**.

This is the more demanding evaluation setting. It is especially informative when the goal
is to study **reaction mechanism proxies**, **reaction-center fidelity**, or the quality
of the extracted transformation rule itself.


Now we consider the **enrichment problem**.

Suppose the reactor generated many candidate reactions, but only **one** corresponds to
the ground truth. We then ask:

> *How concentrated is the correct answer within the generated candidate set?*

In this notebook, enrichment is computed as

$$
\text{Enrichment} = \frac{\#\text{ ground-truth hits}}{\#\text{ unique generated reactions}}
$$

Because recall is binary here, enrichment can be interpreted as the **fraction of the
candidate space occupied by correct solutions**. Smaller candidate sets with the same
recall therefore yield better enrichment.

In practice, recall and enrichment should be read together:

- **high recall, low enrichment** means the rule library finds the answer but generates a
  broad and weakly focused search space;
- **high recall, higher enrichment** means the correct chemistry is recovered in a more
  selective and efficient candidate set.


### 2.2. Without atom map

We now move to a less strict but chemically very important evaluation setting. Instead of
requiring the exact mapped reaction, we only ask whether the prediction recovers the
correct **overall transformation**.

To do this, we:

1. **remove atom mapping** from the generated reactions,
2. **standardize** the resulting representations, and
3. compare them against the standardized ground truth.

This collapses mapping-level variants and focuses on **chemical correctness** rather than
full mechanistic identity. For many practical reaction-prediction use cases, this is the
more directly relevant measure because end users often care first about **what product is
formed**, not whether the internal atom labels are identical.


**Q3 — Evaluate the generated reaction set**

After scaling up forward reaction generation, each test entry may contain a **large number
of candidate reactions**. Because different templates can lead to the **same final
transformation**, evaluation must account for redundancy and the level of equivalence used.

Your task is to:

1. compute **recall** and **enrichment**, and
2. compare performance **with atom mapping** versus **without atom mapping**.

---

<details>
<summary><b>Why this question matters</b></summary>

This exercise forces you to interpret performance scientifically rather than just reading a
single score. A predictor can perform differently depending on whether we care about the
**reaction outcome** or the **full mapped transformation**.

</details>

---

<details>
<summary><b>Hint</b></summary>

You can evaluate the same candidate set under two settings:

- **With atom map**: compare mapped reaction SMARTS after canonicalization.
- **Without atom map**: remove mapping, standardize the reactions, then compare the
  standardized transformations.

The first setting tests **mechanistic consistency**; the second tests **chemical
transformation recovery**.

</details>

---

<details>
<summary><b>Solution</b></summary>

```python
print(f"Recall: {recall(std_fw, std_gt)}")
print(f"Enrichment: {enrichment(std_fw, std_gt):.4f}")
```

**Discussion**

Evaluation **without atom mapping** measures **chemical correctness**: it asks whether the
generated reaction matches the correct net transformation, regardless of atom labels.

Evaluation **with atom mapping** is stricter and measures **mechanistic correctness**:
it asks whether the prediction also preserves the correct reactant-to-product atom
correspondence.

After removing atom maps and standardizing the predictions, the number of unique
candidates often decreases. If recall is preserved while enrichment improves, that means
the generator was already proposing the right chemistry, but part of its apparent
diversity came from **representation-level redundancy** rather than distinct solutions.
</details>


## 3. Backward prediction

In **backward prediction**, we reverse the perspective: instead of mapping reactants to
products, we start from the **product side** and enumerate plausible precursor reactions
by applying the same template library in **inverted mode**.

This is closely related to **single-step retrosynthesis**, but it is usually more
ambiguous than forward prediction. A single product can often be disconnected in multiple
valid ways, and the backward direction is strongly influenced by selectivity, leaving
groups, protecting-group strategy, reagent availability, and reaction conditions. As a
result, backward prediction generally exhibits **larger branching** and **lower
specificity** than the forward task.

Educationally, this section is valuable because it shows that the same reaction rule can
support two distinct tasks depending on the direction in which it is applied.


### 3.1. Inspect one backward prediction example

We first invert a single template application to see how the product-side input expands
into candidate precursor reactions.

This example is useful for building retrosynthetic intuition. In the forward direction,
we usually ask, *“What product can this substrate give?”* In the backward direction, we
instead ask, *“What precursor sets could plausibly lead to this product under this local
transformation?”* That inversion greatly enlarges the space of valid answers.


### 3.2. Scale backward prediction and summarize performance

As in the forward direction, we next apply the full template library to the dataset,
post-process the generated reactions, and compute macro-level summary metrics.

When interpreting the backward results, it is important to keep expectations realistic.
Lower recall or lower enrichment does not necessarily mean the templates are poor. It may
also reflect the intrinsic fact that **many retrosynthetic disconnections are chemically
reasonable**, even if only one corresponds to the recorded reference reaction.


### Backward prediction gallery — ground truth vs candidates

For the same example entry, the inverted template library generates candidate **precursors**. Each candidate is a plausible retrosynthetic disconnection; the ground-truth reaction appears among the candidates when the correct bond-breaking template fires.


### Prediction performance summary

After evaluating all test reactions, the bar chart below aggregates **mean recall** for forward prediction (blue) and backward / retrosynthetic prediction (red). A higher bar means more ground-truth reactions were recovered within the candidate set. Recall without atom mapping is typically higher because it only requires the correct product SMILES, not the exact atom correspondence.


Evaluation **without atom mapping** measures **chemical correctness**: it checks whether
the predicted reaction recovers the correct overall transformation, regardless of how atoms
are labeled.

Evaluation **with atom mapping** is stricter and measures **mechanistic correctness**: it
requires not only the right transformation, but also the right **atom correspondence**
between reactants and products.

This distinction is important. A generated reaction may be **chemically correct** while
still failing to reproduce the expected **atom flow**. Therefore, reaction prediction
should be compared at **both levels**:

- the **chemical level**, to assess recovery of the correct transformation, and
- the **mechanistic level**, to assess whether the reaction is represented with a
  consistent underlying atom mapping.

In practice, the standardized setting often performs slightly better than the mapped
setting, especially in backward prediction, because some failures arise from
**mapping-level variation** rather than incorrect chemistry. That gap is informative:
it reveals where the system already captures the right reaction class, but not yet the
right mechanistic representation.


## 4. Discussion

This notebook illustrates a compact but scientifically meaningful **template-based
one-step prediction pipeline**:

1. extract **reaction-center templates** from a mapped training set,
2. apply the templates to unseen substrates or products,
3. deduplicate the generated candidate set, and
4. evaluate predictions at both the **chemical** and **mechanistic** levels.

Several lessons follow from this workflow.

### 4.1. What the pipeline demonstrates

- **Template extraction compresses reaction knowledge.** Rather than replaying full
  training reactions, the system reuses generalized local transformations.
- **Deduplication is essential.** Raw candidate counts can exaggerate diversity and
  distort enrichment if mapping variants are not collapsed.
- **Chemical and mechanistic evaluation answer different questions.** Transformation-level
  recovery tells us whether the chemistry is right; mapped evaluation tells us whether the
  atom flow is also right.
- **Backward prediction is harder.** Reversing a rule library increases ambiguity because
  many precursor sets can map to the same product.

### 4.2. What the pipeline does *not* yet solve

This notebook is intentionally focused on **generation and interpretation**, not on full
reaction-planning performance. In a production system, one would usually also need:

- **template ranking** or candidate scoring,
- **condition awareness**,
- stronger handling of **stereochemistry** and reagent context,
- filtering by **feasibility** or precursor availability, and
- broader benchmarking across reaction classes.

### 4.3. Scientific takeaway

A self-reproduction experiment of this type is best viewed as a **sanity check** for the
graph-rewriting machinery. If performance is weak, the cause may lie in template quality,
mapping inconsistencies, class imbalance, or overly strict matching constraints—not
necessarily in the idea of template-based prediction itself.

For teaching, however, this framework is especially powerful because every step is
interpretable: we can inspect the template, inspect the match, inspect the generated
reaction, and explain why a prediction succeeds or fails.


## 5. Quiz

Answer the following questions using both **chemical intuition** and **graph-based
reasoning**.

---

a. Why can the number of unique predictions decrease after atom mapping is removed?

Explain how two reactions can be distinct at the **mapped** level but identical at the
**standardized chemical** level.

---

b. What is the difference between chemical correctness and mechanistic correctness?

Why is it useful to report **both** evaluation settings instead of only one?

---

c. Why is backward prediction usually more ambiguous than forward prediction?

Discuss the role of **multiple valid disconnections**, selectivity, leaving groups, and
reaction conditions.

---

d. Why does deduplication matter before computing enrichment?

What would go wrong if repeated or mapping-equivalent candidates were counted separately?

---

e. In a template-based system, what does a recovered hit actually demonstrate?

Does it show exact memorization of a training reaction, recovery of a reaction class, or
successful application of a generalized local transformation? Explain.

---

f. Suppose two template libraries have the same recall, but one has clearly better enrichment.

What does this imply about the **focus** of the generated search space? Why can
enrichment be important even when recall is unchanged?

---

g. Why is atom mapping still valuable, even if many downstream users care mainly about the final product?

Relate your answer to **template extraction**, **mechanistic interpretation**, and
diagnosing errors in reaction generation.


## 6. References

1. **Daylight Theory Manual**. *Reaction SMILES and SMARTS*.  
   Foundational reference for reaction string representations and transform logic.

2. **RDKit Documentation**.  
   Practical reference for molecular representations, canonicalization, and reaction
   handling in cheminformatics workflows.

3. Shervashidze, N.; Schweitzer, P.; van Leeuwen, E. J.; Mehlhorn, K.; Borgwardt, K. M.  
   *Weisfeiler–Lehman Graph Kernels*.  
   **Journal of Machine Learning Research** **2011**, *12*, 2539–2561.  
   Relevant background for WL-based hashing and graph comparison.

4. Coley, C. W.; Green, W. H.; Jensen, K. F.  
   *Machine Learning in Computer-Aided Synthesis Planning*.  
   **Accounts of Chemical Research** **2018**, *51* (5), 1281–1289.  
   Broad overview of one-step prediction and retrosynthetic planning.

5. Segler, M. H. S.; Preuss, M.; Waller, M. P.  
   *Planning Chemical Syntheses with Deep Neural Networks and Symbolic AI*.  
   **Nature** **2018**, *555*, 604–610.  
   Widely cited perspective on one-step disconnection logic within synthesis planning.

6. Coley, C. W.; Thomas, D. A.; Lummiss, J. A. M.; Jaworski, J. N.; Breen, C. P.; Schultz,
   V.; Hart, T.; Fishman, J. S.; Rogers, L.; Gao, H.; et al.  
   *A Robotic Platform for Flow Synthesis of Organic Compounds Informed by AI Planning*.  
   **Science** **2019**, *365*, eaax1566.  
   Useful example of how reaction prediction and planning can connect to practical
   synthesis execution.

7. Coley, C. W.; Green, W. H.; Jensen, K. F.  
   *RDChiral: An RDKit Wrapper for Handling Stereochemistry in Retrosynthetic Template
   Extraction and Application*.  
   **Journal of Chemical Information and Modeling** **2019**, *59* (6), 2529–2537.  
   Highly relevant reference for template extraction and application in reaction
   informatics.

8. Schneider, N.; Stiefl, N.; Landrum, G. A.  
   *What’s What: The (Nearly) Definitive Guide to Reaction Role Assignment*.  
   **Journal of Chemical Information and Modeling** **2016**, *56* (12), 2336–2346.  
   Useful background on reaction representation, roles, and practical cheminformatics.
