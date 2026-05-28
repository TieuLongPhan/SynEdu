# S08: One-Step Reaction Prediction

This talktorial introduces a graph-based, template-driven workflow for single-step reaction prediction and single-step retrosynthesis [\[1\]](#6.-References), [\[2\]](#6.-References). Reaction-center templates are extracted from mapped reactions, applied to unseen substrates, and evaluated after reaction standardization [\[3\]](#6.-References), [\[4\]](#6.-References), [\[5\]](#6.-References).



## Aim of this talktorial

1. Build a compact **forward prediction** pipeline from mapped reactions and reaction-center templates.
2. Evaluate generated candidates with standardized **recall** and **enrichment** metrics.
3. Reuse the same rule library for **backward prediction** and discuss why retrosynthesis is more ambiguous.

---

## Learning outcomes

After completing this talktorial, you will be able to:

- derive a compact template library from mapped reactions using **reaction centers**,
- perform **forward one-step prediction** by applying all templates to a substrate,
- standardize generated reactions before comparing them to a reference,
- evaluate predictions using **recall** and **enrichment**, and
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
  <li><a href="#6.-References">6. References</a></li>
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

We first divide the reactions into a **training split** and a **test split** [\[6\]](#6.-References). The
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

Each training reaction is converted into a **reaction-center graph** [\[7\]](#6.-References), that is, a local
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

For prediction, we standardize the test reactions and keep clean reactant and product
strings. This makes the inference task more realistic: the predictor receives the query
molecule and must generate the correct transformation from the template library.

Standardization is also essential for fair evaluation. The same chemistry may admit
multiple syntactically different reaction strings because of atom ordering, component
ordering, or other representation details. Standardization reduces these superficial
differences so that comparisons focus on the recovered reaction transformation rather
than formatting artifacts.

At this point, the notebook has two clearly separated objects:

- a **template library** extracted from mapped training reactions, and
- an **evaluation set** expressed in a standardized, prediction-oriented form.



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
- Why can one template produce multiple outcomes?
- When are two generated reactions genuinely different, and when are they merely
  different representations of the same chemistry?

Working through one example first makes the later dataset-level evaluation much easier to
interpret.



**Q1 — Inspect reaction uniqueness**

The reactor can return several raw candidate strings for the same substrate-template
application. Some of these candidates are genuinely different transformations, while
others are different representations of the same transformation.

Your task is to:

1. **standardize** each generated reaction,
2. group candidates that collapse to the same standardized representation, and
3. check how many unique standardized reactions remain.

If the number of unique reactions decreases, explain **what type of redundancy** has been
collapsed.

---

<details>
<summary><b>Why this question matters</b></summary>

This exercise teaches an important lesson in reaction informatics: **string diversity is
not the same as chemical diversity**. A candidate list can look large simply because the
same transformation is represented in several equivalent ways.

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

Compare the raw candidates inside each group. What changed in the string, and what stayed
chemically the same?

</details>

---

<details>
<summary><b>Solution</b></summary>

Standardization reduces the candidate list to the set of unique reaction transformations.
The drop occurs because part of the apparent diversity comes from representation-level
redundancy: several raw candidates encode the same standardized chemistry.

This is exactly why candidate post-processing is necessary before evaluation. Otherwise,
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

The reactor returns raw reaction strings, and raw strings are not a reliable evaluation
surface. Equivalent chemistry can appear with different ordering, formatting, or repeated
candidate generation. Before computing metrics, we therefore convert each reaction to a
standardized representation and remove duplicates.

In this talktorial we use one evaluation convention throughout:

1. generate candidate reactions from the template library,
2. standardize each candidate reaction,
3. deduplicate the standardized candidates, and
4. compare them with the standardized reference reaction.

This keeps the metric story focused on the question students usually want to answer
first: did the template system recover the correct reaction transformation?



**Definition (Template-Based One-Step Prediction).**  
Given a template library $\mathcal{T} = \{p_1, \ldots, p_m\}$ of DPO rules and a query molecule $q$, the *one-step forward prediction* generates the candidate set:

$$
\hat{\mathcal{P}}(q) = \bigcup_{p_i \in \mathcal{T}} \{H \mid q \Rightarrow_{p_i} H \text{ via some valid match}\}
$$

**Definition (Standardized Match).**  
Let $S(\cdot)$ denote the reaction standardization procedure. A generated candidate $c$
is counted as a hit for reference reaction $r$ when:

$$
S(c) = S(r).
$$

Standardization makes the comparison less sensitive to representation details and more
focused on the recovered reaction transformation.

**Definition (Recall@K).**  
Let $q^*$ be the ground-truth reaction for query $q$. *Recall at $K$* is:

$$
\mathrm{Recall@}K(q) = \mathbf{1}\!\left[S(q^*) \in S(\hat{\mathcal{P}}_K(q))\right]
$$

where $\hat{\mathcal{P}}_K(q)$ is the top-$K$ candidates ranked by, e.g., template frequency. The *dataset Recall@K* is the mean over all test queries.

**Definition (Enrichment@K).**  
Recall only tells us whether the answer appears. It does **not** tell us how many
candidates a student or algorithm must inspect. We therefore also compute:

$$
\mathrm{Enrichment@}K(q)
= \frac{\#\text{ standardized ground-truth hits in }\hat{\mathcal{P}}_K(q)}
       {\#\text{ unique standardized candidates in }\hat{\mathcal{P}}_K(q)}.
$$

With one ground-truth reaction, this is either $0$ or approximately $1/n$ for a prefix
containing $n$ unique standardized candidates. Higher enrichment means the correct answer
is found in a more focused candidate set.

**Definition (Branching Factor).**  
The *branching factor* of a query $q$ is $|\hat{\mathcal{P}}(q)|$ — the total number of distinct candidates generated. A high branching factor indicates low template specificity; a low branching factor (possibly 0) indicates poor coverage.



### Forward recall@K and enrichment@K trade-off

We now sweep $K$ from 1 to 50 for **forward prediction only**. Recall@K asks whether the ground-truth reaction appears in the top-$K$ candidates. Enrichment@K asks how concentrated that hit is within the inspected candidate prefix.

These curves should be read together: increasing $K$ can improve recall, but it usually lowers enrichment because the search space becomes broader.



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


### 2.1. Standardized evaluation

We now evaluate generated reactions in the same way for every example:

1. **standardize** the generated reactions,
2. **remove duplicate standardized candidates**, and
3. compare them against the standardized ground truth.

This collapses formatting and representation-level variants and focuses the score on the
recovered reaction transformation. For many practical reaction-prediction use cases, this
is the most direct first measure because end users want to know whether the correct
chemistry appears in the generated candidate set.



**Q3 — Evaluate the generated reaction set**

After scaling up forward reaction generation, each test entry may contain a **large
number of candidate reactions**. Because different templates can lead to the **same final
transformation**, evaluation must account for redundancy.

Your task is to:

1. standardize and deduplicate the generated candidates,
2. compute **recall**, and
3. compute **enrichment**.

---

<details>
<summary><b>Why this question matters</b></summary>

This exercise forces you to interpret performance scientifically rather than just reading
a single score. A predictor can recover the correct transformation but still produce a
very broad candidate set. Recall and enrichment separate those two behaviors.

</details>

---

<details>
<summary><b>Hint</b></summary>

Evaluate candidates only after standardization. Then compare the unique standardized
candidate set with the standardized ground truth.

</details>

---

<details>
<summary><b>Solution</b></summary>

```python
print(f"Recall: {recall(std_fw, std_gt)}")
print(f"Enrichment: {enrichment(std_fw, std_gt):.4f}")
```

**Discussion**

Standardized evaluation asks whether the generated reaction matches the correct net
transformation after routine representation cleanup.

After standardizing the predictions, the number of unique candidates often decreases. If
recall is preserved while enrichment improves, that means the generator was already
proposing the right chemistry, but part of its apparent diversity came from redundant
representations rather than distinct solutions.
</details>



## 3. Backward prediction

In **backward prediction**, we reverse the perspective: instead of mapping reactants to
products, we start from the **product side** and enumerate plausible precursor reactions
by applying the same template library in **inverted mode**.

This is closely related to **single-step retrosynthesis** [\[1\]](#6.-References), [\[2\]](#6.-References), but it is usually more
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

After evaluating all test reactions, the bar chart below aggregates standardized **mean
recall** and **mean enrichment** for forward prediction (blue) and backward /
retrosynthetic prediction (gold). Recall shows whether the reference transformation was
recovered at least once; enrichment shows how concentrated that hit was within the unique
candidate set.



## 4. Discussion

This notebook illustrates a compact but scientifically meaningful **template-based
one-step prediction pipeline**:

1. extract **reaction-center templates** from a mapped training set,
2. apply the templates to unseen substrates or products,
3. standardize and deduplicate the generated candidate set, and
4. evaluate predictions with standardized recall and enrichment.

Several lessons follow from this workflow.

### 4.1. What the pipeline demonstrates

- **Template extraction compresses reaction knowledge.** Rather than replaying full
  training reactions, the system reuses generalized local transformations.
- **Deduplication is essential.** Raw candidate counts can exaggerate diversity and
  distort enrichment if repeated representations are not collapsed.
- **Recall and enrichment answer different questions.** Recall asks whether the correct
  transformation appears; enrichment asks how focused the candidate set is.
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

a. Why can the number of unique predictions decrease after standardization?

Explain how two raw reaction strings can describe the same standardized transformation.

---

b. Why is it useful to report both recall and enrichment?

What can happen when recall improves but enrichment decreases?

---

c. Why is backward prediction usually more ambiguous than forward prediction?

Discuss the role of **multiple valid disconnections**, selectivity, leaving groups, and
reaction conditions.

---

d. Why does deduplication matter before computing enrichment?

What would go wrong if repeated candidates were counted separately?

---

e. In a template-based system, what does a recovered hit actually demonstrate?

Does it show exact memorization of a training reaction, recovery of a reaction class, or
successful application of a generalized local transformation? Explain.

---

f. Suppose two template libraries have the same recall, but one has clearly better enrichment.

What does this imply about the **focus** of the generated search space? Why can
enrichment be important even when recall is unchanged?



## 6. References

1. Segler, M. H. S.; Preuss, M.; Waller, M. P. Planning Chemical Syntheses with Deep Neural Networks and Symbolic AI. *Nature* **555**, 604-610 (2018).
2. Coley, C. W. *et al.* A Robotic Platform for Flow Synthesis of Organic Compounds Informed by AI Planning. *Science* **365**, eaax1566 (2019).
3. Coley, C. W.; Green, W. H.; Jensen, K. F. RDChiral: An RDKit Wrapper for Handling Stereochemistry in Retrosynthetic Template Extraction and Application. *Journal of Chemical Information and Modeling* **59**(6), 2529-2537 (2019).
4. Schwaller, P. *et al.* Extraction of organic chemistry grammar from unsupervised learning of chemical reactions. *Science Advances* **7**, eabe4166 (2021). https://doi.org/10.1126/sciadv.abe4166
5. Daylight Theory Manual. *Reaction SMILES and SMARTS*.
6. Coley, C. W.; Green, W. H.; Jensen, K. F. Machine Learning in Computer-Aided Synthesis Planning. *Accounts of Chemical Research* **51**(5), 1281-1289 (2018).
7. Phan, T.-L. *et al.* SynKit: A graph-based framework for rule-based reaction modeling. *Journal of Chemical Information and Modeling* (2025). https://doi.org/10.1021/acs.jcim.5c02123
8. RDKit documentation. https://www.rdkit.org/docs/
9. Shervashidze, N.; Schweitzer, P.; van Leeuwen, E. J.; Mehlhorn, K.; Borgwardt, K. M. Weisfeiler-Lehman Graph Kernels. *Journal of Machine Learning Research* **12**, 2539-2561 (2011).
10. Schneider, N.; Stiefl, N.; Landrum, G. A. What's What: The (Nearly) Definitive Guide to Reaction Role Assignment. *Journal of Chemical Information and Modeling* **56**(12), 2336-2346 (2016).
