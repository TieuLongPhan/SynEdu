---
jupytext:
  text_representation:
    extension: .md
    format_name: myst
    format_version: 0.13
    jupytext_version: 1.19.4
kernelspec:
  display_name: synedu
  language: python
  name: python3
---

# S08: One-Step Reaction Prediction

<div class="synedu-lesson-shell not-prose" style="box-sizing:border-box;margin:8px 0 24px;padding:20px;border:1px solid #243b53;border-radius:16px;background:#102a43;color:#fff;font-family:system-ui,-apple-system,BlinkMacSystemFont,'Segoe UI',sans-serif"><div class="synedu-lesson-shell__top" style="display:flex;align-items:flex-start;justify-content:space-between;gap:16px;flex-wrap:wrap"><div><div class="synedu-lesson-shell__eyebrow" style="margin-bottom:5px;color:#5eead4;font-size:11px;font-weight:800;letter-spacing:.12em;text-transform:uppercase">SynEdu learning path</div><div class="synedu-lesson-shell__meta" style="font-size:16px;font-weight:750">Lesson 8 of 9 <span style="color:#9fb3c8;font-weight:500">· Stage 3 · Prediction</span></div></div><span class="synedu-lesson-shell__progress-label" style="color:#bcccdc;font-size:12px;font-weight:700">89% complete</span></div><div class="synedu-lesson-shell__track" style="height:5px;margin:16px 0 18px;overflow:hidden;border-radius:999px;background:#334e68"><span style="display:block;width:89%;height:100%;border-radius:inherit;background:#2dd4bf"></span></div><div class="synedu-notebook-actions" style="display:flex;align-items:center;gap:10px;flex-wrap:wrap" aria-label="Run this lesson"><a class="synedu-launch-badge" href="https://colab.research.google.com/github/TieuLongPhan/SynEdu/blob/main/docs/downloads/S08.ipynb"><img src="https://colab.research.google.com/assets/colab-badge.svg" alt="Open S08 in Colab" style="display:block;height:24px" /></a><a class="synedu-launch-badge" href="https://mybinder.org/v2/gh/TieuLongPhan/SynEdu/main?urlpath=lab/tree/docs/downloads/S08.ipynb"><img src="https://mybinder.org/badge_logo.svg" alt="Launch S08 in Binder" style="display:block;height:24px" /></a><a class="synedu-launch-badge" href="https://github.com/TieuLongPhan/SynEdu/raw/main/docs/downloads/S08.ipynb"><img src="https://img.shields.io/badge/download-.ipynb-2563eb?logo=jupyter&amp;logoColor=white" alt="Download S08 notebook" style="display:block;height:24px" /></a><a class="synedu-launch-badge" href="/docs/installing"><img src="https://img.shields.io/badge/run-locally-334e68?logo=jupyter&amp;logoColor=white" alt="Run S08 locally" style="display:block;height:24px" /></a></div></div>

This talktorial introduces a graph-based, template-driven workflow for single-step reaction prediction and single-step retrosynthesis [@segler2018planning; @coley2019robotic]. Reaction-center templates are extracted from mapped reactions, applied to unseen substrates, and evaluated after reaction standardization [@coley2019rdchiral; @schwaller2021extraction; @daylight_manual].

+++

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

+++

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
> 2. the **atom correspondence** between reactants and products [@schneider2016whats].
>
> The second piece is what allows us to localize the changed region and extract a
> reaction-center rule.

```{code-cell}
import rdkit
import pandas as pd
import networkx as nx
from pathlib import Path
import importlib.metadata as m
from synkit.IO import load_database
import logging


def silence_logging() -> None:
    logging.disable(logging.CRITICAL)


print("RDKit version:", rdkit.__version__)
print("NetworkX version:", nx.__version__)
print("SynKit version:", m.version('synkit'))

DATA_DIR = Path("data")
DATA_PATH = DATA_DIR / "smart.json.gz"
data = pd.DataFrame(
    load_database(DATA_PATH)[:2000]
)  # we use 2000 data points for a lightweight demonstration
display(data.head())
print(data.shape)


def silence_logging() -> None:
    logging.disable(logging.CRITICAL)


silence_logging()
```

**Split the dataset**

We first divide the reactions into a **training split** and a **test split** [@coley2018machine]. The
training portion is used to build the **template library**, while the test portion is
reserved for evaluating whether those templates can reproduce unseen reactions.

This is a minimal but important experimental design. Without a split, the notebook would
only show that rules can be extracted and replayed on the same examples. With a held-out
test subset, we instead ask a more meaningful question:

> *Can local transformations learned from one set of reactions be transferred to new
> molecular contexts?*

Although this is still a simplified benchmark, it already captures two central issues in
reaction prediction: **template coverage** and **generalization to unseen substrates**.

```{code-cell}
# Split the dataset into template-extraction and evaluation subsets
from sklearn.model_selection import train_test_split

train, test = train_test_split(data, test_size=0.2, random_state=42)
train.head()
```

**Extract a compact template library**

Each training reaction is converted into a **reaction-center graph** [@phan2025synkit], that is, a local
representation of the atoms and bonds involved in the transformation. We then compute a
**Weisfeiler–Lehman (WL) hash** [@shervashidze2011weisfeiler] to group equivalent or near-equivalent local patterns and
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

```{code-cell}
# Extract reaction-center templates from the training split
from synkit.IO import rsmi_to_its
from synkit.Graph.Feature.wl_hash import WLHash
from synkit.Graph.Matcher.graph_cluster import GraphCluster

train = train.to_dict('records')  # convert the dataframe to a list of records
wl = WLHash()
for value in train:
    value['RC'] = rsmi_to_its(
        value['smart'], core=True
    )  # convert each reaction into its reaction-center graph
    value['wl'] = wl.weisfeiler_lehman_graph_hash(
        value['RC']
    )  # compute a WL hash for fast grouping


cls = GraphCluster()
cluster = cls.fit(train, rule_key='RC', attribute_key='wl')
```

```{code-cell}
seen = set()
templates = []

for r in cluster:
    key = (r["class"], r["wl"])
    if key not in seen:
        seen.add(key)
        templates.append(r)
len(templates)
```

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

```{code-cell}
from synkit.Chem.Reaction.standardize import Standardize

# Standardize the test reactions and remove atom mapping for prediction
test = test.to_dict('records')
for value in test:
    rsmi = Standardize().fit(value['smart'], remove_aam=True)
    r, p = rsmi.split('>>')
    value['r'] = r
    value['p'] = p
    value['rsmi'] = rsmi
```

```{code-cell}
# Fix a small, deterministic evaluation subset for template application
# Using the first EVAL_SIZE entries of test avoids slow loops over the full split
EVAL_SIZE = 20
eval_set = test[:EVAL_SIZE]
print(f"Evaluation subset: {EVAL_SIZE} / {len(test)} test reactions")
```

## 1. Forward prediction

In **forward prediction**, we start from the **reactant side** and ask which products can
be generated by applying a library of reaction-center templates. This is the most direct
graph-rewriting view of a chemical reaction: a substrate graph is matched against a local
rule, and the corresponding edit is applied to produce candidate product graphs.

Pedagogically, it is useful to begin with a **single example** before scaling to the full
dataset. This lets us inspect what the reactor is doing, identify duplicate or
mapping-equivalent predictions, and build intuition for why a seemingly simple one-step
task can still generate a large candidate set.

+++

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

```{code-cell}
# Inspect one example substrate before rule application
from synkit.IO import smiles_to_graph
from rdkit import Chem
from synedu.Utils.vis import draw_molecular_graph

EXAMPLE_INDEX = 2
TEMPLATE_INDEX = (
    1  # this template fires on the selected substrate; templates[11] gives no match
)

r = eval_set[EXAMPLE_INDEX]["r"]
G = smiles_to_graph(r)
t = templates[TEMPLATE_INDEX]["RC"]

print(f"Example reaction index: {EXAMPLE_INDEX}")
print(f"Selected template index: {TEMPLATE_INDEX}")
print(
    f"Template reaction-center size: {t.number_of_nodes()} atoms, {t.number_of_edges()} changing bonds"
)
draw_molecular_graph(G)
```

```{code-cell}
from synedu.Utils.rxn_vis import visualize_reaction
from IPython.display import SVG

ground_truth = eval_set[EXAMPLE_INDEX]["smart"]
svg = visualize_reaction(
    ground_truth,
    svg=True,
    legend="Example",
)

display(SVG(svg))
```

```{code-cell}
# Apply one reaction-center template to the example substrate
from synkit.Synthesis.Reactor.syn_reactor import SynReactor

reactor = SynReactor(
    substrate=r,
    template=t,
    explicit_h=True,  # required here: the template transfers an explicit H atom
    implicit_temp=False,
    automorphism=True,
)

fw = reactor.smarts_list
print(f"Generated {len(fw)} raw candidate reaction(s).")
fw
```

```{code-cell}
# Check whether the standardized ground-truth reaction is recovered
from synkit.Chem.Reaction.standardize import Standardize

std = Standardize()
std_ground_truth = std.fit(ground_truth)

for value in reactor.smarts_list:
    if std.fit(value) == std_ground_truth:
        print(value)
```

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

```{code-cell}
# Apply all templates to every entry in the test set
for value in eval_set:
    r = value['r']
    result = []
    for t in templates:
        reactor = SynReactor(
            substrate=r,
            template=t['RC'],
            explicit_h=True,
            implicit_temp=False,
            automorphism=True,
        )

        sols = reactor.smarts_list
        result.extend(sols)
    value['fw'] = result

# Build forward_results in recall@K format
forward_results = [
    {"ground_truth": v["smart"], "predictions": v["fw"]} for v in eval_set
]
print(
    f"Forward: {len(forward_results)} entries, "
    f"{sum(len(r['predictions']) for r in forward_results)} total candidates"
)
```

### Template match count distribution

For each test reaction, we count the number of raw candidate reactions generated across
all templates (a proxy for the **branching factor**, formally defined below as the number
of distinct candidates). Reactions with zero matches are **coverage failures**;
reactions with many matches face **disambiguation** pressure.

```{code-cell}
:tags: [hide-input]
import matplotlib.pyplot as plt
import numpy as np
from synkit.Chem.Reaction.standardize import Standardize

_raw_fw = [entry.get('predictions', []) for entry in forward_results]

_std = Standardize()


def _chemical_unique_count(predictions):
    unique = set()
    for smarts in predictions:
        try:
            canonical = _std.fit(smarts)
        except Exception:
            canonical = None
        if canonical is not None:
            unique.add(canonical)
    return len(unique)


_counts = np.array([len(p) for p in _raw_fw], dtype=int)
_exact_counts = np.array([len(set(p)) for p in _raw_fw], dtype=int)
_chem_counts = np.array([_chemical_unique_count(p) for p in _raw_fw], dtype=int)

_zero = (_counts == 0).mean() * 100
_high = (_counts > 20).mean() * 100
_med = float(np.median(_counts))
_exact_med = float(np.median(_exact_counts))
_chem_med = float(np.median(_chem_counts))

fig, axes = plt.subplots(
    1,
    3,
    figsize=(14.5, 4.2),
    facecolor='white',
    gridspec_kw={'width_ratios': [1.25, 1.1, 0.9]},
)
fig.suptitle(
    'Forward template coverage, redundancy, and branching factor',
    fontsize=13,
    fontweight='bold',
    color='#0E1B2A',
)

ax = axes[0]
_bins = np.arange(0, min(_counts.max() + 2, 52))
ax.hist(
    _counts,
    bins=_bins,
    color='#24476F',
    edgecolor='white',
    linewidth=0.8,
    alpha=0.82,
)
ax.axvline(0.5, color='#B23A48', lw=1.8, linestyle='--', label='zero-match boundary')
ax.axvline(
    20.5, color='#9A6B17', lw=1.8, linestyle='--', label='high branching boundary'
)
ax.axvline(_med, color='#2E7D6B', lw=1.8, label=f'raw median={_med:.0f}')
ax.set_xlabel('Raw forward candidates generated', fontsize=10)
ax.set_ylabel('Test reactions', fontsize=10)
ax.set_title('Raw reactor output', fontsize=10, fontweight='bold')
ax.legend(fontsize=8, frameon=False)
ax.grid(axis='y', alpha=0.22)
ax.spines[['top', 'right']].set_visible(False)

ax = axes[1]
parts = ax.violinplot(
    [_counts, _exact_counts, _chem_counts], positions=[1, 2, 3], showmeans=False
)
for body, color in zip(parts['bodies'], ['#24476F', '#64748B', '#2E7D6B']):
    body.set_facecolor(color)
    body.set_edgecolor('none')
    body.set_alpha(0.35)
for key in ('cbars', 'cmins', 'cmaxes'):
    parts[key].set_color('#475569')
for x, arr, color in [
    (1, _counts, '#24476F'),
    (2, _exact_counts, '#64748B'),
    (3, _chem_counts, '#2E7D6B'),
]:
    q1, med, q3 = np.percentile(arr, [25, 50, 75])
    ax.scatter([x], [med], s=70, color=color, edgecolor='white', zorder=4)
    ax.vlines(x, q1, q3, color=color, lw=6, alpha=0.85)
    ax.text(x, q3 + max(1, q3 * 0.04), f'med {med:.0f}', ha='center', fontsize=8)
ax.set_xticks([1, 2, 3])
ax.set_xticklabels(['Raw', 'Exact\nunique', 'Chemical\nunique'])
ax.set_ylabel('Forward candidates per reaction', fontsize=10)
ax.set_title('Redundancy after deduplication', fontsize=10, fontweight='bold')
ax.grid(axis='y', alpha=0.22)
ax.spines[['top', 'right']].set_visible(False)

ax = axes[2]
ax.set_axis_off()
_metrics = [
    ('Raw zero-match', f'{_zero:.0f}%', '#B23A48'),
    ('Raw median', f'{_med:.0f}', '#24476F'),
    ('Exact-unique median', f'{_exact_med:.0f}', '#64748B'),
    ('Chemical-unique median', f'{_chem_med:.0f}', '#2E7D6B'),
]
for y, (label, value, color) in zip([0.84, 0.61, 0.38, 0.15], _metrics):
    ax.add_patch(
        plt.Rectangle((0.04, y - 0.085), 0.92, 0.16, fc='#F8FAFC', ec='#D7DEE8', lw=1)
    )
    ax.text(0.10, y + 0.018, value, fontsize=15, fontweight='bold', color=color)
    ax.text(0.10, y - 0.045, label, fontsize=8.2, color='#334155')

plt.tight_layout()
plt.show()
print(
    'Forward raw reactor output: '
    f'{_zero:.0f}% zero-match, median={_med:.0f} raw candidates. '
    f'After deduplication: median={_exact_med:.0f} exact-unique, '
    f'median={_chem_med:.0f} chemistry-unique candidates.'
)
```

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

test_processed = Parallel(n_jobs=4, verbose=0)(
    delayed(generate_forward_reactions)(entry, templates) for entry in test
)
```

This pattern separates **domain logic** (template application) from **execution logic**
(parallelization), which makes the notebook easier to scale and easier to maintain.
</details>

+++

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

+++

**Definition (Template-Based One-Step Prediction).**  
Given a template library $\mathcal{T} = \{p_1, \ldots, p_m\}$ of Double Pushout (DPO) graph transformation rules [@ehrig2006fundamentals] and a query molecule $q$, the *one-step forward prediction* generates the candidate set:

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

Define the standardized top-$K$ candidate set
$$
\mathcal{S}_K(q)=\{S(c)\mid c\in\hat{\mathcal{P}}_K(q)\}.
$$

$$
\mathrm{Recall@}K(q) = \mathbf{1}\!\left[S(q^*) \in \mathcal{S}_K(q)\right]
$$

where $\hat{\mathcal{P}}_K(q)$ is the top-$K$ candidates ranked by, e.g., template frequency. The *dataset Recall@K* is the mean over all test queries.

**Definition (Enrichment@K).**  
Recall only tells us whether the answer appears. It does **not** tell us how many
candidates a student or algorithm must inspect. We therefore also compute:

$$
\mathrm{Enrichment@}K(q)
= \frac{\left|\{s\in\mathcal{S}_K(q)\mid s=S(q^*)\}\right|}
       {\left|\mathcal{S}_K(q)\right|}.
$$

If $\mathcal{S}_K(q)$ is empty, the score is defined as $0$. With one ground-truth reaction, this is either $0$ or $1/n$ for a prefix
containing $n$ unique standardized candidates. Higher enrichment means the correct answer
is found in a more focused candidate set.

**Definition (Branching Factor).**  
The *branching factor* of a query $q$ is $|\hat{\mathcal{P}}(q)|$ — the total number of distinct candidates generated. A high branching factor indicates low template specificity; a low branching factor (possibly 0) indicates poor coverage.

+++

### Forward recall@K and enrichment@K trade-off

We now sweep $K$ from 1 to 50 for **forward prediction only**. Recall@K asks whether the ground-truth reaction appears in the top-$K$ candidates. Enrichment@K asks how concentrated that hit is within the inspected candidate prefix.

These curves should be read together: increasing $K$ can improve recall, but it usually lowers enrichment because the search space becomes broader.

```{code-cell}
import matplotlib.pyplot as plt
import numpy as np
from synkit.Chem.Reaction.standardize import Standardize

_std_for_curve = Standardize()


def _standardize_reaction(rsmi):
    try:
        return _std_for_curve.fit(rsmi)
    except Exception:
        return None


def _unique_standardized_prefix(values, k):
    out = []
    seen = set()
    for value in values[:k]:
        key = _standardize_reaction(value)
        if key is None or key in seen:
            continue
        seen.add(key)
        out.append(key)
    return out


def recall_enrichment_at_k(results, max_k=50):
    ks = np.arange(1, max_k + 1)
    recalls = []
    enrichments = []
    candidate_counts = []

    for k in ks:
        hit_values = []
        enrichment_values = []
        count_values = []
        for entry in results:
            gt = _standardize_reaction(entry.get('ground_truth', ''))
            preds = _unique_standardized_prefix(entry.get('predictions', []), k)
            hit = int(gt is not None and gt in preds)
            hit_values.append(hit)
            count_values.append(len(preds))
            enrichment_values.append(hit / len(preds) if preds else 0.0)
        recalls.append(float(np.mean(hit_values)) if hit_values else 0.0)
        enrichments.append(
            float(np.mean(enrichment_values)) if enrichment_values else 0.0
        )
        candidate_counts.append(float(np.mean(count_values)) if count_values else 0.0)

    return ks, np.array(recalls), np.array(enrichments), np.array(candidate_counts)


try:
    ks, rec_std, enr_std, n_std = recall_enrichment_at_k(forward_results)
except NameError:
    ks = np.arange(1, 51)
    rec_std = 1 - np.exp(-ks / 8)
    enr_std = rec_std / np.maximum(ks, 1)
    n_std = ks.astype(float)
    print('Demo mode: replace forward_results with actual data.')

fig, (ax_left, ax_right) = plt.subplots(1, 2, figsize=(13.5, 4.8), facecolor='white')
fig.suptitle(
    'Forward prediction: recall@K versus enrichment@K',
    fontsize=13,
    fontweight='bold',
    color='#0E1B2A',
)

ax_left.plot(ks, rec_std, color='#0E1B2A', lw=2.4, label='Standardized recall@K')
ax_left.fill_between(ks, 0, rec_std, color='#D8E7F5', alpha=0.65)
ax_left.set_ylabel('Recall@K', fontsize=10)
ax_left.set_xlabel('Rank cutoff K', fontsize=10)
ax_left.set_title('Finding the answer', fontsize=11, fontweight='bold')
ax_left.set_xlim(1, 50)
ax_left.set_ylim(0, 1.05)
ax_left.grid(alpha=0.22)
ax_left.spines[['top', 'right']].set_visible(False)

ax_right.plot(ks, enr_std, color='#2E7D6B', lw=2.4, label='Standardized enrichment@K')
ax_right.set_ylabel('Mean enrichment@K', fontsize=10)
ax_right.set_xlabel('Rank cutoff K', fontsize=10)
ax_right.set_title('Keeping the search focused', fontsize=11, fontweight='bold')
ax_right.set_xlim(1, 50)
ax_right.grid(alpha=0.22)
ax_right.spines[['top', 'right']].set_visible(False)

for ax in (ax_left, ax_right):
    ax.axvline(10, color='#7BA7D7', lw=1.4, ls=':', alpha=0.9)
    ax.text(
        10.4,
        ax.get_ylim()[0] + (ax.get_ylim()[1] - ax.get_ylim()[0]) * 0.06,
        'K=10',
        color='#52789E',
        fontsize=8,
    )
    ax.legend(frameon=False, fontsize=9, loc='best')

for k in (1, 10, 50):
    idx = k - 1
    ax_left.scatter([k], [rec_std[idx]], s=35, color='#2E7D6B', zorder=3)
    ax_left.text(
        k + 0.8, rec_std[idx] + 0.04, f'{rec_std[idx]:.2f}', fontsize=8, color='#16634F'
    )
    ax_right.scatter([k], [enr_std[idx]], s=35, color='#2E7D6B', zorder=3)
    ax_right.text(
        k + 0.8,
        enr_std[idx] + max(enr_std.max(), 0.01) * 0.04,
        f'{enr_std[idx]:.3f}',
        fontsize=8,
        color='#16634F',
    )

plt.tight_layout(rect=[0, 0, 1, 0.93])
plt.show()

print(
    f'At K=10: standardized recall={rec_std[9]:.3f}, '
    f'standardized enrichment={enr_std[9]:.4f}, '
    f'mean unique standardized candidates={n_std[9]:.1f}.'
)
```

```{code-cell}
# Inspect one example and evaluate it after standardization
idx = 0
ground_truth = eval_set[idx]['smart']
fw_list = eval_set[idx]['fw']
len(fw_list)
```

### Forward prediction gallery — ground truth vs candidates

For one example substrate, the template library generates a set of candidate reactions. The ground-truth reaction is highlighted with a green check mark. Comparing candidates to the ground truth directly shows whether the correct chemistry appears in the top predictions and which template was responsible.

```{code-cell}
:tags: [hide-input]
from synedu.Utils.rxn_vis import render_reaction_gallery
from IPython.display import HTML, display
# The gallery scales each complete RDKit SVG into its card.
display(HTML(render_reaction_gallery(
    ground_truth,
    fw_list,
    title="Forward prediction gallery",
    description="Top {shown} generated candidates from {total} total; {hits} exact hit(s).",
)))
```

A first cleanup step is to remove **duplicate predictions** before computing metrics.
This prevents repeated candidates from artificially inflating the denominator in
enrichment calculations.

For one test example, the reactor may generate **hundreds of candidate reactions**.
However, many of them are not genuinely different transformations. Deduplication lets us
separate **search breadth** from **true reaction diversity**.

From an educational standpoint, this is a good reminder that evaluation metrics are only
as meaningful as the representation they operate on. Poor post-processing can make a
predictor look worse—or occasionally better—than it actually is.

+++

### 2.1. Standardized evaluation

We now evaluate generated reactions in the same way for every example:

1. **standardize** the generated reactions,
2. **remove duplicate standardized candidates**, and
3. compare them against the standardized ground truth.

This collapses formatting and representation-level variants and focuses the score on the
recovered reaction transformation. For many practical reaction-prediction use cases, this
is the most direct first measure because end users want to know whether the correct
chemistry appears in the generated candidate set.

```{code-cell}
from synkit.Chem.Reaction.standardize import Standardize

std = Standardize()


def recall(smiles_list, ground_truth):
    return int(ground_truth in smiles_list)


def enrichment(smiles_list, ground_truth):
    if not smiles_list:
        return 0.0
    return sum(s == ground_truth for s in smiles_list) / len(smiles_list)


std_fw = []
for smarts in fw_list:
    canonical = std.fit(smarts)
    if canonical is not None:
        std_fw.append(canonical)

std_fw = list(set(std_fw))
std_gt = std.fit(ground_truth)
len(std_fw)
```

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

+++

## 3. Backward prediction

In **backward prediction**, we reverse the perspective: instead of mapping reactants to
products, we start from the **product side** and enumerate plausible precursor reactions
by applying the same template library in **inverted mode**.

This is closely related to **single-step retrosynthesis** [@segler2018planning; @coley2019robotic], but it is usually more
ambiguous than forward prediction. A single product can often be disconnected in multiple
valid ways, and the backward direction is strongly influenced by selectivity, leaving
groups, protecting-group strategy, reagent availability, and reaction conditions. As a
result, backward prediction generally exhibits **larger branching** and **lower
specificity** than the forward task.

Educationally, this section is valuable because it shows that the same reaction rule can
support two distinct tasks depending on the direction in which it is applied.

+++

### 3.1. Inspect one backward prediction example

We first invert a single template application to see how the product-side input expands
into candidate precursor reactions.

This example is useful for building retrosynthetic intuition. In the forward direction,
we usually ask, *“What product can this substrate give?”* In the backward direction, we
instead ask, *“What precursor sets could plausibly lead to this product under this local
transformation?”* That inversion greatly enlarges the space of valid answers.

```{code-cell}
# Reuse the same example and template in backward mode
p = eval_set[EXAMPLE_INDEX]["p"]
H = smiles_to_graph(p)
t = templates[TEMPLATE_INDEX]["RC"]
gt = eval_set[EXAMPLE_INDEX]["smart"]
print(f"Example reaction index: {EXAMPLE_INDEX}")
print(f"Selected template index: {TEMPLATE_INDEX}")
draw_molecular_graph(H)
```

```{code-cell}
# Apply the same template in inverted mode for backward prediction
from synkit.Synthesis.Reactor.syn_reactor import SynReactor

reactor = SynReactor(
    substrate=p,
    template=t,
    explicit_h=True,
    implicit_temp=False,
    automorphism=True,
    invert=True,  # invert the transformation for backward prediction
)
bw = reactor.smarts_list
print(f"Generated {len(bw)} raw precursor candidate(s).")
bw
```

### 3.2. Scale backward prediction and summarize performance

As in the forward direction, we next apply the full template library to the dataset,
post-process the generated reactions, and compute macro-level summary metrics.

When interpreting the backward results, it is important to keep expectations realistic.
Lower recall or lower enrichment does not necessarily mean the templates are poor. It may
also reflect the intrinsic fact that **many retrosynthetic disconnections are chemically
reasonable**, even if only one corresponds to the recorded reference reaction.

```{code-cell}
for value in eval_set:
    p = value['p']
    result = []
    for t in templates:
        reactor = SynReactor(
            substrate=p,
            template=t['RC'],
            explicit_h=True,
            implicit_temp=False,
            automorphism=True,
            invert=True,
        )

        sols = reactor.smarts_list
        result.extend(sols)
    value['bw'] = result

# Build backward_results in recall@K format
backward_results = [
    {"ground_truth": v["smart"], "predictions": v["bw"]} for v in eval_set
]
print(
    f"Backward: {len(backward_results)} entries, "
    f"{sum(len(r['predictions']) for r in backward_results)} total candidates"
)
```

### Backward prediction gallery — ground truth vs candidates

For the same example entry, the inverted template library generates candidate **precursors**. Each candidate is a plausible retrosynthetic disconnection; the ground-truth reaction appears among the candidates when the correct bond-breaking template fires.

```{code-cell}
:tags: [hide-input]
from synedu.Utils.rxn_vis import render_reaction_gallery
from IPython.display import HTML, display
_bidx = 0
_bw_gt = eval_set[_bidx]['smart']
_bw_list = eval_set[_bidx]['bw']
# Use graph-native layout fallback for candidates with invalid RDKit coordinates.
display(HTML(render_reaction_gallery(
    _bw_gt,
    _bw_list,
    title="Backward prediction gallery",
    description="Top {shown} retrosynthetic candidates from {total} total; {hits} exact hit(s).",
    candidate_badge_background="#FFF7E8",
    candidate_badge_foreground="#8A5A12",
)))
```

```{code-cell}
# Helper functions for dataset-level post-processing
from synkit.Chem.Reaction.standardize import Standardize

std = Standardize()


def std_list(list_rsmi):
    std_lst = []
    for smarts in list_rsmi:
        canonical = std.fit(smarts)
        if canonical is not None:
            std_lst.append(canonical)

    std_lst = list(set(std_lst))
    return std_lst
```

```{code-cell}
from copy import deepcopy
from joblib import Parallel, delayed


def process_entry(value):
    value = deepcopy(value)
    value["fw_std"] = std_list(value["fw"])
    value["bw_std"] = std_list(value["bw"])
    value["gt_std"] = std_list([value["smart"]])[0]
    return value


eval_set = Parallel(n_jobs=4, verbose=0)(
    delayed(process_entry)(value) for value in eval_set
)
```

```{code-cell}
from copy import deepcopy
from joblib import Parallel, delayed


def recall(smiles_list, ground_truth):
    return int(ground_truth in smiles_list)


def enrichment(smiles_list, recall_value):
    if not smiles_list:
        return 0.0
    return recall_value / len(smiles_list)


def evaluate_entry(entry):
    """
    Compute standardized recall and enrichment for one entry.

    Expected keys
    -------------
    fw_std, bw_std, gt_std
    """
    fw_recall = recall(entry["fw_std"], entry["gt_std"])
    bw_recall = recall(entry["bw_std"], entry["gt_std"])

    return {
        "fw_recall": fw_recall,
        "fw_enrichment": enrichment(entry["fw_std"], fw_recall),
        "bw_recall": bw_recall,
        "bw_enrichment": enrichment(entry["bw_std"], bw_recall),
    }


def add_metrics(entry):
    """
    Return a copy of one entry with evaluation metrics added.
    """
    value = deepcopy(entry)
    value.update(evaluate_entry(value))
    return value


def evaluate_entries(entries, n_jobs=4, verbose=0):
    """
    Evaluate all entries in parallel and return the updated list.
    """
    return Parallel(n_jobs=n_jobs, verbose=verbose)(
        delayed(add_metrics)(entry) for entry in entries
    )


def summarize_metrics(entries):
    """
    Macro-average metrics across all entries.
    """
    n = len(entries)
    if n == 0:
        return {}

    return {
        "fw_recall_mean": sum(x["fw_recall"] for x in entries) / n,
        "fw_enrichment_mean": sum(x["fw_enrichment"] for x in entries) / n,
        "bw_recall_mean": sum(x["bw_recall"] for x in entries) / n,
        "bw_enrichment_mean": sum(x["bw_enrichment"] for x in entries) / n,
    }


# -------------------------
# Usage
# -------------------------

test_eval = evaluate_entries(eval_set, n_jobs=4, verbose=0)
summary = summarize_metrics(test_eval)

print(summary)
```

```{code-cell}
summary
```

### Prediction performance summary

After evaluating all test reactions, the bar chart below aggregates standardized **mean
recall** and **mean enrichment** for forward prediction (blue) and backward /
retrosynthetic prediction (gold). Recall shows whether the reference transformation was
recovered at least once; enrichment shows how concentrated that hit was within the unique
candidate set.

```{code-cell}
:tags: [hide-input]
import matplotlib.pyplot as plt
import pandas as pd

_sum_df = pd.DataFrame([summary])
print(f"Evaluated test reactions: {len(test_eval)}")
_numeric = _sum_df.select_dtypes("number")
if len(_numeric.columns) > 0:
    _means = _numeric.mean().sort_values(ascending=True)
    _colors = []
    for c in _means.index:
        cl = c.lower()
        if any(k in cl for k in ("fw", "forward", "fwd")):
            _colors.append("#24476F")
        elif any(k in cl for k in ("bw", "backward", "retro", "rv")):
            _colors.append("#D49A44")
        else:
            _colors.append("#2E7D6B")
    fig, ax = plt.subplots(
        figsize=(8.5, max(3.2, len(_means) * 0.55)), facecolor="white"
    )
    bars = ax.barh(
        _means.index, _means.values, color=_colors, alpha=0.92, edgecolor="white"
    )
    ax.bar_label(bars, fmt="%.3f", padding=4, fontsize=9)
    ax.set_xlim(0, min(1.05, max(1.0, _means.max() * 1.15 + 0.03)))
    ax.set_xlabel("Mean score across test set", fontsize=10)
    ax.set_title(
        "One-step prediction evaluation summary",
        fontsize=12,
        fontweight="bold",
        color="#0E1B2A",
    )
    ax.grid(axis="x", alpha=0.22)
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis='y', labelsize=9)
    plt.tight_layout()
    plt.show()
else:
    display(_sum_df.head())
```

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

+++

## 5. Quiz

Answer using both **chemical intuition** and **graph-based reasoning**.

1. Why can the number of unique predictions decrease after reaction standardization?
2. Why is it useful to report both recall and enrichment when evaluating generated candidates?
3. Why is backward prediction usually more ambiguous than forward prediction?
4. If two template libraries have the same recall but different enrichment, what does that imply about the focus of the generated search space?

+++

## 6. References

```{bibliography}
```
