# S08 · Metrics for Rule Application: Coverage, Recall, Branching, and Cost

<div class="alert alert-block alert-info">
<b>Welcome to SynEdu.</b><br>
This talktorial is part of <b>SynEdu</b>, a lightweight, research-oriented teaching series built around the
<b>Syn</b> ecosystem and <b>RDKit</b> for practical, reproducible reaction modeling.
</div>

<div class="alert alert-block alert-success">
<b>What you will gain.</b><br>
Rules are only useful if we can <b>measure</b> how they behave.
In this notebook we define a practical evaluation protocol for one-step rule application and implement
metrics that capture the <b>precision–recall–cost</b> trade-offs that appear in rule-based reaction modeling.
</div>

## Aim of this talktorial

## Learning outcomes
By the end of this notebook, you should be able to:
- Define one-step tasks (forward prediction and retrosynthesis) in terms of candidate generation.
- Implement coverage, top-\(k\) recall, branching factor, and simple cost metrics.
- Build a rule library from mapped reactions and evaluate it on a held-out subset.
- Understand why context size, matching predicates, and deduplication change the metrics.

## Roadmap
- 0. Setup & data
- 1. Theory: what to measure (and why)
- 2. Building a rule library (canonicalized)
- 3. One-step evaluation protocol
- 4. Results & interpretation
- 5. Exercises



## 1. Theory: what to measure

A one-step rule engine is a **candidate generator**.

Given:
- a host \(G\) (reactants, or products for retrosynthesis),
- a rule library \(\mathcal{P}=\{p_1,\dots,p_M\}\),

we generate a set of candidates:
\[
\mathcal{C}(G) = \bigcup_{p\in\mathcal{P}} \{\,\text{apply}(G,p,m)\ :\ m \in \mathrm{Match}(p,G)\,\}.
\]

Let the ground truth target for \(G\) be \(T(G)\) (e.g., the true products of that reaction).

### Metrics

**Coverage**
\[
\mathrm{cov} = \Pr\big(|\mathcal{C}(G)|>0\big)
\]
Fraction of hosts for which at least one rule applies.

**Top-\(k\) recall** (candidate contains truth)
\[
\mathrm{rec@}k = \Pr\big(T(G) \in \mathrm{TopK}(\mathcal{C}(G))\big).
\]
In this notebook we will use an unranked version (set membership), i.e. \(k=\infty\), because ranking is a separate topic.

**Branching factor** (candidate explosion)
\[
b = \mathbb{E}\big[|\mathcal{C}(G)|\big].
\]

**Cost**
- runtime per host,
- number of attempted rule matches,
- number of successful matches.

These metrics typically trade off:
- **more general rules** → higher coverage/recall, higher branching,
- **more specific rules (larger context)** → lower branching, but risk lower coverage.

S08 will make that trade-off systematic by varying context radius.



## 4. Interpretation: what to look for

Typical patterns you will see:

- High **coverage** but low **hit rate** often means rules apply broadly but generate many wrong candidates
  (high branching).
- High **hit rate** with very high branching indicates the truth is present but drowned in many alternatives—
  you need **ranking** (future notebook) or more context (S08).
- Low coverage can mean rules are too specific (context too large, or node/edge matching too strict).

Also note:
- “Hit rate” here is unranked set membership (a generous metric).
- Branching is sensitive to symmetry and deduplication strategy.

In S08 we will vary context radius systematically and plot the trade-off curve.



## 5. Exercises

1. **Indexing.** Improve prefiltering:
   add an edge-count multiset or a degree multiset to the necessary-condition test.
2. **Rule size vs branching.** Plot `nL` (rule size) vs average branching contributed by that rule.
3. **Sanity check.** Reduce `max_matches_per_rule` to 1 and see how branching and hit rate change.
