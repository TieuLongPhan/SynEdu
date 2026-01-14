# S07 · One-Step Rule Application (Forward & Backward): DPO Rules as Reaction Engines

<div class="alert alert-block alert-info">
<b>Welcome to SynEdu.</b><br>
This talktorial is part of <b>SynEdu</b>, a lightweight, research-oriented teaching series built around the
<b>Syn</b> ecosystem and <b>RDKit</b> for practical, reproducible reaction modeling.
</div>

<div class="alert alert-block alert-success">
<b>What you will gain.</b><br>
We turn extracted DPO rules into a practical <b>one-step</b> reaction engine.
You will implement <b>forward synthesis</b> (reactants → products) and <b>retrosynthesis</b>
(products → precursors) using the <b>same</b> rule, via rule inversion.
We also discuss deduplication (symmetry) and basic chemical sanity checks.
</div>

## Aim of this talktorial

## Learning outcomes
By the end of this notebook, you should be able to:
- Recall the DPO span \(L \leftarrow K \rightarrow R\) and what a match \(m: L \hookrightarrow G\) means.
- Implement a practical DPO application that produces a product graph \(G'\) from a host graph \(G\).
- Implement rule inversion \(p^{-1}: R \leftarrow K \rightarrow L\) for one-step retrosynthesis.
- Run a small one-step experiment: reproduce products for a subset of mapped reactions.

## Roadmap
- 0. Setup & data (mapped reactions)
- 1. Theory: DPO application forward/backward
- 2. Implementation: matching, rewriting, and conversion back to RDKit
- 3. Demonstrations (one rule; then many reactions)
- 4. Exercises



## 1. Theory: one-step rewriting in both directions

A DPO rule is a span of typed graphs
\[
p:\quad L \xleftarrow{\,l\,} K \xrightarrow{\,r\,} R
\]
where:
- \(L\) is the **pre-condition** (what must be present to react),
- \(R\) is the **post-condition** (what is produced),
- \(K\) is the **preserved interface** (what stays the same).

Given a host graph \(G\) (a set of reactant molecules as a disjoint union),
a **match** is an injective morphism
\[
m: L \hookrightarrow G
\]
that respects atom and bond labels.

### Forward (synthesis)
The DPO construction removes the part \(L\setminus K\) and glues in \(R\setminus K\), producing \(G'\).

### Backward (retrosynthesis)
The inverse rule is simply
\[
p^{-1}:\quad R \xleftarrow{\,r\,} K \xrightarrow{\,l\,} L
\]
and applying \(p^{-1}\) to a product graph generates plausible precursors.

**Chemistry meaning:** a rule extracted from atom-mapped data can be used in either direction as a
one-step generator; directionality is controlled by which side you treat as “host”.



## 3. Discussion

A self-reproduction study (rule extracted from a reaction, then applied back to its own reactants/products)
is a **sanity check**, not a benchmark:

- If the hit rate is low, it usually indicates one of:
  - parsing/mapping inconsistencies in the dataset,
  - missing atoms (unbalanced reactions),
  - overly strict node/edge matching predicates,
  - stereochemistry or aromaticity handling issues.

In S07 we will move from “self reproduction” to **library-level evaluation**:
apply a *set* of rules to each host and study precision/branching trade-offs.



## 4. Exercises

1. **Strictness toggle.** Modify `node_match()` to include/exclude `formal_charge` or `aromatic`.
   How does forward/backward hit rate change?
2. **Context radius.** Repeat `self_reproduction_study` for `core_radius=0,1,2`.
   Does a larger context reduce spurious candidates?
3. **Match explosion.** Find a reaction where `n_forward_cand` is large.
   Inspect the matches and identify a symmetry source.
