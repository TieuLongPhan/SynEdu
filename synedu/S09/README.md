# S09 · Context Graph Expansion: Systematic Control of Rule Specificity

<div class="alert alert-block alert-info">
<b>Welcome to SynEdu.</b><br>
This talktorial is part of <b>SynEdu</b>, a lightweight, research-oriented teaching series built around the
<b>Syn</b> ecosystem and <b>RDKit</b> for practical, reproducible reaction modeling.
</div>

<div class="alert alert-block alert-success">
<b>What you will gain.</b><br>
We introduce <b>context expansion</b> as a principled knob to control how general or specific a reaction rule is.
Starting from a reaction center, we expand the context by a graph radius \(r\) to obtain a family of rules
\(p_r\). We then measure how \(r\) changes <b>coverage</b>, <b>branching</b>, and <b>hit rate</b>.
</div>

## Aim of this talktorial

## Learning outcomes
By the end of this notebook, you should be able to:
- Define a reaction center \(C\) on the ITS graph and its radius-\(r\) neighborhood \(N_r(C)\).
- Construct a rule family \(p_r\) by extracting \(L_r,K_r,R_r\) on that neighborhood.
- Evaluate rule libraries across radii \(r=0,1,2,\dots\) with the S07 metrics.
- Interpret the specificity–coverage–branching trade-off and choose \(r\) for an application.

## Roadmap
- 0. Setup & data
- 1. Theory: neighborhoods and rule families
- 2. Implementation: extract libraries for multiple radii
- 3. Trade-off curves (coverage vs branching, hit vs branching)
- 4. Exercises



## 1. Theory: neighborhoods and rule families

Let \(G_{\mathrm{ITS}}\) be the ITS graph built from an atom-mapped reaction (S04).

### 1.1 Reaction center
We define the **core center** \(C\subseteq V(G_{\mathrm{ITS}})\) as the set of atoms incident to a changed bond:
\[
C = \{\,v\in V : \exists (u,v)\in E,\ \mathrm{order}_r(u,v)\neq \mathrm{order}_p(u,v)\,\}
\]
(optionally also including formal-charge changes).

### 1.2 Context expansion (radius-\(r\) neighborhood)

For \(r\ge 0\), define the radius-\(r\) neighborhood:
\[
N_r(C) = \{\,v\in V : \mathrm{dist}(v, C)\le r\,\}.
\]

Chemistry meaning:
- \(r=0\): only atoms directly involved in bond changes (most general, most ambiguous).
- \(r=1\): include immediate neighbors (often enough to identify functional groups).
- larger \(r\): include more local environment (more specific, less transferable).

### 1.3 Rule family
From \(N_r(C)\) we extract a DPO rule:
\[
p_r:\quad L_r \leftarrow K_r \rightarrow R_r.
\]
As \(r\) increases, \(L_r\) contains more constraints, so matches become rarer:
\[
\mathrm{Match}(L_{r+1},G)\ \subseteq\ \mathrm{Match}(L_{r},G)
\quad\Rightarrow\quad
\text{coverage and branching typically decrease.}
\]

In the next section we make this trade-off visible with the S07 metrics.



## 3. Choosing a radius in practice

A reasonable workflow:

1. Start with \(r=0\) to identify **core transformation types** (high generality).
2. Increase to \(r=1\) to reduce symmetry-induced ambiguity and “too-many matches”.
3. Use \(r=2\) or \(r=3\) if you need high precision (e.g., retrosynthesis planning),
   but expect lower coverage and more fragmented rule libraries.

If you see:
- **high branching at small \(r\)**: consider increasing \(r\), or tightening node/edge match predicates.
- **low coverage at large \(r\)**: consider decreasing \(r\), or normalizing charges/tautomers upstream.

> Important: context is not the only knob. In real pipelines, ranking (learned or heuristic),
> reagent handling, and stereochemistry decisions often dominate performance.



## 4. Exercises

1. **Backward evaluation.** Repeat the multi-radius study with `direction="backward"`.
   Is the best radius the same as forward?
2. **Center definition.** Modify `reaction_center_core_nodes` to include charge changes or exclude aromatic changes.
   How does the trade-off curve move?
3. **Hybrid context.** Implement a “ring-aware” expansion:
   if any core atom is in a ring, include the full ring as context (even if it exceeds radius).
