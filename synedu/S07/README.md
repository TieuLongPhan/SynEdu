# S06 · Canonicalizing Atom-Mapped Reactions and Rules: Equivalence, Symmetry, Determinism

<div class="alert alert-block alert-info">
<b>Welcome to SynEdu.</b><br>
This talktorial is part of <b>SynEdu</b>, a lightweight, research-oriented teaching series built around the
<b>Syn</b> ecosystem and <b>RDKit</b> for practical, reproducible reaction modeling.
</div>

<div class="alert alert-block alert-success">
<b>What you will gain.</b><br>
Atom-mapped reactions are <b>not unique</b>: the same chemistry may admit many valid maps and many map-numberings.
In this notebook we turn mapped reactions into <b>deterministic, comparable objects</b> by canonicalizing (i) the <b>reaction string</b>
(order of components) and (ii) the <b>atom-map numbering</b>. We then show how canonicalization stabilizes
<b>reaction-center clustering</b> and <b>DPO rule extraction</b> downstream.
</div>

## Aim of this talktorial

## Learning outcomes
By the end of this notebook, you should be able to:
- Explain why atom maps are non-unique (symmetry) and why map numbers should be canonicalized.
- Implement a deterministic atom-map reindexing \(\pi: \mathbb{N} \to \mathbb{N}\) based on map-invariant structural ranks.
- Canonicalize a mapped reaction SMILES (molecule order + map ids) without changing chemistry.
- Quantify how canonicalization affects the number of unique centers/rules (hash + isomorphism).

## Roadmap
- 0. Setup & data (mapped reactions from `./data/smart.json.gz`)
- 1. Theory: equivalence of maps and why “canonical” matters
- 2. Canonicalizing mapped reaction SMILES
- 3. From canonical maps to canonical DPO rules
- 4. Mini-study: duplicates before/after
- 5. Exercises

<div class="alert alert-block alert-warning">
<b>Note.</b> This notebook assumes your dataset already contains <b>atom-mapped</b> reaction SMILES
(e.g., from RXNMapper). If the file is not present, the code will raise a clear error.
</div>



## 0. Setup & data

We load a JSON (optionally gzipped) file that contains **atom-mapped reaction SMILES**.
Throughout SynEdu we will treat the mapped reaction string as *raw evidence* from which we derive
graph objects (ITS, centers, rules).

We will do two things immediately:

1. Normalize the JSON into a list of records (`records`).
2. Guess which field contains the mapped reaction SMILES (`mapped_key`), since datasets vary.

> In many datasets the mapped reaction may appear under keys like `mapped_rxn`, `rxn_smiles_mapped`, `rxn`, etc.



## 1. Theory (formal, but chemist-friendly)

### 1.1 Atom maps as witnesses, not identities

For a reaction \(r\), let \(\mathcal{M}(r)\) denote the set of **valid atom mappings**.
A mapping \(m\in\mathcal{M}(r)\) is best viewed as a **witness** that reactant atoms correspond to product atoms.

In practice, \(|\mathcal{M}(r)|\) is often \(>1\) because of **molecular symmetry**:
if a molecule has automorphisms (e.g., benzene), then multiple atom correspondences are equally valid.

### 1.2 Two sources of non-uniqueness

Even if you fix the chemistry, two mapped reaction strings may differ by:

1. **Component order** (which molecule appears first on each side).
2. **Map-number permutation** (renaming atom map ids).

Formally, a **map renaming** is a bijection
\[
\pi:\{1,\dots,n\}\to\{1,\dots,n\},
\]
and it induces an equivalent mapped reaction by replacing every \([*:i]\) with \([*:\pi(i)]\).

### 1.3 Canonicalization as a quotient

We define an equivalence relation on mapped reactions:
\[
\rho_1 \sim \rho_2
\quad\Longleftrightarrow\quad
\rho_2 = \pi(\rho_1) \text{ for some map renaming }\pi \text{ (and possibly permuted components).}
\]

A **canonicalization** procedure chooses a representative
\[
\operatorname{can}(\rho)\in[\rho]
\]
so that equivalent inputs yield the same output.

**Chemistry meaning:** we do not change the reaction—only its *encoding*.

### 1.4 Why this matters downstream

If you extract a rule (or a reaction center) from a mapped reaction, map-number permutations
can create many syntactically different but chemically identical rules.

Canonicalization makes:
- **rule libraries reproducible** (stable hashing and clustering),
- **statistics meaningful** (counts of “unique rules” stop drifting),
- **benchmarking fair** (different mappers can be compared under the same quotient).



## 2. Canonicalizing mapped reaction SMILES

We canonicalize a mapped reaction in two stages:

1. **Canonicalize molecule order** on each side by sorting molecules by their **unmapped canonical SMILES**.
2. **Canonicalize map ids** by assigning new ids \(1,2,\dots\) in a deterministic order derived from
   **map-invariant RDKit canonical atom ranks** on each molecule.

### 2.1 Deterministic reindexing rule (practical)

For every atom-map id \(i\), we build a key

\[
\kappa(i) = \big(\text{side},\ \text{mol\_key},\ \text{rank},\ \text{symbol},\ \text{charge},\ \text{aromatic}\big),
\]

where:
- `side = 0` if the mapped atom occurs on the **reactant** side, else `1` (product-only atoms, if any);
- `mol_key` is the **unmapped canonical SMILES** of the parent molecule;
- `rank` is the RDKit **canonical rank** of the atom in that unmapped molecule.

Then we sort map ids by \(\kappa(i)\) and assign new ids \(1,2,\dots\).

This does **not** solve symmetry in a philosophically perfect way (symmetry cannot always be broken uniquely),
but RDKit's canonicalization gives a consistent tie-breaker that is independent of the original atom numbering.



## 3. From canonical maps to canonical DPO rules

Once a mapped reaction is canonicalized, every downstream graph we build becomes more stable:

- ITS graphs use **map ids as node identifiers**.
- Reaction centers are extracted as subgraphs of the ITS.
- DPO rules \(L \leftarrow K \rightarrow R\) reuse those ids.

So if we canonicalize map ids up front, then:
- identical chemistry yields identical ITS node ids,
- “unique rule counts” stop depending on incidental numbering.

Below we reuse the same ITS + center + rule extraction ideas as S04, but apply them to:
- raw mapped reactions, and
- canonicalized mapped reactions,
to quantify the difference.



## 4. Discussion

If canonicalization works as intended, then:

- Many reactions that only differ by *incidental map numbering* collapse to the same rule signature.
- The fraction `same?` (raw vs canonicalized) may be less than 1 because **raw datasets can contain**
  differences beyond simple renumbering (e.g., different mapping conventions, incomplete mapping, or agents handling).
  Canonicalization does **not** fix those—by design.

In S06–S08 we will assume rules were extracted from **canonicalized mapped reactions** to get:
- stable identifiers for rules,
- stable clustering of reaction centers,
- comparable evaluation metrics.



## 5. Exercises

1. **Alternative tie-breaker.** Modify `canonical_map_renaming` to use product-side ranks first.  
   Compare how many signatures change.
2. **Radius sensitivity.** Repeat the mini-study for `core_radius=0` and `core_radius=2`.  
   Does canonicalization reduce duplicates more at small or large context?
3. **Sanity test.** Pick a reaction, permute map ids with different seeds, and verify canonicalization collapses them.
