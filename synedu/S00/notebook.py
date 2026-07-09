# ---
# jupyter:
#   jupytext:
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.19.4
#   kernelspec:
#     display_name: synedu
#     language: python
#     name: python3
# ---

# %% [markdown]
# # S00 : [Short Descriptive Title Here]

# %% [markdown]
# ## Aim of this talktorial
#
# <!-- One paragraph: what problem does this talktorial solve?
#      Reference the previous talktorial and state what gap this one fills. -->
#
# Building on **[SXX]**, this talktorial (**S00**) introduces **[CORE CONCEPT]**
# for [CHEMICAL / GRAPH APPLICATION].
#
# Concretely, we focus on:
#
# 1. **[Topic A]**  
#    [One-sentence description of what the learner will do / build.]
#
# 2. **[Topic B]**  
#    [One-sentence description.]
#
# 3. **[Topic C]** *(optional)*  
#    [One-sentence description.]
#
#
# ## Learning outcomes
#
# After completing S00, you can:
#
# 1. [Verb] [concept] — e.g., *Explain why atom maps are non-unique.*
# 2. [Verb] [concept]
# 3. [Verb] [concept]
# 4. [Verb] [concept]
#
#
# ## Outline
#
# <ul class="synedu-outline">
#   <li><a href="#0-setup--data">0. Setup &amp; data</a></li>
#   <li><a href="#1-section-one-title">1. Section One Title</a></li>
#   <li><a href="#2-section-two-title">2. Section Two Title</a></li>
#   <li><a href="#3-section-three-title">3. Section Three Title</a></li>
#   <li><a href="#4-quiz">4. Quiz</a></li>
#   <li><a href="#5-discussion">5. Discussion</a></li>
#   <li><a href="#6-references">6. References</a></li>
# </ul>

# %% [markdown]
# ### Pipeline overview
#
# <!-- Optional: edit figure/pipeline.tex to match your talktorial workflow,
#      then compile it with `latexmk -pdf pipeline.tex`. Delete this cell if
#      no static pipeline figure is needed. -->

# %% [markdown]
# <figure style="text-align:center; margin: 1.5rem 0;">
#   <object data="figure/pipeline.pdf" type="application/pdf" width="100%" height="360">
#     <a href="figure/pipeline.pdf">Open the S00 pipeline overview PDF</a>
#   </object>
#   <figcaption style="font-size:0.88em; color:#555; margin-top:0.4rem;">
#     <b>Figure&nbsp;1.</b> S00 pipeline overview. Edit <code>figure/pipeline.tex</code>
#     and recompile it to update this PDF figure.
#   </figcaption>
# </figure>

# %% [markdown]
# <a id="0-setup--data"></a>
#
# ## 0. Setup & data

# %% [markdown]
# <!-- State which package(s) this notebook depends on and how to install them.
#      Use <a href="#ref-N">[N]</a> for inline citations. -->
#
# This notebook uses **[Package]** <a href="#ref-1">[1]</a> and **[Package2]** <a href="#ref-2">[2]</a>.
# Install via:
# ```bash
# pip install [package] [package2]
# ```

# %%
# %matplotlib inline
import rdkit
import networkx as nx
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import importlib.metadata as m

# Synedu utilities — extend this list as needed
from synedu.Utils import draw_molecular_graph, draw_rxn_graph

print("RDKit version:", rdkit.__version__)
print("NetworkX version:", nx.__version__)

# ── Load dataset ─────────────────────────────────────────────────────────
# Replace with your actual data path and format.
DATA_DIR = Path("data")
DATA_PATH = DATA_DIR / "your_dataset.json.gz"

# from synkit.IO import load_database
# data = pd.DataFrame(load_database(DATA_PATH)[:10000])
# display(data.head())
# print(data.shape)

# %% [markdown]
# <a id="1-section-one-title"></a>
#
# ## 1. Section One Title
#
# <!-- Introductory paragraph: what problem does this section address?
#      Keep it to 2-4 sentences. -->

# %% [markdown]
# **Definition ([Concept Name]).**  
# <!-- One formal definition in math or precise prose.
#      Use $...$ for inline math, $$...$$ for display math.
#      Follow the pattern used in S02–S06. -->
#
# A *[concept]* is a [formal description]:
#
# $$
# \mathcal{X} = \{ x \in V \mid P(x) \}
# $$
#
# where $V$ is [set description] and $P$ is [predicate description].

# %% [markdown]
# <div style="padding: 1rem; border-left: 6px solid #4D96FF; background: #F3F8FF; border-radius: 10px;">
#
# **Key idea** — [one-sentence summary of the core insight of this section.]  
# This matters because [why it connects to the broader workflow].
#
# </div>

# %%
# ── [Descriptive comment for what this cell demonstrates] ────────────────
from rdkit import Chem
from synedu.Utils import smiles_to_graph

example_smi = 'c1ccccc1O'  # phenol — replace with your molecule
mol = Chem.MolFromSmiles(example_smi)
G = smiles_to_graph(example_smi)

draw_molecular_graph(G, title='[Your molecule]', show_indices=True, label_mode='all')

# %% [markdown]
# **Q1 — [Question topic]**
#
# [One clear question that asks the learner to apply or extend the concept just introduced.]
#
# ---
#
# <details>
# <summary><b>Solution:</b></summary>
#
# ```python
# # Solution code here
# def your_function(arg):
#     ...
#     return result
# ```
#
# </details>

# %% [markdown]
# ### 1.1 [Sub-section title]
#
# <!-- Sub-section narrative: 1-2 sentences introducing the specific topic. -->

# %%
# ── [What this cell shows] ────────────────────────────────────────────────

# %% [markdown]
# ### 1.2 [Sub-section with figure]
#
# <!-- Use this pattern when you have a static figure compiled from TikZ/LaTeX
#      or generated externally. Store editable sources and PDFs in figure/.
#      Delete this cell if no static figure is needed. -->

# %% [markdown]
# <figure style="text-align:center; margin: 1.5rem 0;">
#   <object data="figure/your_figure.pdf" type="application/pdf" width="100%" height="360">
#     <a href="figure/your_figure.pdf">Open the figure PDF</a>
#   </object>
#   <figcaption style="font-size:0.88em; color:#555; margin-top:0.4rem;">
#     <b>Figure&nbsp;1.</b> [Caption: what the figure shows and why it matters.
#     Use HTML entities for math: &rarr; for &rarr;, &larr; for &larr;, <i>K</i> for italics.]
#   </figcaption>
# </figure>

# %% [markdown]
# <a id="2-section-two-title"></a>
#
# ## 2. Section Two Title
#
# <!-- What does this section build on top of Section 1? -->

# %% [markdown]
# **Definition ([Second Concept]).**  
#
# Given [context], the *[concept]* is defined as:
#
# $$
# f(x) = \sum_{i} w_i \, g(x_i)
# $$

# %%
# ── [What this cell computes] ─────────────────────────────────────────────

# %% [markdown]
# <!-- Narrative between code cells: explain what the output above means
#      and what the learner should notice. -->

# %%
# ── [Visualization / result display] ─────────────────────────────────────

# %% [markdown]
# **Q2 — [Question topic]**
#
# [Question.]
#
# ---
#
# <details>
# <summary><b>Solution:</b></summary>
#
# ```python
# # Solution
# ```
#
# </details>

# %% [markdown]
# <a id="3-section-three-title"></a>
#
# ## 3. Section Three Title *(optional — delete if not needed)*
#
# <!-- Typically: bringing everything together into one end-to-end example,
#      or introducing a more advanced / applied topic. -->

# %%
# ── End-to-end example / advanced application ─────────────────────────────

# %% [markdown]
# **Q3 — [Question topic]**
#
# [Question — typically asks learner to integrate concepts from all sections.]
#
# ---
#
# <details>
# <summary><b>Solution:</b></summary>
#
# ```python
# # Solution
# ```
#
# </details>

# %% [markdown]
# <a id="4-quiz"></a>
#
# ## 4. Quiz
#
# <!-- Short conceptual questions (no code required).
#      3-5 questions, each testing one learning outcome.
#      Use bullet points for multi-part questions. -->
#
# 1. **[Topic A]**  
#    [Question text.]
#
# 2. **[Topic B]**  
#    [Question text.]
#
# 3. **[Topic C]**  
#    [Question text.]
#
# 4. **Challenge**  
#    [Open-ended or harder question that extends the core material.]

# %%
# ── Quiz helpers ──────────────────────────────────────────────────────────
# Provide any pre-built objects the learner needs to answer the quiz questions.
# Keep this cell minimal — it should scaffold, not solve.

# %% [markdown]
# <a id="5-discussion"></a>
#
# ## 5. Discussion
#
# <!-- 4-6 bullet points summarising the key takeaways.
#      Each bullet should be one concrete, memorable insight.
#      Mirror the learning outcomes, but phrased as conclusions. -->
#
# - **[Takeaway A]**: [One sentence explaining what the learner should remember.]
# - **[Takeaway B]**: [One sentence.]
# - **[Takeaway C]**: [One sentence.]
# - **[Takeaway D]**: [One sentence — ideally a bridge to the next talktorial.]
#
# ### What comes next
#
# **[SXX+1]** will extend this to [next topic] by [brief description of the next step].

# %% [raw] raw_mimetype="text/x-rst"
# .. _s00-references:
#
# .. rubric:: 6. References
#
# .. _ref-1:
#
# **[1]** RDKit documentation: https://www.rdkit.org/docs/
#
# .. _ref-2:
#
# **[2]** NetworkX documentation: https://networkx.org/documentation/stable/
#
# .. _ref-3:
#
# **[3]** Phan, T.-L. *et al.*
# SynKit: A graph-based framework for rule-based reaction modeling.
# *Journal of Chemical Information and Modeling* (2025).
# https://doi.org/10.1021/acs.jcim.5c02123
#
# .. _ref-4:
#
# **[4]** [Author(s)]. [Title].
# *[Journal]* **[Volume]** ([Year]).
# https://doi.org/...
#
# .. _ref-5:
#
# **[5]** [Author(s)]. [Title].
# *[Journal]* **[Volume]** ([Year]).
# https://doi.org/...
#
