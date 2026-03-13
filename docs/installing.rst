Installation
============

.. raw:: html

   <div class="synedu-home-hero synedu-home-hero--compact synedu-home-hero--install">
     <div class="synedu-hero-main">
       <div class="synedu-home-kicker"><span class="synedu-dot"></span> Local setup</div>
       <div class="synedu-home-title">Run SynEdu end-to-end on your own machine.</div>
       <div class="synedu-home-subtitle">
         Install the package, open the talktorials, and start exploring reaction informatics
         workflows locally with a clean Python environment.
       </div>
     </div>
   </div>

.. raw:: html

   <div class="synedu-contact-highlights synedu-install-highlights">
     <span class="synedu-contact-pill">Python ≥ 3.11</span>
     <span class="synedu-contact-pill">RDKit via conda/mamba</span>
     <span class="synedu-contact-pill">JupyterLab ready</span>
     <span class="synedu-contact-pill">PyPI + editable install</span>
   </div>

.. raw:: html

   <div class="synedu-install-grid">
     <div class="synedu-install-card">
       <div class="synedu-install-card__icon">⚙️</div>
       <div class="synedu-install-card__title">Environment</div>
       <div class="synedu-install-card__desc">
         Start from a fresh Python environment with RDKit and JupyterLab installed.
       </div>
     </div>
     <div class="synedu-install-card">
       <div class="synedu-install-card__icon">📦</div>
       <div class="synedu-install-card__title">Package</div>
       <div class="synedu-install-card__desc">
         Install SynEdu from PyPI for normal use or in editable mode for development.
       </div>
     </div>
     <div class="synedu-install-card">
       <div class="synedu-install-card__icon">📘</div>
       <div class="synedu-install-card__title">Talktorials</div>
       <div class="synedu-install-card__desc">
         Launch the notebook-based teaching material directly in JupyterLab.
       </div>
     </div>
     <div class="synedu-install-card">
       <div class="synedu-install-card__icon">🚀</div>
       <div class="synedu-install-card__title">Quick start</div>
       <div class="synedu-install-card__desc">
         Follow a minimal setup path and begin exploring the local workflows immediately.
       </div>
     </div>
   </div>

Overview
--------

SynEdu is designed to run locally from start to finish, including the
interactive talktorials and supporting workflows.

A typical local setup includes a modern Python environment, an RDKit
installation, the ``synedu`` package itself, and JupyterLab for opening the
notebooks.

Requirements
------------

Before installing SynEdu, make sure you have:

- **Python 3.11 or newer**
- **RDKit** installed in the active environment
- **pip** available for Python package installation

.. note::

   RDKit is most reliably installed with **conda** or **mamba**, especially on
   Linux and macOS.

Recommended environment setup
-----------------------------

For the smoothest local experience, start from a clean environment:

.. code-block:: bash

   conda create -n synedu python=3.11
   conda activate synedu
   conda install -c conda-forge rdkit jupyterlab

You can then install SynEdu with ``pip`` inside that environment.

Install SynEdu
--------------

For a standard installation from PyPI:

.. code-block:: bash

   pip install synedu

For local development with an editable install:

.. code-block:: bash

   pip install -e ".[dev]"

.. note::

   In ``zsh``, extras should usually be quoted, for example ``".[dev]"``, to
   avoid shell glob expansion.

Launch the talktorials
----------------------

Once installed, open the SynEdu talktorial workspace:

.. code-block:: bash

   jupyter lab synedu

This launches the notebook-based learning material in JupyterLab.

Optional CLI workflow
---------------------

If you prefer the command-line interface:

.. code-block:: bash

   synedu start .

This creates a local workspace directory containing the talktorial notebooks and
associated files.

Quick start
-----------

.. raw:: html

   <div class="synedu-install-quickstart">
     <div class="synedu-install-quickstart__title">Minimal setup path</div>
     <div class="synedu-install-quickstart__desc">
       Create a fresh environment, install the package, and launch the talktorials in a few commands.
     </div>
   </div>

.. code-block:: bash

   conda create -n synedu python=3.11
   conda activate synedu
   conda install -c conda-forge rdkit jupyterlab
   pip install synedu
   jupyter lab synedu

Troubleshooting
---------------

**RDKit import issues**
   Make sure RDKit is installed in the same environment where ``synedu`` is installed.

**Command not found: jupyter**
   Install JupyterLab in the active environment:

   .. code-block:: bash

      pip install jupyterlab

**Command not found: synedu**
   Ensure the package installed successfully and that the correct environment is activated.

Next step
---------

After installation, continue to the talktorials to explore the main SynEdu
workflows interactively.