:orphan:

Installation
============

SynEdu is designed to run locally from start to finish. A normal setup needs a
clean Python environment, RDKit, the ``synedu`` package, and JupyterLab for the
notebook talktorials.

Setup checklist
---------------

.. rst-class:: synedu-info-table

.. list-table::
   :widths: 24 76
   :header-rows: 1

   * - Item
     - Notes
   * - Python
     - Use Python 3.11 or newer.
   * - RDKit
     - Install with conda or mamba for the most reliable chemistry stack.
   * - JupyterLab
     - Needed for opening and running the talktorial notebooks.
   * - SynEdu
     - Install from PyPI for normal use or editable mode for development.

Recommended environment setup
-----------------------------

Start from a fresh environment:

.. code-block:: bash

   conda create -n synedu python=3.11
   conda activate synedu
   conda install -c conda-forge rdkit jupyterlab

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

Optional CLI workflow
---------------------

If you prefer the command-line interface:

.. code-block:: bash

   synedu start .

This creates a local workspace directory containing the talktorial notebooks and
associated files.

Minimal setup path
------------------

.. code-block:: bash

   conda create -n synedu python=3.11
   conda activate synedu
   conda install -c conda-forge rdkit jupyterlab
   pip install synedu
   jupyter lab synedu

Troubleshooting
---------------

**RDKit import issues**
   Make sure RDKit is installed in the same environment where ``synedu`` is
   installed.

**Command not found: jupyter**
   Install JupyterLab in the active environment:

   .. code-block:: bash

      pip install jupyterlab

**Command not found: synedu**
   Ensure the package installed successfully and that the correct environment is
   activated.

Next step
---------

After installation, continue to the :doc:`talktorials/index` to explore the main
SynEdu workflows interactively.
