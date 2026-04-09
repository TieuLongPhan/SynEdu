Installation
============

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
