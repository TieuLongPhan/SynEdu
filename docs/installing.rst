Run locally
==========

SynEdu is designed to be runnable end-to-end on your machine.

Requirements
------------

- Python >= 3.11
- A working RDKit installation (via conda/mamba is recommended)

Install
-------

If you installed SynEdu from PyPI:

.. code-block:: bash

   pip install synedu

For development (editable) installs:

.. code-block:: bash

   pip install -e ".[dev]"

.. note::

   In ``zsh``, you typically must quote extras (e.g. ``".[dev]"``) to avoid
   glob expansion.

Open the talktorials
--------------------

From the repository root:

.. code-block:: bash

   jupyter lab synedu

If you use the CLI (optional):

.. code-block:: bash

   synedu start .

This creates a workspace directory with the talktorial notebooks.
