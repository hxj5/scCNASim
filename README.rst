scCNVSim - Simulation of somatic CNVs in single cells
=====================================================

scCNVSim is a python package for simulation of somatic copy number variations
(CNVs) in single cells.
It mainly takes existing alignment file(s) and a clonal CNV profile as input,
and outputs new alignments with designated signals of CNVs and clonal 
structure.


News
----
Release notes are at `docs/release.rst <./docs/release.rst>`_.


Installation
------------
.. code-block:: bash

   pip install -U git+https://github.com/hxj5/scCNVSim


Manual
------
The full manual is at `docs/manual.rst <./docs/manual.rst>`_.


FAQ and feedback
----------------
For troubleshooting, please have a look of `docs/FAQ.rst <./docs/FAQ.rst>`_,
and we welcome reporting any issue_ for bugs, questions and 
new feature requests.


Acknowledgement
---------------
The simulator has a precursor named scCNASimulator_, which has been used in
XClone_ to demonstrate its robustness to detect allele-specific CNVs.

scCNASimulator implements a naive strategy for CNV simulation, which 
multiplies the UMI/read counts directly by copy number fold to generate the
new counts of CNV features, whereas this new simulator models the counts
with certain probability distribution and encodes the CN fold in the updated
distribution parameters before generating new simulated counts.


.. _issue: https://github.com/hxj5/scCNVSim/issues
.. _scCNASimulator: https://github.com/hxj5/scCNASimulator
.. _XClone: https://github.com/single-cell-genetics/XClone
