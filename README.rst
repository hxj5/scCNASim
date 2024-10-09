scCNASim - Simulation of somatic CNAs in single cells
=====================================================

scCNASim is a python package for simulation of somatic copy number alterations
(CNAs) in single cells.
It mainly takes existing alignment file(s) and a clonal CNA profile as input,
and outputs new alignments with designated signals of CNAs and clonal 
structure.


News
----
Release notes are at `docs/release.rst <./docs/release.rst>`_.


Installation
------------
.. code-block:: bash

   pip install -U git+https://github.com/hxj5/scCNASim


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
XClone_ to demonstrate its robustness to detect allele-specific CNAs.

scCNASimulator implements a naive strategy for CNA simulation, which 
multiplies the UMI/read counts directly by copy number fold to generate the
new counts of CNA features, whereas this new simulator models the counts
with certain probability distribution and encodes the CN fold in the updated
distribution parameters before generating new simulated counts.


.. _issue: https://github.com/hxj5/scCNASim/issues
.. _scCNASimulator: https://github.com/hxj5/scCNASimulator
.. _XClone: https://github.com/single-cell-genetics/XClone
