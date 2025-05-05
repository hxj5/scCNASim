scCNASim
========

scCNASim - Haplotype-aware simulation of somatic CNAs from single-cell and spatial transcriptomics
--------------------------------------------------------------------------------------------------

scCNASim is a python package designed for simulation of allele-specific 
somatic copy number alterations (CNAs) from single-cell and spatial 
transcriptomics.
It mainly takes existing alignment file, phased SNPs, and a clonal CNA profile
as input, and outputs new alignments with designated signals of CNAs and 
clonal structure. 
The core idea involves processing haplotype-specific reads separately, 
including fitting and simulating haplotype-specific gene expression counts, 
followed by UMI (read) sampling.



News
----
Release notes are at `docs/release.rst <./docs/release.rst>`_.



Installation
------------
Currently, only Python 3.11 (compatible) and 3.7 (not compatible) were tested.
Therefore, we strongly recommend to install the package with Python >= 3.11.


Dependency
~~~~~~~~~~
* Python >= 3.11

.. code-block:: bash

   pip install -U git+https://github.com/hxj5/scCNASim


Potential Issues
~~~~~~~~~~~~~~~~
If you encounter an error
``"configure: error: liblzma development files not found"``
when installing scCNASim, it is actually an installation issue of pysam.

You can fix the error easily by installing pysam via conda, if you are
installing scCNASim in an conda env, i.e., run

.. code-block:: bash

   conda config --add channels bioconda
   conda config --add channels conda-forge
   conda install pysam

and then re-install scCNASim.
See `Issue 3 <https://github.com/hxj5/scCNASim/issues/3>`_ for details.



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
