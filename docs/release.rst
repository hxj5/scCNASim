
..
   History
   =======


Release v0.1.2 (11/02/2025)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
This version mainly aims to reduce the variation in the simulated BAF signals
of normal features/regions.

* afc: set default min_count=20, min_maf=0.1.
  It may filter some input phased SNPs whose expression levels in the seed
  data are low.
  Motivation: if one gene contains mainly lowly-expressed SNPs, then its
  haplotype-specific counts (Hap-A and Hap-B) will be small, and its simulated
  Hap-A and Hap-B counts are probably also small, hence AF may be biased
  towards 0 or 1.
* cs: use Poisson distribution when fitting NB failed.
  Previously, empirical parameters of NB were used when fitting NB failed.
  We expect Poisson to produce lower variation level in simulated counts, 
  compared to NB, especially for lowly-expressed features.
* pp: set "quantile2" as default option of ``merge_features_how``.
  Both "quantile2" and "quantile2_union" strategies can remove features that
  overlap large number of other features, while "quantile2" seem to produce
  stronger CNA signals.
* main: logging APP and VERSION.


Release v0.1.1 (03/02/2025)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
* pp: add ``merge_features_how`` - How to merge overlapping features.
* Support both INT and FLOAT as value of ``--minINCLUDE``.
  If float between (0, 1), it is the minimum fraction of included length.
* Set default value of ``--minINCLUDE`` or ``min_include`` as 0.9.
* docs: add TODO.


Release v0.1.0 (06/12/2024)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Add ``--minINCLUDE`` option for read filtering.

* ``--minINCLUDE`` is the minimum length of included part within specific
  feature. 
* For example, if the genomic range of a feature is chr1:1000-3000, and one
  fetched read (100bp) aligned to two locus, chr1:601-660 (60bp) and 
  chr1:3801-3840 (40bp), then no any part of the read is actually included 
  within the feature, hence it will be filtered by ``--minINCLUDE=30``, 
  whereas older versions of scCNASim may keep the read.
  Note, when features are processed independently, one read filtered by
  --minINCLUDE in one feature may still be fetched and counted by other 
  features.
* Previously, there is noise present in inferCNV heatmap that both signals 
  of duplication and deletion present in a strip of genes, even in the
  reference cells.
  By using ``--minINCLUDE`` (default 30), the noise is largely removed.
  
Others

* rs: do not output sampled reads of multi-feature UMIs for non-overlapping
  features.
  If one multi-read UMI is sampled by specific feature (in rs module), and
  some of its reads are not included within the feature (``--minINCLUDE``),
  then those reads will not be outputted to BAM for this feature.
  Without this step, there will be inflation of UMI counts in rs BAM, compared
  to the simulated counts in cs module, considering the non-included reads may
  be counted by other features.
* rs: output sampled UMIs aligned to distinct alleles in different features.
  Assume there is a multi-feature UMI (due to error in UMI collapse?) 
  aligned to distinct alleles in different features, e.g., Hap-B in one 
  feature and Hap-U in another feature.
  If the UMI is sampled by both features, then the UMI is outputted for both
  features, while mimicking the real scRNA-seq BAM (error in UMI collapse?).
  Previously, this UMI is only outputted once for one (first iterated) 
  feature, which may result in the decrease of UMI counts in rs BAM, compared
  to the simulated counts in cs module.
* pp: filter features by chromosomes.
  Filter features whose chromosomes are not in the input chrom list.
* convert column chrom astype str in anndata.
  Previously, the chrom column will be of int dtype if all chromosome names are
  numeric strings, e.g., "1", "2", etc.
* init setting random seed.
  Currently the whole simulation results are not reproducible with a seed,
  possibly due to the parallel computing.
* cs: also output the counts into sparse matrices, in addition to the
  ``h5ad`` file.
* pp and afc: rename ``utils`` to ``io``.


Bug fix:

* utils: fix bug in ``xbarcode.str2int()``.


Release v0.0.2 (12/10/2024)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
* rename CNV to CNA.
* allow input empty CNA profile file.
* require Python>=3.11.
* fix typos.


Release v0.0.1 (17/09/2024)
~~~~~~~~~~~~~~~~~~~~~~~~~~~
Implement a pipeline wrapping four modules:

#. ``pp``: preprocessing.
#. ``afc``: allele-specific feature counting.
#. ``cs``: count simulation.
#. ``rs``: read simulation.
