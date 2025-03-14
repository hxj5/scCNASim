
Manual
======

.. contents:: Contents
   :depth: 2
   :local:


Quick Usage
-----------
First, please look at `Input`_ to prepare the input data.

Then call the ``main_wrapper()`` function to run the simulation pipeline.

An example is:

.. code-block:: python

   main_wrapper(
       sam_fn = "{sample}.bam",
       cell_anno_fn = "cell_anno.tsv", 
       feature_fn = "hg38.features.tsv",
       phased_snp_fn = "phased.snp.vcf.gz",
       clone_meta_fn = "clone_anno.tsv",
       cna_profile_fn = "cna_profile.tsv", 
       refseq_fn = "hg38.fa",
       out_dir = "./simu_result",
       cell_tag = "CB", 
       umi_tag = "UB",
       ncores = 10, 
       verbose = False
   )

The full parameters can be found at ``main_wrapper()`` in the ``API`` page.

You may also run each step (module) explicitly by calling corresponding 
wrapper functions (see ``Tutorial`` page).

See `Implementation`_ for details of the four modules.


Input
-----
The inputs to the simulator include:

* Alignment file of seed data (BAM file).
* Cell annotations of seed data (TSV file).
* Feature annotations of seed data (TSV file).
* Phased SNPs of seed data (TSV or VCF file).
* Reference genome sequence of seed data (FASTA file).
* Clone annotations of simulated data (TSV file).
* Clonal CNA profiles of simulated data  (TSV file).


Alignment file of seed data (BAM file)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The aligned reads stored in either one single BAM file (from droplet-based 
sequencing platform) or a list of BAM files (from well-based sequencing 
platform).


Cell annotations of seed data (TSV file)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The cell annotation stored in a header-free TSV file.
Its first two columns are ``cell`` and ``cell_type``, where

cell : str
    Cell barcodes (droplet-based data) or sample ID (well-based data).

cell_type : str
    Cell type.

An example is as follows:

.. code-block::

   AAAGATGGTCCGAAGA-1    immune
   AACCATGTCTCGTATT-1    immune
   AACGTTGTCTCTTGAT-1    epithelial
   AACTCAGAGCCTATGT-1    immune
   AAGACCTAGATGTAAC-1    epithelial
   AAGCCGCTCCTCAATT-1    epithelial


Feature annotations of seed data (TSV file)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The feature annotation stored in a header-free TSV file.
Its first five columns are ``chrom``, ``start``, ``end``, ``feature``,
and ``strand``, where

chrom : str
    Chromosome name of the feature.

start : int
    Start genomic position of the feature, 1-based and inclusive.

end : int
    End genomic position of the feature, 1-based and inclusive.

feature : str
    Feature name.
    
strand : str
    DNA strand orientation of the feature, "+" (positive) or "-" (negative).

An example is as follows:

.. code-block::

   chr1       29554   31109   MIR1302-2HG     +
   chr1       34554   36081   FAM138A -
   chr1       65419   71585   OR4F5   +
   chr2       38814   46870   FAM110C -
   chr2       197569  202605  AC079779.1      +
   chr3       23757   24501   LINC01986       +


Phased SNPs of seed data (TSV or VCF file)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The phased SNPs stored in either a TSV file or a VCF file.

Phased SNPs in TSV format
+++++++++++++++++++++++++
If it is in a TSV file, it should be header-free and its first 6 columns
should be ``chrom``, ``pos``, ``ref``, ``alt``, ``ref_hap``, and 
``alt_hap``, where

chrom : str
    The chromosome name of the SNP.

pos : int
    The genomic position of the SNP, 1-based.

ref : str
    The reference (REF) allele of the SNP, one of ``{'A', 'C', 'G', 'T'}``.

alt : str
    The alternative (ALT) allele of the SNP, one of ``{'A', 'C', 'G', 'T'}``.

ref_hap : int
    The haplotype index of ``ref``, one of ``{0, 1}``.

alt_hap : int
    The haplotype index of ``alt``, one of ``{1, 0}``.
 
An example is as follows:

.. code-block::

   chr1    986336   C       A   0   1
   chr1    1007256  G       A   1   0
   chr1    1163041  C       T   1   0
   chr2    264895   G       C   0   1
   chr2    277003   A       G   0   1
   chr2    3388055  C       T   1   0


Phased SNPs in VCF format
+++++++++++++++++++++++++
If it is in VCF format, the file should contain the ``GT`` in its
``FORMAT`` field (i.e., the 9th column).
The corresponding phased genotype could be delimited by either ``'/'`` or
``'|'``, e.g., "0/1", or "0|1".

.. note::
   * As reference phasing, e.g., with Eagle2, is not perfect, one UMI may 
     cover two SNPs with conflicting haplotype states.
   * Reference phasing tends to have higher rate in longer distance.
     Therefore, further local phasing (e.g., in gene level) and global phasing
     (e.g., in bin level) could be used to reduce error rate, e.g., with the
     3-step phasing used by CHISEL_ in scDNA-seq data and XClone_ in scRNA-seq
     data.
     

Reference genome sequence of seed data (FASTA file)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The sequence of reference genome, e.g., the human genome version hg38, 
should be stored in a FASTA file.
Its version should match the one used for generating the alignment (BAM)
file of seed data.


Clone annotations of simulated data (TSV file)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
Clone annotation stored in a header-free TSV file.
Its first 3 columns should be ``clone``, ``source_cell_type``, and ``n_cell``,
where

clone : str
    The clone ID.

source_cell_type : str
    The source cell type of ``clone``.

n_cell : int
    Number of cells in the ``clone``.
    If negative, then it will be set as the number of cells in 
    ``source_cell_type``.
 
An example is as follows:

.. code-block::

   clone1_normal    immune  -1
   clone2_normal    epithelial  -1
   clone3_cancer    epithelial  -1
   clone4_cancer    epithelial  -1
   clone5_cancer    epithelial  -1

.. note::
   The simulator is designed for diploid genome.
   Generally, it is recommended to use normal cells as ``source_cell_type``
   for simulation of somatic CNAs.


Clonal CNA profiles of simulated data (TSV file)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The clonal CNA profile stored in a header-free TSV file.
Its first 6 columns should be ``chrom``, ``start``, ``end``,
``clone``, ``cn_ale0``, and ``cn_ale1``, where

chrom : str
    The chromosome name of the CNA region.

start : int
    The start genomic position of the CNA region, 1-based and inclusive.

end : int or "Inf"
    The end genomic position of the CNA region, 1-based and inclusive.
    To specify the end of the whole chromosome, you can use either the actual
    genomic position or simply ``Inf``.

clone : str
    The clone ID.

cn_ale0 : int
    The copy number of the first allele (haplotype).

cn_ale1 : int
    The copy number of the second allele (haplotype).
 
One clone-specific CNA per line.
An example is as follows:

.. code-block::

   chr8 1   Inf clone3_cancer   1   2
   chr6 1   Inf clone4_cancer   0   1
   chr8 1   Inf clone4_cancer   1   2
   chr6 1   Inf clone5_cancer   1   0
   chr8 1   Inf clone5_cancer   1   2
   chr11    1   Inf clone5_cancer   2   0


**Support all three major CNA types**

By specifying different values for ``cn_ale0`` and ``cn_ale1``, you may
specify various CNA types, including copy gain (e.g., setting ``1, 2``), 
copy loss (e.g., setting ``0, 1``), LOH (e.g., setting ``2, 0``).

**Support allele-specific CNA**

This format fully supports allele-specific CNAs.
For instance, to simulate the scenario that two subclones have copy loss in
the same region while on distinct alleles, setting ``cn_ale0, cn_ale1``
to ``0, 1`` and ``1, 0`` in two subclones, respectively, as the example of
copy loss in chr6.

**Support whole genome duplication (WGD)**

It also supports whole genome duplication (WGD), e.g., by setting 
``cn_ale0, cn_ale1`` of all chromosomes to ``2, 2``.
Generally, detecting WGD from scRNA-seq data is challenging, as it is hard
to distinguish WGD from high library size.
One scenario eaiser to detect WGD is that a balanced copy loss occurred 
after WGD, e.g., setting ``cn_ale0, cn_ale1`` of chr3 to ``1, 1``, while
``2, 2`` for all other chromosomes.
In this case, chr3 may have signals of balanced BAF while copy-loss RDR,
which should not happen on normal diploid genome.

**Notes**

* All CNA clones ``clone`` in this file must be in the clone annotation file.
* Only the CNA clones are needed to be listed in this file. Do not list normal
  clones in this file.


Output
------
The final output is available at folder ``{out_dir}/4_rs``.
It contains

* Simulated alignment file (BAM file).
* Simulated cell annotation (TSV file).


Simulated alignment file (BAM file)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The simulated reads stored in either one single BAM file (from droplet-based
sequencing platform) or a list of BAM files (from well-based sequencing 
platform).
The BAM file(s) are available at folder ``{out_dir}/4_rs/bam``.


Simulated cell annotation (TSV file)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The simulated cell annotation stored in a header-free TSV file, located at
``{out_dir}/4_rs/rs.cell_anno.tsv``.
It has two columns ``cell`` and ``clone``, where

cell : str
    The cell barcode (droplet-based data) or sample ID (well-based).

clone : str
    The clone ID.

Note that there is a one-column TSV file storing ``cell`` (cell barcodes or
sample ID) only, located at ``{out_dir}/4_rs/rs.samples.tsv``.


Implementation
--------------
The simulator outputs simulated haplotype-aware alignments for clonal single 
cells based on user-specified CNA profiles, by training on input BAM files.

It mainly includes four modules:

#. ``pp``: preprocessing.
#. ``afc``: allele-specific feature counting.
#. ``cs``: count simulation.
#. ``rs``: read simulation.


The ``pp`` (preprocessing) module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This module is implemented in the function ``pp.main.pp_wrapper()``.
The results of this module are stored in the folder ``{out_dir}/1_pp``.

It preprocesses the inputs, including:

* Check and merge overlapping features in the input feature annotation file.
* Check and merge overlapping CNA profiles in the input clonal CNA profile 
  file.


The ``afc`` (allele-specific feature counting) module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This module extracts and counts allele-specific UMIs/reads in single cells.

It is implemented in the function ``afc.main.afc_wrapper()``.
The results of this module are stored in the folder ``{out_dir}/2_afc``.

To speedup, features are splitted into batches for multi-processing.
In one feature, the haplotype state of each UMI/read is inferred by
integrating haplotype information from all SNPs covered by the UMI/read.

Haplotype state of UMI/read pair in SNP scale
+++++++++++++++++++++++++++++++++++++++++++++
Firstly, SNP pileup is performmed to fetch the reads covering the SNPs that
are located within the feature.
For each SNP, the haplotype state of its fetched UMI/read pair is inferred by
comparing the fetched SNP allele to the phased ones, e.g., 
if the fetched allele is 'A' in one UMI/read pair, and the phased REF and ALT
alleles are 'A' and 'C', respectively, then the UMI/read pair would be 
inferred as from the REF haplotype given the SNP.

All available haplotype state of UMI/read pair in SNP scale are listed below:

.. list-table:: SNP-scale haplotype state
   :align: center
   :widths: 15 30 55
   :header-rows: 1

   * - Index
     - String
     - Brief Description
   * - 0
     - ref (reference)
     - The fetched SNP allele is on the reference haplotype.
   * - 1
     - alt (alternative)
     - The fetched SNP allele is on the alternative haplotype.
   * - -1
     - oth (others)
     - Some allele is fetched but is on neither the reference nor 
       alternative haplotype.
   * - -2 
     - unknown
     - No allele is fetched (the value is None).


Haplotype state of UMI/read pair in feature scale
+++++++++++++++++++++++++++++++++++++++++++++++++
Secondly, the final haplotype state of one UMI/read pair (in feature scale)
is inferred by integrating information from all SNPs covered by the 
UMI/read pair.

All available haplotype state of UMI/read pair in feature scale are listed
below:

.. list-table:: Feature-scale haplotype state
   :align: center
   :widths: 15 30 55
   :header-rows: 1

   * - Index
     - String
     - Brief Description
   * - 0
     - ref (reference)
     - Reference haplotype has supporting SNPs but alternative haplotype does
       not.
   * - 1
     - alt (alternative)
     - Alternative haplotype has supporting SNPs but reference haplotype does
       not.
   * - 2 
     - both
     - Both reference and alternative haplotypes have supporting SNPs.
   * - -1
     - oth (others)
     - Neither reference nor alternative haplotype has supporting SNPs, but
       other alleles (bases) in SNP level are fetched.
   * - -2
     - unknown
     - The UMI/read pair is fetched by some SNPs, but no any alleles (bases)
       are fetched.


Final haplotype state of UMI/read pair
++++++++++++++++++++++++++++++++++++++
Lastly, all reads of the feature will be iterated, including both fetched 
reads of given SNPs and other reads covering no SNPs.
The haplotype state of each iterated UMI/read pair is determined based on 
previous step.

All final haplotype state of UMI/read pair in feature scale are listed
below:

.. list-table:: Final haplotype state
   :widths: 15 25 60
   :header-rows: 1

   * - Index
     - String
     - Brief Description
   * - 0
     - A (Haplotype-A; ref)
     - Haplotype A has supporting SNPs but haplotype B does not.
   * - 1
     - B (Haplotype-B; alt)
     - Haplotype B has supporting SNPs but haplotype A does not.
   * - 2 
     - D (Duplicate; both)
     - Both haplotype A and B have supporting SNPs.
   * - -1
     - O (Others)
     - Neither haplotype A nor B has supporting SNPs, but other alleles 
       (bases) in SNP level are fetched.
   * - -2
     - U (Unknown)
     - The UMI/read pair is fetched by some SNPs, but no any alleles (bases)
       are fetched.
   * - -3
     - U (Unknown)
     - The UMI/read pair is not fetched by any SNPs.

The output allele-specific *feature x cell* count matrices are at folder 
``{out_dir}/2_afc/counts``.

Additionally, all the count matrices are also saved into one anndata ".h5ad"
file, ``{out_dir}/2_afc/afc.counts.cell_anno.h5ad``, which will be used by 
downstream ``cs`` module.


The allele-specific UMIs
++++++++++++++++++++++++
Although there are in total 6 haplotype states available, only 3 of them,
including "A", "B", "U" (merged from haplotype index -2 and -3), will be used
for downstream ``rs`` module.
Reads of the other 3 haplotype states are excluded from downstream analysis.

The extracted allele-specific UMIs (CUMIs) are stored in each corresponding
header-free TSV file, i.e., 
``{out_dir}/2_afc/{batch}/{feature}/{feature}.{haplotype}.aln.afc.tsv``.
These files contain 2 columns ``cell``, ``UMI``, where

cell : str
    The cell barcode (droplet-based data) or sample ID (well-based data).

UMI : str
    The UMI barcode (droplet-based data) or query name (well-based data).

An example is as follows:

.. code-block::

   ACCCACTCAGTTTACG-1      AGCAGATCAG
   ACGATGTTCACCTCGT-1      AATTTACGCA
   AGCGTATAGCCGCCTA-1      GGTCTCAGCT


The ``cs`` (count simulation) module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This module simulates new allele-specific *cell x feature* count matrices
based on existing matrices.

It is implemented in the function ``cs.main.cs_wrapper()``.
The results of this module are stored in the folder ``{out_dir}/3_cs``.

This module processes the count matrices of haplotypes "A", "B", "U",
separately, mainly following three steps:

#. Fit feature-specific counts with a specific distribution.
#. Update the fitted feature-specific parameters based on the CNA profile.
#. Generate new feature-specific counts based on the updated parameters.


Fit input feature-specific counts
+++++++++++++++++++++++++++++++++
For each haplotype-specific *cell x feature* count matrix, features are 
processed separately within each cell type using multi-processing.
For feature counts in a specific cell type, the counts are modelled with one
of the four distribution: "poi" (Poisson), "nb" (Negative Binomial), "zip" 
(Zero-Inflated Poisson), and "zinb" (Zero-Inflated Negative Binomial), either
speficied by users or using a data-driven auto-selected strategy.


Update parameters based on CNA profile
++++++++++++++++++++++++++++++++++++++
The fitted feature-specific parameters are updated, multiplying a coefficient
of copy number fold based on the CNA profile.
For example, if one feature overlaps a copy loss region (e.g., either 1,0 or
0,1) in certain CNA clone, then the CN fold of this feature in this clone
would be less than 1.0 (e.g., 0.5).
If the feature overlaps a copy gain region (e.g., 1,2 or 2,1), then the CN
fold is larger than 1.0 (e.g., 1.5).
If the feature overlaps a LOH region (e.g., either 0,2 or 2,0), then the CN
fold is 1.0.


Generate new feature-specific counts
++++++++++++++++++++++++++++++++++++
The updated parameters are used for generation of new haplotype-specific
*cell x feature* count matrices.
All haplotype-specific count matrices will then be merged to construct a 
anndata ".h5ad" file, located at ``{out_dir}/3_cs/cs.counts.h5ad``.
Additionally, the parameters are also outputted, to a python pickle file
``{out_dir}/3_cs/cs.params.pickle``.


The ``rs`` (read simulation) module
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
This module simulates new reads for new clonal single cells by sampling reads
from the input BAM file(s) according to the simulated counts.

It is implemented in the function ``rs.main.rs_wrapper()``.
The results of this module are stored in the folder ``{out_dir}/4_rs``.

Specifically, it includes following steps:

#. Sample *cell x feature* CUMIs based on simulated counts.
#. Extract output reads according to the sampled CUMIs.


Sample *cell x feature* CUMIs
+++++++++++++++++++++++++++++
To generate new reads, CUMIs are first sampled instead of reads themselves to
improve computational efficiency.
CUMIs are representative strings of reads sharing the same cell and UMI
barcodes (droplet-based data) or sample ID and read query name (well-based
data).
The haplotype-specific CUMIs have been extracted for each feature and
stored in TSV files, e.g.,
``{out_dir}/2_afc/{batch}/{feature}/{feature}.{haplotype}.aln.afc.tsv``.
For each feature, CUMIs from all input cells are sampled with replacement in
a pseudo-bulk manner, based on the simulated counts, to generate a list of 
sampled CUMIs (and corresponding new CUMIs) for each new single cell.


Extract output reads
++++++++++++++++++++
The simulator extracts the reads from the input BAM file(s) based on the
sampled CUMIs by matching the CUMIs of reads with all sampled CUMIs.
To speedup, reads of different chromosomes are iterated in parallel using
multi-processing.
For one iterated read, if its CUMI matches some sampled CUMIs, the read
will be extracted and modified before being outputted:

* assigned with corresponding new CUMI(s), i.e., new cell and UMI tags
  (droplet-based data) or new sample ID and query name (well-based data).
* its query name is added with a unique suffix.
* other information of the read, such as FLAG, CIGAR, SEQ, and QUAL, is not 
  changed.

Note that

* one source read could be sampled and outputted multiple times
  (e.g., sampled multiple times by one new cell or sampled by more than one 
  new cells).
  The combination of new CUMI and query name ensures that each read in the 
  output BAM file(s) is unique. 
* sampling is performmed in units of UMI group (droplet-based data) or
  read pair (well-based data).
  The output reads will be assigned the same new CUMI, i.e., in the same
  output UMI group or read pair, each time their source group of reads are 
  sampled.

The output reads of all chromosomes will be merged into new BAM file(s) and
stored in folder ``{out_dir}/4_rs/bam``.


.. _CHISEL: https://www.nature.com/articles/s41587-020-0661-6
.. _XClone: https://www.biorxiv.org/content/10.1101/2023.04.03.535352v2
