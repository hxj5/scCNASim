# config.py - global configuration.


import sys
from .afc.config import DefaultConfig as AFC_DefConf


class Config:
    """Configuration.

    Attributes
    ----------
    sam_fn : str or None
        Comma separated indexed BAM file.
        Note that one and only one of `sam_fn` and `sam_list_fn` should be
        specified.
    cell_anno_fn : str
        The cell annotation file. 
        It is header-free and its first two columns are:
        - "cell" (str): cell barcodes.
        - "cell_type" (str): cell type.
    feature_fn : str
        A TSV file listing target features. 
        It is header-free and its first 4 columns shoud be: 
        - "chrom" (str): chromosome name of the feature.
        - "start" (int): start genomic position of the feature, 1-based
          and inclusive.
        - "end" (int): end genomic position of the feature, 1-based and
          inclusive.
        - "feature" (str): feature name.
    snp_fn : str
        A TSV or VCF file listing phased SNPs.
        If TSV, it is a header-free file containing SNP annotations, whose
        first six columns should be:
        - "chrom" (str): chromosome name of the SNP.
        - "pos" (int): genomic position of the SNP, 1-based.
        - "ref" (str): the reference allele of the SNP.
        - "alt" (str): the alternative allele of the SNP.
        - "ref_hap" (int): the haplotype index of `ref`, one of {0, 1}.
        - "alt_hap" (int): the haplotype index of `alt`, one of {1, 0}.
        If VCF, it should contain "GT" in its "FORMAT" field.
    clone_meta_fn : str
        A TSV file listing clonal meta information.
        It is header-free and its first 3 columns are:
        - "clone" (str): clone ID.
        - "source_cell_type" (str): the reference cell type of `clone`.
        - "n_cell" (int): number of cells in the `clone`. If negative, 
          then it will be set as the number of cells in `source_cell_type`.
    cna_profile_fn : str
        A TSV file listing clonal CNA profiles. 
        It is header-free and its first 7 columns are:
        - "chrom" (str): chromosome name of the CNA region.
        - "start" (int): start genomic position of the CNA region, 1-based
          and inclusive.
        - "end" (int): end genomic position of the CNA region, 1-based and
          inclusive.
        - "region" (str): ID of the CNA region.
        - "clone" (str): clone ID.
        - "cn_ale0" (int): copy number of the first allele.
        - "cn_ale1" (int): copy number of the second allele.
    refseq_fn : str
        A FASTA file storing reference genome sequence.
    out_dir : str
        The output folder.
    sam_list_fn : str or None, default None
        A file listing indexed BAM files, each per line.
    sample_ids : str or None, default None
        Comma separated sample IDs.
        It should be specified for well-based or bulk data.
        When `barcode_fn` is not specified, the default value will be
        "SampleX", where "X" is the 0-based index of the BAM file(s).
        Note that `sample_ids` and `sample_id_fn` should not be specified
        at the same time.
    sample_id_fn : str or None, default None
        A file listing sample IDs, each per line.
    size_factor : str or None, default "libsize"
        The type of size factor.
        Currently, only support "libsize" (library size).
        Set to `None` if do not use size factors for model fitting.
    marginal : {"auto", "poisson", "nb", "zinb"}
        Type of marginal distribution.
        One of
        - "auto" (auto select).
        - "poisson" (Poisson).
        - "nb" (Negative Binomial).
        - "zinb" (Zero-Inflated Negative Binomial).
    kwargs_fit_sf : dict
        The additional kwargs passed to function 
        :func:`~marginal.fit_libsize_wrapper` for fitting size factors.
        The available arguments are:
        - dist : {"normal", "t"}
            Type of distribution.
    kwargs_fit_rd : dict
        The additional kwargs passed to function 
        :func:`~marginal.fit_RD_wrapper` for fitting read depth.
        The available arguments are:
        - min_nonzero_num : int, default 3
            The minimum number of cells that have non-zeros for one feature.
            If smaller than the cutoff, then the feature will not be fitted
            (i.e., its mean will be directly treated as 0).
        - max_iter : int, default 1000
            Number of maximum iterations in model fitting.
        - pval_cutoff : float, default 0.05
            The p-value cutoff for model selection with GLR test.
    chroms : str, default "1,2,...22"
        Comma separated chromosome names.
        Reads in other chromosomes will not be used for sampling and hence
        will not be present in the output BAM file(s).
    cell_tag : str or None, default "CB"
        Tag for cell barcodes, set to None when using sample IDs.
    umi_tag : str or None, default "UB"
        Tag for UMI, set to None when reads only.
    umi_len : int, default 10
        Length of output UMI barcode.
    ncores : int, default 1
        Number of cores.
    verbose : bool, default False
        Whether to show detailed logging information.
    min_mapq : int, default 20
        Minimum MAPQ for read filtering.
    min_len : int, default 30
        Minimum mapped length for read filtering.
    min_include : int, default 30
        Minimum length of included part within specific feature.
    incl_flag : int, default 0
        Required flags: skip reads with all mask bits unset.
    excl_flag : int, default -1
        Filter flags: skip reads with any mask bits set.
        Value -1 means setting it to 772 when using UMI, or 1796 otherwise.
    no_orphan : bool, default True
        If `False`, do not skip anomalous read pairs.
    """
    def __init__(self):
        self.afc_def_conf = AFC_DefConf()

        # input and output files.
        self.sam_fn = None
        self.cell_anno_fn = None
        self.feature_fn = None
        self.snp_fn = None
        self.clone_meta_fn = None
        self.cna_profile_fn = None
        self.refseq_fn = None
        self.out_dir = None
        self.sam_list_fn = None
        self.sample_ids = None
        self.sample_id_fn = None

        # count simulation.
        self.size_factor = "libsize"
        self.marginal = "auto"
        self.kwargs_fit_sf = dict()
        self.kwargs_fit_rd = dict()

        # optional arguments.
        self.chroms = ",".join([str(c) for c in range(1, 23)])
        self.cell_tag = self.afc_def_conf.CELL_TAG
        self.umi_tag = self.afc_def_conf.UMI_TAG
        self.umi_len = 10
        self.ncores = self.afc_def_conf.NCORES
        self.verbose = False

        # read filtering.
        self.min_mapq = self.afc_def_conf.MIN_MAPQ
        self.min_len = self.afc_def_conf.MIN_LEN
        self.min_include = self.afc_def_conf.MIN_INCLUDE
        self.incl_flag = self.afc_def_conf.INCL_FLAG
        self.excl_flag = -1
        self.no_orphan = self.afc_def_conf.NO_ORPHAN

    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stdout

        s =  "%s\n" % prefix
        s += "%ssam_file = %s\n" % (prefix, self.sam_fn)
        s += "%scell_anno_file = %s\n" % (prefix, self.cell_anno_fn)
        s += "%sfeature_file = %s\n" % (prefix, self.feature_fn)
        s += "%sphased_snp_file = %s\n" % (prefix, self.snp_fn)
        s += "%sclone_meta_file = %s\n" % (prefix, self.clone_meta_fn)
        s += "%scna_profile_file = %s\n" % (prefix, self.cna_profile_fn)
        s += "%srefseq_file = %s\n" % (prefix, self.refseq_fn)
        s += "%sout_dir = %s\n" % (prefix, self.out_dir)
        s += "%ssam_list_file = %s\n" % (prefix, self.sam_list_fn)
        s += "%ssample_ids = %s\n" % (prefix, self.sample_ids)
        s += "%ssample_id_file = %s\n" % (prefix, self.sample_id_fn)
        s += "%s\n" % prefix

        s += "%ssize_factor = %s\n" % (prefix, self.size_factor)
        s += "%smarginal = %s\n" % (prefix, self.marginal)
        s += "%skwargs_fit_sf = %s\n" % (prefix, str(self.kwargs_fit_sf))
        s += "%skwargs_fit_rd = %s\n" % (prefix, str(self.kwargs_fit_rd))
        s += "%s\n" % prefix

        s += "%schroms = %s\n" % (prefix, self.chroms)
        s += "%scell_tag = %s\n" % (prefix, self.cell_tag)
        s += "%sumi_tag = %s\n" % (prefix, self.umi_tag)
        s += "%sumi_len = %d\n" % (prefix, self.umi_len)
        s += "%snumber_of_cores = %d\n" % (prefix, self.ncores)
        s += "%sverbose = %s\n" % (prefix, self.verbose)
        s += "%s\n" % prefix

        s += "%smin_mapq = %d\n" % (prefix, self.min_mapq)
        s += "%smin_len = %d\n" % (prefix, self.min_len)
        s += "%smin_include = %d\n" % (prefix, self.min_include)
        s += "%sinclude_flag = %d\n" % (prefix, self.incl_flag)
        s += "%sexclude_flag = %d\n" % (prefix, self.excl_flag)
        s += "%sno_orphan = %s\n" % (prefix, self.no_orphan)
        s += "%s\n" % prefix

        fp.write(s)

    def use_barcodes(self):
        return self.cell_tag is not None

    def use_umi(self):
        return self.umi_tag is not None
