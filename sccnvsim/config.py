# config.py - global configuration.


import sys

from .afc.config import Config as AFC_Conf
from .afc.config import DefaultConfig as AFC_DefConf
from .cs.config import Config as CS_Conf
from .pp.config import Config as PP_Conf


class Config:
    def __init__(self):
        self.g = GlobalConfig()
        self.pp = PPConfig()
        self.afc = AFCConfig()
        self.cs = CSConfig()
        self.rs = None

    def show(self, fp = None):
        if fp is None:
            fp = sys.stdout
        
        fp.write("Global Config:\n")
        self.g.show(fp = fp, prefix = "\t")

        fp.write("Preprocessing Config:\n")
        self.pp.show(fp = fp, prefix = "\t")

        fp.write("Allele-Specific Feature Counting Config:\n")
        self.afc.show(fp = fp, prefix = "\t")

        fp.write("Count Simulation Config:\n")
        self.cs.show(fp = fp, prefix = "\t")

        fp.write("Read Simulation Config:\n")
        #self.rs.show(fp = fp, prefix = "\t")   # TODO: implement


class GlobalConfig:
    """Global configuration

    Attributes
    ----------
    out_dir : str
        The output folder.
    cell_tag : str
        Tag for cell barcodes, set to `None` when using sample IDs.
        Default is `"CB"`.
    umi_tag : str
        Tag for UMI, set to `None` when reads only.
        Default is `"UB"`.
    ncores : int
        Number of processes.
        Default is `1`.
    verbose : bool
        Whether to show detailed logging information.
        Default is `False`.
    """
    def __init__(self):
        self._afc_def_conf = AFC_DefConf()

        self.out_dir = None
        self.cell_tag = self._afc_def_conf.CELL_TAG
        self.umi_tag = self._afc_def_conf.UMI_TAG
        self.ncores = self._afc_def_conf.NPROC
        self.verbose = False

    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stdout
        
        s =  "%s\n" % prefix
        s += "%sout_dir = %s\n" % (prefix, self.out_dir)
        s += "%scell_tag = %s\n" % (prefix, self.cell_tag)
        s += "%sumi_tag = %s\n" % (prefix, self.umi_tag)
        s += "%sncores = %s\n" % (prefix, self.ncores)
        s += "%sverbose = %s\n" % (prefix, self.verbose)
        s += "%s\n" % prefix

        fp.write(s)


class PPConfig(PP_Conf):
    """Configuration of preprocessing.

    Attributes
    ----------
    cell_anno_fn : str
        The cell annotation file. It is header-free and its first two columns
        are `cell` and `cell_type`.
    feature_fn : str
        A TSV file listing target features. It is header-free and its first 
        4 columns shoud be: `chrom`, `start`, `end` (both start and end are
        1-based and inclusive), and `feature_name`.
    snp_fn : str
        A TSV or VCF file listing phased SNPs (i.e., containing phased GT).
    cnv_profile_fn : str
        A TSV file listing clonal CNV profiles. It is header-free and its first
        7 columns are "chrom" (str), "start" (int), "end" (int), 
        "reg_id" (str), "clone_id" (str), "cn_ale0" (int), "cn_ale1" (int).
        Note that both "start" and "end" are 1-based and inclusive.
    clone_meta_fn : str
        A TSV file listing clonal meta information. It is header-free and its
        first 3 columns are "clone_id" (str), "ref_cell_type" (str),
        "n_cells" (int). If "n_cells" is negative, then it will be set as
        the number of cells in "ref_cell_type".
    """
    def __init__(self):
        super().__init__()

    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stdout
        
        s =  "%s\n" % prefix
        s += "%scell_anno_fn = %s\n" % (prefix, self.cell_anno_fn)
        s += "%sfeature_fn = %s\n" % (prefix, self.feature_fn)
        s += "%ssnp_fn = %s\n" % (prefix, self.snp_fn)
        s += "%scnv_profile_fn = %s\n" % (prefix, self.cnv_profile_fn)
        s += "%sclone_meta_fn = %s\n" % (prefix, self.clone_meta_fn)
        s += "%s\n" % prefix

        fp.write(s)


# only list selected attributes here.
class AFCConfig(AFC_Conf):
    """Configuration of allele-specific feature counting.

    Attributes
    ----------
    sam_fn : str
        Comma separated indexed sam/bam/cram file.
    sam_list_fn : str
        A list file containing bam files, each per line.
    sample_id_str : str
        Comma separated sample IDs. Typically used in well-based data (e.g.,
        SMART-seq2).
    sample_id_fn : str
        A list file containing sample IDs, each per line.
        Typically used in well-based data (e.g., SMART-seq2).
    min_count : int
        Mininum aggragated count for SNP.
        Default is `1`.
    min_maf : int
        Mininum minor allele fraction for SNP.
        Default is `0`.
    min_mapq : int
        Minimum MAPQ for read filtering.
        Default is `20`.
    min_len : int
        Minimum mapped length for read filtering.
        Default is `30`.
    incl_flag : int
        Required flags: skip reads with all mask bits unset.
        Default is `0`.
    excl_flag : int
        Filter flags: skip reads with any mask bits set.
        Default is `772` (when use UMI) or `1796` (otherwise).
    no_orphan : bool
        If `True`, skip anomalous read pairs.
        Default is `True`.
    """
    def __init__(self):
        super().__init__()

    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stdout

        s =  "%s\n" % prefix
        s += "%ssam_file = %s\n" % (prefix, self.sam_fn)
        s += "%ssam_list_file = %s\n" % (prefix, self.sam_list_fn)
        s += "%ssample_id_str = %s\n" % (prefix, self.sample_id_str)
        s += "%ssample_id_file = %s\n" % (prefix, self.sample_id_fn)

        s += "%smin_count = %d\n" % (prefix, self.min_count)
        s += "%smin_maf = %f\n" % (prefix, self.min_maf)

        s += "%smin_mapq = %d\n" % (prefix, self.min_mapq)
        s += "%smin_len = %d\n" % (prefix, self.min_len)
        s += "%sinclude_flag = %d\n" % (prefix, self.incl_flag)
        s += "%sexclude_flag = %d\n" % (prefix, self.excl_flag)
        s += "%sno_orphan = %s\n" % (prefix, self.no_orphan)
        s += "%s\n" % prefix

        fp.write(s)


class CSConfig(CS_Conf):
    """Configuration of count simulation.

    Attributes
    ----------
    size_factor : str
        The type of size factor.
        Currently, only support "libsize" (library size).
        Set to `None` if do not use size factors for model fitting.
        Default is `"libsize"`.
    marginal : str
        Type of marginal distribution.
        One of "auto" (auto select), "poisson" (Poisson), 
        "nb" (Negative Binomial),
        and "zinb" (Zero-Inflated Negative Binomial).
        Default is `"auto"`.
    kwargs_fit_sf : dict
        The additional kwargs passed to function 
        :func:`~marginal.fit_libsize_wrapper` for fitting size factors.
        Default is `{}`.
        The available arguments are:
        dist : str
            Type of distribution. One of "normal" (normal) and "t" (t).
            Default is `"normal"`.
    kwargs_fit_rd : dcit
        The additional kwargs passed to function 
        :func:`~marginal.fit_RD_wrapper` for fitting read depth.
        Default is `{}`.
        The available arguments are:
        min_nonzero_num : int
            The minimum number of cells that have non-zeros for one feature.
            If smaller than the cutoff, then the feature will not be fitted
            (i.e., its mean will be directly treated as 0).
            Default is `3`.
        max_iter : int
            Number of maximum iterations in model fitting.
            Default is `1000`.
        pval_cutoff : float
            The p-value cutoff for model selection with GLR test.
            Default is `0.05`.
    """
    def __init__(self):
        super().__init__()

    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stdout
        
        s =  "%s\n" % prefix
        s += "%ssize_factor = %s\n" % (prefix, self.size_factor)
        s += "%smarginal = %s\n" % (prefix, self.marginal)
        s += "%skwargs_fit_sf = %s\n" % (prefix, self.kwargs_fit_sf)
        s += "%skwargs_fit_rd = %s\n" % (prefix, self.kwargs_fit_rd)
        s += "%s\n" % prefix

        fp.write(s)