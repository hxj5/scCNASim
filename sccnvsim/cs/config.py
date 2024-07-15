# config.py

import sys


COMMAND = "cs"


class Config:
    """Configuration of count simulation.

    Attributes
    ----------
    count_fn : str
        A h5ad file storing the *cell x feature* count matrices for allele
        A, B, U in three layers "A", "B", "U", respectively.
        It should contain columns "cell", "cell_type" in its `.obs`; and 
        contain columns "chrom", "start", "end", "feature" in its `.var`.
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
    out_dir : str
        The output folder.
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
    ncores : int
        The number of cores/sub-processes.
        Default is `1`.
    verbose : bool
        Whether to show detailed logging information.
        Default is `False`.
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
        self.count_fn = None
        self.cnv_profile_fn = None
        self.clone_meta_fn = None
        self.out_dir = None

        self.size_factor = "libsize"
        self.marginal = "auto"
        self.ncores = 1
        self.verbose = False

        self.kwargs_fit_sf = dict()
        self.kwargs_fit_rd = dict()

        self.adata = None
        self.cnv_profile = None
        self.clone_meta = None

        self.out_prefix = COMMAND + "."

    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stdout
        
        s =  "%s\n" % prefix
        s += "%scount_fn = %s\n" % (prefix, self.count_fn)
        s += "%scnv_profile_fn = %s\n" % (prefix, self.cnv_profile_fn)
        s += "%sclone_meta_fn = %s\n" % (prefix, self.clone_meta_fn)
        s += "%sout_dir = %s\n" % (prefix, self.out_dir)

        s += "%ssize_factor = %s\n" % (prefix, self.size_factor)
        s += "%smarginal = %s\n" % (prefix, self.marginal)
        s += "%sncores = %s\n" % (prefix, self.ncores)
        s += "%sverbose = %s\n" % (prefix, self.verbose)

        s += "%skwargs_fit_sf = %s\n" % (prefix, self.kwargs_fit_sf)
        s += "%skwargs_fit_rd = %s\n" % (prefix, self.kwargs_fit_rd)

        s += "%sout_prefix = %s\n" % (prefix, self.out_prefix)
        s += "%s\n" % prefix

        fp.write(s)