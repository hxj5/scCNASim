# config.py

import sys


COMMAND = "cs"


class Config:
    """Configuration of the `cs` (count simulation) module.

    Attributes
    ----------
    count_fn : str
        A h5ad file storing the *cell x feature* count matrices for allele
        A, B, U in three layers "A", "B", "U", respectively.
        Its `.obs` should contain columns:
        - "cell" (str): cell barcodes.
        - "cell_type" (str): cell type.
        Its `.var` should contain columns:
        - "chrom" (str): chromosome name of the feature.
        - "start" (int): start genomic position of the feature, 1-based
          and inclusive.
        - "end" (int): end genomic position of the feature, 1-based and
          inclusive.
        - "feature" (str): feature name.
    clone_meta_fn : str
        A TSV file listing clonal meta information.
        It is header-free and its first 3 columns are:
        - "clone" (str): clone ID.
        - "source_cell_type" (str): the source cell type of `clone`.
        - "n_cell" (int): number of cells in the `clone`. If negative, 
          then it will be set as the number of cells in `source_cell_type`.
    cnv_profile_fn : str
        A TSV file listing clonal CNV profiles.
        It is header-free and its first 7 columns are:
        - "chrom" (str): chromosome name of the CNV region.
        - "start" (int): start genomic position of the CNV region, 1-based
          and inclusive.
        - "end" (int): end genomic position of the CNV region, 1-based and
          inclusive.
        - "region" (str): ID of the CNV region.
        - "clone" (str): clone ID.
        - "cn_ale0" (int): copy number of the first allele.
        - "cn_ale1" (int): copy number of the second allele.
    out_dir : str
        The output folder.
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
    ncores : int, default 1
        The number of cores/sub-processes.
    verbose : bool, default False
        Whether to show detailed logging information.
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
    """
    def __init__(self):
        # command-line arguments/parameters.
        self.count_fn = None
        self.clone_meta_fn = None
        self.cnv_profile_fn = None
        self.out_dir = None

        self.size_factor = "libsize"
        self.marginal = "auto"
        self.ncores = 1
        self.verbose = False

        self.kwargs_fit_sf = dict()
        self.kwargs_fit_rd = dict()

        # derived parameters.

        # adata : anndata.AnnData
        #   The loaded allele-specific count matrices.
        self.adata = None

        # clone_meta : pandas.DataFrame
        #   The loaded clone annotations.
        self.clone_meta = None

        # cnv_profile : pandas.DataFrame
        #   The loaded clonal CNV profile.
        self.cnv_profile = None

        # out_prefix : str
        #   Prefix to the output files.
        self.out_prefix = COMMAND + "."

    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stdout
        
        # command-line arguments/parameters.
        s =  "%s\n" % prefix
        s += "%scount_fn = %s\n" % (prefix, self.count_fn)
        s += "%sclone_meta_fn = %s\n" % (prefix, self.clone_meta_fn)
        s += "%scnv_profile_fn = %s\n" % (prefix, self.cnv_profile_fn)
        s += "%sout_dir = %s\n" % (prefix, self.out_dir)

        s += "%ssize_factor = %s\n" % (prefix, self.size_factor)
        s += "%smarginal = %s\n" % (prefix, self.marginal)
        s += "%sncores = %s\n" % (prefix, self.ncores)
        s += "%sverbose = %s\n" % (prefix, self.verbose)

        s += "%skwargs_fit_sf = %s\n" % (prefix, self.kwargs_fit_sf)
        s += "%skwargs_fit_rd = %s\n" % (prefix, self.kwargs_fit_rd)

        # derived parameters.
        s += "%sout_prefix = %s\n" % (prefix, self.out_prefix)
        s += "%s\n" % prefix

        fp.write(s)
