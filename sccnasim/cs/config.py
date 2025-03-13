# config.py

import sys


COMMAND = "cs"


class Config:
    """Configuration of the `cs` (count simulation) module.

    Attributes
    ----------
    See `cs::main::cs_wrapper()`.
    """
    def __init__(self):
        # command-line arguments/parameters.
        self.count_fn = None
        self.clone_meta_fn = None
        self.cna_profile_fn = None
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

        # cna_profile : pandas.DataFrame
        #   The loaded clonal CNA profile.
        self.cna_profile = None

        # out_prefix : str
        #   Prefix to the output files.
        self.out_prefix = COMMAND + "."


        # internal parameters.
        
        # def_kwargs_fit_sf : dict
        #   Default settings passing to `kwargs_fit_sf`.
        self.def_kwargs_fit_sf = {
            "dist": "lognormal"
        }
        
        # def_kwargs_fit_rd : dict
        #   Default settings passing to `kwargs_fit_rd`.
        self.def_kwargs_fit_rd = {
            "min_nonzero_num": 5,
            "max_iter": 1000,
            "pval_cutoff": 0.05
        }
        
        # qc_min_library_size : int
        #   Minimum library size of one cell.
        self.qc_min_library_size = 1000
        
        # qc_max_library_size : int or None
        #   Maximum library size of one cell.
        #   If ``None``, there is no limit.
        self.qc_max_library_size = None
        
        # qc_min_features : int or float
        #   Minimum expressed features in one cell.
        #   If ``float``, it is the fraction of all input features.
        self.qc_min_features = 0.01
        
        # qc_cw_low_quantile : float
        #   The lower quantile of cell-wise statistics.
        self.qc_cw_low_quantile = 0.01
        
        # qc_cw_up_quantile : float
        #   The upper quantile of cell-wise statistics.
        self.qc_cw_up_quantile = 0.99
        

    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stdout
        
        # command-line arguments/parameters.
        s =  "%s\n" % prefix
        s += "%scount_fn = %s\n" % (prefix, self.count_fn)
        s += "%sclone_meta_fn = %s\n" % (prefix, self.clone_meta_fn)
        s += "%scna_profile_fn = %s\n" % (prefix, self.cna_profile_fn)
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
        
        # internal parameters.
        s += "%sqc_min_library_size = %s\n" % (prefix, self.qc_min_library_size)
        s += "%sqc_max_library_size = %s\n" % (prefix, self.qc_max_library_size)
        s += "%sqc_min_features = %s\n" % (prefix, self.qc_min_features)
        s += "%sqc_cw_low_quantile = %s\n" % (prefix, self.qc_cw_low_quantile)
        s += "%sqc_cw_up_quantile = %s\n" % (prefix, self.qc_cw_up_quantile)
        s += "%s\n" % prefix        

        fp.write(s)
