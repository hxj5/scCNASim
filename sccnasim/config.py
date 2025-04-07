# config.py - global configuration.


import sys


class Config:
    """Configuration.

    Attributes
    ----------
    See :func:`~.main.main_wrapper()`.
    """
    def __init__(self):
        self.defaults = Defaults()

        # input and output files.
        self.sam_fn = None
        self.cell_anno_fn = None
        self.feature_fn = None
        self.snp_fn = None
        self.clone_anno_fn = None
        self.cna_profile_fn = None
        self.refseq_fn = None
        self.out_dir = None
        self.sam_list_fn = None
        self.sample_ids = None
        self.sample_id_fn = None
        
        # preprocessing.
        self.merge_features_how = "quantile"

        # count simulation.
        self.size_factor = "libsize"
        self.marginal = "auto"
        self.kwargs_fit_sf = dict()
        self.kwargs_fit_rd = dict()

        # optional arguments.
        self.chroms = ",".join([str(c) for c in range(1, 23)])
        self.cell_tag = self.defaults.CELL_TAG
        self.umi_tag = self.defaults.UMI_TAG
        self.umi_len = 10
        self.ncores = self.defaults.NCORES
        self.seed = 123
        self.verbose = False

        # snp filtering
        self.min_count = 1
        self.min_maf = 0
        
        # read assignment
        self.strandness = self.defaults.STRANDNESS
        self.min_include = self.defaults.MIN_INCLUDE

        # read filtering.
        self.min_mapq = self.defaults.MIN_MAPQ
        self.min_len = self.defaults.MIN_LEN
        self.incl_flag = self.defaults.INCL_FLAG
        self.excl_flag = -1
        self.no_orphan = self.defaults.NO_ORPHAN

        # others
        self.debug_level = 0

    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stdout

        s =  "%s\n" % prefix
        s += "%ssam_file = %s\n" % (prefix, self.sam_fn)
        s += "%scell_anno_file = %s\n" % (prefix, self.cell_anno_fn)
        s += "%sfeature_file = %s\n" % (prefix, self.feature_fn)
        s += "%sphased_snp_file = %s\n" % (prefix, self.snp_fn)
        s += "%sclone_anno_file = %s\n" % (prefix, self.clone_anno_fn)
        s += "%scna_profile_file = %s\n" % (prefix, self.cna_profile_fn)
        s += "%srefseq_file = %s\n" % (prefix, self.refseq_fn)
        s += "%sout_dir = %s\n" % (prefix, self.out_dir)
        s += "%ssam_list_file = %s\n" % (prefix, self.sam_list_fn)
        s += "%ssample_ids = %s\n" % (prefix, self.sample_ids)
        s += "%ssample_id_file = %s\n" % (prefix, self.sample_id_fn)
        s += "%s\n" % prefix
        
        s += "%smerge_features_how = %s\n" % (prefix, str(self.merge_features_how))
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
        s += "%sseed = %s\n" % (prefix, str(self.seed))
        s += "%sverbose = %s\n" % (prefix, self.verbose)
        s += "%s\n" % prefix

        s += "%smin_count = %d\n" % (prefix, self.min_count)
        s += "%smin_maf = %f\n" % (prefix, self.min_maf)
        s += "%s\n" % prefix

        s += "%sstrandness = %s\n" % (prefix, self.strandness)
        s += "%smin_include = %f\n" % (prefix, self.min_include)
        s += "%s\n" % prefix

        s += "%smin_mapq = %d\n" % (prefix, self.min_mapq)
        s += "%smin_len = %d\n" % (prefix, self.min_len)
        s += "%sinclude_flag = %d\n" % (prefix, self.incl_flag)
        s += "%sexclude_flag = %d\n" % (prefix, self.excl_flag)
        s += "%sno_orphan = %s\n" % (prefix, self.no_orphan)
        s += "%s\n" % prefix

        s += "%sdebug_level = %d\n" % (prefix, self.debug_level)
        s += "%s\n" % prefix

        fp.write(s)


    def is_stranded(self):
        return self.strandness in ("forward", "reverse")

    def use_barcodes(self):
        return self.cell_tag is not None

    def use_umi(self):
        return self.umi_tag is not None

    

class Defaults:
    def __init__(self):
        self.DEBUG = 0
        self.NCORES = 1
        
        self.CELL_TAG = "CB"
        self.UMI_TAG = "UB"
        self.UMI_TAG_BC = "UB"    # the default umi tag for 10x data.

        self.MIN_COUNT = 20
        self.MIN_MAF = 0.1
        
        self.STRANDNESS = "forward"
        self.MIN_INCLUDE = 0.9

        self.MIN_MAPQ = 20
        self.MIN_LEN = 30
        self.INCL_FLAG = 0
        self.EXCL_FLAG_UMI = 772
        self.EXCL_FLAG_XUMI = 1796
        self.NO_ORPHAN = True
