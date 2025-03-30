# config.py - configuration


import sys

COMMAND = "afc"

class Config:
    """Configuration of the `afc` (allele-specific counting) module.

    Attributes
    ----------
    See `afc::main::afc_wrapper()`.
    """
    def __init__(self):
        # defaults : DefaultConfig
        #   The default values of parameters.
        self.defaults = DefaultConfig()

        # argv : list of str or None, default None
        #   A list of command line arguments, typically from sys.argv.
        self.argv = None

        # command-line arguments/parameters.
        self.sam_fn = None
        self.sam_list_fn = None
        self.barcode_fn = None
        self.sample_ids = None
        self.sample_id_fn = None
        self.feature_fn = None
        self.snp_fn = None
        self.out_dir = None
        self.debug = self.defaults.DEBUG

        self.cell_tag = self.defaults.CELL_TAG
        self.umi_tag = self.defaults.UMI_TAG
        self.ncores = self.defaults.NCORES
        
        # snp filtering.
        self.min_count = self.defaults.MIN_COUNT
        self.min_maf = self.defaults.MIN_MAF
                
        # read assignment.
        self.strandness = self.defaults.STRANDNESS
        self.min_include = self.defaults.MIN_INCLUDE
        
        # read filtering.
        self.min_mapq = self.defaults.MIN_MAPQ
        self.min_len = self.defaults.MIN_LEN
        self.incl_flag = self.defaults.INCL_FLAG
        self.excl_flag = -1
        self.no_orphan = self.defaults.NO_ORPHAN
        
        # out_feature_dirs : list of str
        #   A list of output folders for feature-specific results.
        self.out_feature_dirs = None


        # derived parameters.

        # barcodes : list of str or None
        #   A list of cell barcodes.
        #   None if sample IDs are used.
        self.barcodes = None

        # sample_ids : list of str or None
        #   A list of sample IDs.
        #   None if cell barcodes are used.
        self.sample_ids = None

        # reg_list : list of utils.gfeature.Feature
        #   A list of features.
        self.reg_list = None

        # snp_set : utils.gfeature.SNPSet
        #   The object storing a set of SNPs.
        self.snp_set = None

        # sam_fn_list : list of str
        #   A list of input SAM/BAM files. 
        self.sam_fn_list = None

        # samples : list of str
        #   A list of cell barcodes (droplet-based data) or sample IDs (
        #   well-based data).
        #   It will be used as output IDs of each cell.        
        self.samples = None

        # count_dir : str
        #   The output folder for count matrices.
        self.count_dir = None

        # feature_obj_fn : str
        #   Path to a python pickle file storing the `reg_list`.
        self.feature_obj_fn = None

        # alleles : tuple of str
        #   All alleles.
        self.alleles = ("A", "B", "D", "O", "U")

        # cumi_alleles : tuple of str
        #   Alleles whose CUMIs will be outputed for read sampling.
        self.cumi_alleles = ("A", "B", "U")
        
        # hap_idx_tag : int
        #   Tag for haplotype index.
        self.hap_idx_tag = "HI"

        # out_prefix : str
        #   The prefix of the output files.
        self.out_prefix = COMMAND + "."

        # out_feature_fn : str
        #   Path to the output feature file.
        self.out_feature_fn = None

        # out_sample_fn : str
        #   Path to the output sample file.
        self.out_sample_fn = None

        # out_ale_fns : dict of {str : str}
        #   The allele-specific *feature x cell* count matrices.
        #   Keys are allele names and values are pathes to the count matrix
        #   files.
        self.out_ale_fns = {ale:None for ale in self.alleles}

        # out_adata_fn : str
        #   Path to a ".adata" file storing a :class:`~anndata.Anndata`
        #   object, which contains all allele-specific *feature x cell* count
        #   matrices.
        self.out_adata_fn = None

    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stdout

        s =  "%s\n" % prefix
        s += "%ssam_file = %s\n" % (prefix, self.sam_fn)
        s += "%ssam_list_file = %s\n" % (prefix, self.sam_list_fn)
        s += "%sbarcode_file = %s\n" % (prefix, self.barcode_fn)
        s += "%ssample_ids = %s\n" % (prefix, self.sample_ids)
        s += "%ssample_id_file = %s\n" % (prefix, self.sample_id_fn)
        s += "%sfeature_file = %s\n" % (prefix, self.feature_fn)
        s += "%ssnp_file = %s\n" % (prefix, self.snp_fn)
        s += "%sout_dir = %s\n" % (prefix, self.out_dir)
        s += "%sdebug_level = %d\n" % (prefix, self.debug)
        s += "%s\n" % prefix

        s += "%scell_tag = %s\n" % (prefix, self.cell_tag)
        s += "%sumi_tag = %s\n" % (prefix, self.umi_tag)
        s += "%snumber_of_cores = %d\n" % (prefix, self.ncores)
        s += "%s\n" % prefix
        
        # snp filtering.
        s += "%smin_count = %d\n" % (prefix, self.min_count)
        s += "%smin_maf = %f\n" % (prefix, self.min_maf)
        s += "%s\n" % prefix

        # read assignment.
        s += "%sstrandness = %s\n" % (prefix, self.strandness)
        s += "%smin_include = %f\n" % (prefix, self.min_include)
        s += "%s\n" % prefix

        # read filtering.
        s += "%smin_mapq = %d\n" % (prefix, self.min_mapq)
        s += "%smin_len = %d\n" % (prefix, self.min_len)
        s += "%sinclude_flag = %d\n" % (prefix, self.incl_flag)
        s += "%sexclude_flag = %d\n" % (prefix, self.excl_flag)
        s += "%sno_orphan = %s\n" % (prefix, self.no_orphan)
        s += "%s\n" % prefix
        
        s += "%slen(out_feature_dirs) = %d\n" % (prefix, \
            len(self.out_feature_dirs) if self.out_feature_dirs else 0)
        s += "%s\n" % prefix

        
        # derived parameters.
        
        s += "%snumber_of_BAMs = %d\n" % (prefix, len(self.sam_fn_list) if \
                self.sam_fn_list is not None else -1)
        s += "%snumber_of_barcodes = %d\n" % (prefix, len(self.barcodes) if \
                self.barcodes is not None else -1)
        s += "%snumber_of_sample_IDs = %d\n" % (prefix, len(self.sample_ids) \
                if self.sample_ids is not None else -1)
        s += "%snumber_of_features = %d\n" % (prefix, len(self.reg_list) if \
                self.reg_list is not None else -1)
        s += "%snumber_of_snps = %d\n" % (prefix, self.snp_set.get_n() if \
                self.snp_set is not None else -1)
        s += "%s\n" % prefix

        s += "%scount_dir = %s\n" % (prefix, self.count_dir)
        s += "%sfeature_obj_fn = %s\n" % (prefix, self.feature_obj_fn)
        s += "%s\n" % prefix

        s += "%salleles = %s\n" % (prefix, str(self.alleles))
        s += "%scumi_alleles = %s\n" % (prefix, str(self.cumi_alleles))
        s += "%shap_idx_tag = %s\n" % (prefix, self.hap_idx_tag)
        s += "%s\n" % prefix

        s += "%soutput_feature_file = %s\n" % (prefix, self.out_feature_fn)
        s += "%soutput_sample_file = %s\n" % (prefix, self.out_sample_fn)
        for ale, fn in self.out_ale_fns.items():
            s += "%soutput_ale_%s_file = %s\n" % (prefix, ale, fn)
        s += "%s\n" % prefix

        s += "%sout_adata_fn = %s\n" % (prefix, self.out_adata_fn)
        s += "%s\n" % prefix

        fp.write(s)
        
        
    def is_stranded(self):
        return self.strandness in ("forward", "reverse")

    def use_barcodes(self):
        return self.cell_tag is not None

    def use_umi(self):
        return self.umi_tag is not None


class DefaultConfig:
    def __init__(self):
        self.DEBUG = 0
        self.CELL_TAG = "CB"
        self.UMI_TAG = "UB"
        self.UMI_TAG_BC = "UB"    # the default umi tag for 10x data.
        self.NCORES = 1
        
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


if __name__ == "__main__":
    conf = Config()
    conf.show()
