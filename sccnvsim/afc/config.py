# config.py - configuration

import sys

COMMAND = "afc"


class Config:
    def __init__(self):
        self.defaults = DefaultConfig()
        self.argv = None

        self.sam_fn = None
        self.sam_list_fn = None
        self.barcode_fn = None
        self.sample_id_str = None
        self.sample_id_fn = None
        self.feature_fn = None
        self.snp_fn = None
        self.out_dir = None
        self.debug = self.defaults.DEBUG

        self.cell_tag = self.defaults.CELL_TAG
        self.umi_tag = self.defaults.UMI_TAG
        self.nproc = self.defaults.NPROC
        self.min_count = self.defaults.MIN_COUNT
        self.min_maf = self.defaults.MIN_MAF

        self.min_mapq = self.defaults.MIN_MAPQ
        self.min_len = self.defaults.MIN_LEN
        self.incl_flag = self.defaults.INCL_FLAG
        self.excl_flag = -1
        self.no_orphan = self.defaults.NO_ORPHAN

        self.barcodes = None     # list of barcode strings.
        self.sample_ids = None
        self.reg_list = None     # list of gene/block objects.
        self.snp_set = None      # set of SNPs.

        self.sam_fn_list = None
        self.samples = None

        self.aln_dir = None
        self.count_dir = None

        self.out_prefix = COMMAND + "."
        self.out_feature_fn = None
        self.out_sample_fn = None
        self.out_ale_fns = {ale:None for ale in ("A", "B", "D", "O", "U")}

        # `out_feature_meta_fn`: a python pickle file storing the 
        # `self.reg_list`.
        # It will be saved after extracting reads from input BAM(s),
        # and be re-loaded for read sampling.
        self.out_feature_meta_fn = None
        self.out_adata_fn = None

        # Unique UMI tag (typically for cell barcode+UMI).
        self.uumi_tag = "UU"

    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stdout

        s =  "%s\n" % prefix
        s += "%ssam_file = %s\n" % (prefix, self.sam_fn)
        s += "%ssam_list_file = %s\n" % (prefix, self.sam_list_fn)
        s += "%sbarcode_file = %s\n" % (prefix, self.barcode_fn)
        s += "%ssample_id_str = %s\n" % (prefix, self.sample_id_str)
        s += "%ssample_id_file = %s\n" % (prefix, self.sample_id_fn)
        s += "%sfeature_file = %s\n" % (prefix, self.feature_fn)
        s += "%ssnp_file = %s\n" % (prefix, self.snp_fn)
        s += "%sout_dir = %s\n" % (prefix, self.out_dir)
        s += "%sdebug_level = %d\n" % (prefix, self.debug)
        s += "%s\n" % prefix

        s += "%scell_tag = %s\n" % (prefix, self.cell_tag)
        s += "%sumi_tag = %s\n" % (prefix, self.umi_tag)
        s += "%snumber_of_processes = %d\n" % (prefix, self.nproc)
        s += "%smin_count = %d\n" % (prefix, self.min_count)
        s += "%smin_maf = %f\n" % (prefix, self.min_maf)
        s += "%s\n" % prefix

        s += "%smin_mapq = %d\n" % (prefix, self.min_mapq)
        s += "%smin_len = %d\n" % (prefix, self.min_len)
        s += "%sinclude_flag = %d\n" % (prefix, self.incl_flag)
        s += "%sexclude_flag = %d\n" % (prefix, self.excl_flag)
        s += "%sno_orphan = %s\n" % (prefix, self.no_orphan)
        s += "%s\n" % prefix

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

        s += "%saln_dir = %s\n" % (prefix, self.aln_dir)
        s += "%scount_dir = %s\n" % (prefix, self.count_dir)
        s += "%s\n" % prefix

        s += "%soutput_feature_file = %s\n" % (prefix, self.out_feature_fn)
        s += "%soutput_sample_file = %s\n" % (prefix, self.out_sample_fn)
        for ale, fn in self.out_ale_fns.items():
            s += "%soutput_ale_%s_file = %s\n" % (prefix, ale, fn)
        s += "%s\n" % prefix

        s += "%sout_feature_meta_fn = %s\n" % (prefix, self.out_feature_meta_fn)
        s += "%sout_adata_fn = %s\n" % (prefix, self.out_adata_fn)
        s += "%s\n" % prefix

        s += "%suumi_tag = %s\n" % (prefix, self.uumi_tag)
        s += "%s\n" % prefix

        fp.write(s)

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
        self.NPROC = 1
        self.MIN_COUNT = 1 
        self.MIN_MAF = 0

        self.MIN_MAPQ = 20
        self.MIN_LEN = 30
        self.INCL_FLAG = 0
        self.EXCL_FLAG_UMI = 772
        self.EXCL_FLAG_XUMI = 1796
        self.NO_ORPHAN = True


if __name__ == "__main__":
    conf = Config()
    conf.show()