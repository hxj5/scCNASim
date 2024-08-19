# config.py - configuration

import sys
from ..afc.config import DefaultConfig as AFC_Def_Conf

COMMAND = "rs"


class Config:
    def __init__(self):
        self.defaults = DefaultConfig()
        self.argv = None

        self.sam_fn = None
        self.sam_list_fn = None
        self.barcode_fn = None
        self.sample_id_str = None
        self.sample_id_fn = None
        self.count_fn = None
        self.feature_fn = None
        self.refseq_fn = None
        self.out_dir = None
        self.debug = self.defaults.DEBUG

        self.chroms = self.defaults.CHROMS
        self.cell_tag = self.defaults.CELL_TAG
        self.umi_tag = self.defaults.UMI_TAG
        self.umi_len = self.defaults.UMI_LEN
        self.nproc = self.defaults.NPROC

        self.min_mapq = self.defaults.MIN_MAPQ
        self.min_len = self.defaults.MIN_LEN
        self.incl_flag = self.defaults.INCL_FLAG
        self.excl_flag = -1
        self.no_orphan = self.defaults.NO_ORPHAN

        self.chrom_list = None
        self.barcodes = None     # list of barcode strings.
        self.sample_ids = None
        self.adata = None        # the adata containing count matrices.
        self.reg_list = None     # list of gene/block objects.

        self.sam_fn_list = None
        self.samples = None

        self.out_prefix = COMMAND + "."

        self.hap_tag = "HT"       # tag for haplotype in BAM file.
        self.alleles = ("A", "B", "U")

        #self.cumi_max_pool = (1000, 1000, 10000)  # for allele A,B,U
        self.cumi_max_pool = (0, 0, 0)     # 0 means ulimited.

        self.out_sam_dir = None
        self.out_step_dir = None


    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stdout

        s =  "%s\n" % prefix
        s += "%ssam_file = %s\n" % (prefix, self.sam_fn)
        s += "%ssam_list_file = %s\n" % (prefix, self.sam_list_fn)
        s += "%sbarcode_file = %s\n" % (prefix, self.barcode_fn)
        s += "%ssample_id_str = %s\n" % (prefix, self.sample_id_str)
        s += "%ssample_id_file = %s\n" % (prefix, self.sample_id_fn)
        s += "%scount_file = %s\n" % (prefix, self.count_fn)
        s += "%sfeature_file = %s\n" % (prefix, self.feature_fn)
        s += "%srefseq_file = %s\n" % (prefix, self.refseq_fn)
        s += "%sout_dir = %s\n" % (prefix, self.out_dir)
        s += "%sdebug_level = %d\n" % (prefix, self.debug)
        s += "%s\n" % prefix

        s += "%schroms = %s\n" % (prefix, self.chroms)
        s += "%scell_tag = %s\n" % (prefix, self.cell_tag)
        s += "%sumi_tag = %s\n" % (prefix, self.umi_tag)
        s += "%sumi_len = %s\n" % (prefix, self.umi_len)
        s += "%snumber_of_processes = %d\n" % (prefix, self.nproc)
        s += "%s\n" % prefix

        s += "%smin_mapq = %d\n" % (prefix, self.min_mapq)
        s += "%smin_len = %d\n" % (prefix, self.min_len)
        s += "%sinclude_flag = %d\n" % (prefix, self.incl_flag)
        s += "%sexclude_flag = %d\n" % (prefix, self.excl_flag)
        s += "%sno_orphan = %s\n" % (prefix, self.no_orphan)
        s += "%s\n" % prefix

        s += "%schrom_list = %s\n" % (prefix, str(self.chrom_list))
        s += "%snumber_of_BAMs = %d\n" % (prefix, len(self.sam_fn_list) if \
                self.sam_fn_list is not None else -1)
        s += "%snumber_of_barcodes = %d\n" % (prefix, len(self.barcodes) if \
                self.barcodes is not None else -1)
        s += "%snumber_of_sample_IDs = %d\n" % (prefix, len(self.sample_ids) \
                if self.sample_ids is not None else -1)
        s += "%sshape_of_adata = %s\n" % (prefix, self.adata.shape if \
                self.adata is not None else "None")
        s += "%snumber_of_features = %d\n" % (prefix, len(self.reg_list) if \
                self.reg_list is not None else -1)
        s += "%s\n" % prefix

        s += "%shap_tag = %s\n" % (prefix, self.hap_tag)
        s += "%salleles = %s\n" % (prefix, str(self.alleles))
        s += "%scumi_max_pool = %s\n" % (prefix, str(self.cumi_max_pool))
        s += "%s\n" % prefix

        s += "%sout_sam_dir = %s\n" % (prefix, self.out_sam_dir)
        s += "%sout_step_dir = %s\n" % (prefix, self.out_step_dir)
        s += "%s\n" % prefix

        fp.write(s)

    def use_barcodes(self):
        return self.cell_tag is not None

    def use_umi(self):
        return self.umi_tag is not None


class DefaultConfig(AFC_Def_Conf):
    def __init__(self):
        super().__init__()
        self.UMI_LEN = 10
        self.CHROMS = ",".join([str(i) for i in range(1, 23)] + ["X", "Y"])


if __name__ == "__main__":
    conf = Config()
    conf.show()