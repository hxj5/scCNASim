# config.py - configuration

import sys
from ..afc.config import DefaultConfig as AFC_Def_Conf

COMMAND = "rs"


class Config:
    """Configuration of the `rs` (read simulation) module.

    Attributes
    ----------
    See `rs::main::rs_wrapper()`.
    """
    def __init__(self):
        # defaults : DefaultConfig
        #   The default values of parameters.
        self.defaults = DefaultConfig()

        # argv : list of str or None, default None
        #   A list of command line arguments, typically from sys.argv.
        self.argv = None

        # command-line arguments/parameters.
        self.count_fn = None
        self.feature_fn = None
        self.refseq_fn = None
        self.out_dir = None
        self.debug = self.defaults.DEBUG

        self.cell_tag = self.defaults.CELL_TAG
        self.umi_tag = self.defaults.UMI_TAG
        self.umi_len = self.defaults.UMI_LEN
        self.ncores = self.defaults.NCORES
        
        # derived variables

        # adata : anndata.Anndata
        #   The object storing count matrices.
        self.adata = None

        # reg_list : list of utils.gfeature.Feature
        #   A list of features.
        self.reg_list = None

        # out_prefix : str
        #   The prefix of the output files.
        self.out_prefix = COMMAND + "."

        # internal parameters.

        # hap_tag : str
        #   Tag for haplotype in the output BAM file.
        self.hap_tag = "HT"
        
        # cell_raw_tag : str or None
        #   Tag for uncorrected raw cell tag in seed and simulated BAM.
        #   Set to None if do not use it.
        self.cell_raw_tag = "CR"
        
        # backup_cell_raw_tag : str or None
        #   Tag for backup uncorrected raw cell barcode (from seed BAM) 
        #   in simulated BAM. Set to None if do not use it.
        self.backup_cell_raw_tag = "RC"
        
        # umi_raw_tag : str or None
        #   Tag for uncorrected raw umi tag in seed and simulated BAM.
        #   Set to None if do not use it.
        self.umi_raw_tag = "UR"
        
        # backup_umi_raw_tag : str or None
        #   Tag for backup uncorrected raw umi barcode (from seed BAM) 
        #   in simulated BAM. Set to None if do not use it.
        self.backup_umi_raw_tag = "RU"
        
        # alleles : tuple of str
        #   All alleles.
        self.alleles = ("A", "B", "U")

        # out_sam_dir : str
        #   Output folder for SAM/BAM file(s).
        self.out_sam_dir = None
        
        # out_step_dir : str
        #   Output folder for step-wise results.
        self.out_step_dir = None


    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stdout

        s =  "%s\n" % prefix
        s += "%scount_file = %s\n" % (prefix, self.count_fn)
        s += "%sfeature_file = %s\n" % (prefix, self.feature_fn)
        s += "%srefseq_file = %s\n" % (prefix, self.refseq_fn)
        s += "%sout_dir = %s\n" % (prefix, self.out_dir)
        s += "%sdebug_level = %d\n" % (prefix, self.debug)
        s += "%s\n" % prefix

        s += "%scell_tag = %s\n" % (prefix, self.cell_tag)
        s += "%sumi_tag = %s\n" % (prefix, self.umi_tag)
        s += "%sumi_len = %s\n" % (prefix, self.umi_len)
        s += "%snumber_of_processes = %d\n" % (prefix, self.ncores)
        s += "%s\n" % prefix

        # derived variables

        s += "%sshape_of_adata = %s\n" % (prefix, self.adata.shape if \
                self.adata is not None else "None")
        s += "%snumber_of_features = %d\n" % (prefix, len(self.reg_list) if \
                self.reg_list is not None else -1)
        s += "%s\n" % prefix
        
        # internal parameters.

        s += "%shap_tag = %s\n" % (prefix, self.hap_tag)
        s += "%scell_raw_tag = %s\n" % (prefix, self.cell_raw_tag)
        s += "%sbackup_cell_raw_tag = %s\n" % (prefix, self.backup_cell_raw_tag)
        s += "%sumi_raw_tag = %s\n" % (prefix, self.umi_raw_tag)
        s += "%sbackup_umi_raw_tag = %s\n" % (prefix, self.backup_umi_raw_tag)
        s += "%salleles = %s\n" % (prefix, str(self.alleles))
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


if __name__ == "__main__":
    conf = Config()
    conf.show()
