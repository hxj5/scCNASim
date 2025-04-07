# config.py - configuration


import sys
from ..config import Defaults as MainDefaults

COMMAND = "rs"


class Config:
    """Configuration of the `rs` (read simulation) module.

    Attributes
    ----------
    See :func:`~.main.rs_wrapper()`.
    """
    def __init__(self):
        self.defaults = Defaults()

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

        # internal parameters.
        
        # out_prefix : str
        #   The prefix of the output files.
        self.out_prefix = COMMAND
        
        # cell_raw_tag : str
        #   Tag for uncorrected raw cell tag in seed and simulated BAM.
        self.cell_raw_tag = "CR"
        
        # umi_raw_tag : str or None
        #   Tag for uncorrected raw umi tag in seed and simulated BAM.
        #   Set to None if do not use it.
        self.umi_raw_tag = "UR"
        
        # backup_cell_tag : str
        #   Tag for backuping corrected cell barcode (from seed BAM) 
        #   in simulated BAM.
        self.backup_cell_tag = "KC"
        
        # backup_umi_tag : str or None
        #   Tag for backuping corrected umi barcode (from seed BAM) 
        #   in simulated BAM.
        #   Set to None if do not use it.
        self.backup_umi_tag = "KU"
        
        # alleles : tuple of str
        #   All alleles.
        self.alleles = ("A", "B", "U")


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
        
        # internal parameters.

        s += "%sout_prefix = %s\n" % (prefix, self.out_prefix)
        s += "%scell_raw_tag = %s\n" % (prefix, self.cell_raw_tag)
        s += "%sumi_raw_tag = %s\n" % (prefix, self.umi_raw_tag)
        s += "%sbackup_cell_tag = %s\n" % (prefix, self.backup_cell_tag)
        s += "%sbackup_umi_tag = %s\n" % (prefix, self.backup_umi_tag)
        s += "%salleles = %s\n" % (prefix, str(self.alleles))
        s += "%s\n" % prefix

        fp.write(s)


    def use_barcodes(self):
        return self.cell_tag is not None

    def use_umi(self):
        return self.umi_tag is not None



class Defaults(MainDefaults):
    def __init__(self):
        super().__init__()
        self.UMI_LEN = 10



if __name__ == "__main__":
    conf = Config()
    conf.show()
