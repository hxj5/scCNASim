# config.py

import sys

class Config:
    """Configuration of preprocessing

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
    out_dir : str
        The output folder.
    verbose : bool
        Whether show detailed logging information.
    """
    def __init__(self):
        self.cell_anno_fn = None
        self.feature_fn = None
        self.snp_fn = None
        self.cnv_profile_fn = None
        self.clone_meta_fn = None
        self.out_dir = None
        self.verbose = True

    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stdout
        
        s =  "%s\n" % prefix
        s += "%scell_anno_fn = %s\n" % (prefix, self.cell_anno_fn)
        s += "%sfeature_fn = %s\n" % (prefix, self.feature_fn)
        s += "%ssnp_fn = %s\n" % (prefix, self.snp_fn)
        s += "%scnv_profile_fn = %s\n" % (prefix, self.cnv_profile_fn)
        s += "%sclone_meta_fn = %s\n" % (prefix, self.clone_meta_fn)
        s += "%sout_dir = %s\n" % (prefix, self.out_dir)
        s += "%scell_tag = %s\n" % (prefix, self.cell_tag)
        s += "%sumi_tag = %s\n" % (prefix, self.umi_tag)
        s += "%snproc = %s\n" % (prefix, self.nproc)
        s += "%sverbose = %s\n" % (prefix, self.verbose)