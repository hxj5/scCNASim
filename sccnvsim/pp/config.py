# config.py

import sys


class Config:
    """Configuration of preprocessing.

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
    """
    def __init__(self):
        self.cell_anno_fn = None
        self.feature_fn = None
        self.snp_fn = None
        self.cnv_profile_fn = None
        self.clone_meta_fn = None
        self.out_dir = None

        self.out_prefix_raw = "raw."
        self.out_prefix_pp = "pp."

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
        s += "%s\n" % prefix

        s += "%sout_prefix_raw = %s\n" % (prefix, self.out_prefix_raw)
        s += "%sout_prefix_pp = %s\n" % (prefix, self.out_prefix_pp)
        s += "%s\n" % prefix

        fp.write(s)