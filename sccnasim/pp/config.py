# config.py


import sys



COMMAND = "pp"


class Config:
    """Configuration of the `pp` (preprocessing) module.
    
    Attributes
    ----------
    See :func:`~.main.pp_wrapper()`.
    """
    def __init__(self):
        # command-line arguments/parameters.
        self.cell_anno_fn = None
        self.feature_fn = None
        self.snp_fn = None
        self.clone_anno_fn = None
        self.cna_profile_fn = None
        self.out_dir = None
        self.chroms = ",".join([str(c) for c in range(1, 23)])
        self.strandness = "forward"
        self.overlap_features_how = "raw"

        # internal parameters.
        
        # chrom_list : list of str or None
        #   A list of chromosome names.
        #   It is used when pileup whole chromosomes.
        self.chrom_list = None

        # out_prefix_raw : str
        #   Prefix to the output raw files.
        self.out_prefix_raw = "raw"

        # out_prefix_pp : str
        #   Prefix to the output preprocess-ed files.
        self.out_prefix_pp = COMMAND

    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stdout
        
        # command-line arguments/parameters.
        s =  "%s\n" % prefix
        s += "%scell_anno_fn = %s\n" % (prefix, self.cell_anno_fn)
        s += "%sfeature_fn = %s\n" % (prefix, self.feature_fn)
        s += "%ssnp_fn = %s\n" % (prefix, self.snp_fn)
        s += "%sclone_anno_fn = %s\n" % (prefix, self.clone_anno_fn)
        s += "%scna_profile_fn = %s\n" % (prefix, self.cna_profile_fn)
        s += "%sout_dir = %s\n" % (prefix, self.out_dir)
        s += "%schroms = %s\n" % (prefix, self.chroms)
        s += "%sstrandness = %s\n" % (prefix, self.strandness)
        s += "%soverlap_features_how = %s\n" % (prefix, str(self.overlap_features_how))
        s += "%s\n" % prefix

        # internal parameters.
        s += "%sout_prefix_raw = %s\n" % (prefix, self.out_prefix_raw)
        s += "%sout_prefix_pp = %s\n" % (prefix, self.out_prefix_pp)
        s += "%s\n" % prefix

        fp.write(s)
        
    def is_stranded(self):
        return self.strandness in ("forward", "reverse")
