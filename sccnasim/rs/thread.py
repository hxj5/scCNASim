# thread.py - wrapper of data for each thread.


import sys

class BatchData:
    """Batch data."""

    def __init__(self, 
        idx,
        conf,
        reg_obj_fn,
        reg_idx_b,
        reg_idx_e,
        alleles,
        refseq_fn,
        tmp_dir
    ):
        """
        Parameters
        ----------
        idx : int
            The 0-based index of batch.
        conf : rs.config.Config
            The configuration object.
        reg_obj_fn : str
            Path to the file storing a list of feature objects, i.e., the
            :class:`~utils.gfeature.Feature` objects.
        reg_idx_b : int
            The index (within transcriptomics scale) of the first feature
            in this batch. 0-based and inclusive.
        reg_idx_e : int
            The index (within transcriptomics scale) of the last feature
            in this batch. 0-based and exclusive.
        alleles : list of str
            A list of alleles.
        refseq_fn : str
            A FASTA file storing reference genome sequence.
        tmp_dir : str
            Path to the folder storing temporary files.
        """
        self.idx = idx
        self.conf = conf

        self.reg_obj_fn = reg_obj_fn
        self.reg_idx_b = reg_idx_b
        self.reg_idx_e = reg_idx_e
        
        self.alleles = alleles
        
        self.refseq_fn = refseq_fn

        self.tmp_dir = tmp_dir
        
        # ret : int
        #   Return code of this batch. 0 if success, negative otherwise.
        self.ret = -1
        
        # reg_sam_fns : list of str
        #   A list of simulated feature-specific BAM files.
        self.reg_sam_fns = []

    def destroy(self):
        pass

    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stdout
        
        s = "%s\n" % prefix
        s += "%sindex = %d\n" % (prefix, self.idx)

        s += "%sreg_obj filename = %s\n" % (prefix, self.reg_obj_fn)
        s += "%sreg_idx_b = %d\n" % (prefix, self.reg_idx_b)
        s += "%sreg_idx_e = %d\n" % (prefix, self.reg_idx_e)
        
        s += "%salleles = %s\n" % (prefix, self.alleles)
        
        s += "%srefseq_fn = %s\n" % (prefix, self.refseq_fn)

        s += "%stmp_dir = %s\n" % (prefix, self.tmp_dir)

        s += "%sreturn code = %d\n" % (prefix, self.ret)
        s += "%s#reg_sam_fns = %d\n" % (prefix, len(self.reg_sam_fns))
        s += "%s\n" % prefix

        fp.write(s)
