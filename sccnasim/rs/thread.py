# thread.py - wrapper of data for each thread.

import sys

class ThreadData:
    """Thread Data for processing chrom-specific reads."""

    def __init__(self, 
        idx,
        conf,
        reg_obj_fn,
        out_samples,
        out_sam_fn
    ):
        """
        Parameters
        ----------
        idx : int
            The 0-based index of thread.
        conf : rs.config.Config
            The configuration object.
        reg_obj_fn : str
            Path to the file storing a list of feature objects, i.e., the
            :class:`~utils.gfeature.Feature` objects.
        out_samples : list of str
            Output cell barcodes (droplet-based platform) or sample IDs (
            well-based platform).
        out_sam_fn : str
            Output SAM/BAM file.
        """
        self.idx = idx
        self.conf = conf

        self.reg_obj_fn = reg_obj_fn

        self.out_samples = out_samples
        self.out_sam_fn = out_sam_fn
        
        # ret : int
        #   Return code of this thread. 0 if success, negative otherwise.
        self.ret = -1

    def destroy(self):
        pass

    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stdout
        
        s = "%s\n" % prefix
        s += "%sindex = %d\n" % (prefix, self.idx)

        s += "%sreg_obj filename = %s\n" % (prefix, self.reg_obj_fn)

        s += "%s#out_samples = %d\n" % (prefix, len(self.out_samples))
        s += "%sout_sam_fn = %s\n" % (prefix, self.out_sam_fn)

        s += "%sreturn code = %d\n" % (prefix, self.ret)
        s += "%s\n" % prefix

        fp.write(s)
