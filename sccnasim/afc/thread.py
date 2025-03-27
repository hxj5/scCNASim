# thread.py - wrapper of data for each thread.


import sys

class ThreadData:
    """Thread Data"""
    def __init__(self, 
        idx, conf, 
        reg_obj_fn,
        out_feature_fn,
        out_ale_fns
    ):
        """
        Parameters
        ----------
        idx : int
            The 0-based index of thread.
        conf : afc.config.Config
            The global configuration.
        reg_obj_fn : str
            Path to the python pickle file storing a list of features (
            :class:`~utils.gfeature.Feature` objects).
        out_feature_fn : str
            Path to the output feature TSV file in this thread.
        out_ale_fns : dict of {str : str}
            Path to the allele-specific count matrix file in this read.
            Keys are allele names, values are pathes to the count matrix file.
        """
        self.idx = idx
        self.conf = conf

        self.reg_obj_fn = reg_obj_fn

        self.out_feature_fn = out_feature_fn
        self.out_ale_fns = out_ale_fns

        # nr_reg : int
        #   Number of unique features (row indexes) in the output count matrix
        #   file.
        self.nr_reg = 0

        # nr_ale : dict of {str : int}
        #   Number of allele-specific records in the output count matrix
        #   file(s).
        #   Keys are allele names, values are number of records.
        self.nr_ale = {ale:0 for ale in out_ale_fns.keys()}
        
        # ret : int
        #   Return code. 0 if success, negative otherwise.
        self.ret = -1

    def destroy(self):
        pass

    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stdout
        
        s = "%s\n" % prefix
        s += "%sindex = %d\n" % (prefix, self.idx)

        s += "%sreg_obj filename = %s\n" % (prefix, self.reg_obj_fn)

        s += "%sout_feature_fn = %s\n" % (prefix, self.out_feature_fn)
        for ale, fn in self.out_ale_fns.items():
            s += "%sout_ale_%s_fn = %s\n" % (prefix, ale, fn)
        
        s += "%snum_record_feature = %d\n" % (prefix, self.nr_reg)
        for ale, nr in self.nr_ale.items():
            s += "%snum_record_%s = %d\n" % (prefix, ale, nr)

        s += "%sreturn code = %d\n" % (prefix, self.ret)
        s += "%s\n" % prefix

        fp.write(s)
