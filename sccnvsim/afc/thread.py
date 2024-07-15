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
        self.idx = idx
        self.conf = conf

        self.reg_obj_fn = reg_obj_fn

        self.out_feature_fn = out_feature_fn
        self.out_ale_fns = out_ale_fns

        self.nr_reg = 0
        self.nr_ale = {ale:0 for ale in out_ale_fns.keys()}
        
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
