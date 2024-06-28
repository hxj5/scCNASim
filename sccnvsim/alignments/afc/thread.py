# thread.py - wrapper of data for each thread.

import sys

class ThreadData:
    """Thread Data"""
    def __init__(self, 
        idx, conf, 
        reg_obj_fn,
        out_feature_fn,
        out_ale_a_fn, out_ale_b_fn, out_ale_o_fn, out_ale_u_fn,
        out_fn = None
    ):
        self.idx = idx
        self.conf = conf

        self.reg_obj_fn = reg_obj_fn

        self.out_feature_fn = out_feature_fn
        self.out_ale_a_fn = out_ale_a_fn
        self.out_ale_b_fn = out_ale_b_fn
        self.out_ale_o_fn = out_ale_o_fn
        self.out_ale_u_fn = out_ale_u_fn

        self.out_fn = out_fn

        self.nr_reg = 0
        self.nr_a = 0
        self.nr_b = 0
        self.nr_o = 0
        self.nr_u = 0
        
        self.ret = -1

    def destroy(self):
        pass

    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stderr
        
        s = "%s\n" % prefix
        s += "%sindex = %d\n" % (prefix, self.idx)

        s += "%sreg_obj filename = %s\n" % (prefix, self.reg_obj_fn)

        s += "%sout_feature_fn = %s\n" % (prefix, self.out_feature_fn)
        s += "%sout_ale_a_fn = %s\n" % (prefix, self.out_ale_a_fn)
        s += "%sout_ale_b_fn = %s\n" % (prefix, self.out_ale_b_fn)
        s += "%sout_ale_o_fn = %s\n" % (prefix, self.out_ale_o_fn)
        s += "%sout_ale_u_fn = %s\n" % (prefix, self.out_ale_u_fn)

        s += "%sout_fn = %s\n" % (prefix, self.out_fn)
        
        s += "%snum_record_feature = %d\n" % (prefix, self.nr_reg)
        s += "%snum_record_a = %d\n" % (prefix, self.nr_a)
        s += "%snum_record_b = %d\n" % (prefix, self.nr_b)
        s += "%snum_record_o = %d\n" % (prefix, self.nr_o)
        s += "%snum_record_u = %d\n" % (prefix, self.nr_u)

        s += "%sreturn code = %d\n" % (prefix, self.ret)
        s += "%s\n" % prefix

        fp.write(s)
