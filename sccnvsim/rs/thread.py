# thread.py - wrapper of data for each thread.

import sys

class ThreadData:
    """Thread Data"""
    def __init__(self, 
        idx, conf,
        chrom,
        reg_obj_fn,
        reg_idx0,
        cumi_fn,
        out_samples,
        out_sam_fn_list
    ):
        self.idx = idx
        self.conf = conf

        self.chrom = chrom
        self.reg_obj_fn = reg_obj_fn
        self.reg_idx0 = reg_idx0
        self.cumi_fn = cumi_fn

        self.out_samples = out_samples
        self.out_sam_fn_list = out_sam_fn_list
        
        self.ret = -1

    def destroy(self):
        pass

    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stdout
        
        s = "%s\n" % prefix
        s += "%sindex = %d\n" % (prefix, self.idx)

        s += "%schrom = %s\n" % (prefix, self.chrom)
        s += "%sreg_obj filename = %s\n" % (prefix, self.reg_obj_fn)
        s += "%sreg_idx0 = %s\n" % (prefix, self.reg_idx0)
        s += "%sCUMI filename = %s\n" % (prefix, self.cumi_fn)

        s += "%s#out_samples = %d\n" % (prefix, len(self.out_samples))
        s += "%s#out_sam_fn_list = %d\n" % (prefix, len(self.out_sam_fn_list))

        s += "%sreturn code = %d\n" % (prefix, self.ret)
        s += "%s\n" % prefix

        fp.write(s)
