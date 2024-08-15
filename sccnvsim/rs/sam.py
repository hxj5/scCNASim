# sam.py

import pysam
from ..utils.grange import format_chrom
from ..utils.sam import sam_fetch, sam_merge, \
    BAM_FPAIRED, BAM_FPROPER_PAIR
from ..utils.xbarcode import Barcode


class SAMInput:
    def __init__(
        self, 
        sams, n_sam, samples, chrom,
        cell_tag, umi_tag,
        min_mapq = 20, min_len = 30,
        incl_flag = 0, excl_flag = None,
        no_orphan = True
    ):
        self.sams = sams
        if n_sam == 1:
            if not isinstance(sams, list):
                self.sams = [sams]
        else:
            assert len(sams) == n_sam
        self.n_sam = n_sam
        self.samples = samples

        self.chrom = format_chrom(chrom)

        self.cell_tag = cell_tag
        self.umi_tag = umi_tag

        self.min_mapq = min_mapq
        self.min_len = min_len
        self.incl_flag = incl_flag
        self.excl_flag = excl_flag
        self.no_orphan = no_orphan

        if not self.use_barcodes():
            assert len(self.samples) == n_sam

        self.idx = 0
        self.fp = pysam.AlignmentFile(self.sams[self.idx], "r")
        self.iter = sam_fetch(self.fp, self.chrom, None, None)

    def __fetch_read(self):
        """Fetch one read."""
        read = next(self.iter, None)
        if read is None:
            self.idx += 1
            if self.idx >= self.n_sam:
                return(None)
            self.fp.close()
            self.fp = pysam.AlignmentFile(self.sams[self.idx], "r")
            self.iter = sam_fetch(self.fp, self.chrom, None, None)
            return(self.__fetch_read())
        return(read)
    
    def check_read1(self, read):
        # partial filtering to speed up.
        if read.mapq < self.min_mapq:
            return(-2)
        if self.cell_tag and not read.has_tag(self.cell_tag):
            return(-11)
        if self.umi_tag and not read.has_tag(self.umi_tag):
            return(-12)        
        return(0)
    
    def check_read2(self, read):
        # partial filtering to speed up.
        if self.excl_flag and read.flag & self.excl_flag:
            return(-3)
        if self.incl_flag and not read.flag & self.incl_flag:
            return(-4)
        if self.no_orphan and read.flag & BAM_FPAIRED and not \
            read.flag & BAM_FPROPER_PAIR:
            return(-5)
        if len(read.positions) < self.min_len:
            return(-21)
        return(0)

    def close(self):
        if self.fp:
            self.fp.close()
        self.fp = None
        self.iter = None
    
    def fetch(self):
        while True:
            read = self.__fetch_read()
            if read is None:
                return(None)
            if self.check_read1(read) == 0:
                break
        cell = umi = None
        if self.use_barcodes():
            cell = read.get_tag(self.cell_tag)
        else:
            cell = self.samples[self.idx]
        if self.use_umi():
            umi = read.get_tag(self.umi_tag)
        else:
            umi = read.query_name
        return((read, cell, umi))

    def use_barcodes(self):
        return self.cell_tag is not None

    def use_umi(self):
        return self.umi_tag is not None


class SAMOutput:
    def __init__(
        self,
        sams, n_sam, samples, ref_sam,
        cell_tag, umi_tag, umi_len
    ):
        self.sams = sams
        if n_sam == 1:
            if not isinstance(sams, list):
                self.sams = [sams]
        else:
            assert len(sams) == n_sam
        self.n_sam = n_sam
        self.samples = samples

        self.cell_tag = cell_tag
        self.umi_tag = umi_tag
        self.umi_len = umi_len

        if not self.use_barcodes():
            assert len(self.samples) == n_sam

        in_sam = pysam.AlignmentFile(ref_sam, "r")
        self.fps = [pysam.AlignmentFile(fn, "wb", template = in_sam)  \
                    for fn in self.sams]
        in_sam.close()

        self.b = Barcode(self.umi_len)

    def close(self):
        if self.fps:
            for fp in self.fps:
                fp.close()
        self.fps = None

    def write(self, read, cell_idx, umi_int, reg_idx, qname = None):
        fp = None
        if self.use_barcodes():
            fp = self.fps[0]
            read.set_tag(self.cell_tag, self.samples[cell_idx])
        else:
            fp = self.fps[cell_idx]
        if self.use_umi():
            umi = self.b.int2str(umi_int)
            read.set_tag(self.umi_tag, umi)
        else:
            suffix = "_%d_%d" % (reg_idx, umi_int)
            if qname is None:
                read.query_name += suffix
            else:
                read.query_name = qname + suffix
        fp.write(read)

    def use_barcodes(self):
        return self.cell_tag is not None

    def use_umi(self):
        return self.umi_tag is not None

 
def sam_merge_and_index(in_fn_list, out_fn):
    sam_merge(in_fn_list, out_fn)
    pysam.index(out_fn)
