# sam.py

import pysam
from ..utils.grange import format_chrom
from ..utils.sam import sam_fetch, sam_merge, \
    BAM_FPAIRED, BAM_FPROPER_PAIR
from ..utils.xbarcode import Barcode


class SAMInput:
    """Input SAM/BAM file object."""

    def __init__(
        self, 
        sams, n_sam, samples, chrom,
        cell_tag, umi_tag,
        min_mapq = 20, min_len = 30,
        incl_flag = 0, excl_flag = None,
        no_orphan = True
    ):
        """
        Parameters
        ----------
        sams : str or list of str
            The input SAM/BAM file(s).
        n_sam : int
            Number of input SAM/BAM file(s).
        samples : list of str
            A list of cell barcodes (droplet-based platform) or sample IDs
            (well-based platform).
        chrom : str
            The chromosome name.
        cell_tag : str or None
            Tag for cell barcodes, set to None when using sample IDs.
        umi_tag : str or None
            Tag for UMI, set to None when reads only.
        min_mapq : int, default 20
            Minimum MAPQ for read filtering.
        min_len : int, default 30
            Minimum mapped length for read filtering.
        incl_flag : int, default 0
            Required flags: skip reads with all mask bits unset.
        excl_flag : int
            Filter flags: skip reads with any mask bits set.
        no_orphan : bool, default True
            If `False`, do not skip anomalous read pairs.
        """
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

        assert excl_flag is not None
        self.min_mapq = min_mapq
        self.min_len = min_len
        self.incl_flag = incl_flag
        self.excl_flag = excl_flag
        self.no_orphan = no_orphan

        if not self.use_barcodes():
            assert len(self.samples) == n_sam

        # idx : int
        #   The index of currently iterated BAM file(s), 0-based.
        self.idx = 0

        # fp : pysam.AlignmentFile
        #   The SAM/BAM file object.
        self.fp = pysam.AlignmentFile(
            self.sams[self.idx], "r", require_index = True)

        # iter : pysam.IteratorRow
        #   An iterator over a collection of chrom-specific reads.
        self.iter = sam_fetch(self.fp, self.chrom, None, None)

    def __fetch_read(self):
        """Fetch one read."""
        read = next(self.iter, None)
        if read is None:
            self.idx += 1
            if self.idx >= self.n_sam:
                return(None)
            self.fp.close()
            self.fp = pysam.AlignmentFile(
                self.sams[self.idx], "r", require_index = True)
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
        """Fetch one valid read.
        
        This function iterately fetch reads until some valid read is fetched,
        i.e., this read passes all filtering.

        Returns
        -------
        tuple of (pysam.AlignedSegment, str, str) or None
            The fetched valid read data, including
            - (pysam.AlignedSegment) the fetched valid read.
            - (str) the cell barcode (droplet) or sample ID (well) of the read.
            - (str) the UMI barcode (droplet) or query name (well) of the read.
            None if reaching the end of all input SAM/BAM file(s).
        """
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
    """Output SAM object."""

    def __init__(
        self,
        sams, n_sam, samples, ref_sam,
        cell_tag, umi_tag, umi_len
    ):
        """
        Parameters
        ----------
        sams : list of str
            The output SAM/BAM file(s).
        n_sam : int
            Number of output SAM/BAM file(s).
        samples : list of str
            Output cell barcodes (droplet-based platform) or sample IDs (
            well-based platform).
        ref_sam : str
            The reference SAM/BAM file used as a template for output file(s).
        cell_tag : str or None
            Tag for cell barcodes, set to None when using sample IDs.
        umi_tag : str or None
            Tag for UMI, set to None when reads only.
        umi_len : int
            Length of the UMI barcode.
        """
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

        # fps : list of pysam.AlignmentFile
        #   A list of file objects for output SAM/BAM files.
        in_sam = pysam.AlignmentFile(ref_sam, "r")
        self.fps = [pysam.AlignmentFile(fn, "wb", template = in_sam)  \
                    for fn in self.sams]
        in_sam.close()

        # b : utils.xbarcode.Barcode
        #   The object used for transformming the integer UMIs into string
        #   format.
        self.b = Barcode(self.umi_len)

    def close(self):
        if self.fps:
            for fp in self.fps:
                fp.close()
        self.fps = None

    def write(self, read, cell_idx, umi_int, reg_idx, idx, qname = None):
        """Write one read into output SAM/BAM file.

        Parameters
        ----------
        read : pysam.AlignedSegment
            The read to be write.
        cell_idx : int
            The 0-based index of the new cell/sample within `samples`.
            The cell barcode or sample ID is part of the new CUMI.
        umi_int : int
            The int format of the new UMI barcode, which is part of the 
            new CUMI.
        reg_idx : int
            The 0-based index of the feature within transcriptomics-scale.
        idx : int
            0-based index of the new CUMI within all new CUMIs assigned to
            one old sampled CUMI.
        qname : str
            The query name of the read.

        Returns
        -------
        Void.
        """
        fp = None
        if self.use_barcodes():
            fp = self.fps[0]
            read.set_tag(self.cell_tag, self.samples[cell_idx])
        else:
            fp = self.fps[cell_idx]
        if self.use_umi():
            umi = self.b.int2str(umi_int)
            read.set_tag(self.umi_tag, umi)
        # suffix = "_%d_%d" % (idx, umi_int)
        suffix = "_%d" % (idx, )
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
