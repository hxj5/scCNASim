# mcount_feature_all.py - counting machine for features (all allele types).


class SCount:
    """Counting for single sample

    Attributes
    ----------
    mcnt : MCount object
        A MCount object that the SCount object belongs to.
    conf : config::Config object
        Configuration.
    prev : mcount_feature::SCount object
        It contains the previous allele-specific feature counting results.
    allele_cnt : dict
        The key is the allele index in UMI level, i.e., 0 (ref), 1 (alt),
        2 (both), -1 (oth), -2 (unknown; covered SNPs without fetched alleles),
        -3 (unknown; no any SNPs covered), the value is the number of UMIs.
    umi_cnt : dict
        HashMap of <int:str> for allele_idx:UMI_set pairs.
    is_reset : bool
        Has this object been reset.
    """
    def __init__(self, mcnt, conf):
        self.mcnt = mcnt
        self.conf = conf
        self.prev = None

        self.allele_cnt = {0:0, 1:0, 2:0, -1:0, -2:0, -3:0}
        self.umi_cnt = {0:set(), 1:set(), 2:set(), -1:set(), -2:set(), -3:set()}
        self.is_reset = False

    def add_prev(self, prev):
        self.prev = prev
        self.mark_reset_false()

    def mark_reset_false(self):
        self.is_reset = False

    def mark_reset_true(self):
        self.is_reset = True

    def push_read(self, read):
        """
        Returns
        -------
        int
            0 if success, 1 if invalid UMI/read barcode, -1 error.
        str
            UMI or read barcode. Could be empty or `None`.
        int
            The allele index. One of -3, -2, -1, 0, 1, 2.
        """
        conf = self.conf
        umi = None
        ale_idx = None
        if conf.use_umi():
            umi = read.get_tag(conf.umi_tag)
        else:
            umi = read.query_name
        if not umi:
            return((1, umi, ale_idx))
        if umi in self.prev.umi_cnt:
            ale_idx = self.prev.umi_cnt[umi].allele_idx
        else:
            ale_idx = -3
        self.umi_cnt[ale_idx].add(umi)
        return((0, umi, ale_idx))

    def reset(self):
        if self.is_reset:
            return
        for ale_idx in self.allele_cnt.keys():
            self.allele_cnt[ale_idx] = 0
        for ale_idx in self.umi_cnt.keys():
            self.umi_cnt[ale_idx].clear()
        self.mark_reset_true()

    def stat(self):
        for ale_idx in self.umi_cnt.keys():
            self.allele_cnt[ale_idx] = len(self.umi_cnt[ale_idx])
        return(0)


class MCount:
    """Counting for multiple samples

    Attributes
    ----------
    samples : list
        A list of cell barcodes or sample IDs [list of str].
    conf : config::Config object
        Configuration
    reg : gfeature::BlockRegion object
        The region in which feature counting is done.
    prev : mcount_feature::MCount object
        It contains the previous allele-specific feature counting results.
    cell_cnt : dict
        HashMap of <str, SCount> for sample:SCount pair.
    is_reset : bool
        Has this object been reset.
    """
    def __init__(self, samples, conf):
        self.samples = samples
        self.conf = conf

        self.reg = None
        self.prev = None
        self.cell_cnt = {}
        for smp in self.samples:
            if smp in self.cell_cnt:    # duplicate samples
                return(-2)
            self.cell_cnt[smp] = SCount(self, self.conf)
        self.is_reset = False

    def add_feature(self, reg, prev):
        self.reset()
        self.reg = reg
        self.prev = prev
        for smp in self.samples:
            scnt_prev = prev.cell_cnt[smp]
            self.cell_cnt[smp].add_prev(scnt_prev)
        self.mark_reset_false(recursive = False)
        return(0)
    
    def mark_reset_false(self, recursive = True):
        self.is_reset = False
        if recursive:
            for scnt in self.cell_cnt.values():
                scnt.mark_reset_false()

    def mark_reset_true(self, recursive = False):
        self.is_reset = True
        if recursive:
            for scnt in self.cell_cnt.values():
                scnt.mark_reset_true()

    def push_read(self, read, sid = None):
        """Push one read into this counting machine.

        Parameters
        ----------
        read : pysam::AlignedSegment object
            A BAM read to be counted.
        sid : str
            The ID of the sample that the read belongs to. 
            Set to `None` if cell barcodes are used.

        Returns
        -------
        int
            0 if success, -1 error, 1 if read filtered.
        str
            Cell barcode or sample name. Could be empty or `None`.
        str
            UMI or read barcode. Could be empty or `None`.
        int
            The allele index. One of -3, -2, -1, 0, 1, 2.
        """
        conf = self.conf
        smp = None
        ale_idx = None
        if conf.use_barcodes():
            smp = read.get_tag(conf.cell_tag)
        else:
            smp = sid
        scnt = None
        if smp in self.cell_cnt:
            scnt = self.cell_cnt[smp]
        else:
            return((1, smp, None, ale_idx))

        ret, umi, ale_idx = scnt.push_read(read)
        if ret < 0: 
            return((-1, smp, umi, ale_idx))
        elif ret > 0:
            return((1, smp, umi, ale_idx))
        return((0, smp, umi, ale_idx))

    def reset(self):
        if self.is_reset:
            return
        self.reg = None
        self.prev = None
        if self.cell_cnt:
            for scnt in self.cell_cnt.values():
                scnt.reset()
        self.mark_reset_true()

    def stat(self):
        for smp in self.samples:
            scnt = self.cell_cnt[smp]
            if scnt.stat() < 0:
                return(-1)
        return(0)