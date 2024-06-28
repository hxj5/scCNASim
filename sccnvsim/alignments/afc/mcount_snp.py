# mcount_snp.py - counting machine for SNPs.

from ...utils.sam import get_query_bases


class UCount:
    """Counting Unit.

    TODO: UMI/read collapsing.

    Attributes
    ----------
    scnt : SCount object
        A SCount object that the UCount object belongs to.
    conf : config::Config object
        Configuration.
    allele : str
        The allele for the query SNP in this UMI.
    allele_idx : int
        The index of phased allele. 
        0 (ref), 1 (alt), -1 (oth).
    """
    def __init__(self, scnt, conf):
        self.scnt = scnt
        self.conf = conf
        self.allele = None
        self.allele_idx = -2

    def push_read(self, read):
        """Push one read into this count machine.
        
        Parameters
        ----------
        read : pysam::AlignedSegment object
            A BAM read to be counted.

        Returns
        -------
        int
            0 if success, -1 otherwise.
        """
        snp = self.scnt.mcnt.snp
        try:
            idx = read.positions.index(snp.pos - 1)
        except:
            self.allele = None
            self.allele_idx = -2
        else:
            bases = get_query_bases(read, full_length = False)
            self.allele = bases[idx].upper()
            self.allele_idx = snp.get_feature_allele_index(self.allele)
        return(0)

    def stat(self):
        return(0)


class SCount:
    """Counting for single sample

    Attributes
    ----------
    mcnt : MCount object
        A MCount object that the SCount object belongs to.
    conf : config::Config object
        Configuration.
    tcount : list
        Total read / UMI counts for A/C/G/T/N bases, only for this 
        sample [list of int; 5 elements].
    umi_cnt : dict
        HashMap of <str:UCount> for umi:UCount pair, mainly for 10x data.
    is_reset : bool
        Has this object been reset.
    """
    def __init__(self, mcnt, conf):
        self.mcnt = mcnt
        self.conf = conf

        self.tcount = [0] * 5
        self.umi_cnt = {}
        self.is_reset = False

    def mark_reset_false(self):
        self.is_reset = False

    def mark_reset_true(self):
        self.is_reset = True

    def push_read(self, read):
        conf = self.conf
        umi = None
        if conf.use_umi():
            umi = read.get_tag(conf.umi_tag)
        else:
            umi = read.query_name
        if not umi:
            return(0)
        if umi in self.umi_cnt:
            return(0)
        else:
            ucnt = UCount(self, conf)
            self.umi_cnt[umi] = ucnt
            ret = ucnt.push_read(read)
            if ret < 0:
                return(-1)
            return(0)

    def reset(self):
        if self.is_reset:
            return
        for i in range(len(self.tcount)):
            self.tcount[i] = 0
        self.umi_cnt.clear()
        self.mark_reset_true()

    def stat(self):
        for ucnt in self.umi_cnt.values():
            if ucnt.stat() < 0:
                return(-1)
            allele = ucnt.allele
            if allele:
                idx = self.mcnt.base_idx["N"] 
                if allele in self.mcnt.base_idx:
                    idx = self.mcnt.base_idx[allele]
                self.tcount[idx] += 1
        return(0)


class MCount:
    """Counting for multiple samples

    Attributes
    ----------
    samples : list
        A list of cell barcodes or sample IDs [list of str].
    conf : config::Config object
        Configuration
    snp : gfeature::SNP object
        The SNP being pileuped.
    tcount : list
        Total read / UMI counts for A/C/G/T/N bases, aggregated for all 
        samples [list of int; 5 elements].
    base_idx : dict
        The mapping from base (str) to index (int) for `tcount`.
    cell_cnt : dict
        HashMap of <str, SCount> for sample:SCount pair.
    is_reset : bool
        Has this object been reset.
    """
    def __init__(self, samples, conf):
        self.samples = samples
        self.conf = conf

        self.snp = None
        self.tcount = [0] * 5
        self.base_idx = {"A":0, "C":1, "G":2, "T":3, "N":4}
        self.cell_cnt = {}
        for smp in self.samples:
            if smp in self.cell_cnt:    # duplicate samples
                return(-2)
            self.cell_cnt[smp] = SCount(self, self.conf)
        self.is_reset = False

    def add_snp(self, snp):
        self.reset()
        self.snp = snp
        self.mark_reset_false()
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
            0 if success, -1 error, -2 read filtered.
        """
        conf = self.conf
        if conf.use_barcodes():
            smp = read.get_tag(conf.cell_tag)
        else:
            smp = sid
        scnt = None
        if smp in self.cell_cnt:
            scnt = self.cell_cnt[smp]
        else:
            return(-2)

        ret = scnt.push_read(read)
        if ret < 0: 
            return(-1)
        return(0)

    def reset(self):
        if self.is_reset:
            return
        self.snp = None
        if self.tcount:
            for i in range(len(self.tcount)):
                self.tcount[i] = 0
        if self.cell_cnt:
            for scnt in self.cell_cnt.values():
                scnt.reset()
        self.mark_reset_true()

    def stat(self):
        for smp in self.samples:
            scnt = self.cell_cnt[smp]
            if scnt.stat() < 0:
                return(-1)
            for j in range(len(self.tcount)):
                self.tcount[j] += scnt.tcount[j]
        return(0)