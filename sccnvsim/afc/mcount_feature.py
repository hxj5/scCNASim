# mcount_feature.py - counting machine for features (all allele types).


class SCount:
    """Counting for single cell.

    This class counts the UMIs/reads of all haplotype states for specific 
    feature in one cell, leveraging previous counting results of mainly
    haplotype A and B.

    The haplotype state of one UMI/read pair is listed below:
    - A (Haplotype-A; internal index: 0)
        Haplotype A has supporting SNPs but haplotype B does not.
    - B (Haplotype-B; internal index: 1)
        Haplotype B has supporting SNPs but haplotype A does not.
    - D (Duplicate; internal index: 2)
        Both haplotype A and B have supporting SNPs.
    - O (Others; internal index: -1)
        Neither haplotype A nor B has supporting SNPs, but other alleles 
        (bases) in SNP level are fetched.
    - U (Unknown; internal index: -2)
        The UMI/read pair is fetched by some SNPs, but no any alleles are
        fetched.
    - U (Unknown; internal index: -3)
        The UMI/read pair is not fetched by any SNPs.
    """
    def __init__(self, mcnt, conf):
        """
        Parameters
        ----------
        mcnt : afc.mcount_feature.MCount
            A MCount object (multiple cells) that the SCount object belongs to.
        conf : afc.config.Config
            Global configuration object.
        """
        self.mcnt = mcnt
        self.conf = conf

        # ab : afc.mcount_ab.SCount
        #   The object containing the feature counting results (mainly of
        #   haplotype A and B) in this cell.
        self.ab = None

        # hap_cnt : dict of {int : int}
        #   The haplotype-specific UMI/read counts in feature level.
        #   The key is the haplotype index in UMI level, 
        #   the value is the number of UMIs/reads.
        self.hap_cnt = {0:0, 1:0, 2:0, -1:0, -2:0, -3:0}

        # umi_cnt : dict of {int : set of str}
        #   The cell-wise supporting UMIs/read pairs for each haplotype state.
        #   Keys are haplotype indexes, values are set of UMI barcodes or read
        #   query names.
        self.umi_cnt = {
            0:set(), 1:set(), 2:set(), 
            -1:set(), -2:set(), -3:set()
        }

        # is_reset : bool
        #   Whether this object has been reset?
        self.is_reset = False

    def add_ab(self, ab):
        """Add feature counting results of haplotype A and B.

        Parameters
        ----------
        ab : afc.mcount_ab.SCount
            The object containing cell-wise feature counting results of mainly
            haplotype A and B.
        
        Returns
        -------
        Void.        
        """
        self.ab = ab
        self.mark_reset_false()

    def mark_reset_false(self):
        self.is_reset = False

    def mark_reset_true(self):
        self.is_reset = True

    def push_read(self, read):
        """Count one BAM read.

        Parameters
        ----------
        read : pysam.AlignedSegment
            The BAM read to be counted.

        Returns
        -------
        int
            Return code. 0 if success, 1 if read filtered, -1 error.
        str
            UMI barcode (droplet-based) or read query name (well-based) of
            the read.
        int
            The haplotype index of this read.
        """
        conf = self.conf
        umi = None
        hap_idx = None
        if conf.use_umi():
            umi = read.get_tag(conf.umi_tag)
        else:
            umi = read.query_name
        if not umi:
            return((1, umi, hap_idx))
        if umi in self.ab.umi_cnt:
            hap_idx = self.ab.umi_cnt[umi].hap_idx
        else:
            hap_idx = -3
        self.umi_cnt[hap_idx].add(umi)
        return((0, umi, hap_idx))

    def reset(self):
        if self.is_reset:
            return
        for hap_idx in self.hap_cnt.keys():
            self.hap_cnt[hap_idx] = 0
        for hap_idx in self.umi_cnt.keys():
            self.umi_cnt[hap_idx].clear()
        self.mark_reset_true()

    def stat(self):
        for hap_idx in self.umi_cnt.keys():
            self.hap_cnt[hap_idx] = len(self.umi_cnt[hap_idx])
        return(0)


class MCount:
    """Counting for multiple cells.
    
    This class counts the UMIs/reads of all haplotype states for specific 
    feature in individual cells, leveraging previous counting results of 
    mainly haplotype A and B.

    Use the :func:`~afc.mcount_feature.MCount.add_feature` function to add
    the feature to be counted.
    """
    def __init__(self, samples, conf):
        """
        Parameters
        ----------
        samples : list of str
            A list of cell barcodes (droplet-based) or sample IDs (well-based).
        conf : afc.config.Config
            Global configuration object.
        """
        self.samples = samples
        self.conf = conf

        # reg : afc.gfeature.BlockRegion
        #   The feature to be counted.
        self.reg = None

        # ab : afc.mcount_ab.MCount
        #   The object containing (multi-cell) feature counting results of 
        #   mainly haplotype A and B.
        self.ab = None

        # cell_cnt : dict of {str : afc.mcount_feature.SCount}
        #   The cell-specific counting data.
        #   Keys are cell barcodes (droplet-based) or sample IDs (well-based),
        #   values are the associated :class:`~afc.mcount_feature.SCount`
        #   objects.
        self.cell_cnt = {}
        for smp in self.samples:
            if smp in self.cell_cnt:    # duplicate samples
                return(-2)
            self.cell_cnt[smp] = SCount(self, self.conf)

        # is_reset : bool
        #   Whether this object has been reset.
        self.is_reset = False

    def add_feature(self, reg, ab):
        """Add one feature to be counted.

        Parameters
        ----------
        reg : afc.gfeature.BlockRegion
            A feature to be counted.
        ab : afc.mcount_ab.MCount
            The object containing (multi-cell) feature counting results of 
            mainly haplotype A and B.
        
        Returns
        -------
        int
            Return code. 0 if success, negative otherwise.
        """
        self.reset()
        self.reg = reg
        self.ab = ab
        for smp in self.samples:
            scnt_ab = ab.cell_cnt[smp]
            self.cell_cnt[smp].add_ab(scnt_ab)
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
        """Count one BAM read.

        Parameters
        ----------
        read : pysam::AlignedSegment
            A BAM read to be counted.
        sid : str
            The ID of the sample that the read belongs to. 
            Set to `None` if cell barcodes are used.

        Returns
        -------
        int
            Return code. 0 if success, -1 error, 1 if read filtered.
        str
            Cell barcode (droplet-based) or sample name (well-based) of the
            read.
        str or None
            UMI barcode (droplet-based) or read query name (well-based) of
            the read.
            None if read is filtered.
        int
            The haplotype index of this read.
        """
        conf = self.conf
        smp = None
        hap_idx = None
        if conf.use_barcodes():
            smp = read.get_tag(conf.cell_tag)
        else:
            smp = sid
        scnt = None
        if smp in self.cell_cnt:
            scnt = self.cell_cnt[smp]
        else:
            return((1, smp, None, hap_idx))

        ret, umi, hap_idx = scnt.push_read(read)
        if ret < 0: 
            return((-1, smp, umi, hap_idx))
        elif ret > 0:
            return((1, smp, umi, hap_idx))
        return((0, smp, umi, hap_idx))

    def reset(self):
        if self.is_reset:
            return
        self.reg = None
        self.ab = None
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
