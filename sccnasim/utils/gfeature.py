# gfeature.py - genomic features, supporting interval query.

import os
from .grange import Region, RegionSet


class SNP(Region):
    """Phased SNP."""

    def __init__(self, chrom, pos, ref, alt, ref_hap, alt_hap):
        """
        Parameters
        ----------
        chrom : str
            Chromosome name.
        pos : int
            1-based genomic position.
        ref : str
            The reference (REF) base of the SNP.
        alt : str
            The alternative (ALT) base of the SNP.
        ref_hap : {0, 1}
            The haplotype index for the `ref` base.
        alt_hap : {1, 0}
            The haplotype index for the `alt` base.
        """
        super().__init__(chrom, pos, pos + 1)
        ref = ref.upper()
        alt = alt.upper()

        self.pos = pos
        self.ref = ref
        self.alt = alt
        self.ref_hap = ref_hap
        self.alt_hap = alt_hap
        self.gt = {ref:ref_hap, alt:alt_hap}
        self.hap = {ref_hap:ref, alt_hap:alt}

    def get_id(self):
        return "%s_%d" % (self.chrom, self.pos)

    def get_hap_idx(self, base):
        """Get the haplotype index of the SNP give base.
        
        Parameters
        ----------
        base : str
            The base, one of {"A", "C", "G", "T"}.

        Returns
        -------
        int
            The haplotype index.
            0 or 1 if the `base` if one of the "REF" or "ALT" base of the SNP,
            -1 otherwise.
        """
        base = base.upper()
        return self.gt[base] if base in self.gt else -1

    def get_hap_base(self, hap):
        """Get the base given the haplotype index.
        
        Parameters
        ----------
        hap : int
            The haplotype index.
            
        Returns
        -------
        str or None
            The base for the `hap`. None if `hap` is not 0 and 1.
        """
        if hap not in self.hap:
            return(None)
        return(self.hap[hap])



class SNPSet(RegionSet):
    """A set of phased SNPs."""
    def __init__(self, is_uniq = False):
        super().__init__(is_uniq)

    def add(self, snp):
        return super().add(snp)



class AlleleData(object):
    """Internal wrapper of all allelic data."""
    
    def __init__(self, allele, feature, res_dir):
        """
        Parameters
        ----------
        allele : str
            Allele name.
        feature : str
            Feature name.
        res_dir : str
            Path to the folder storing the results.
        """
        self.allele = allele
        self.feature = feature
        self.res_dir = res_dir

        # derived parameters
        
        # sam_fn : str
        #   The output alignment file.
        self.sam_fn = None
        
        # cumi_fn : str
        #   The output file containing allele-specific cell-umi barcodes.
        self.cumi_fn = None
        
        # internal parameters
        
        # out_prefix : str
        #   The prefix to the output files.
        self.out_prefix = allele + "."



class Feature(Region):
    """Feature with covered SNPs."""

    def __init__(self, chrom, start, end, name, strand, 
                 snp_list = None, res_dir = None):
        """
        Parameters
        ----------
        chrom : str
            Chromosome name.
        start : int
            1-based genomic start position of the feature, inclusive.
        end : int
            1-based genomic end position of the feature, exclusive.
        name : str
            Name of the feature.
        strand : str
            DNA strand orientation of the feature, "+" (positive) or 
            "-" (negative).
        snp_list : list of utils.gfeature.SNP
            A list of SNPs located within the feature.
        res_dir : str
            Path to the folder storing the results of this feature.
        """
        super().__init__(chrom, start, end)
        self.name = name
        self.strand = strand
        self.snp_list = snp_list
        self.res_dir = res_dir

        # derived parameters
        
        # allele_data : dict of {str : str}
        #   Keys are the alleles, values are the `AlleleData` objects.
        self.allele_data = None
        
        
    def init_allele_data(self, alleles):
        """Initialize the allele-specific data.
        
        Parameters
        ----------
        alleles : list of str
            A list of alleles.
            
        Returns
        -------
        Void.
        """
        assert self.allele_data is None
        self.allele_data = {}
        for ale in alleles:
            ale_dir = os.path.join(self.res_dir, ale)
            os.makedirs(ale_dir, exist_ok = True)
            ale_data = AlleleData(
                allele = ale,
                feature = self.name,
                res_dir = ale_dir
            )
            self.allele_data[ale] = ale_data
            ale_data.sam_fn = os.path.join(ale_data.res_dir, \
                ale_data.out_prefix + "bam")
            ale_data.cumi_fn = os.path.join(ale_data.res_dir, \
                ale_data.out_prefix + "cumi.tsv")
            


def assign_feature_batch(feature_names, root_dir, batch_size = 1000):
    """Assign features into several batches.

    This function assign features into several batches of result folders, 
    to avoid exceeding the maximum number of files/sub-folders in one folder.
    
    Parameters
    ----------
    feature_names : list of str
        A list of feature names.
    root_dir : str
        The root folder, under which batch subfolders will be created.
    batch_size : int, default 1000
        Number of features in each batch.
    
    Returns
    -------
    list of str
        A list of paths to result dir of every feature.
    """
    batch_idx = -1
    batch_dir = None
    feature_dirs = []
    for fet_idx, fet_name in enumerate(feature_names):
        if fet_idx % batch_size == 0:
            batch_idx += 1
            batch_dir = os.path.join(root_dir, "b%d" % batch_idx)
            os.makedirs(batch_dir, exist_ok = True)
        fet_dir = os.path.join(
            batch_dir, "%d_%s" % (fet_idx, fet_name))
        os.makedirs(fet_dir, exist_ok = True)
        feature_dirs.append(fet_dir)
    return(feature_dirs)
