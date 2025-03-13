# config.py - configuration

import sys
from ..afc.config import DefaultConfig as AFC_Def_Conf

COMMAND = "rs"


class Config:
    """Configuration of the `rs` (read simulation) module.

    Attributes
    ----------
    sam_fn : str or None, default None
        Comma separated indexed BAM file.
        Note that one and only one of `sam_fn` and `sam_list_fn` should be
        specified.
    sam_list_fn : str or None, default None
        A file listing indexed BAM files, each per line.
    barcode_fn : str or None, default None
        A plain file listing all effective cell barcode.
        It should be specified for droplet-based data.
    sample_ids : str or None, default None
        Comma separated sample IDs.
        It should be specified for well-based or bulk data.
        When `barcode_fn` is not specified, the default value will be
        "SampleX", where "X" is the 0-based index of the BAM file(s).
        Note that `sample_ids` and `sample_id_fn` should not be specified
        at the same time.
    sample_id_fn : str or None, default None
        A file listing sample IDs, each per line.
    count_fn : str
        An ".adata" file storing count matrices.
        Typically it is returned by the `cs` module.
        This file should contain several layers for the allele-specific count 
        matrices.
        Its ".obs" should contain two columns "cell" and "cell_type".
        Its ".var" should contain two columns "feature" and "chrom".
    feature_fn : str
        A pickle object file storing target features.
        Typically it is returned by the `afc` module.
        This file contains a list of :class:`~afc.gfeature.BlockRegion`
        objects, whose order should be the same with the 
        ".var["feature"]" in `count_fn`.
    refseq_fn : str
        A FASTA file storing reference genome sequence.
    out_dir : str
        Output directory.
    chroms : str, default "1,2,...22"
        Comma separated chromosome names.
        Reads in other chromosomes will not be used for sampling and hence
        will not be present in the output BAM file(s).
    cell_tag : str or None, default "CB"
        Tag for cell barcodes, set to None when using sample IDs.
    umi_tag : str or None, default "UB"
        Tag for UMI, set to None when reads only.
    umi_len : int, default 10
        Length of output UMI barcode.
    ncores : int, default 1
        Number of cores.
    min_mapq : int, default 20
        Minimum MAPQ for read filtering.
    min_len : int, default 30
        Minimum mapped length for read filtering.
    min_include : int or float, default 0.9
        Minimum length of included part within specific feature.
        If float between (0, 1), it is the minimum fraction of included length.
    incl_flag : int, default 0
        Required flags: skip reads with all mask bits unset.
    excl_flag : int, default -1
        Filter flags: skip reads with any mask bits set.
        Value -1 means setting it to 772 when using UMI, or 1796 otherwise.
    no_orphan : bool, default True
        If `False`, do not skip anomalous read pairs.
    """
    def __init__(self):
        # defaults : DefaultConfig
        #   The default values of parameters.
        self.defaults = DefaultConfig()

        # argv : list of str or None, default None
        #   A list of command line arguments, typically from sys.argv.
        self.argv = None

        # command-line arguments/parameters.
        self.sam_fn = None
        self.sam_list_fn = None
        self.barcode_fn = None
        self.sample_ids = None
        self.sample_id_fn = None
        self.count_fn = None
        self.feature_fn = None
        self.refseq_fn = None
        self.out_dir = None
        self.debug = self.defaults.DEBUG

        self.chroms = self.defaults.CHROMS
        self.cell_tag = self.defaults.CELL_TAG
        self.umi_tag = self.defaults.UMI_TAG
        self.umi_len = self.defaults.UMI_LEN
        self.ncores = self.defaults.NCORES

        self.min_mapq = self.defaults.MIN_MAPQ
        self.min_len = self.defaults.MIN_LEN
        self.min_include = self.defaults.MIN_INCLUDE
        self.incl_flag = self.defaults.INCL_FLAG
        self.excl_flag = -1
        self.no_orphan = self.defaults.NO_ORPHAN

        # derived variables

        # chrom_list : list of str or None
        #   A list of chromosome names.
        #   It is used when pileup whole chromosomes.
        self.chrom_list = None

        # barcodes : list of str or None
        #   A list of cell barcodes.
        #   None if sample IDs are used.
        self.barcodes = None

        # sample_ids : list of str or None
        #   A list of sample IDs.
        #   None if cell barcodes are used.
        self.sample_ids = None

        # adata : anndata.Anndata
        #   The object storing count matrices.
        self.adata = None

        # reg_list : list of afc.gfeature.BlockRegion
        #   A list of features.
        self.reg_list = None

        # sam_fn_list : list of str
        #   A list of input SAM/BAM files.
        self.sam_fn_list = None

        # samples : list of str
        #   A list of cell barcodes (droplet-based data) or sample IDs (
        #   well-based data).
        #   It will be used as output IDs of each cell.
        self.samples = None

        # out_prefix : str
        #   The prefix of the output files.
        self.out_prefix = COMMAND + "."
        
        
        # internal parameters.

        # hap_tag : str
        #   Tag for haplotype in the output BAM file.
        self.hap_tag = "HT"
        
        # cell_raw_tag : str or None
        #   Tag for uncorrected raw cell tag in seed and simulated BAM.
        #   Set to None if do not use it.
        self.cell_raw_tag = "CR"
        
        # backup_cell_raw_tag : str or None
        #   Tag for backup uncorrected raw cell barcode (from seed BAM) 
        #   in simulated BAM. Set to None if do not use it.
        self.backup_cell_raw_tag = "RC"
        
        # umi_raw_tag : str or None
        #   Tag for uncorrected raw umi tag in seed and simulated BAM.
        #   Set to None if do not use it.
        self.umi_raw_tag = "UR"
        
        # backup_umi_raw_tag : str or None
        #   Tag for backup uncorrected raw umi barcode (from seed BAM) 
        #   in simulated BAM. Set to None if do not use it.
        self.backup_umi_raw_tag = "RU"
        
        # alleles : tuple of str
        #   All alleles.
        self.alleles = ("A", "B", "U")

        # cumi_max_pool : tuple of int
        #   The maximum size of sampling pool for each allele-specific CUMIs.
        #   Its length and order should match `alleles`.
        #   Element value 0 means no limit.
        #   One example setting for allele A,B,U is:
        #       self.cumi_max_pool = (1000, 1000, 10000)
        self.cumi_max_pool = (0, 0, 0)

        # out_sam_dir : str
        #   Output folder for SAM/BAM file(s).
        self.out_sam_dir = None

        # out_step_dir : str
        #   Output folder for step-wise results.
        self.out_step_dir = None


    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stdout

        s =  "%s\n" % prefix
        s += "%ssam_file = %s\n" % (prefix, self.sam_fn)
        s += "%ssam_list_file = %s\n" % (prefix, self.sam_list_fn)
        s += "%sbarcode_file = %s\n" % (prefix, self.barcode_fn)
        s += "%ssample_ids = %s\n" % (prefix, self.sample_ids)
        s += "%ssample_id_file = %s\n" % (prefix, self.sample_id_fn)
        s += "%scount_file = %s\n" % (prefix, self.count_fn)
        s += "%sfeature_file = %s\n" % (prefix, self.feature_fn)
        s += "%srefseq_file = %s\n" % (prefix, self.refseq_fn)
        s += "%sout_dir = %s\n" % (prefix, self.out_dir)
        s += "%sdebug_level = %d\n" % (prefix, self.debug)
        s += "%s\n" % prefix

        s += "%schroms = %s\n" % (prefix, self.chroms)
        s += "%scell_tag = %s\n" % (prefix, self.cell_tag)
        s += "%sumi_tag = %s\n" % (prefix, self.umi_tag)
        s += "%sumi_len = %s\n" % (prefix, self.umi_len)
        s += "%snumber_of_processes = %d\n" % (prefix, self.ncores)
        s += "%s\n" % prefix

        s += "%smin_mapq = %d\n" % (prefix, self.min_mapq)
        s += "%smin_len = %d\n" % (prefix, self.min_len)
        s += "%smin_include = %f\n" % (prefix, self.min_include)
        s += "%sinclude_flag = %d\n" % (prefix, self.incl_flag)
        s += "%sexclude_flag = %d\n" % (prefix, self.excl_flag)
        s += "%sno_orphan = %s\n" % (prefix, self.no_orphan)
        s += "%s\n" % prefix

        # derived variables

        s += "%schrom_list = %s\n" % (prefix, str(self.chrom_list))
        s += "%snumber_of_BAMs = %d\n" % (prefix, len(self.sam_fn_list) if \
                self.sam_fn_list is not None else -1)
        s += "%snumber_of_barcodes = %d\n" % (prefix, len(self.barcodes) if \
                self.barcodes is not None else -1)
        s += "%snumber_of_sample_IDs = %d\n" % (prefix, len(self.sample_ids) \
                if self.sample_ids is not None else -1)
        s += "%sshape_of_adata = %s\n" % (prefix, self.adata.shape if \
                self.adata is not None else "None")
        s += "%snumber_of_features = %d\n" % (prefix, len(self.reg_list) if \
                self.reg_list is not None else -1)
        s += "%s\n" % prefix
        
        # internal parameters.

        s += "%shap_tag = %s\n" % (prefix, self.hap_tag)
        s += "%scell_raw_tag = %s\n" % (prefix, self.cell_raw_tag)
        s += "%sbackup_cell_raw_tag = %s\n" % (prefix, self.backup_cell_raw_tag)
        s += "%sumi_raw_tag = %s\n" % (prefix, self.umi_raw_tag)
        s += "%sbackup_umi_raw_tag = %s\n" % (prefix, self.backup_umi_raw_tag)
        s += "%salleles = %s\n" % (prefix, str(self.alleles))
        s += "%scumi_max_pool = %s\n" % (prefix, str(self.cumi_max_pool))
        s += "%s\n" % prefix

        s += "%sout_sam_dir = %s\n" % (prefix, self.out_sam_dir)
        s += "%sout_step_dir = %s\n" % (prefix, self.out_step_dir)
        s += "%s\n" % prefix

        fp.write(s)

    def use_barcodes(self):
        return self.cell_tag is not None

    def use_umi(self):
        return self.umi_tag is not None


class DefaultConfig(AFC_Def_Conf):
    def __init__(self):
        super().__init__()
        self.UMI_LEN = 10
        self.CHROMS = ",".join([str(i) for i in range(1, 23)])


if __name__ == "__main__":
    conf = Config()
    conf.show()
