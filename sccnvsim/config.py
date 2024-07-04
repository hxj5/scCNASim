# config.py - global configuration.

from .aln.afc.config import Config as Raw_AFCConf


class Config:
    def __init__(self):
        self.g = GlobalConfig()
        self.afc = AFCConfig()
        self.cs = None
        self.rs = None


# only list selected attributes here.
class AFCConfig(Raw_AFCConf):
    """Configuration of allele-specific feature counting.

    Attributes
    ----------
    sam_fn : str
        Comma separated indexed sam/bam/cram file.
    sam_list_fn : str
        A list file containing bam files, each per line.
    feature_fn : str
        A TSV file listing target features. The first 4 columns shoud be:
        chrom, start, end (both 1-based and inclusive), name.
    ...
    """
    def __init__(self):
        super().__init__()


class GlobalConfig:
    def __init__(self):
        self.out_dir = None