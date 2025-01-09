# config.py

import sys


# TODO:
# 1. Add new options for `merge_features_how`
#    "smallest" or "largest" - only keep the smallest/largest feature.
#        Only keep the smallest/largest of a group of consecutively 
#        overlapping features.
class Config:
    """Configuration of the `pp` (preprocessing) module.
    
    Attributes
    ----------
    cell_anno_fn : str
        The cell annotation file. 
        It is header-free and its first two columns are:
        - "cell" (str): cell barcodes.
        - "cell_type" (str): cell type.
    feature_fn : str
        A TSV file listing target features. 
        It is header-free and its first 4 columns shoud be: 
        - "chrom" (str): chromosome name of the feature.
        - "start" (int): start genomic position of the feature, 1-based
          and inclusive.
        - "end" (int): end genomic position of the feature, 1-based and
          inclusive.
        - "feature" (str): feature name.
    snp_fn : str
        A TSV or VCF file listing phased SNPs.
        If TSV, it is a header-free file containing SNP annotations, whose
        first six columns should be:
        - "chrom" (str): chromosome name of the SNP.
        - "pos" (int): genomic position of the SNP, 1-based.
        - "ref" (str): the reference allele of the SNP.
        - "alt" (str): the alternative allele of the SNP.
        - "ref_hap" (int): the haplotype index of `ref`, one of {0, 1}.
        - "alt_hap" (int): the haplotype index of `alt`, one of {1, 0}.
        If VCF, it should contain "GT" in its "FORMAT" field.
    clone_meta_fn : str
        A TSV file listing clonal meta information.
        It is header-free and its first 3 columns are:
        - "clone" (str): clone ID.
        - "source_cell_type" (str): the source cell type of `clone`.
        - "n_cell" (int): number of cells in the `clone`. If negative, 
          then it will be set as the number of cells in `source_cell_type`.
    cna_profile_fn : str
        A TSV file listing clonal CNA profiles. 
        It is header-free and its first 7 columns are:
        - "chrom" (str): chromosome name of the CNA region.
        - "start" (int): start genomic position of the CNA region, 1-based
          and inclusive.
        - "end" (int): end genomic position of the CNA region, 1-based and
          inclusive.
        - "region" (str): ID of the CNA region.
        - "clone" (str): clone ID.
        - "cn_ale0" (int): copy number of the first allele.
        - "cn_ale1" (int): copy number of the second allele.
    out_dir : str
        The output folder.
    chroms : str, default "1,2,...22"
        Comma separated chromosome names.
    merge_features_how : str, default "none"
        How to merge overlapping features.
        "none" - do not merge overlapping features.
        "bidel" - remove overlapping bi-features.
        "first1" - only keep the first feature.
            Only keep the first of the consecutively overlapping features.
        "first2" - only keep the first feature.
            Keep the first feature and remove features overlapping with it.
        "quantile1" - remove outliers of bi-features.
            Remove outliers given specific quantile among all features.
        "quantile1_union" - "quantile1" followed by "union".
        "quantile2" - remove outliers of bi-features.
            Remove outliers given specific quantile among all features 
            overlapping with at least one features.
        "quantile2_union" - "quantile2" followed by "union".
        "union" - keep the union range.
            Keep the union genomic range of a group of consecutively
            overlapping features.
    """
    def __init__(self):
        # command-line arguments/parameters.
        self.cell_anno_fn = None
        self.feature_fn = None
        self.snp_fn = None
        self.clone_meta_fn = None
        self.cna_profile_fn = None
        self.out_dir = None
        self.chroms = ",".join([str(c) for c in range(1, 23)])
        self.merge_features_how = "none"

        # derived parameters.
        
        # chrom_list : list of str or None
        #   A list of chromosome names.
        #   It is used when pileup whole chromosomes.
        self.chrom_list = None

        # out_prefix_raw : str
        #   Prefix to the output raw files.
        self.out_prefix_raw = "raw."

        # out_prefix_pp : str
        #   Prefix to the output preprocess-ed files.
        self.out_prefix_pp = "pp."

    def show(self, fp = None, prefix = ""):
        if fp is None:
            fp = sys.stdout
        
        # command-line arguments/parameters.
        s =  "%s\n" % prefix
        s += "%scell_anno_fn = %s\n" % (prefix, self.cell_anno_fn)
        s += "%sfeature_fn = %s\n" % (prefix, self.feature_fn)
        s += "%ssnp_fn = %s\n" % (prefix, self.snp_fn)
        s += "%sclone_meta_fn = %s\n" % (prefix, self.clone_meta_fn)
        s += "%scna_profile_fn = %s\n" % (prefix, self.cna_profile_fn)
        s += "%sout_dir = %s\n" % (prefix, self.out_dir)
        s += "%schroms = %s\n" % (prefix, self.chroms)
        s += "%s\n" % prefix

        # derived parameters.
        s += "%schrom_list = %s\n" % (prefix, str(self.chrom_list))
        s += "%sout_prefix_raw = %s\n" % (prefix, self.out_prefix_raw)
        s += "%sout_prefix_pp = %s\n" % (prefix, self.out_prefix_pp)
        s += "%smerge_features_how = %s\n" % (prefix, str(self.merge_features_how))
        s += "%s\n" % prefix

        fp.write(s)
