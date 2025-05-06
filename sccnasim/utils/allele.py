# allele.py - genomic alleles.



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
        
        # seed_sam_fn : str
        #   The alignment file of seed data.
        self.seed_sam_fn = None
        
        # seed_cumi_fn : str
        #   The file containing allele-specific cell-umi barcodes of 
        #   seed data.
        self.seed_cumi_fn = None
        
        
        # seed_smpl_cumi_fn : str
        #   The file containing sampled seed CUMIs.
        #   Its length and order would match `simu_cumi_fn`.
        self.seed_smpl_cumi_fn = None
        

        # simu_sam_fn : str
        #   The alignment file of simulated data.
        self.simu_sam_fn = None
        
        # simu_cumi_fn : str
        #   The file containing allele-specific cell-umi barcodes of 
        #   simulated data.
        self.simu_cumi_fn = None
        
        
        # internal parameters
        
        # out_prefix : str
        #   The prefix to the output files.
        self.out_prefix = allele + "."
