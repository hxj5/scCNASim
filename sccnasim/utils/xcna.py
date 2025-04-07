# xcna.py - clonal CNA profile.


from logging import error
from .grange import Region, RegionSet
from .zfile import zopen
from ..io.base import load_cnas



class CNARegCN(Region):
    """Allele-specific copy numbers of CNA region."""

    def __init__(self, chrom, start, end, name, cn_ale0, cn_ale1):
        """
        Parameters
        ----------
        chrom : str
            Chromosome name.
        start : int
            The 1-based genomic start pos, inclusive.
        end : int
            The 1-based genomic end pos, exclusive.
        name : str
            The ID of the region.
        cn_ale0 : int
            Copy Number of the first allele.
        cn_ale1 : int
            Copy Number of the second allele.
        """
        super().__init__(chrom, start, end, name)
        self.name = name
        self.cn_ale0 = cn_ale0
        self.cn_ale1 = cn_ale1


        
class CNAProfile:
    """CNA profile of one clone."""

    def __init__(self):
        # rs : afc.grange.RegionSet
        #   The CNA regions associated with this clone.
        self.rs = RegionSet()

    def add_cna(self, chrom, start, end, name, cn_ale0, cn_ale1):
        """Add a new record of CNA profile.

        Parameters
        ----------
        chrom : str
            Chromosome name.
        start : int
            The 1-based genomic start pos, inclusive.
        end : int
            The 1-based genomic end pos, exclusive.
        name : str
            The ID of the region.
        cn_ale0 : int
            Copy Number of the first allele.
        cn_ale1 : int
            Copy Number of the second allele.

        Returns
        -------
        int
            Return code.
            0 if success, 1 discarded as duplicate region, -1 error.
        """
        reg = CNARegCN(chrom, start, end, name, cn_ale0, cn_ale1)
        ret = self.rs.add(reg)
        return(ret)

    def fetch(self, chrom, start, end):
        """Get the CNA profile records for the query region.

        This function returns the CNA profile records of the overlapping
        regions in this clone to the query region.

        Parameters
        ----------
        chrom : str
            Chromosome name.
        start : int
            The 1-based genomic start pos, inclusive.
        end : int
            The 1-based genomic end pos, exclusive.

        Returns
        -------
        int
            Number of overlapping regions, -1 if error.
        list
            A list of tuples of CNA profiles. Each tuple contains the 
            region-specific CNA profile:
            - int
                Copy numbers of the first allele.
            - int
                Copy numbers  of the first allele.
            - str
                The region ID.
        """
        hits = self.rs.fetch(chrom, start, end)
        res = [(reg.cn_ale0, reg.cn_ale1, reg.name) for reg in hits]
        return((len(res), res))

    def get_all(self):
        """Returns the whole CNA profile of this clone."""

        reg_list = self.rs.get_regions(sort = True)
        dat = {
                "chrom":[],
                "start":[],
                "end":[],
                "name":[],
                "cn_ale0":[],
                "cn_ale1":[]
            }
        for reg in reg_list:
            dat["chrom"].append(reg.chrom)
            dat["start"].append(reg.start)
            dat["end"].append(reg.end - 1)
            dat["name"].append(reg.name)
            dat["cn_ale0"].append(reg.cn_ale0)
            dat["cn_ale1"].append(reg.cn_ale1)
        return(dat)     
        
    def query(self, name):
        """Return the CNA profile records given the ID of one query region.
        
        Parameters
        ----------
        name : str
            The ID of the query region.

        Returns
        -------
        int
            Number of regions whose ID is `name`.
        list
            A list of tuples of CNA profiles. Each tuple contains the 
            region-specific CNA profile:
            - int
                Copy numbers of the first allele.
            - int
                Copy numbers  of the first allele.
            - str
                The region ID.
        """
        hits = self.rs.query(name)
        res = [(reg.cn_ale0, reg.cn_ale1, reg.name) for reg in hits]
        return((len(res), res))


    
class CloneCNAProfile:
    """CNA profiles of all clones."""

    def __init__(self):
        # dat : dict of {str : utils.xcna.CNAProfile}
        #   The CNA profiles of all clones.
        #   Keys are clone IDs, values are clone-specific CNA profiles.
        self.dat = {}

    def add_cna(self, chrom, start, end, name, cn_ale0, cn_ale1, clone_id):
        """Add a new record of CNA profile.

        Parameters
        ----------
        chrom : str
            Chromosome name.
        start : int
            The 1-based genomic start pos, inclusive.
        end : int
            The 1-based genomic end pos, exclusive.
        name : str
            The ID of the region.
        cn_ale0 : int
            Copy Number of the first allele.
        cn_ale1 : int
            Copy Number of the second allele.
        clone_id : str
            Clone ID.

        Returns
        -------
        int
            Return code.
            0 if success, 1 discarded as duplicate region, -1 error.
        """
        if clone_id not in self.dat:
            self.dat[clone_id] = CNAProfile()
        cp = self.dat[clone_id]
        ret = cp.add_cna(chrom, start, end, name, cn_ale0, cn_ale1)
        return(ret)

    def fetch(self, chrom, start, end, clone_id):
        """Get the CNA profile records for the query region and clone.

        This function returns the CNA profile records of the overlapping
        regions in specific clone to the query region.

        Parameters
        ----------
        chrom : str
            Chromosome name.
        start : int
            The 1-based genomic start pos, inclusive.
        end : int
            The 1-based genomic end pos, exclusive.
        clone_id : str
            The ID of the CNA clone.

        Returns
        -------
        int
            Number of overlapping regions, -1 if error.
        list
            A list of tuples of CNA profiles. Each tuple contains the 
            region-specific CNA profile:
            - int
                Copy numbers of the first allele.
            - int
                Copy numbers  of the first allele.
            - str
                The region ID.
        """
        if clone_id in self.dat:
            cp = self.dat[clone_id]    # cna profile
            ret, hits = cp.fetch(chrom, start, end)
            return((ret, hits))
        else:
            return((0, []))

    def get_all(self):
        """Returns the whole CNA profile of all clones."""
        dat_list = {
                "clone":[],
                "chrom":[],
                "start":[],
                "end":[],
                "name":[],
                "cn_ale0":[],
                "cn_ale1":[]
            }
        for clone_id in sorted(self.dat.keys()):
            cp = self.dat[clone_id]
            cp_dat = cp.get_all()
            n = len(cp_dat["chrom"])
            dat_list["clone"].extend([clone_id] - n)
            dat_list["chrom"].extend(cp_dat["chrom"])
            dat_list["start"].extend(cp_dat["start"])
            dat_list["end"].extend(cp_dat["end"])
            dat_list["name"].extend(cp_dat["name"])
            dat_list["cn_ale0"].extend(cp_dat["cn_ale0"])
            dat_list["cn_ale1"].extend(cp_dat["cn_ale1"])
        return(dat_list)

    def get_clones(self):
        """Get all clone IDs."""
        clones = sorted([c for c in self.dat.keys()])
        return(clones)

    def query(self, name, clone_id):
        """Return the CNA profile records given the query region ID and 
        clone ID.

        This function returns regions whose ID is `name` in clone `clone_id`.
        
        Parameters
        ----------
        name : str
            The ID of the query region.
        clone_id : str
            The ID of the query clone.

        Returns
        -------
        int
            Number of regions whose ID is `name` in `clone_id`.
        list
            A list of tuples of CNA profiles. Each tuple contains the 
            region-specific CNA profile:
            - int
                Copy numbers of the first allele.
            - int
                Copy numbers  of the first allele.
            - str
                The region ID.        
        """
        if clone_id in self.dat:
            cp = self.dat[clone_id]    # cna profile
            ret, hits = cp.query(name)
            return((ret, hits))
        else:
            return((0, []))


        
def load_cna_profile(fn, sep = "\t"):
    """Load CNA profiles from file.

    Parameters
    ----------
    fn : str
        Path to the input file storing clonal CNA profile.
    sep : str, default "\t"
        The file delimiter.

    Returns
    -------
    utils.xcna.CloneCNAProfile or None
        The object of loaded clonal CNA profile. `None` if error.
    """
    try:
        df = load_cnas(fn, sep = sep)
    except:
        error("load CNA profile file failed.")
        return(None)

    dat = CloneCNAProfile()    
    for i in range(df.shape[0]):
        rec = df.iloc[i, ]
        dat.add_cna(
            rec["chrom"],
            rec["start"],
            rec["end"] + 1,
            rec["region"],
            rec["cn_ale0"],
            rec["cn_ale1"],
            rec["clone"]
        )
    return(dat)
                        

    
def save_cna_profile(dat, fn):
    """Save CNA profile to file.
    
    Parameters
    ----------
    dat : utils.xcna.CloneCNAProfile
        The object of clonal CNA profile to be saved.
    fn : str
        Path to the output file.

    Returns
    -------
    Void.
    """
    fp = zopen(fn, "wt")
    cp = dat.get_all()
    for i in range(len(cp["chrom"])):
        s = "\t".join([
                cp["chrom"][i],
                str(cp["start"][i]),
                str(cp["end"][i]),
                cp["clone"][i],
                str(cp["cn_ale0"][i]),
                str(cp["cn_ale1"][i])
            ]) + "\n"
        fp.write(s)
    fp.close()
