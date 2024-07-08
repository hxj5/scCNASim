# xcnv.py - clonal CNV profile.


from logging import error
from .grange import Region, RegionSet
from .zfile import zopen
from ..io.base import load_cnvs


class CNVRegCN(Region):
    """Allele-specific copy numbers of CNV region.

    Attributes
    ----------
    chrom : str
        Chromosome name.
    start : int
        The 1-based start pos, inclusive.
    end : int
        The 1-based end pos, exclusive.
    name : str
        The ID of the region.
    cn_ale0 : int
        Copy Number of the first allele.
    cn_ale1 : int
        Copy Number of the second allele.
    """
    def __init__(self, chrom, start, end, name, cn_ale0, cn_ale1):
        super().__init__(chrom, start, end, name)
        self.name = name
        self.cn_ale0 = cn_ale0
        self.cn_ale1 = cn_ale1


class CNVProfile:
    def __init__(self):
        self.rs = RegionSet()

    def add_cnv(self, chrom, start, end, name, cn_ale0, cn_ale1):
        """Add a new CNV profile.

        Parameters
        ----------
        chrom : str
            Chromosome name.
        start : int
            The 1-based start pos, inclusive.
        end : int
            The 1-based end pos, exclusive.
        name : str
            The ID of the region.
        cn_ale0 : int
            Copy Number of the first allele.
        cn_ale1 : int
            Copy Number of the second allele.

        Returns
        -------
        int
            0 success, 1 discarded as duplicate, -1 error.
        """
        reg = CNVRegCN(chrom, start, end, name, cn_ale0, cn_ale1)
        ret = self.rs.add(reg)
        return(ret)

    def fetch(self, chrom, start, end):
        """Get the CNV profile for the query region.

        Parameters
        ----------
        chrom : str
            Chromosome name.
        start : int
            The 1-based start pos, inclusive.
        end : int
            The 1-based end pos, exclusive.

        Returns
        -------
        n : int
            Number of overlapping regions, -1 if error.
        profile : list
            A list of tuples of CNV profiles; `None` if error:
            cn_ale0 : int
                Copy numbers of the first allele.
            cn_ale1 : int
                Copy numbers  of the first allele.
            reg_id : str
                The region ID.
        """
        hits = self.rs.fetch(chrom, start, end)
        res = [(reg.cn_ale0, reg.cn_ale1, reg.name) for reg in hits]
        return((len(res), res))

    def get_all(self):
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
        """Query CNV profile for the given region and clone."""
        hits = self.rs.query(name)
        res = [(reg.cn_ale0, reg.cn_ale1, reg.name) for reg in hits]
        return((len(res), res))


class CloneCNVProfile:
    def __init__(self):
        self.dat = {}

    def add_cnv(self, chrom, start, end, name, cn_ale0, cn_ale1, clone_id):
        if clone_id not in self.dat:
            self.dat[clone_id] = CNVProfile()
        cp = self.dat[clone_id]
        ret = cp.add_cnv(chrom, start, end, name, cn_ale0, cn_ale1)
        return(ret)

    def fetch(self, chrom, start, end, clone_id):
        """Get the CNV profile for the query region and cell.

        Parameters
        ----------
        chrom : str
            Chromosome name.
        start : int
            The 1-based start pos, inclusive.
        end : int
            The 1-based end pos, exclusive.
        clone_id : str
            The ID of the CNV clone.

        Returns
        -------
        See @return of :func:`CNVProfile.fetch`.
        """
        if clone_id in self.dat:
            cp = self.dat[clone_id]    # cnv profile
            ret, hits = cp.fetch(chrom, start, end)
            return((ret, hits))
        else:
            return((0, []))

    def get_all(self):
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
            dat_list["clone"].extend([clone_id] * n)
            dat_list["chrom"].extend(cp_dat["chrom"])
            dat_list["start"].extend(cp_dat["start"])
            dat_list["end"].extend(cp_dat["end"])
            dat_list["name"].extend(cp_dat["name"])
            dat_list["cn_ale0"].extend(cp_dat["cn_ale0"])
            dat_list["cn_ale1"].extend(cp_dat["cn_ale1"])
        return(dat_list)

    def get_clones(self):
        clones = sorted([c for c in self.dat.keys()])
        return(clones)

    def query(self, name, clone_id):
        """Query CNV profile for the given region and clone.

        Returns
        -------
        See @return of :func:`CNVProfile.fetch`.
        """
        if clone_id in self.dat:
            cp = self.dat[clone_id]    # cnv profile
            ret, hits = cp.query(name)
            return((ret, hits))
        else:
            return((0, []))


def load_cnv_profile(fn, sep = "\t"):
    """Load CNV profiles from file.

    Returns
    -------
    CloneCNVProfile object.
        `None` if error.
    """
    try:
        df = load_cnvs(fn, sep = sep)
    except:
        error("load CNV profile file failed.")
        return(None)

    dat = CloneCNVProfile()    
    for i in range(df.shape[0]):
        rec = df.loc[i, ]
        dat.add_cnv(
            rec["chrom"],
            rec["start"],
            rec["end"] + 1,
            rec["region"],
            rec["cn_ale0"],
            rec["cn_ale1"],
            rec["clone"]
        )
    return(dat)
                        

def save_cnv_profile(dat, fn):
    """Save CNV profile to file."""
    fp = zopen(fn, "wt")
    cp = dat.get_all()
    for i in range(len(cp["chrom"])):
        s = "\t".join([
                cp["chrom"][i],
                str(cp["start"][i]),
                str(cp["end"][i]),
                cp["name"][i],
                cp["clone"][i],
                str(cp["cn_ale0"][i]),
                str(cp["cn_ale1"][i])
            ]) + "\n"
        fp.write(s)
    fp.close()
