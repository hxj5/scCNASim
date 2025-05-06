# gcna.py - clonal CNA profile.


import functools
import numpy as np
import pandas as pd
from logging import error
from logging import warning as warn
from .grange import cmp_two_intervals
from ..xlib.xbase import is_file_empty
from ..xlib.xfile import zopen
from ..xlib.xrange import Region, RegionSet,  \
    format_chrom, format_start, format_end, reg2str



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
        
        

def load_cnas(fn, sep = "\t", cna_mode = "hap-aware"):
    """Load clonal CNA profile from a header-free file.
    
    Parameters
    ----------
    fn : str
        Path to a a header-free file containing clonal CNA profile, whose
        first several columns should be:
        - "chrom" (str): chromosome name of the CNA region.
        - "start" (int): start genomic position of the CNA region, 1-based
          and inclusive.
        - "end" (int): end genomic position of the CNA region, 1-based and
          inclusive.
        - "clone" (str): clone ID.
        if `cna_mode` is "hap-aware":
        - "cn_ale0" (int): copy number of the first allele.
        - "cn_ale1" (int): copy number of the second allele.
        otherwise:
        - "cn" (int): copy number of both alleles.
    sep : str, default "\t"
        File delimiter.
    cna_mode : {"hap-aware", "hap-unknown"}
        The mode of CNA profiles.
        - "hap-aware": haplotype/allele aware.
        - "hap-unknown": haplotype/allele unknown.

    Returns
    -------
    pandas.DataFrame
        The loaded clonal CNA profile, whose first several columns are
        "chrom", "start", "end", "clone", 
        and
        - if `cna_mode` is "hap-aware": "cn_ale0", and "cn_ale1", "region";
        - otherwise: "cn", "region".
        Note that "region" is a formatted string combining "chrom", "start",
        and "end".
    """
    df = None
    if is_file_empty(fn):
        if cna_mode == "hap-aware":
            df = pd.DataFrame(columns = ["chrom", "start", "end", "clone", \
                        "cn_ale0", "cn_ale1", "region"])
        else:
            df = pd.DataFrame(columns = ["chrom", "start", "end", "clone", \
                        "cn", "region"])
        return(df)
    
    df = pd.read_csv(fn, sep = sep, header = None, dtype = {0: str})
    df.columns = df.columns.astype(str)
    
    if cna_mode == "hap-aware":
        assert df.shape[1] >= 6
        df.columns.values[:6] = [
            "chrom", "start", "end", "clone", "cn_ale0", "cn_ale1"]
    else:
        assert df.shape[1] >= 5
        df.columns.values[:5] = [
            "chrom", "start", "end", "clone", "cn"]        

    df["chrom"] = df["chrom"].map(format_chrom)
    df["start"] = df["start"].map(format_start)
    df["end"] = df["end"].map(format_end)
    df["region"] = [reg2str(
                        df["chrom"].iloc[i], 
                        df["start"].iloc[i], 
                        df["end"].iloc[i]
                    ) for i in range(df.shape[0])]
    return(df)



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
    
    

def check_dup_cna(fn):
    """Check duplicated records in the CNA profile file.
    
    Parameters
    ----------
    fn : str
        Path to the CNA profile file.
    
    Returns
    -------
    int
        Number of duplicates in the file.
    """
    df = load_cnas(fn, sep = "\t")
    bool_dup = df.duplicated()
    n_dup = np.sum(bool_dup)
    if n_dup > 0:
        warn("%d/%d duplicates in CNA profiles." % (n_dup, df.shape[0]))
    return(n_dup)
    

    
def merge_cna_profile(in_fn, out_fn, max_gap = 1):
    """Merge adjacent regions with the same CNA profiles.

    Merge adjacent regions with the same allele-specific copy number
    profile in each CNA clone.

    Parameters
    ----------
    in_fn : str
        Path to input file.
    out_fn : str
        Path to output file.
    max_gap : int, default 1
        The maximum gap length that is allowed between two adjacent regions.
        `1` for strict adjacence.

    Returns
    -------
    int
        The return code. 0 if success, negative if error.
    int
        Number of records before merging.
    int
        Number of records after merging.
    """
    sep = "\t"
    n_old, n_new = -1, -1

    # load data
    try:
        df = load_cnas(in_fn, sep = sep)
    except Exception as e:
        error("load CNA profile file failed '%s'." % str(e))
        return((-3, n_old, n_new))
    n_old = df.shape[0]

    dat = {}
    for i in range(df.shape[0]):
        rec = df.iloc[i, ]
        chrom = rec["chrom"]
        clone = rec["clone"]
        if clone not in dat:
            dat[clone] = {}
        if chrom not in dat[clone]:
            dat[clone][chrom] = {}
        cns = "%d_%d" % (rec["cn_ale0"], rec["cn_ale1"])
        if cns not in dat[clone][chrom]:
            dat[clone][chrom][cns] = []
        dat[clone][chrom][cns].append((rec["start"], rec["end"]))


    # merge (clone-specific) adjacent CNAs.
    for clone, cl_dat in dat.items():
        for chrom, ch_dat in cl_dat.items():
            for cns in ch_dat.keys():
                iv_list = sorted(
                    ch_dat[cns], 
                    key = functools.cmp_to_key(cmp_two_intervals)
                )
                s1, e1 = iv_list[0]
                new_list = []
                for s2, e2 in iv_list[1:]:
                    if s2 <= e1 + max_gap:    # overlap adjacent region
                        e1 = max(e1, e2)
                    else:                     # otherwise
                        new_list.append((s1, e1))
                        s1, e1 = s2, e2
                new_list.append((s1, e1))
                ch_dat[cns] = new_list


    # check whether there are (strictly) overlapping regions with 
    # distinct profiles.
    for clone, cl_dat in dat.items():
        for chrom, ch_dat in cl_dat.items():
            iv_list = []
            for cns in ch_dat.keys():
                cn_ale0, cn_ale1 = [int(x) for x in cns.split("_")]
                iv_list.extend(
                    [(s, e, cn_ale0, cn_ale1) for s, e in ch_dat[cns]])
            iv_list = sorted(
                iv_list, 
                key = functools.cmp_to_key(cmp_two_intervals)
            )
            s1, e1 = iv_list[0][:2]
            for iv in iv_list[1:]:
                s2, e2 = iv[:2]
                if s2 <= e1:    # overlap adjacent region
                    error("distinct CNA profiles '%s', (%d, %d) and (%d, %d)." % 
                        (chrom, s1, e1, s2, e2))
                    return((-5, n_old, n_new))
            cl_dat[chrom] = iv_list


    # save profile
    n_new = 0
    fp = open(out_fn, "w")
    for clone in sorted(dat.keys()):
        cl_dat = dat[clone]
        for chrom in sorted(cl_dat.keys()):
            ch_dat = cl_dat[chrom]
            for s, e, cn_ale0, cn_ale1 in ch_dat:
                fp.write("\t".join([chrom, str(s), str(e), \
                    clone, str(cn_ale0), str(cn_ale1)]) + "\n")
                n_new += 1
    fp.close()
    return((0, n_old, n_new))
