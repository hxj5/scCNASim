# grange.py - genomic range/interval routine.


import numpy as np
from functools import cmp_to_key
from intervaltree import IntervalTree



class RegPos(int):
    """Genomic position compatible with "infinite" value.
    
    This class inherits from the `int` class, but is compatible with 
    "infinite" value:
    - Any value is capped as `REG_MAX_POS`.
    - The infinite value will be set as `REG_MAX_POS`.
    - In type conversion, the string starts with "inf" will be recognized
      as infinite.
    - The string format of value no less than `REG_MAX_POS` is "Inf".
    """

    # note that ``int`` is immutable, hence __new__() should be used
    # instead of __init__().

    def __new__(cls, x, *args, **kwargs):
        if isinstance(x, str):
            if x.lower().startswith("inf"):
                x = REG_MAX_POS
        elif np.isinf(x):
            x = REG_MAX_POS
        return super(RegPos, cls).__new__(cls, x)

    def __str__(self):
        if self >= REG_MAX_POS:
            return("Inf")
        else:
            return("%d" % int(self))

    def __add__(self, x):
        if self > REG_MAX_POS - x:
            return(REG_MAX_POS)
        else:
            return(int(self) + x)


        
class Region:
    """Genomic region."""

    def __init__(self, chrom, start, end, rid = None):
        """
        Parameters
        ----------
        chrom : str
            Chromosome name.
        start : int
            1-based start pos, inclusive.
        end : int
            1-based end pos, exclusive.
        rid : str or None, default None
            Region ID.
            If `None`, it will be set as "{chrom}_{start}_{end}".
        """
        self.chrom = format_chrom(chrom)
        self.start = start
        self.end = end
        self._rid = rid

        # len : int
        #   Length of the region.
        self.len = self.end - self.start

    def compare(self, region):
        """Compare with another region.

        Parameters
        ----------
        region : utils.grange.Region
            The region to be compared with `self`.
        
        Returns
        -------
        int
            Comparison result.
            - negative integer if `self` is smaller; 
            - 0 if equal;
            - positive integer if `self` is bigger.
        """
        if self.chrom == region.chrom:
            if self.start == region.start:
                return(self.end - region.end)
            else:
                return(self.start - region.start)
        elif self.chrom < region.chrom:
            return(-1)
        else:
            return(1)

    def get_id(self):
        if self._rid is None:
            self._rid = "%s_%d_%d" % (self.chrom, self.start, self.end)
        return self._rid

    def get_len(self):
        return self.len


    
class RegionSet:
    """Region set with payload.

    This class enables efficiently querying overlapping regions (with payload)
    for a specific genomic range, by utilizing the data structure of 
    "interval tree".

    The region objects could be simply of class :class:`~utils.grange.Region`
    or its subclasses with payloads (i.e., additional attributes/data).
    """
    def __init__(self, is_uniq = False):
        """
        Parameters
        ----------
        is_uniq : bool, default False
            Whether payloads of duplicate regions should be discarded.
        """
        self.is_uniq = is_uniq

        # creg : dict of {str : list of utils.grange.Region}
        #   Regions in each chromosome.
        #   Keys are chromosome names, values are lists of regions.
        self.creg = {}

        # is_sorted : dict of {str : bool}
        #   Whether the regions in each chromosome have been sorted.
        #   Keys are chromosome names, values are bool values.
        self.is_sorted = {}

        # ctree : dict of {str : intervaltree.IntervalTree}
        #   Intervaltree for each chromosome.
        self.ctree = {}

        # cid : dict of {str : list of utils.grange.Region}
        #   Mapping from region ID to a list of regions.
        #   Keys are region IDs, values are a list of regions with that ID.
        self.cid = {}

        # n : int
        #   Total number of regions.
        self.n = 0

    def __sort_chrom_regions(self, chrom):
        if chrom in self.creg:
            if not self.is_sorted[chrom]:
                self.creg[chrom] = self.__sort_regions(self.creg[chrom])
                self.is_sorted[chrom] = True            

    def __sort_regions(self, reg_list):
        cmp_reg = lambda r1, r2: r1.compare(r2)
        return sorted(reg_list, key = cmp_to_key(cmp_reg))
    
    def add(self, region):
        """Add a new region into the set.

        Parameters
        ----------
        region : utils.grange.Region
            The region to be added.
        
        Returns
        -------
        int
            Return code. 0 success, 1 discarded as duplicate, -1 error.
        """
        chrom = format_chrom(region.chrom)
        if chrom not in self.creg:
            self.creg[chrom] = list()
            self.ctree[chrom] = IntervalTree()

        item = None
        rid = region.get_id()
        if rid in self.cid:
            if self.is_uniq:
                return(1)
            else:
                item = region
        else:
            self.cid[rid] = list()
            item = region
        if item is not None:
            self.creg[chrom].append(item)
            self.is_sorted[chrom] = False
            self.ctree[chrom][region.start:region.end] = item
            self.cid[rid].append(item)
            self.n += 1
        return(0)

    def destroy(self):
        self.reset()

    def fetch(self, chrom, start, end):
        """Fetch overlapping regions from the set.

        This function fetches all regions in the set that overlap with the
        query range "{chrom}:{start}-{end}".

        Parameters
        ----------
        chrom : str
            Chromosome name.
        start : int
            1-based start pos, inclusive.
        end : int
            1-based end pos, exclusive.

        Returns
        -------
        list of utils.grange.Region
            All overlapping regions.
        """
        chrom = format_chrom(chrom)
        if chrom not in self.ctree:
            return([])
        tree = self.ctree[chrom]
        hits = [region for begin, end, region in tree[start:end]]
        return(hits)

    def get_n(self):
        return(self.n)

    def get_regions(self, chrom = None, sort = False):
        """Get regions in specific chromosome(s).

        Parameters
        ----------
        chrom : str or None, default None
            Chromosome name;
            set to `None` to use all chromosomes.
        sort : bool, default False
            Whether to sort the returned regions within each chromosome.

        Returns
        -------
        list of utils.grange.Region
            A list of regions.
        """
        ch_list = []
        if chrom is None:
            ch_list = self.creg.keys()
        else:
            chrom = format_chrom(chrom)
            if chrom in self.creg:
                ch_list = [chrom]
            else:
                return([])

        lst = []        
        for ch in ch_list:
            if sort:
                self.__sort_chrom_regions(ch)
            lst.extend(self.creg[ch])
        return(lst)

    def merge(self, rs):
        """Merge another region set.

        Parameters
        ----------
        rs : utils.grange.RegionSet
            The set of regions to be merged.
        
        Returns
        -------
        int
            Number of regions merged if success, -1 if error.
        """
        k = 0
        reg_list = rs.get_regions()
        for region in reg_list:
            ret = self.add(region)
            if ret != 0:
                if ret < 0:
                    return(-1)
            else:
                k += 1
        return(k)

    def query(self, rid):
        """Query region(s) given its ID.

        Parameters
        ----------        
        rid : str
            Region ID.

        Returns
        -------
        list of utils.grange.Region
            A list of regions whose ID are `rid`.
        """
        if rid in self.cid:
            return(self.cid[rid])
        return([])

    def reset(self):
        for chrom, reg_list in self.creg.items():
            reg_list.clear()
        self.creg.clear()
        self.n = 0
        self.is_sorted.clear()
        for chrom, tree in self.ctree.items():
            tree.clear()
        self.ctree.clear()
        for chrom, id_set in self.cid.items():
            id_set.clear()
        self.cid.clear()

    def sort(self):
        """Sort the regions within each chromosome."""
        for chrom in self.citem:
            if not self.is_sorted[chrom]:
                self.citem[chrom] = self.__sort_items(self.citem[chrom])
                self.is_sorted[chrom] = True

                

def format_chrom(chrom):
    """Format chromosome name, removing the "chr" prefix (if available)."""
    return chrom[3:] if chrom.lower().startswith("chr") else chrom



def format_start(x, base = 1):
    """Format start genomic position of one region.
    
    Parameters
    ----------
    x : int or str
        The start genomic position of one region.
    base : int, default 1
        The index of the first base in the reference genome.

    Returns
    -------
    utils.grange.RegPos
        The formatted start genomic position of one region.
    """
    x = RegPos(x)
    if x < base:
        x = base
    return(x)



def format_end(x, base = 1):
    """Format end genomic position of one region.
    
    Parameters
    ----------
    x : int or str
        The end genomic position of one region.
    base : int, default 1
        The index of the first base in the reference genome.

    Returns
    -------
    utils.grange.RegPos
        The formatted end genomic position of one region.
    """
    x = RegPos(x)
    if x < base:
        x = base
    return(x)



def reg2str(chrom, start, end, base = 1):
    """concatenate chrom, start, and end into a string.
    
    Parameters
    ----------
    chrom : str
        Chromosome name.
    start : int or str
        The start genomic position of one region.
    end : int or str
        The end genomic position of one region.
    base : int, default 1
        The index of the first base in the reference genome.

    Returns
    -------
    str
        The concatenated string for the region.    
    """
    chrom = format_chrom(chrom)
    start = format_start(start, base = base)
    end = format_end(end, base = base)
    s = None
    if end >= REG_MAX_POS:
        s = "%s:%s-" % (chrom, start)
    else:
        s = "%s:%s-%s" % (chrom, start, end)
    return(s)



def str2tuple(s):
    """Convert a string of genomic region into 3-element tuple.
    
    Parameters
    ----------
    s : str
        The string of genomic region.
    
    Returns
    -------
    tuple of (str, int, int) or None
        A tuple of 3 elements: chrom, start, and end.
        `None` if the input `s` is invalid.
    """
    if s is None or not isinstance(s, str):
        return(None)
    if len(s) <= 0:
        return(None)
    if ":" in s:
        if s.count(":") != 1:
            return(None)
        chrom, coord = s.split(":")
        if len(chrom) <= 0:
            return(None)
        if len(coord) <= 0:
            return((chrom, None, None))
        if "-" in coord:
            if coord.count("-") != 1:
                return(None)
            start, end = coord.split("-")
            if len(start) <= 0:
                start = None
            else:
                try:
                    start = int(start)
                except:
                    return(None)
            if len(end) <= 0:
                end = None
            else:
                try:
                    end = int(end)
                except:
                    return(None)
            return((chrom, start, end))
        else:
            try:
                start = int(coord)
            except:
                return(None)
            else:
                return((chrom, start, None))
    else:
        chrom = s
        if "-" in chrom:
            return(None)
        else:
            return((chrom, None, None))


        
REG_MAX_POS = 0x7fffffff    # same with setting of pysam
