# xbarcode.py - barcodes generation and processing.


import numpy as np


class Barcode:
    """Barcode.
    
    This class processes barcodes (e.g., cell barcodes and UMI barcodes),
    implementing multiple operations, such as sampling.
    
    Internally, it uses one integer (e.g., int64) to store barcode string.
    Each of the four bases ('A'/'C'/'G'/'T') occupies two bits, hence one
    barcode string will take up "2*m" bits.
    """
    def __init__(self, m):
        """
        Parameters
        ----------
        m : int
            Length of the barcode string.
        """
        self.m = m

        # d : str
        #   The string of the four bases whose order encodes their bit format,
        #   i.e., 00, 01, 02, 03 for the four bases in order.
        self.d = "ACGT"

        assert m <= 31

    def int2str(self, i):
        """Convert one integer into barcode string.
        
        Parameters
        ----------
        i : int
            One integer.
        
        Returns
        -------
        str
            The converted barcode string.
        """
        x = [self.d[(i >> 2*j) & 3] for j in range(self.m)]
        x.reverse()
        s = "".join(x)
        return(s)
    
    def str2int(self, s):
        """Convert one barcode string into its integer format.
        
        Parameters
        ----------
        s : str
            The barcode string.
        
        Returns
        -------
        int
            The integer format of the barcode `s`.
        """
        assert len(s) == self.m
        i = 0
        for j in range(self.m):
            i |= self.d.index(s[j]) << 2*j
        return(i)
    
    def __randint(self, m, n, b, e):
        """Generate a random sample of barcodes within specific range.
        
        This function generate a random sample of barcodes, in integer format,
        within specific sampling space/range.
        As the whole size of sampling space may be huge, it uses a strategy of
        divide-and-conquer that recursively split the sampling space into four
        subspaces with equal size.
        
        Parameters
        ----------
        m : int
            Length of barcode string.
        n : int
            Size of the sample.
        b : int
            The start index of the sampling space/range.
        e : int
            The end index of the sampling space/range.

        Returns
        -------
        numpy.ndarray
            The sampled barcodes in integer format.
        """
        assert (e - b) % 4 == 0
        if m <= 6:
            x = np.random.choice(range(b, e), size = n, replace = False)
            return(x)
        k1, k2, k3, k4 = np.random.multinomial(n, [0.25] * 4)
        u = (e - b) // 4
        x1 = self.__randint(m - 1, k1, b, b + u)
        x2 = self.__randint(m - 1, k2, b + u, b + 2*u)
        x3 = self.__randint(m - 1, k3, b + 2*u, b + 3*u)
        x4 = self.__randint(m - 1, k4, b + 3*u, b + 4*u)
        return(np.concatenate([x1, x2, x3, x4]))
    
    def sample_int(self, n, sort = True):
        """Generate a random sample of barcodes in integer format.
        
        Parameters
        ----------
        n : int
            The sample size.
        sort : bool, default True
            Whether to sort the sampled barcodes.
        
        Returns
        -------
        numpy.ndarray
            The sampled barcodes in integer format.
        """
        assert n <= 4**self.m
        x = self.__randint(self.m, n, 0, 4**self.m)
        if sort:
            x = np.sort(x)     # as self.d is in ascending order.
        assert len(x) == len(np.unique(x))
        return(x)

    def sample(self, n, sort = True):
        """Generate a random sample of barcode strings.
        
        Parameters
        ----------
        n : int
            The sample size.
        sort : bool, default True
            Whether to sort the sampled barcodes.
        
        Returns
        -------
        list of str
            The sampled barcode strings.
        """
        x = self.sample_int(n, sort = sort)
        s = [self.int2str(i) for i in x]
        return(s)


def rand_cell_barcodes(m, n, suffix = "-1", sort = True):
    """Generate (unique) random cell barcodes.
    
    Parameters
    ----------
    m : int
        Length of cell barcodes.
    n : int
        Number of barcodes to be generated.
    suffix : str, default "-1"
        Suffix appended to the barcodes.
    sort : bool, default True
        Whether to sort the returned barcodes.
    
    Returns
    -------
    list of str
        A list of generated cell barcodes.
    """
    b = Barcode(m)
    s = b.sample(n, sort)
    if suffix:
        s = [x + suffix for x in s]
    return(s)


def rand_umi(m, n, sort = True):
    """Generate (unique) random UMI barcodes.
    
    Parameters
    ----------
    m : int
        Length of UMI barcode.
    n : int
        Number of UMIs to be generated.
    sort : bool, default True
        Whether to sort the returned UMIs.
    
    Returns
    -------
    list of str
        A list of generated UMIs.
    """
    b = Barcode(m)
    s = b.sample(n, sort)
    return(s)
