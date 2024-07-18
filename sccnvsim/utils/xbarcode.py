# xbarcode.py - barcodes generation and processing.


import numpy as np


class Barcode:
    def __init__(self, m):
        self.m = m
        self.d = "ACGT"    # in ascending order.

        assert m <= 31

    def int2str(self, i):
        x = [self.d[(i >> 2*j) & 3] for j in range(self.m)]
        x.reverse()
        s = "".join(x)
        return(s)
    
    def str2int(self, s):
        assert len(s) == self.m
        i = 0
        for j in range(self.m):
            i |= self.d.index(s[j]) << 2*j
        return(i)
    
    def __randint(self, m, n, b, e):
        """Generate a random sample """
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
        assert n <= 4**self.m
        x = self.__randint(self.m, n, 0, 4**self.m)
        if sort:
            x = np.sort(x)     # as self.d is in ascending order.
        assert len(x) == len(np.unique(x))
        return(x)

    def sample(self, n, sort = True):
        x = self.sample_int(n, sort = sort)
        s = [self.int2str(i) for i in x]
        return(s)


def rand_cell_barcodes(m, n, suffix = "-1", sort = True):
    """Generate random unique cell barcodes.
    
    Parameters
    ----------
    m : int
        Length of cell barcodes.
    n : int
        Number of barcodes to be generated.
    suffix : str
        Suffix appended to the barcodes. Default is "-1".
    sort : bool
        Whether to sort the returned barcodes.
    
    Returns
    -------
    list
        A list of generated cell barcodes.
    """
    b = Barcode(m)
    s = b.sample(n, sort)
    if suffix:
        s = [x + suffix for x in s]
    return(s)


def rand_umi(m, n, sort = True):
    """Generate random unique UMIs.
    
    Parameters
    ----------
    m : int
        Length of UMI.
    n : int
        Number of UMIs to be generated.
    sort : bool
        Whether to sort the returned UMIs.
    
    Returns
    -------
    list
        A list of generated UMIs.
    """
    b = Barcode(m)
    s = b.sample(n, sort)
    return(s)