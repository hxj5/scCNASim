# uumi.py - unique UMI.


import numpy as np
from ..utils.xbarcode import Barcode


class UMIGenerator:
    """To generate cell x feature UMIs to be used in new BAM.
    
    Attributes
    ----------
    X: np.array
        The new *cell x feature* matrix of simulated UMI counts.
    m : int
        Length of one new UMI barcode.
    b : xbarcode::Barcode object.
        The object to sample UMIs.
    dat : list
        The object storing the UMIs, *cell x feature x UMIs(or UMI_ints)*.
        Its length equals to number of rows in `X`. Each element is a list of
        arrays, while each array stores the UMIs for one feature.
    """
    def __init__(self, X, m):
        self.X = X
        self.m = m
        assert m <= 31

        self.b = Barcode(m)
        self.dat = None

        self.__sample_umi()

    def __sample_umi(self):
        assert self.dat is None
        self.dat = []
        n, p = self.X.shape
        for i in range(n):
            x = self.b.sample_int(n = np.sum(self.X[i, :]), sort = False)
            res = []
            k = 0
            for j in range(p):
                res.append(x[k:(k+self.X[i, j])])
                k += self.X[i, j]
            self.dat.append(res)

    def get_feature_umi(self, feature_idx):
        """Return the UMIs of each cell for this feature."""
        res = []
        for i in range(self.X.shape[0]):
            res.append(self.dat[i][feature_idx])
        return(res)

    def int2str(self, i):
        return self.b.int2str(i)


class UUMISampler:
    """A sampler to sample cell x Feature UUMIs.

    Here UUMI is short for unique UMI barcode, which is typically set as a
    combination of cell and UMI barcodes.
    This class is used for sampling old UUMIs from input `cells` and `umis`
    barcodes based on simulated UMI counts in `X`, and then assign new UUMI,
    i.e., cell (+UMI, if `use_umi` is `True`) barcode, to each sampled one.

    Note that one old UUMI can be assigned multiple new UUMIs, if it is 
    sampled more than once.
    
    Attributes
    ----------
    X: np.array
        The new *cell x feature* matrix of simulated UMI counts.
    m : int
        Length of one new UMI barcode.
    use_umi : int
        Whether use UMI in UUMI.
        If False, only cell barcode is included in UUMI.
    max_pool : int
        Maximum pool size, i.e., max number of old UUMIs to be used for
        sampling for one feature. Set to `0` if unlimited.
    """
    def __init__(self, X, m = 10, use_umi = True, max_pool = 0):
        self.X = X
        self.m = m
        self.use_umi = use_umi
        self.max_pool = max_pool

        self.umi_gen = None
        if use_umi:
            self.umi_gen = UMIGenerator(X, m)
        self.dat = {}

    def int2str(self, i):
        return self.umi_gen.int2str(i)

    def query(self, cell, umi):
        """Query new UUMI.

        Parameters
        ----------
        cell : str
            The cell barcode from old UUMI.
        umi : str
            The UMI barcode from old UUMI.
        
        Returns
        -------
        list
            A list of new UUMIs assigned to the old one. Each element is a
            tuple, containing
            int
                The index of new cell, 0-based;
            int
                The integer format of new UMI, can be transformed to string
                format with `int2str()`. `None` if `use_umi` is `False`.
            int
                The index of region/feature, 0-based.
            `None` if the query UUMI (`cell`+`umi`) is not sampled.
        """
        if cell not in self.dat or umi not in self.dat[cell]:
            return(None)
        return(self.dat[cell][umi])
    
    def sample(self, cells, umis, reg_idx):
        """Sample UUMIs.

        This function will sample UUMIs from input `cells` and `umis` barcodes
        based on the simulated UMI counts in `X`.
        If `use_umi` is `True`, then each sampled UUMI will also be assigned
        a newly generated UMI barcode (note new cell barcodes have been
        generated elsewhere beforehand).

        Note that one old UUMI can be assigned multiple new UUMIs, if it is 
        sampled more than once.
        
        Parameters
        ----------
        cells : list-like
            The cell barcodes (str), as part of UUMIs, to be sampled from.
        umis : list-like
            The UMI barcodes (str), as part of UUMIs, to be sampled from.
            Its length and order should match `cells`.
        reg_idx : int
            The index (0-based) of the region/feature.
        
        Returns
        -------
        Void.
        """
        n, p = self.X.shape
        m = len(cells)
        assert len(umis) == m

        if self.max_pool > 0 and m > self.max_pool:
            idx = np.random.choice(
                range(m), 
                size = self.max_pool,
                replace = False
            )
            idx = sorted(list(idx))
            cells = [cells[i] for i in idx]
            umis = [umis[i] for i in idx]
            m = self.max_pool
        
        n_list = self.X[:, reg_idx]
        uumi_idx_list = sample_uumi_idx(m, n_list)

        new_umi_list = None
        if self.use_umi:
            new_umi_list = self.umi_gen.get_feature_umi(reg_idx)
        else:
            new_umi_list = []
            for i in range(n):
                new_umi_list.append([None] * self.X[i, reg_idx])

        assert len(uumi_idx_list) == n
        assert len(new_umi_list) == n
        for i in range(n):
            uumi_idxes = uumi_idx_list[i]
            new_umis = new_umi_list[i]
            assert len(uumi_idxes) == self.X[i, reg_idx]
            assert len(new_umis) == self.X[i, reg_idx]
            
            for uumi_idx, new_umi in zip(uumi_idxes, new_umis):
                cell, umi = cells[uumi_idx], umis[uumi_idx]
                if cell not in self.dat:
                    self.dat = {}
                if umi not in self.dat[cell]:
                    self.dat[cell][umi] = []
                self.dat[cell][umi].append((i, new_umi, reg_idx))


def sample_uumi_idx(m, n_list):
    """Sampling UUMI indexes.

    Parameters
    ----------
    m : int
        The number of UUMIs to sample from.
    n_list : list-like
        Its elements are the numbers of UUMIs to be sampled in each 
        new sample/cell.

    Returns
    -------
    list
        A list of arrays. Its length and order match `n_list`.
        Each element is an array storing the original indexes (0-based) 
        of sampled UUMIs.
    """
    n_total = np.sum(n_list)
    idx = None
    if n_total <= m:
        idx = np.random.choice(
            range(m), 
            size = n_total,
            replace = False
        )
    else:
        idx = np.random.choice(
            range(m),
            size = n_total,
            replace = True
        )
    res = []
    k = 0
    for n in n_list:
        res.append(np.sort(idx[k:(k+n)]))
        k += n
    return(res)