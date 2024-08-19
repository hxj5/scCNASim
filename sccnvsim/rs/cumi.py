# cumi.py - cell-specific UMI.


import anndata as ad
import multiprocessing
import numpy as np
import os
import pandas as pd
import pickle
from ..utils.base import is_file_empty
from ..utils.xbarcode import Barcode


def gen_umis(
    xdata, 
    chrom_list, chrom_reg_idx_range_list,
    alleles, m, out_dir, ncores = 1,
):
    """To generate cell x feature UMIs to be used in new BAM.
    
    Parameters
    ----------
    xdata : adata object
        The adata object containing *cell x feature* matrices of simulated 
        UMI counts.
    alleles : list
        A list of alleles, e.g., "A", "B", "U".
    m : int
        Length of one new UMI barcode.
    ncores : int
        Number of cores.
    out_dir : str
        Output folder.
    """
    assert m <= 31
    assert len(chrom_list) == len(chrom_reg_idx_range_list)
    n, p = xdata.shape

    ncores = min(n, ncores)
    k = None
    if n % ncores == 0:
        k = n // ncores
    else:
        k = n // ncores + 1

    tmp_dir = os.path.join(out_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok = True)

    fn_list = []
    for idx, i in enumerate(range(0, n, k)):
        fn = os.path.join(tmp_dir, "adata_batch%d.h5ad" % idx)
        xdata_batch = xdata[i:(i+k), :]
        xdata_batch.write_h5ad(fn)
        fn_list.append(fn)
    del xdata

    mp_res = []
    pool = multiprocessing.Pool(processes = ncores)
    for idx, fn in enumerate(fn_list):
        mp_res.append(pool.apply_async(
            func = gen_umis_thread,
            args = (idx, fn, chrom_list, chrom_reg_idx_range_list, alleles, m),
            callback = None
        ))
    pool.close()
    pool.join()
    mp_res = [res.get() for res in mp_res]

    dat = []     # chrom x allele
    for _ in chrom_list:
        chrom_dat = [[] for __ in range(len(alleles))]
        dat.append(chrom_dat)

    for idx_res, res in enumerate(mp_res):
        for c in range(len(chrom_list)):
            for k in range(len(alleles)):
                dat[c][k].extend(res[c][k])

    # save chrom-specific UMI data into file.
    # TODO: this step is time-consuming, use multiprocessing.
    out_fn_list = []
    for c, chrom in enumerate(chrom_list):
        fn = os.path.join(
            out_dir, "chrom%s.allele_x_cell_x_feature.umis.pickle" % chrom)
        with open(fn, "wb") as fp:
            pickle.dump(dat[c], fp)
        out_fn_list.append(fn)
    return(out_fn_list)


def gen_umis_thread(idx, fn, chrom_list, chrom_reg_idx_range_list, alleles, m):
    xdata = ad.read_h5ad(fn)
    RD = None
    for i, ale in enumerate(alleles):
        assert ale in xdata.layers
        if i == 0:
            RD = xdata.layers[ale].copy()
        else:
            RD += xdata.layers[ale]
    n, p = RD.shape

    dat = []     # chrom x allele x cell
    for _ in chrom_list:
        chrom_dat = []
        for __ in range(len(alleles)):
            ale_dat = [[] for ___ in range(n)]
            chrom_dat.append(ale_dat)
        dat.append(chrom_dat)

    b = Barcode(m)
    for i in range(n):
        x = b.sample_int(n = np.sum(RD[i, :]), sort = False)
        r = 0
        for k, ale in enumerate(alleles):
            X = xdata.layers[ale]
            for c, (s, e) in enumerate(chrom_reg_idx_range_list):
                if s is None or e is None:
                    continue
                for j in range(s, e):
                    dat[c][k][i].append(x[r:(r+X[i, j])])
                    r += X[i, j]

    return(dat)     # chrom x allele x cell x feature


def gen_cumis(
    xdata_fn_list, reg_fn_list, reg_idx_range_list,
    umi_fn_list, chrom_list,
    alleles, out_dir,
    use_umi = True,
    max_pool = None, ncores = 1
):
    assert len(xdata_fn_list) == len(chrom_list)
    assert len(reg_fn_list) == len(chrom_list)
    assert len(reg_idx_range_list) == len(chrom_list)
    assert len(umi_fn_list) == len(chrom_list)
    if max_pool is None:
        max_pool = [0] * len(alleles)
    else:
        assert len(max_pool) == len(alleles)

    mp_res = []
    pool = multiprocessing.Pool(processes = ncores)
    for idx, chrom in enumerate(chrom_list):
        mp_res.append(pool.apply_async(
            func = gen_cumis_chrom,
            kwds = dict(
                xdata_fn = xdata_fn_list[idx],
                reg_fn = reg_fn_list[idx],
                reg_idx_range = reg_idx_range_list[idx], 
                umi_fn = umi_fn_list[idx],
                alleles = alleles,
                out_fn = os.path.join(
                    out_dir, "chrom%s.cumi.msampler.pickle" % chrom),
                use_umi = use_umi,
                max_pool = max_pool
            ),
            callback = None
        ))
    pool.close()
    pool.join()
    mp_res = [res.get() for res in mp_res]

    out_fn_list = mp_res
    return(out_fn_list)
    

def gen_cumis_chrom(
    xdata_fn, reg_fn, reg_idx_range, 
    umi_fn, alleles, out_fn,
    use_umi = True, max_pool = None
):
    xdata = ad.read_h5ad(xdata_fn)
    with open(reg_fn, "rb") as fp:
        reg_list = pickle.load(fp)
    with open(umi_fn, "rb") as fp:
        umi_list = pickle.load(fp)

    s, e = reg_idx_range[:2]
    reg_idx_list = None
    if s is None or e is None:
        reg_idx_list = []
    else:
        reg_idx_list = range(s, e)

    ale_samplers = []
    for idx, ale in enumerate(alleles):
        sampler = CUMISampler(
            X = xdata.layers[ale],
            reg_idx_list = reg_idx_list,
            allele_fn_list = [reg.aln_fns[ale] for reg in reg_list],
            use_umi = use_umi,
            umis = umi_list[idx], 
            max_pool = max_pool[idx]
        )
        sampler.sample()
        ale_samplers.append(sampler)

    ms = MergedSampler(
        {ale:sampler for ale, sampler in zip(alleles, ale_samplers)})
    with open(out_fn, "wb") as fp:
        pickle.dump(ms, fp)
    return(out_fn)


class CUMISampler:
    """A sampler to sample cell x feature CUMIs.

    Here CUMI is short for unique UMI barcode, which is typically set as a
    combination of cell and UMI barcodes (or read query name).
    This class is used for sampling old CUMIs from input `cells` and `umis`
    barcodes based on simulated UMI counts in `X`, and then assign new CUMI,
    i.e., cell (+UMI, if `use_umi` is `True`) barcode, to each sampled one.

    Note that one old CUMI can be assigned multiple new CUMIs, if it is 
    sampled more than once.
    
    Attributes
    ----------
    X: np.array
        The new *cell x feature* matrix of simulated UMI counts.
    reg_list : list
        A list of region objects.
    use_umi : int
        Whether use UMI in CUMI.
        If False, only cell barcode is included in CUMI.
    umis : object
        The *cell x feature* new UMIs returned by :func:`gen_umis`.
    max_pool : int
        Maximum pool size, i.e., max number of old CUMIs to be used for
        sampling for one feature. Set to `0` if unlimited. 
        This option is designed to speed up by reducing the reads to be masked,
        since the read mask is time consuming.
    """
    def __init__(
        self, X, reg_idx_list, allele_fn_list,
        use_umi = True, umis = None, 
        max_pool = 0
    ):
        self.X = X
        self.reg_idx_list = reg_idx_list
        self.allele_fn_list = allele_fn_list
        self.use_umi = use_umi
        self.umis = umis
        self.max_pool = max_pool

        assert len(reg_idx_list) == X.shape[1]
        assert len(reg_idx_list) == len(allele_fn_list)
        if umis is not None:
            assert len(umis) == X.shape[0]
            for dat in umis:
                assert len(dat) == X.shape[1]

        self.dat = {}
    
    def __sample_for_feature(self, cells, umis, reg_idx, reg_idx_whole = None):
        """Sample CUMIs for one feature.

        This function will sample CUMIs from input `cells` and `umis` barcodes
        based on the simulated UMI counts in `X`.
        If `use_umi` is `True`, then each sampled CUMI will also be assigned
        a newly generated UMI barcode (note new cell barcodes have been
        generated elsewhere beforehand).

        Note that one old CUMI can be assigned multiple new CUMIs, if it is 
        sampled more than once.
        
        Parameters
        ----------
        cells : list-like
            The cell barcodes (str), as part of CUMIs, to be sampled from.
        umis : list-like
            The UMI barcodes (str), as part of CUMIs, to be sampled from.
            Its length and order should match `cells`.
        reg_idx : int
            The index (0-based) of the region/feature.
        
        Returns
        -------
        Void.
        """
        n, p = self.X.shape
        m = len(cells)

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
        old_cumi_idx_list = cumi_sample_for_cells(m, n_list)

        new_umi_list = []
        for i in range(n):
            new_umi_list.append(self.umis[i][reg_idx])

        assert len(old_cumi_idx_list) == n
        assert len(new_umi_list) == n
        for i in range(n):
            old_cumi_idxes = old_cumi_idx_list[i]
            new_umis = new_umi_list[i]
            assert len(old_cumi_idxes) == self.X[i, reg_idx]
            if len(new_umis) != self.X[i, reg_idx]:
                import logging
                logging.info("n=%d; p=%d; i=%d; reg_idx=%d; reg_idx_whole=%d" % (
                    n, p, i, reg_idx, reg_idx_whole))
                logging.info("len(new_umi_list)=%d;" % len(new_umi_list))
                logging.info("len(new_umis)=%d; Xij=%d" % (len(new_umis), self.X[i, reg_idx]))
                logging.info("new_umis=%s" % str(new_umis))
            assert len(new_umis) == self.X[i, reg_idx]
            
            for old_cumi_idx, new_umi in zip(old_cumi_idxes, new_umis):
                cell, umi = cells[old_cumi_idx], umis[old_cumi_idx]
                if cell not in self.dat:
                    self.dat[cell] = {}
                if umi not in self.dat[cell]:
                    self.dat[cell][umi] = []
                self.dat[cell][umi].append((i, new_umi, reg_idx_whole))

    def query(self, cell, umi):
        """Query new CUMI.

        Parameters
        ----------
        cell : str
            The cell barcode from old CUMI.
        umi : str
            The UMI barcode from old CUMI.
        
        Returns
        -------
        list
            A list of new CUMIs assigned to the old one. Each element is a
            tuple, containing
            int
                The index of new cell, 0-based;
            int
                The integer format of new UMI, can be transformed to string
                format with `int2str()`.
            int
                The index of region/feature, 0-based.
            `None` if the query CUMI (`cell`+`umi`) is not sampled.
        """
        if cell not in self.dat or umi not in self.dat[cell]:
            return(None)
        return(self.dat[cell][umi])

    def sample(self):
        for reg_idx, reg_idx_whole in enumerate(self.reg_idx_list):
            fn = self.allele_fn_list[reg_idx]
            if is_file_empty(fn):
                continue
            dat = load_cumi(fn, sep = "\t")
            self.__sample_for_feature(
                dat["cell"], dat["umi"], reg_idx, reg_idx_whole)


class MergedSampler:
    """Merged object from all allele-specific CUMI samplers.

    Attributes
    ----------
    samplers : dict
        Allele-specific CUMI samplers (CUMISampler object). Keys are the
        alleles (str) and values are samplers.
    """
    def __init__(self, samplers):
        self.samplers = samplers

    def query(self, cell, umi):
        """Query CUMI given cell and umi barcodes.
        
        Parameters
        ----------
        cell : str
            The cell ID. The cell barcode (10x) or sample ID (SMART-seq).
        umi : str
            The read ID. The UMI barcode (10x) or read query name (SMART-seq).
        
        Returns
        -------
        str
            The allele where the query CUMI comes from. `None` if the CUMI
            is not from any sampler.
        list
            The meta data assigned to this CUMI. See the returned value of 
            :func:`CUMISampler.query`. 
            `None` if the CUMI is not from any sampler.
        """
        for allele, sampler in self.samplers.items():
            res = sampler.query(cell, umi)
            if res is None:
                continue
            return((allele, res))
        return((None, None))


def cumi_sample_for_cells(m, n_list):
    """Sampling CUMI indexes.

    Parameters
    ----------
    m : int
        The number of CUMIs to sample from.
    n_list : list-like
        Its elements are the numbers of CUMIs to be sampled in each 
        new sample/cell.

    Returns
    -------
    list
        A list of arrays. Its length and order match `n_list`.
        Each element is an array storing the original indexes (0-based) 
        of sampled CUMIs.
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


def load_cumi(fn, sep = "\t"):
    dat = pd.read_csv(fn, sep = sep, header = None)
    dat.columns = ["cell", "umi"]
    return(dat)
