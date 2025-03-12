# cumi.py - cell-specific UMI.


# TODO:
# 1. consider generaing CUMIs (UMIs) by iterating integers from 1 to N, where
#    N is the total number of CUMIs in the simulated count matrix.
#    Current strategy of randomly sampling CUMIs is time and memory consuming,
#    becuase it requires storing all CUMIs, wheras proposed strategy almost
#    stores only one CUMI and should be much more efficient.


import anndata as ad
import multiprocessing
import numpy as np
import os
import pandas as pd
import pickle
from ..io.base import load_h5ad, save_h5ad
from ..utils.base import is_file_empty
from ..utils.xbarcode import Barcode


def gen_umis(
    xdata, 
    chrom_list, chrom_reg_idx_range_list,
    alleles, m, out_dir, ncores = 1
):
    """Generate UMIs to be used in new BAM.

    This function generates chrom-specific *allele x cell x feature* UMIs
    based on the input count matrices in `xdata`. 
    
    Parameters
    ----------
    xdata : anndata.Anndata
        The ".adata" object containing the allele-specific *cell x feature*
        matrices of simulated UMI counts.
        It should contain three layers "A", "B", "U".
    chrom_list : list of str
        A list of chromosome names.
    chrom_reg_idx_range_list : list of tuple
        The range of chrom-specific feature indexes.
        Each element in the list is a tuple of (int, int) that are the 0-based
        start (inclusive) and end (exclusive) indexes of features.
        The index is within transcriptomics-scale instead of chrom-scale.
        The tuple would be (None, None) if one chromosome covers no features.
        Note that the order of the elements should match `chrom_list`. 
    alleles : list of str
        A list of alleles, e.g., ["A", "B", "U"].
    m : int
        Length of one new UMI barcode.
    out_dir : str
        Output folder.
    ncores : int, default 1
        Number of cores.

    Returns
    -------
    list of str
        A list of files, each is a pickle object file storing a "multi-layer" 
        list of chrom-specific *allele x cell x feature* UMIs.
        The UMIs can be accessed by, e.g., 
        dat[allele_idx][cell_idx][feature_idx] (which is a list of UMIs),
        where the 
        - "allele_idx" 0-based, matching the order of `alleles`;
        - "cell_idx" 0-based, matching the row order of ".obs" in `xdata`;
        - "feature_idx" 0-based, matching the row order of ".var"  in `xdata`.
        Note that the returned UMIs are stored as `int` type, which can be
        converted to `str` type by calling 
        :func:`~utils.xbarcode.Barcode.int2str`.
    """
    assert m <= 31
    assert len(chrom_list) == len(chrom_reg_idx_range_list)
    n, p = xdata.shape

    # split cells for multi-processing
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
        save_h5ad(xdata_batch, fn)
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

    # merge chrom-specific UMIs
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


def gen_umis_thread(
    idx, fn, 
    chrom_list, chrom_reg_idx_range_list, 
    alleles, m
):
    xdata = load_h5ad(fn)
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

    # UMIs are generated in a cell-specific manner, mimicking real data that
    # the UMI of each transcript should be unique within one cell.
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


def sample_cumis(
    xdata_fn_list, reg_fn_list, reg_idx_range_list,
    umi_fn_list, chrom_list,
    alleles, out_dir,
    use_umi = True,
    max_pool = None, ncores = 1
):
    """Sampling cell-specific UMIs.

    This function samples cell-specific UMIs (CUMIs) from the previously
    extracted CUMIs for each feature, and then assign each sampled CUMI with a
    new CUMI barcode.

    Note that each CUMI represents a group of reads sharing the same cell+UMI
    barcode (droplet-based platforms) or cell+query_name (well-based 
    platforms).
    
    Parameters
    ----------
    xdata_fn_list : list of str
        A list of ".adata" filenames, each stores chrom-specific count matrices
        for every allele.
        Note that `xdata_fn_list`, `reg_fn_list`, `reg_idx_range_list`,
        `umi_fn_list`, and `chrom_list` should have the same length and order.
    reg_fn_list : list of afc.gfeature.BlockRegion
        A list of :class:`~afc.gfeature.BlockRegion` objects, each stores the
        allele-specific old CUMIs to be sampled from in its `.aln_fns` 
        attribute.
    reg_idx_range_list : list of tuple
        The range of chrom-specific feature indexes.
        Each element in the list is a tuple of (int, int) that are the 0-based
        start (inclusive) and end (exclusive) indexes of features.
        The index is within transcriptomics-scale instead of chrom-scale.
        The tuple would be (None, None) if one chromosome covers no features.
        Note that the order of the elements should match `chrom_list`.        
    umi_fn_list : list of str
        A list of files, each stores chrom-specific new UMIs to be used as
        part of new CUMIs.
    chrom_list : list of str
        A list of chromosome names.
    alleles : list of str
        A list of alleles, e.g., ["A", "B", "U"].
    out_dir : str
        Output folder.
    use_umi : bool, default True
        Whether the sequencing platform uses UMIs.
    max_pool : list of int or None, default None, meaning no limit
        A list of maximum size of sampling pool of old CUMIs for each allele,
        0 means no limit for specific allele.
        If None, every allele has no limit.
    ncores : int, default 1
        Number of cores.

    Returns
    -------
    list of str
        A list of pickle files, each storing a chrom-specific `MergedSampler`
        object that has `query` method for accessing sampled CUMIs.
        Note that its length and order match `chrom_list`.
    """
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
            func = sample_cumis_chrom,
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
    

def sample_cumis_chrom(
    xdata_fn, reg_fn, reg_idx_range, 
    umi_fn, alleles, out_fn,
    use_umi = True, max_pool = None
):
    xdata = load_h5ad(xdata_fn)
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
            umis = umi_list[idx],
            use_umi = use_umi,
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
    """A CUMI sampler.

    This class is used for sampling old CUMIs from input cell and umi barcodes
    based on simulated UMI counts in `X`, and then assign new CUMI to each 
    sampled one.
    
    Note that one old CUMI can be assigned multiple new CUMIs, if it is 
    sampled more than once.
    """
    def __init__(
        self, X, reg_idx_list, allele_fn_list,
        umis, use_umi = True, 
        max_pool = 0
    ):
        """
        Parameters
        ----------
        X : numpy.ndarray
            The allele-specific *cell x feature* matrix of simulated UMI/CUMI 
            counts.
        reg_idx_list : list of int
            A list of chrom-specific 0-based feature index.
            The index is for all features of all chromosomes, which means the
            index is not always starting from 0 in each chromosome.
        allele_fn_list : list of str
            Allele-specific list of files, each contains the feature-specific
            old CUMIs to be sampled from.
        umis : list or None
            The *cell x feature* new UMIs.
            It contains *n* cell-specific (sub-)lists, each (sub)list contains
            *p* lists of feature-specific UMIs, where (n, p) is the shape
            of `X`.
            None means do not use it as part of new CUMIs.
        use_umi : bool, default True
            Whether the sequencing platform uses UMIs.
        max_pool : int, default 0
            Maximum size of sampling pool of old CUMIs, 0 means no limit.
            This option is designed to speed up by reducing the number of reads
            to be masked, since the read mask is time consuming.
        """
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

        # dat : dict
        #   The data structure stores the mapping between sampled old CUMIs
        #   to the new assigned CUMIs.
        #   It is a two-layer dict, with "cell (cell barcode/sample ID, str)"
        #   and "UMI (UMI barcode/query name, str)" of sampled old CUMIs as
        #   their keys, respectively, and list of new assigned CUMIs as values.
        #   Each element of the list is a tuple(int, int, int), containing
        #   - (int) The index of new cell, 0-based.
        #   - (int) The integer format of new UMI barcode, can be transformed
        #     to string format with `int2str()`.
        #   - (int) The 0-based index of feature within transcriptomics-scale.
        self.dat = {}
    
    def __sample_for_feature(self, cells, umis, reg_idx, reg_idx_whole = None):
        """Sample CUMIs for one feature.
        
        Parameters
        ----------
        cells : list of str
            The feature-specific cell barcodes, as part of CUMIs, to be sampled
            from.
        umis : list of str
            The feature-specific UMI barcodes, as part of CUMIs, to be sampled
            from.
            Its length and order should match `cells`.
        reg_idx : int
            The index (0-based) of the feature within chrom-scale.
        reg_idx_whole : int
            The index (0-based) of the feature within transcriptomics-scale.

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
        
        # sample old CUMIs.
        n_list = self.X[:, reg_idx]
        old_cumi_idx_list = cumi_sample_for_cells(m, n_list)

        new_umi_list = []
        for i in range(n):
            new_umi_list.append(self.umis[i][reg_idx])

        # assign new CUMIs to sampled old UMIs.
        assert len(old_cumi_idx_list) == n
        assert len(new_umi_list) == n
        for i in range(n):
            old_cumi_idxes = old_cumi_idx_list[i]
            new_umis = new_umi_list[i]
            assert len(old_cumi_idxes) == self.X[i, reg_idx]
            assert len(new_umis) == self.X[i, reg_idx]
            for old_cumi_idx, new_umi in zip(old_cumi_idxes, new_umis):
                cell, umi = cells[old_cumi_idx], umis[old_cumi_idx]
                if cell not in self.dat:
                    self.dat[cell] = {}
                if umi not in self.dat[cell]:
                    self.dat[cell][umi] = []
                self.dat[cell][umi].append((i, new_umi, reg_idx_whole))

    def query(self, cell, umi):
        """Query new CUMI based on old cell and UMI barcodes.

        Parameters
        ----------
        cell : str
            The cell barcode from old CUMI.
        umi : str
            The UMI barcode from old CUMI.
        
        Returns
        -------
        list or None
            A list of new CUMIs assigned to the old one. 
            Each element of the list is a tuple(int, int, int), containing
            - (int) The index of new cell, 0-based.
            - (int) The integer format of new UMI barcode, can be transformed
              to string format with `int2str()`.
            - (int) The 0-based index of feature within transcriptomics-scale.
            `None` if the query CUMI (`cell`+`umi`) is not sampled.
        """
        if cell not in self.dat or umi not in self.dat[cell]:
            return(None)
        return(self.dat[cell][umi])

    def sample(self):
        """Sample CUMIs."""
        for reg_idx, reg_idx_whole in enumerate(self.reg_idx_list):
            fn = self.allele_fn_list[reg_idx]
            if is_file_empty(fn):
                continue
            dat = load_cumi(fn, sep = "\t")
            self.__sample_for_feature(
                dat["cell"], dat["umi"], reg_idx, reg_idx_whole)


class MergedSampler:
    """Merged object from all allele-specific CUMI samplers."""

    def __init__(self, samplers):
        """
        Parameters
        ----------
        samplers : dict of {str : CUMISampler}
            Allele-specific CUMI samplers (CUMISampler object).
            Keys are the alleles (str) and values are samplers.
        """
        self.samplers = samplers

    def query(self, cell, umi):
        """Query new CUMI(s) given old cell and UMI IDs.
        
        Parameters
        ----------
        cell : str
            The cell ID. The cell barcode (10x) or sample ID (SMART-seq).
        umi : str
            The read ID. The UMI barcode (10x) or read query name (SMART-seq).
        
        Returns
        -------
        list
            A list of hits. Each element in the list is a tuple:
            str
                The allele where the query CUMI comes from. 
            list
                A list of new CUMI(s) assigned to this query CUMI.
                See the returned value of :func:`~CUMISampler.query()`.
        """
        hits = []
        for allele, sampler in self.samplers.items():
            res = sampler.query(cell, umi)
            if res is None:
                continue
            hits.append((allele, res))
        return(hits)


def cumi_sample_for_cells(m, n_list):
    """Sample CUMIs for a list of samples/cells.

    Parameters
    ----------
    m : int
        The number of CUMIs to sample from.
    n_list : list of int
        Its elements are the numbers of CUMIs to be sampled in each 
        new sample/cell.

    Returns
    -------
    list of numpy.ndarray
        The 0-based indexes of sampled CUMIs for each sample/cell.
        Its length and order match `n_list`.
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
    """Load old CUMIs from file."""
    dat = pd.read_csv(fn, sep = sep, header = None)
    dat.columns = ["cell", "umi"]
    return(dat)
