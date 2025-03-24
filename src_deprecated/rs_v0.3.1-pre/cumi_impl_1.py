# cumi.py - cell-umi barcodes.


# Note, this is the first version of re-implementation of gen_cumi().
# - the previous version (v0.3.0 and before) stores all CUMIs of all cells
#   and features in memory, hence taking huge memory, and is also time
#   consuming.
# - while this version uses much less memory, it takes even much longer time
#   than previous version, largely due to the extreamly high frequency of
#   IO operations, e.g., open and close lots of files for multiple times in
#   gen_cumi_thread().


# Discussion:
# 1. consider generaing CUMIs (UMIs) by iterating integers from 1 to N, where
#    N is the total number of CUMIs in the simulated count matrix.
#    Pros:
#      Current strategy of randomly sampling CUMIs is time and memory
#      consuming, becuase it requires storing all CUMIs, wheras proposed 
#      strategy almost stores only one CUMI and should be much more efficient.
#    Cons:
#      The generated UMIs are not randomly distributed, different from real 
#      data.


import anndata as ad
import multiprocessing
import numpy as np
import os
import pandas as pd
import pickle
from ..io.base import load_h5ad, save_h5ad
from ..utils.base import is_file_empty
from ..utils.xbarcode import Barcode
from ..utils.xthread import split_n2m



def __merge_small_files(in_fn_list, out_fn, inplace = False):
    """Merge small files into one."""
    out_fp = open(out_fn, "w")
    for in_fn in in_fn_list:
        with open(in_fn, "r") as in_fp:
            s = in_fp.read()
            out_fp.write(s)
    out_fp.close()
    
    if inplace:
        for fn in in_fn_list:
            os.remove(fn)


def gen_cumi(
    count_fn,
    alleles,
    umi_len,
    out_files,
    tmp_dir,
    ncores = 1
):
    """Generate CUMIs to be used in new BAM.

    This function generates new CUMIs for each feature based on the input
    count matrices in `xdata`.
    
    Parameters
    ----------
    count_fn : str
        Path to `anndata.Anndata` object containing the allele-specific 
        *cell x feature* matrices of simulated UMI counts.
        The alleles are specified in `alleles`.
    alleles : list of str
        A list of alleles, e.g., ["A", "B", "U"].
    umi_len : int
        Length of one new UMI barcode.
    out_files : list of list of str
        Two layers of lists specifying output files storing generated CUMIs.
        The first layer of list matches features in `count_fn` while the 
        second layer of list matches `alleles`.
    tmp_dir : str
        Path to folder storing temporary data.
    ncores : int, default 1
        Number of cores.

    Returns
    -------
    Void.
    """
    # check args.
    xdata = load_h5ad(count_fn)
    n, p = xdata.shape
    
    assert "cell" in xdata.obs
        
    for ale in alleles:
        assert ale in xdata.layers
        
    assert umi_len <= 31
    
    assert len(out_files) == p
    for out_fns in out_files:
        assert len(out_fns) == len(alleles)
        
    os.makedirs(tmp_dir, exist_ok = True)


    # split cells for multi-processing
    td_m, td_n, td_cell_indices = split_n2m(n, ncores)
    count_fn_list = []
    for idx, (b, e) in enumerate(td_cell_indices):
        fn = os.path.join(tmp_dir, "%d.adata.h5ad" % idx)
        xdata_batch = xdata[b:e, :]
        save_h5ad(xdata_batch, fn)
        count_fn_list.append(fn)
    del xdata
    
    # prepare CUMI files for each thread.
    td_out_fn_list = []     # three layers: thread - feature - allele
    for idx in range(td_m):
        td_fn_list = []
        for fet_fn_list in out_files:
            td_fn_list.append(["%s.%d" % (ale_fn, idx)  \
                for ale_fn in fet_fn_list])
        td_out_fn_list.append(td_fn_list)

    # multi-processing
    mp_res = []
    pool = multiprocessing.Pool(processes = td_m)
    for idx in range(td_m):
        mp_res.append(pool.apply_async(
            func = gen_cumi_thread,
            kwds = dict(
                idx = idx,
                count_fn = count_fn_list[idx],
                alleles = alleles,
                umi_len = umi_len,
                out_files = td_out_fn_list[idx]
            ),
            callback = None
        ))
    pool.close()
    pool.join()

    # merge feature-specific UMIs
    for j in range(p):
        for k in range(len(alleles)):
            fn = out_files[j][k]
            __merge_small_files(
                in_fn_list = [td_files[j][p] for td_files in td_out_fn_list],
                out_fn = fn,
                inplace = False
            )


def gen_cumi_thread(
    idx,
    count_fn, 
    alleles,
    umi_len,
    out_files
):
    """Generate CUMIs to be used in new BAM for a batch of cells.
    
    Parameters
    ----------
    count_fn : str
        Path to `anndata.Anndata` object containing the allele-specific 
        *cell x feature* matrices of simulated UMI counts.
        The alleles are specified in `alleles`.
    alleles : list of str
        A list of alleles, e.g., ["A", "B", "U"].
    umi_len : int
        Length of one new UMI barcode.
    out_files : list of list of str
        Two layers of lists specifying output files storing generated CUMIs.
        The first layer of list matches features in `count_fn` while the 
        second layer of list matches `alleles`.

    Returns
    -------
    Void.
    """
    # check args.
    xdata = load_h5ad(count_fn)
    n, p = xdata.shape
    
    assert "cell" in xdata.obs
    cell_list = xdata.obs["cell"].values
    
    assert len(out_files) == p
    for out_fns in out_files:
        assert len(out_fns) == len(alleles)
    
    RD = None
    for i, ale in enumerate(alleles):
        assert ale in xdata.layers
        if i == 0:
            RD = xdata.layers[ale].copy()
        else:
            RD += xdata.layers[ale]

    # UMIs are generated in a cell-specific manner, mimicking real data that
    # the UMI of each transcript should be unique within one cell.
    b = Barcode(umi_len)
    for i in range(n):
        cell = cell_list[i]
        uint_list = b.sample_int(n = np.sum(RD[i, :]), sort = False)
        r = 0
        fmode = "w" if i == 0 else "a"
        for j in range(p):
            ale_fns = out_files[j]
            for k, ale in enumerate(alleles):
                x = xdata.layers[ale][i, j]
                ale_uint_list = uint_list[r:(r+x)]
                r += x
                
                fn = ale_fns[k]
                with open(fn, fmode) as fp:
                    s = ""
                    for uint in ale_uint_list:
                        umi = b.int2str(uint)
                        s += "%s\t%s\n" % (cell, umi)
                    fp.write(s)
