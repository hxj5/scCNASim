# cumi.py - cell-umi barcodes.


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
from logging import info, error
from .io import merge_tsv
from ..io.base import load_h5ad, save_h5ad
from ..utils.base import is_file_empty
from ..utils.xbarcode import Barcode
from ..utils.xthread import split_n2m
from ..utils.zfile import zopen, ZF_F_PLAIN



def __split_features(
    in_fn,
    out_dir,
    out_prefix,
    p,
    ncores,
    max_per_batch
):
    """Split CUMIs of all features into batches.
    
    Parameters
    ----------
    in_fn : str
        Path to file storing CUMIs of all features.
    out_dir : str
        Path to splitted files of CUMIs.
    out_prefix : str
        Prefix to output batch files.
    p : int
        Number of features.
    ncores : int
        Number of cores.
    max_per_batch : int
        Maximum features per batch.
        
    Returns
    -------
    list of tuple
        Information of splitted batches:
        int
            transcriptomics-scale index (0-based, inclusive) of the first
            feature within this batch.
        int
            transcriptomics-scale index (0-based, inclusive) of the last
            feature within this batch.
        str
            Path to the file storing CUMIs of features in this batch.
    """
    n_batch = p // max_per_batch + 1 - (p % max_per_batch == 0)
    n_batch = max(n_batch, ncores)
    
    td_m, td_n, td_indices = split_n2m(p, n_batch)
    td_fp_list = []
    idx_map = {}
    res = []
    for idx, (b, e) in enumerate(td_indices):
        td_fn = os.path.join(out_dir, "%s%d.cumi.tsv" % (out_prefix, idx))
        td_fp = zopen(td_fn, "w", ZF_F_PLAIN)
        td_fp_list.append(td_fp)
        for i in range(b, e):
            assert i not in idx_map
            idx_map[i] = td_fp
        res.append((b, e - 1, td_fn))
    
    in_fp = open(in_fn, "r")
    for line in in_fp:
        fet, _, cumi = line.partition("\t")
        fet = int(fet)
        assert fet in idx_map
        fp = idx_map[fet]
        fp.write(line)
    in_fp.close()
    
    for fp in td_fp_list:
        fp.close()
        
    return(res)



def extract_simu_cumi(
    in_fn,
    out_files,
    b,
    e
):
    """Extract CUMIs for each feature.
    
    Parameters
    ----------
    in_fn : str
        Path to the file storing CUMIs of a batch of features.
    out_files : list of str
        A list of feature-specific CUMI files.
    b : int
        The transcriptomics-scale index of the first feature in this batch.
        0-based, inclusive.
    e : int
        The transcriptomics-scale index of the last feature in this batch.
        0-based, inclusive.
        
    Returns
    -------
    int
        Return code. 0 if success, negative otherwise.
    """
    ret = -1
    
    # check args.
    assert len(out_files) == e - b + 1

    fet_fp_list = [zopen(out_fn, "w", ZF_F_PLAIN) for out_fn in out_files]
    idx_map = {}
    for i, fp in zip(range(b, e + 1), fet_fp_list):
        idx_map[i] = fp

    in_fp = open(in_fn, "r")
    for line in in_fp:
        fet, _, cumi = line.partition("\t")
        fet = int(fet)
        assert fet in idx_map
        fp = idx_map[fet]
        fp.write(cumi)
    in_fp.close()
    
    for fp in fet_fp_list:
        fp.close()
        
    ret = 0
    return(ret)



def gen_simu_cumi(
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
        The first layer of list matches `alleles` while the second layer of
        list matches features in `count_fn`.
    tmp_dir : str
        Path to folder storing temporary data.
    ncores : int, default 1
        Number of cores.

    Returns
    -------
    int
        Return code. 0 if success, negative otherwise.
    """
    ret = -1
    
    # check args.
    xdata = load_h5ad(count_fn)
    n, p = xdata.shape
    
    assert "cell" in xdata.obs
        
    for ale in alleles:
        assert ale in xdata.layers
        
    assert umi_len <= 31
    
    assert len(out_files) == len(alleles)
    for out_fns in out_files:
        assert len(out_fns) == p
        
    os.makedirs(tmp_dir, exist_ok = True)


    # split cells for multi-processing
    td_m, td_n, td_cell_indices = split_n2m(n, ncores)
    cell_count_fn_list = []
    for idx, (b, e) in enumerate(td_cell_indices):
        fn = os.path.join(tmp_dir, "cell.b%d.adata.h5ad" % idx)
        xdata_batch = xdata[b:e, :]
        save_h5ad(xdata_batch, fn)
        cell_count_fn_list.append(fn)
    del xdata
    
    
    # prepare CUMI files for each thread.
    cell_cumi_fn_list = []      # two layers: thread - allele
    for idx in range(td_m):
        fn_list = []
        for ale in alleles:
            fn = os.path.join(tmp_dir, "cell.b%d.%s.cumi.tsv" % (idx, ale))
            fn_list.append(fn)
        cell_cumi_fn_list.append(fn_list)


    # multi-processing for generating CUMIs.
    mp_res = []
    pool = multiprocessing.Pool(processes = td_m)
    for idx in range(td_m):
        mp_res.append(pool.apply_async(
            func = gen_simu_cumi_thread,
            kwds = dict(
                idx = idx,
                count_fn = cell_count_fn_list[idx],
                alleles = alleles,
                umi_len = umi_len,
                out_files = cell_cumi_fn_list[idx]
            ),
            callback = None
        ))
    pool.close()
    pool.join()
    
    ret_list = [res.get() for res in mp_res]
    for r in ret_list:
        if r != 0:
            return(-3)
    
    
    # merge CUMI files.
    ale_cumi_fn_list = []
    for idx, ale in enumerate(alleles):
        fn = os.path.join(tmp_dir, "%s.cumi.tsv" % ale)
        merge_tsv(
            in_fn_list = [fn_list[idx] for fn_list in cell_cumi_fn_list],
            out_fn = fn,
            remove = True
        )
        ale_cumi_fn_list.append(fn)


    # clean tmp files.
    for fn in cell_count_fn_list:
        os.remove(fn)
        
        
    # split features for multi-processing.
    fet_cumi_dir = os.path.join(tmp_dir, "tmp_fet_cumi")
    os.makedirs(fet_cumi_dir, exist_ok = True)
    
    fet_cumi_fn_list = []
    fet_ale_cumi_dirs = []
    for idx, ale in enumerate(alleles):
        ale_cumi_dir = os.path.join(fet_cumi_dir, ale)
        os.makedirs(ale_cumi_dir, exist_ok = True)
        fet_ale_cumi_dirs.append(ale_cumi_dir)
        
        fn_list_batch = __split_features(
            in_fn = ale_cumi_fn_list[idx],
            out_dir = ale_cumi_dir,
            out_prefix = "%s.fet.b" % ale,
            p = p,
            ncores = ncores,
            max_per_batch = 500       # to avoid "max open files" issue.
        )
        fet_cumi_fn_list.append(fn_list_batch)


    # multi-processing to extract feature-specific CUMIs.
    mp_res = []
    pool = multiprocessing.Pool(processes = ncores)
    for idx, ale in enumerate(alleles):
        fn_list_batch = fet_cumi_fn_list[idx]
        for fn_batch in fn_list_batch:
            b, e, cumi_fn = fn_batch
            mp_res.append(pool.apply_async(
                func = extract_simu_cumi,
                kwds = dict(
                    in_fn = cumi_fn,
                    out_files = out_files[idx][b:(e+1)],
                    b = b,
                    e = e
                ),
                callback = None
            ))
    pool.close()
    pool.join()
    
    ret_list = [res.get() for res in mp_res]
    for r in ret_list:
        if r != 0:
            return(-5)

    
    # clean tmp files.
    # here we do not use methods like shutil.rmtree() to remove entire folder,
    # to avoid removing some files unexpectedly.
    for fn_list in fet_cumi_fn_list:
        for b, e, fn in fn_list:
            os.remove(fn)
    for ale_dir in fet_ale_cumi_dirs:
        os.rmdir(ale_dir)
    os.rmdir(fet_cumi_dir)

    ret = 0
    return(ret)



def gen_simu_cumi_thread(
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
    out_files : list of str
        A list of allele-specific files storing generated CUMIs.

    Returns
    -------
    int
        Return code. 0 if success, negative otherwise.
    """
    ret = -1
    
    # check args.
    xdata = load_h5ad(count_fn)
    n, p = xdata.shape
    
    assert "cell" in xdata.obs
    cell_list = xdata.obs["cell"].values
    
    assert len(out_files) == len(alleles) 
    
    RD = None
    for i, ale in enumerate(alleles):
        assert ale in xdata.layers
        if i == 0:
            RD = xdata.layers[ale].copy()
        else:
            RD += xdata.layers[ale]

    # UMIs are generated in a cell-specific manner, mimicking real data that
    # the UMI of each transcript should be unique within one cell.
    fp_list = [zopen(fn, "w") for fn in out_files]
    b = Barcode(umi_len)
    for i in range(n):
        cell = cell_list[i]
        uint_list = b.sample_int(n = np.sum(RD[i, :]), sort = False)
        r = 0
        for j in range(p):
            for k, (ale, fp) in enumerate(zip(alleles, fp_list)):
                x = xdata.layers[ale][i, j]
                ale_uint_list = uint_list[r:(r+x)]
                r += x
                s = ""
                for uint in ale_uint_list:
                    umi = b.int2str(uint)
                    s += "%d\t%s\t%s\n" % (j, cell, umi)
                fp.write(s)
                
    for fp in fp_list:
        fp.close()
        
    ret = 0
    return(ret)


        
def load_cumi(fn, sep = "\t"):
    """Load CUMIs from file."""
    if is_file_empty(fn):
        df = pd.DataFrame(columns = ["cell", "umi"])
        return(df)
    dat = pd.read_csv(fn, sep = sep, header = None)
    dat.columns = ["cell", "umi"]
    return(dat)



def smpl_seed_cumi(
    count_fn,
    feature_fn,
    alleles
):
    """Sampling feature-specific CUMIs of seed data.
    
    Parameters
    ----------
    count_fn : str
        Path to ".h5ad" file storing allele-specifc *cell x feature* count
        matrices.
    feature_fn : str
        Path to python pickle file storing a list of features, i.e., the
        `~utils.gfeature.Feature` objects.
    alleles : list of str
        A list of alleles.

    Returns
    -------
    int
        Return code. 0 if success, negative otherwise.
    """
    ret = -1
    
    # check args.
    adata = load_h5ad(count_fn)
    
    with open(feature_fn, "rb") as fp:
        reg_list = pickle.load(fp)
        
    assert "feature" in adata.var
    assert np.all(adata.var["feature"].values == \
                  [reg.name for reg in reg_list])
    
    for ale in alleles:
        assert ale in adata.layers
        
    # sampling CUMIs.
    for idx, reg in enumerate(reg_list):
        for ale in alleles:
            in_fn = reg.allele_data[ale].seed_cumi_fn
            out_fn = reg.allele_data[ale].seed_smpl_cumi_fn
            cumis = load_cumi(in_fn)
            x = adata.layers[ale][:, idx]
            if cumis.shape[0] == 0:
                assert np.sum(x) == 0
                with open(out_fn, "w") as fp:
                    pass
                continue
            smpl_cumi(
                cells = cumis["cell"].values,
                umis = cumis["umi"].values,
                counts = x,
                out_fn = out_fn
            )
            
    ret = 0
    return(ret)



def smpl_cumi(
    cells,
    umis,
    counts,
    out_fn
):
    """Sampling feature-specific CUMIs for simulated cells and output to file.
    
    Parameters
    ----------
    cells : numpy.ndarray
        Cell barcodes of CUMIs to be sampled from.
    umis : numpy.ndarray
        UMI barcodes of CUMIs to be sampled from.
    counts : numpy.ndarray
        Number of CUMIs to sample for each newly simulated cell.
    out_fn : str
        Path to the output CUMI file.
        
    Returns
    -------
    Void.
    """
    assert len(cells) == len(umis)
    
    m = len(cells)
    n = np.sum(counts)

    indices = np.random.choice(
        range(m), 
        size = n,
        replace = n > m
    )

    fp = zopen(out_fn, "w", ZF_F_PLAIN)
    for idx in indices:
        fp.write("%s\t%s\n" % (cells[idx], umis[idx]))
    fp.close()
