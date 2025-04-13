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


import gc
import multiprocessing
import numpy as np
import os
import pandas as pd
import shutil

from logging import info
from .io import merge_tsv
from ..io.base import load_h5ad, save_h5ad, load_feature_objects
from ..utils.base import is_file_empty
from ..utils.xbarcode import Barcode
from ..utils.xdata import sum_layers, array_to_sparse
from ..utils.xthread import split_n2batch, mp_error_handler
from ..utils.zfile import zopen, ZF_F_PLAIN



def cumi_simu_main(
    count_fn,
    n,
    p,
    alleles,
    umi_len,
    out_files,
    tmp_dir,
    ncores = 1
):
    """Main function of simulating *cell x feature* CUMIs for every allele.

    This function generates simulated *cell x feature* CUMIs for every allele
    based on the input count matrices stored in `count_fn`.
    
    Parameters
    ----------
    count_fn : str
        Path to `anndata.Anndata` object containing the allele-specific 
        *cell x feature* UMI count matrices.
        The alleles are specified in `alleles`.
    n : int
        Number of cells in `count_fn`.
    p : int
        Number of features in `count_fn`.
    alleles : list of str
        A list of alleles, e.g., ["A", "B", "U"].
    umi_len : int
        Length of UMI barcodes.
    out_files : list of list of str
        Output files storing generated CUMIs, specified in two layers of lists
        that match `alleles` and the features in `count_fn`, respectively.
        Each file contains simulated allele- and feature-specific CUMIs from
        all cells.
    tmp_dir : str
        Path to folder storing temporary data.
    ncores : int, default 1
        Number of cores.

    Returns
    -------
    int
        Return code. 0 if success, negative otherwise.
    """
    info("start ...")

    # check args.
    assert umi_len <= 31
    
    assert len(out_files) == len(alleles)
    for fns in out_files:
        assert len(fns) == p
        
    os.makedirs(tmp_dir, exist_ok = True)
    
    
    # generate *cell x feature* CUMIs for every allele.
    info("generate cell x feature CUMIs for every allele ...")
    
    ale_cumi_fn_list = []
    for idx, ale in enumerate(alleles):
        fn = os.path.join(tmp_dir, "%s.cumi.tsv" % ale)
        ale_cumi_fn_list.append(fn)
    
    cs_dir = os.path.join(tmp_dir, "tmp_cs")
    os.makedirs(cs_dir, exist_ok = True)
    
    if cumi_simu_cs_main(
        count_fn = count_fn,
        alleles = alleles,
        umi_len = umi_len,
        out_files = ale_cumi_fn_list,
        tmp_dir = cs_dir,
        ncores = ncores
    ) < 0:
        return(-3)
    
    shutil.rmtree(cs_dir)
    
    
    # extract feature-specific CUMIs for every allele.
    info("extract feature-specific CUMIs for every allele ...")
    
    for ale, fn, out_fns in zip(alleles, ale_cumi_fn_list, out_files):
        info("extract for allele '%s' ..." % ale)

        fs_dir = os.path.join(tmp_dir, "tmp_fs_%s" % ale)
        os.makedirs(fs_dir, exist_ok = True)
        cumi_extract_fs_main(
            in_fn = fn,
            out_files = out_fns,
            tmp_dir = fs_dir,
            ncores = ncores
        )
        shutil.rmtree(fs_dir)
        
    return(0)



def cumi_simu_cs_main(
    count_fn,
    alleles,
    umi_len,
    out_files,
    tmp_dir,
    ncores = 1
):
    """Main function for simulating CUMIs in cell-specific manner for
    every allele.
    
    Parameters
    ----------
    count_fn : str
        Path to `anndata.Anndata` object containing the allele-specific 
        *cell x feature* UMI count matrices.
        The alleles are specified in `alleles`.
    alleles : list of str
        A list of alleles, e.g., ["A", "B", "U"].
    umi_len : int
        Length of UMI barcodes.
    out_files : list of str
        A list of allele-specific files storing generated CUMIs of all cells
        and features.
        Its length and order match `alleles`.
    tmp_dir : str
        Path to folder storing temporary data.
    ncores : int, default 1
        Number of cores.

    Returns
    -------
    int
        Return code. 0 if success, negative otherwise.
    """
    # check args.
    adata = load_h5ad(count_fn)
    n, p = adata.shape
    
    assert "cell" in adata.obs
    assert len(alleles) == len(out_files)
    for ale in alleles:
        assert ale in adata.layers
        
        
    # here use "csr" to make row (cell) slicing efficient.
    adata = array_to_sparse(adata, which = "csr", layers = alleles)

    
    # split cells for multi-processing
    # Note, here
    # - max_n_batch: to account for the max allowed files and subfolders in
    #   one folder.
    #   Currently, 4 files output in each batch.
    bd_m, bd_n, bd_cell_indices = split_n2batch(
        n, ncores, max_n_batch = 7000)
    cell_count_fn_list = []
    for idx, (b, e) in enumerate(bd_cell_indices):
        fn = os.path.join(tmp_dir, "cell.b%d.adata.h5ad" % idx)
        adata_batch = adata[b:e, :]
        save_h5ad(adata_batch, fn)
        cell_count_fn_list.append(fn)
    del adata
    
    
    # prepare CUMI files for each batch.
    cell_cumi_fn_list = []      # two layers: batch - allele
    for idx in range(bd_m):
        fn_list = []
        for ale in alleles:
            fn = os.path.join(tmp_dir, "cell.b%d.%s.cumi.tsv" % (idx, ale))
            fn_list.append(fn)
        cell_cumi_fn_list.append(fn_list)


    # multi-processing for generating CUMIs in each batch.
    mp_res = []
    pool = multiprocessing.Pool(processes = min(ncores, bd_m))
    for idx in range(bd_m):
        mp_res.append(pool.apply_async(
            func = cumi_simu_cs,
            kwds = dict(
                count_fn = cell_count_fn_list[idx],
                alleles = alleles,
                umi_len = umi_len,
                out_files = cell_cumi_fn_list[idx]
            ),
            callback = None,
            error_callback = mp_error_handler
        ))
    pool.close()
    pool.join()
    
    ret_list = [res.get() for res in mp_res]
    for ret in ret_list:
        if ret != 0:
            return(-3)
    
    
    # merge CUMI files.
    for idx, ale in enumerate(alleles):
        merge_tsv(
            in_fn_list = [fn_list[idx] for fn_list in cell_cumi_fn_list],
            out_fn = out_files[idx],
            remove = True
        )


    # clean tmp files.
    for fn in cell_count_fn_list:
        os.remove(fn)
        
    return(0)



def cumi_simu_cs(
    count_fn, 
    alleles,
    umi_len,
    out_files
):
    """Simulate CUMIs in cell-specific manner for every allele.
    
    Parameters
    ----------
    count_fn : str
        Path to `anndata.Anndata` object containing the allele-specific 
        *cell x feature* UMI count matrices.
        The alleles are specified in `alleles`.
    alleles : list of str
        A list of alleles, e.g., ["A", "B", "U"].
    umi_len : int
        Length of UMI barcodes.
    out_files : list of str
        A list of allele-specific files storing generated CUMIs of all cells
        and features.
        Its length and order match `alleles`.

    Returns
    -------
    int
        Return code. 0 if success, negative otherwise.
    """
    # check args.
    adata = load_h5ad(count_fn)
    n, p = adata.shape
    
    assert "cell" in adata.obs
    cell_list = adata.obs["cell"].values
    
    assert len(out_files) == len(alleles) 
    
    RD = sum_layers(adata, layers = alleles)

    # UMIs are generated in a cell-specific manner, mimicking real data that
    # the UMI of each transcript should be unique within one cell.
    fp_list = [zopen(fn, "w") for fn in out_files]
    b = Barcode(umi_len)
    for i in range(n):
        cell = cell_list[i]
        uint_list = b.sample_int(n = RD[i, :].sum(), sort = False)
        r = 0
        for j in range(p):
            for k, (ale, fp) in enumerate(zip(alleles, fp_list)):
                x = adata.layers[ale][i, j]
                ale_uint_list = uint_list[r:(r+x)]
                r += x
                s = ""
                for uint in ale_uint_list:
                    umi = b.int2str(uint)
                    s += "%d\t%s\t%s\n" % (j, cell, umi)
                fp.write(s)
                
    for fp in fp_list:
        fp.close()

    del adata
    del RD
    gc.collect()

    return(0)



def cumi_extract_fs_main(
    in_fn,
    out_files,
    tmp_dir,
    ncores
):
    """Main function for extracting feature-specific CUMIs from the combined
    file (of specific allele).
    
    Parameters
    ----------
    in_fn : str
        File storing CUMIs of combined features.
    out_files : list of str
        A list of output files storing feature-specific CUMIs from all cells.
    tmp_dir : str
        Path to folder storing temporary data.
    ncores : int, default 1
        Number of cores.

    Returns
    -------
    int
        Return code. 0 if success, negative otherwise.
    """
    ret = __cumi_extract_fs_batch(
        in_fn = in_fn,
        b0 = 0,
        e0 = len(out_files) - 1,
        out_files = out_files,
        tmp_dir = tmp_dir,
        ncores = ncores,
        max_per_batch = 300,
        depth = 0
    )
    return(ret)



def __cumi_extract_fs_batch(
    in_fn,
    b0,
    e0,
    out_files,
    tmp_dir,
    ncores,
    max_per_batch,
    depth
):
    """Recursive function for `cumi_extract_fs_main()`.
    
    To avoid the issue of `max open files` in cumi_extract_fs(), this function
    recursively splits large combined file into smaller batches, until the 
    batch size is small than given `max_per_batch`.
    
    Parameters
    ----------
    in_fn : str
        File storing CUMIs of combined features.
    b0 : int
        The transcriptomics-scale index of the first feature in this batch.
        0-based, inclusive.
    e0 : int
        The transcriptomics-scale index of the last feature in this batch.
        0-based, inclusive.    
    out_files : list of str
        A list of output files storing feature-specific CUMIs from all cells.
    tmp_dir : str
        Path to folder storing temporary data.
    ncores : int, default 1
        Number of cores.
    max_per_batch : int
        Maximum number of features allowed to be processed simultaneously.
    depth : int
        Depth index, 0-based.

    Returns
    -------
    int
        Return code. 0 if success, negative otherwise.
    """
    p = len(out_files)
    assert p == e0 - b0 + 1
    
    if p <= max_per_batch:
        ret = cumi_extract_fs(
            in_fn = in_fn,
            out_files = out_files,
            b = b0,
            e = e0
        )
        return(ret)


    # split the input CUMI file into smaller batches.
    # Note, here
    # - max_n_batch: to account for the issue of "max open files" when
    #   splitting the large combined file into smaller batches.
    #   It will open every batch-specific splitted file simultaneously, 
    #   in total `n_batch` files.
    bd_m, bd_n, bd_indices = split_n2batch(
        p, ncores, min_n_batch = 50, max_n_batch = 300)
    
    bd_fp_list = []
    idx_map = {}
    bd_batches = []
    for idx, (b, e) in enumerate(bd_indices):
        b += b0
        e += b0
        fn = os.path.join(tmp_dir, "%d_%d.cumi.tsv" % (depth, idx))
        fp = zopen(fn, "w", ZF_F_PLAIN)
        bd_fp_list.append(fp)
        for reg_idx in range(b, e):
            assert reg_idx not in idx_map
            idx_map[reg_idx] = fp
        bd_batches.append((b, e - 1, fn))
    
    in_fp = open(in_fn, "r")
    for line in in_fp:
        reg_idx, _, cumi = line.partition("\t")
        reg_idx = int(reg_idx)
        assert reg_idx in idx_map
        fp = idx_map[reg_idx]
        fp.write(line)
    in_fp.close()
    
    for fp in bd_fp_list:
        fp.close()
        

    # next round of extracting and splitting.
    if ncores <= 1:
        for idx, (b, e, fn) in enumerate(bd_batches):
            res_dir = os.path.join(tmp_dir, "%d_%d" % (depth, idx))
            os.makedirs(res_dir, exist_ok = True)
            __cumi_extract_fs_batch(
                in_fn = fn,
                b0 = b,
                e0 = e,
                out_files = out_files[(b-b0):(e-b0+1)],
                tmp_dir = res_dir,
                ncores = 1,
                max_per_batch = max_per_batch,
                depth = depth + 1
            )
    else:
        mp_res = []
        pool = multiprocessing.Pool(processes = min(ncores, bd_m))
        for idx, (b, e, fn) in enumerate(bd_batches):
            res_dir = os.path.join(tmp_dir, "%d_%d" % (depth, idx))
            os.makedirs(res_dir, exist_ok = True)
            mp_res.append(pool.apply_async(
                func = __cumi_extract_fs_batch,
                kwds = dict(
                    in_fn = fn,
                    b0 = b,
                    e0 = e,
                    out_files = out_files[(b-b0):(e-b0+1)],
                    tmp_dir = res_dir,
                    ncores = 1,
                    max_per_batch = max_per_batch,
                    depth = depth + 1
                ),
                callback = None,
                error_callback = mp_error_handler
            ))
        pool.close()
        pool.join()
        
        
    del bd_fp_list
    del idx_map
    del bd_batches
    gc.collect()

    return(0)



def cumi_extract_fs(
    in_fn,
    out_files,
    b,
    e
):
    """Extract feature-specific CUMIs from combined file.
    
    Parameters
    ----------
    in_fn : str
        Path to the file storing CUMIs of combined features.
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
    # check args.
    assert len(out_files) == e - b + 1

    reg_fp_list = [zopen(fn, "w", ZF_F_PLAIN) for fn in out_files]
    idx_map = {}
    for reg_idx, reg_fp in zip(range(b, e + 1), reg_fp_list):
        idx_map[reg_idx] = reg_fp

    in_fp = open(in_fn, "r")
    for line in in_fp:
        reg_idx, _, cumi = line.partition("\t")
        reg_idx = int(reg_idx)
        assert reg_idx in idx_map
        reg_fp = idx_map[reg_idx]
        reg_fp.write(cumi)
    in_fp.close()
    
    for fp in reg_fp_list:
        fp.close()

    return(0)



def load_cumi(fn, sep = "\t"):
    """Load CUMIs from file."""
    if is_file_empty(fn):
        df = pd.DataFrame(columns = ["cell", "umi"])
        return(df)
    dat = pd.read_csv(fn, sep = sep, header = None)
    dat.columns = ["cell", "umi"]
    return(dat)



def cumi_sample_seed_main(
    count_fn,
    feature_fn,
    alleles,
    index = None
):
    """Main function for sampling feature-specific CUMIs of seed data for
    every allele.
    
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
    index : int or None
        The index of the batch.
        None means there is no batch information.

    Returns
    -------
    int
        Return code. 0 if success, negative otherwise.
    """
    if index is not None:
        info("[Batch-%d] start ..." % index)

    # check args.
    adata = load_h5ad(count_fn)
    
    # here use "csr" to make column (feature) slicing efficient.
    adata = array_to_sparse(adata, which = "csc", layers = alleles)
    
    reg_list = load_feature_objects(feature_fn)
        
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
                assert x.sum() == 0
                with open(out_fn, "w") as fp:
                    pass
                continue
            cumi_sample(
                cells = cumis["cell"].values,
                umis = cumis["umi"].values,
                counts = x,
                out_fn = out_fn
            )
            
    del adata
    del reg_list
    gc.collect()
    
    if index is not None:
        info("[Batch-%d] done!" % index)

    return(0)



def cumi_sample(
    cells,
    umis,
    counts,
    out_fn
):
    """Sampling CUMIs based on given counts.
    
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
    n = counts.sum()

    indices = np.random.choice(
        range(m), 
        size = n,
        replace = n > m
    )

    fp = zopen(out_fn, "w", ZF_F_PLAIN)
    for idx in indices:
        fp.write("%s\t%s\n" % (cells[idx], umis[idx]))
    fp.close()
