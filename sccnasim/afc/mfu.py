# mfu.py - multi-feature UMIs.


import gc
import multiprocessing
import os
import pandas as pd
from logging import info
from ..utils.cumi import load_cumi
from ..utils.gfeature import load_feature_objects, load_features
from ..utils.io import load_samples
from ..xlib.xbase import is_file_empty
from ..xlib.xfile import zopen, zconcat, ZF_F_GZIP, ZF_F_PLAIN
from ..xlib.xio import load_pickle, save_pickle
from ..xlib.xthread import split_n2batch, mp_error_handler



def mfu_main(
    alleles,
    multi_mapper_how,
    fet_obj_fn,
    sample_fn,
    feature_fn,
    count_dir,
    tmp_dir,
    out_prefix,
    ncores
):
    """Main function for processing multi-feature UMIs.
    
    This function aims to process multi-feature UMIs, specifically,
    (1) find the multi-feature UMIs from combined allele-specific UMIs of all
        features and process them according to `multi_mapper_how`.
    (2) re-calculate the allele-specific cell x gene count matrices, and
        update the corresponding file paths in each feature object.
    
    Parameters
    ----------
    alleles : list of str
        A list of alleles.
    multi_mapper_how : {"discard", "dup"}
        How to process the multi-feature UMIs (reads).
        - "discard": discard the UMI.
        - "dup": count the UMI for every mapped gene.
    fet_obj_fn : str
        The python pickle file storing `~..utils.gfeature.Feature` objects.
    sample_fn : str
        File storing cell IDs or barcodes.
        Its order should match that in count matrices.
    feature_fn : str
        File storing feature annotations.
    count_dir : str
        The output folder to store the count matrices.
    tmp_dir : str
        Folder to store temporary data.
    out_prefix : str
        The prefix to output files.
    ncores : int
        Number of cores.

    Returns
    -------
    Dict
        The results.
    """
    # check args.
    reg_list = load_feature_objects(fet_obj_fn)
    samples = load_samples(sample_fn)         # ndarray
    features = load_features(feature_fn)      # DataFrame
    features = features["feature"].to_numpy()
    
    assert len(reg_list) == features.shape[0]
    for i, reg in enumerate(reg_list):
        assert reg.name == features[i]
    
    os.makedirs(count_dir, exist_ok = True)
    os.makedirs(tmp_dir, exist_ok = True)
    
    del reg_list
    gc.collect()
    
    step = 1
    
    
    # merge allele-specific CUMI files of all features.
    info("merge allele-specific CUMI files of all features ...")
    
    cumi_merge_fn = os.path.join(count_dir, "%s.cumi.merged.tsv" % out_prefix)
    res_dir = os.path.join(tmp_dir, "%d_cumi_merge" % step)
    os.makedirs(res_dir, exist_ok = True)
    merge_cumis(
        alleles = alleles,
        fet_obj_fn = fet_obj_fn,
        out_fn = cumi_merge_fn,
        tmp_dir = res_dir,
        ncores = ncores
    )
    step += 1
    
    
    # process cell-specific multi-feature CUMIs.
    info("process cell-specific multi-feature CUMIs ...")
    
    cumi_cs_fn = os.path.join(count_dir, "%s.cumi.merged.cs.tsv" % out_prefix)
    res_dir = os.path.join(tmp_dir, "%d_cumi_cs" % step)
    os.makedirs(res_dir, exist_ok = True)
    mfu_cs_main(
        in_fn = cumi_merge_fn,
        out_fn = cumi_cs_fn,
        samples = samples,
        multi_mapper_how = multi_mapper_how,
        tmp_dir = res_dir,
        ncores = ncores
    )
    step += 1
    
    #samples : dict of {str : int}
    #Each key is a cell barcode, value is 0-based index of the cell.
    
    
    
def mfu_cs_main(
    in_fn,
    out_fn,
    samples,
    multi_mapper_how,
    tmp_dir,
    ncores
):
    """Main function of processing cell-specific multi-feature CUMIs.
    
    Parameters
    ----------
    in_fn : str
        Path to the input 4-column CUMI file.
    out_fn : str
        Path to the output 4-column CUMI file.
    samples : list of str
        A list of cell barcodes.
        Its order should match that in count matrices.
    multi_mapper_how : {"discard", "dup"}
        How to process the multi-feature UMIs (reads).
        - "discard": discard the UMI.
        - "dup": count the UMI for every mapped gene.
    tmp_dir : str
        Folder to store temporary data.
    ncores : int
        Number of cores.

    Returns
    -------
    Dict
        The results.
    """
    # check args.
    assert multi_mapper_how in ("discard", "dup")
    os.makedirs(tmp_dir, exist_ok = True)
    
    
    __mfu_cs_batch(
        in_fn = in_fn,
        out_fn = out_fn,
        samples = samples,
        multi_mapper_how = multi_mapper_how,
        tmp_dir = tmp_dir,
        ncores = ncores,
        max_per_batch = 500,
        depth = 0
    )
    

    
def __mfu_cs_batch(
    in_fn,
    out_fn,
    samples,
    multi_mapper_how,
    tmp_dir,
    ncores,
    max_per_batch,
    depth
):
    """Recursive function for `mfu_cs_main()`.
    
    To avoid too many cells (and hence CUMIs) in one batch, this function
    recursively splits large combined file into smaller batches, until the 
    batch size is small than given `max_per_batch`.
    
    Parameters
    ----------
    in_fn : str
        Path to the input 4-column CUMI file.
    out_fn : str
        Path to the output 4-column CUMI file.
    samples : list of str
        A list of cell barcodes.
    multi_mapper_how : {"discard", "dup"}
        How to process the multi-feature UMIs (reads).
        - "discard": discard the UMI.
        - "dup": count the UMI for every mapped gene.
    tmp_dir : str
        Path to folder storing temporary data.
    ncores : int, default 1
        Number of cores.
    max_per_batch : int
        Maximum number of `samples` allowed to be processed simultaneously.
    depth : int
        Depth index, 0-based.

    Returns
    -------
    int
        Return code. 0 if success, negative otherwise.
    """
    n = len(samples)   
    
    if n <= max_per_batch:
        ret = mfu_cs(
            in_fn = in_fn,
            out_fn = out_fn,
            multi_mapper_how = multi_mapper_how
        )
        return(ret)
    
    os.makedirs(tmp_dir, exist_ok = True)


    # split the input CUMI file into smaller batches.
    # Note, here
    # - max_n_batch: to account for the issue of "max open files" when
    #   splitting the large combined file into smaller batches.
    #   It will open every batch-specific splitted file simultaneously, 
    #   in total `n_batch` files.
    bd_m, bd_n, bd_indices = split_n2batch(
        n, ncores, min_per_batch = 200, 
        max_per_batch = max(200, min(500, max_per_batch)), max_n_batch = 300)
    
    fp_list = []
    idx_map = {}
    batches = []
    for idx, (b, e) in enumerate(bd_indices):
        bd_in_fn = os.path.join(tmp_dir, "%d_%d.in.cumi.tsv" % (depth, idx))
        bd_out_fn = os.path.join(tmp_dir, "%d_%d.out.cumi.tsv" % (depth, idx))
        fp = zopen(bd_in_fn, "w", ZF_F_PLAIN)
        fp_list.append(fp)
        for i in range(b, e):
            assert samples[i] not in idx_map
            idx_map[samples[i]] = fp
        batches.append((b, e, bd_in_fn, bd_out_fn))
    
    in_fp = open(in_fn, "r")
    for line in in_fp:
        cell, _, s = line.partition("\t")
        assert cell in idx_map
        fp = idx_map[cell]
        fp.write(line)
    in_fp.close()
    
    for fp in fp_list:
        fp.close()
        

    # next round of extracting and splitting.
    if ncores <= 1:
        for idx, (b, e, bd_in_fn, bd_out_fn) in enumerate(batches):
            res_dir = os.path.join(tmp_dir, "%d_%d" % (depth, idx))
            os.makedirs(res_dir, exist_ok = True)
            __mfu_cs_batch(
                in_fn = bd_in_fn,
                out_fn = bd_out_fn,
                samples = samples[b:e],
                multi_mapper_how = multi_mapper_how,
                tmp_dir = res_dir,
                ncores = 1,
                max_per_batch = max_per_batch,
                depth = depth + 1
            )
    else:
        mp_res = []
        pool = multiprocessing.Pool(processes = min(ncores, bd_m))
        for idx, (b, e, bd_in_fn, bd_out_fn) in enumerate(batches):
            res_dir = os.path.join(tmp_dir, "%d_%d" % (depth, idx))
            os.makedirs(res_dir, exist_ok = True)
            mp_res.append(pool.apply_async(
                func = __mfu_cs_batch,
                kwds = dict(
                    in_fn = bd_in_fn,
                    out_fn = bd_out_fn,
                    samples = samples[b:e],
                    multi_mapper_how = multi_mapper_how,
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
        
        
    # merge batch-specific CUMI files.
    zconcat(
        in_fn_list = [x[3] for x in batches],
        in_format = ZF_F_PLAIN,
        out_fn = out_fn, 
        out_fmode = "w", 
        out_format = ZF_F_PLAIN, 
        remove = False
    )

    del fp_list
    del idx_map
    del batches
    gc.collect()

    return(0)



def mfu_cs(
    in_fn,
    out_fn,
    multi_mapper_how
):
    """Process cell-specific multi-feature CUMIs for one batch of cells.
    
    Parameters
    ----------
    in_fn : str
        Path to the input 4-column CUMI file.
    out_fn : str
        Path to the output 4-column CUMI file.
    multi_mapper_how : {"discard", "dup"}
        How to process the multi-feature UMIs (reads).
        - "discard": discard the UMI.
        - "dup": count the UMI for every mapped gene.
        
    Returns
    -------
    int
        Return code. 0 if success, negative otherwise.
    """
    df = mfu_load_cumi(in_fn)
    if multi_mapper_how == "discard":
        df = df.drop_duplicates(["cell", "umi"], keep = False, 
                                ignore_index = True)
    else:
        pass
    mfu_save_cumi(df, out_fn)
    return(0)

    
    
def merge_cumis(
    alleles,
    fet_obj_fn,
    out_fn,
    tmp_dir,
    ncores
):
    """Merge all allele-specific CUMI files of all features.
    
    Parameters
    ----------
    alleles : list of str
        A list of alleles.
    fet_obj_fn : str
        The python pickle file storing `~..utils.gfeature.Feature` objects.
    out_fn : str
        Output file storing all CUMIs of all features.
    tmp_dir : str
        Folder to store temporary data.
    ncores : int
        Number of cores.

    Returns
    -------
    Void.
    """
    os.makedirs(tmp_dir, exist_ok = True)
    
    
    # split features into batches.
    reg_list = load_feature_objects(fet_obj_fn)
    p = len(reg_list)
    bd_m, bd_n, bd_indices = split_n2batch(
        p, ncores, min_per_batch = 200, max_per_batch = 500)
    
    batches = []
    for idx, (b, e) in enumerate(bd_indices):
        reg_fn = os.path.join(tmp_dir, "%d.features.pickle" % idx)
        save_pickle(reg_list[b:e], reg_fn)
        cumi_fn = os.path.join(tmp_dir, "%d.cumis.tsv" % idx)
        batches.append((reg_fn, cumi_fn))
    del reg_list
    gc.collect()
    
    
    # merge CUMIs in each batch.
    if ncores <= 1:
        for idx, (reg_fn, cumi_fn) in enumerate(batches):
            merge_cumis_batch(
                fet_obj_fn = reg_fn,
                alleles = alleles,
                out_fn = cumi_fn
            )
    else:
        mp_res = []
        pool = multiprocessing.Pool(processes = min(ncores, bd_m))
        for idx, (reg_fn, cumi_fn) in enumerate(batches):
            mp_res.append(pool.apply_async(
                func = merge_cumis_batch,
                kwds = dict(
                    fet_obj_fn = reg_fn,
                    alleles = alleles,
                    out_fn = cumi_fn
                ),
                callback = None,
                error_callback = mp_error_handler
            ))
        pool.close()
        pool.join()
        
        
    # merge batch-specific CUMI files.
    zconcat(
        in_fn_list = [x[1] for x in batches],
        in_format = ZF_F_PLAIN, 
        out_fn = out_fn, 
        out_fmode = "w", 
        out_format = ZF_F_PLAIN, 
        remove = False
    )
    
    
    
def merge_cumis_batch(
    fet_obj_fn,
    alleles,
    out_fn
):
    """Merge CUMIs for a batch of features.
    
    Parameters
    ----------
    fet_obj_fn : str
        The python pickle file storing `~..utils.gfeature.Feature` objects.
    alleles : list of str
        A list of alleles.
    out_fn : str
        Path to output file.

    Returns
    -------
    Void.
    """
    reg_list = load_feature_objects(fet_obj_fn)
    dat = []
    for idx, reg in enumerate(reg_list):
        df = merge_cumis_feature(
            fn_list = [reg.allele_data[ale].seed_cumi_fn for ale in alleles],
            alleles = alleles,
            name = reg.name
        )
        dat.append(df)
    df = pd.concat(dat, ignore_index = True)
    mfu_save_cumi(df, out_fn)



def merge_cumis_feature(
    fn_list,
    alleles,
    name
):
    """Merge CUMI files for one feature.
    
    Parameters
    ----------
    fn_list : list of str
        A list of allele-specific CUMI files.
    alleles : list of str
        A list of alleles.
    name : str
        Feature name.

    Returns
    -------
    pandas.DataFrame
        The merged CUMI data, containing four columns:
        - "cell" (str): cell barcode;
        - "umi" (str): UMI barcode;
        - "feature" (str): feature name;
        - "allele" (str): allele
    """
    # check args.
    assert len(fn_list) == len(alleles)
    
    dat = []
    for ale, fn in zip(alleles, fn_list):
        df = load_cumi(fn)
        if df.shape[0] == 0:
            df = pd.DataFrame(columns = ["cell", "umi", "feature", "allele"],
                             dtype = "string")
        else:
            df["feature"] = name
            df["allele"] = ale
        dat.append(df)
    df = pd.concat(dat, ignore_index = True)
    return(df)



def mfu_load_cumi(fn):
    """Load 4-column CUMIs from file."""
    df = None
    if is_file_empty(fn):
        df = pd.DataFrame(columns = ["cell", "umi", "feature", "allele"],
                        dtype = "string")
    else:
        df = pd.read_table(fn, header = None)
        df.columns = ["cell", "umi", "feature", "allele"]
    return(df)


def mfu_save_cumi(df, fn):
    df.to_csv(fn, sep = "\t", header = False, index = False)
