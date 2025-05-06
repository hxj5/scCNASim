# mfu.py - multi-feature UMIs.


import gc
import multiprocessing
import os
import pandas as pd
from logging import info
from ..utils.cumi import load_cumi
from ..utils.gfeature import load_feature_objects, load_features
from ..utils.io import load_samples
from ..xlib.xfile import zconcat, ZF_F_GZIP, ZF_F_PLAIN
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
    
    cumi_fn = os.path.join(count_dir, "%s.merged.cumi.tsv" % out_prefix)
        
    res_dir = os.path.join(tmp_dir, "%d_merge_cumi" % step)
    os.makedirs(res_dir, exist_ok = True)
    merge_cumis(
        alleles = alleles,
        fet_obj_fn = fet_obj_fn,
        out_fn = cumi_fn,
        tmp_dir = res_dir,
        ncores = ncores
    )
    step += 1
    
    #samples : dict of {str : int}
    #Each key is a cell barcode, value is 0-based index of the cell. 

    
    
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
        batches.append((b, reg_fn, cumi_fn))
    del reg_list
    gc.collect()
    
    
    # merge CUMIs in each batch.
    if ncores <= 1:
        for idx, (b, reg_fn, cumi_fn) in enumerate(batches):
            merge_cumis_batch(
                fet_obj_fn = reg_fn,
                b = b,
                alleles = alleles,
                out_fn = cumi_fn
            )
    else:
        mp_res = []
        pool = multiprocessing.Pool(processes = min(ncores, bd_m))
        for idx, (b, reg_fn, cumi_fn) in enumerate(batches):
            mp_res.append(pool.apply_async(
                func = merge_cumis_batch,
                kwds = dict(
                    fet_obj_fn = reg_fn,
                    b = b,
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
        in_fn_list = [x[2] for x in batches],
        in_format = ZF_F_PLAIN, 
        out_fn = out_fn, 
        out_fmode = "w", 
        out_format = ZF_F_PLAIN, 
        remove = False
    )      
    
    
    
def merge_cumis_batch(
    fet_obj_fn,
    b,
    alleles,
    out_fn
):
    """Merge CUMIs for a batch of features.
    
    Parameters
    ----------
    fet_obj_fn : str
        The python pickle file storing `~..utils.gfeature.Feature` objects.
    b : int
        The transcriptomics-scale 0-based index of the first feature in this 
        batch.
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
    df.to_csv(out_fn, sep = "\t", header = False, index = False)



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
