# xsam.py - sam alignment processing.


import gc
import multiprocessing
import os
import pysam
import shutil
import subprocess

from logging import info, error
from .xbase import is_vector
from .xio import file2list, list2file
from .xthread import split_n2batch, mp_error_handler



def sam_cat(in_fns, out_fn, ncores = 1):
    """Concatenate a list of BAM files.

    Parameters
    ----------
    in_fns : str or list of str
        If str, path to the file listing BAM files to be concatenated;
        If list of str, a list of BAM files to be concatenated.
    out_fn : str
        Path to the output concatenated BAM file.
    ncores : int
        Number of cores.
    
    Returns
    -------
    Void.
    """
    is_input_list = is_vector(in_fns)
    
    list_fn = None
    if is_input_list:
        list_fn = out_fn + ".lst"
        list2file(in_fns, list_fn)
    else:
        list_fn = in_fns
        
    sam_cat_from_file(list_fn, out_fn, ncores = ncores)
    
    if is_input_list:
        os.remove(list_fn)

        
        
def sam_cat_from_file(list_fn, out_fn, ncores = 1):
    # pysam v0.23.0 has a bug that it does not recognize the "-@/--threads"
    # option.
    try:
        pysam.cat(
            "-b", list_fn,
            "-o", out_fn,
            "--no-PG", 
            "--threads", str(ncores-1)
        )
    except pysam.utils.SamtoolsError as e:
        error("catch pysam.utils.SamtoolsError:")
        error(str(e))
        info("try using samtools instead of pysam ...")
        proc = subprocess.Popen(
            args = "samtools cat -b %s -o %s --no-PG --threads %d" % \
                        (list_fn, out_fn, ncores - 1),
            shell = True,
            executable = "/bin/bash",
            stdout = subprocess.PIPE,
            stderr = subprocess.PIPE
        )
        outs, errs = proc.communicate()
        ret = proc.returncode
        if ret != 0:
            raise RuntimeError(str(errs.decode()))



def sam_fetch(sam, chrom, start = None, end = None):
    """Wrapper for pysam.fetch method that could automatically handle 
    chromosome names with or without the "chr" prefix.

    Parameters
    ----------
    sam : pysam.AlignmentFile
        The BAM file object.
    chrom : str
        Chromosome name.
    start : int or None, default None
        1-based, inclusive. `None` means minimum possible value.
    end : int or None, default None
        1-based, inclusive. `None` means maximum possible value.

    Returns
    -------
    Iterator
        Iterator of reads/alignments if success, `None` otherwise.
    """
    # sam.fetch(): start and stop denote 0-based, half-open intervals.
    if start is not None:
        start = start - 1

    try:   
        itr = sam.fetch(chrom, start, end)
    except:
        pass
    else:
        if itr:
            return itr
    chrom = chrom[3:] if chrom.startswith("chr") else "chr" + chrom
    try:
        itr = sam.fetch(chrom, start, end)
    except:
        return None
    else:
        return itr if itr else None
    

    
def sam_index(sam_fn_list, ncores = 1):
    """Index BAM file(s).
    
    Parameters
    ----------
    sam_fn_list : list of str
        A list of BAM files to be indexed.
    ncores : int, default 1
        Number of cores.
    
    Returns
    -------
    int
        Return code. 0 if success, negative otherwise.
    """
    if ncores == 1:
        for sam_fn in sam_fn_list:
            pysam.index(sam_fn)
    else:
        pool = multiprocessing.Pool(processes = ncores)
        mp_res = []
        for sam_fn in sam_fn_list:
            mp_res.append(pool.apply_async(
                func = pysam.index,
                args = (sam_fn, ),
                callback = None,
                error_callback = mp_error_handler
            ))
        pool.close()
        pool.join()
    return(0)

    

def sam_merge(in_fns, out_fn, ncores = 1):
    """Merge a list of BAM files.

    Parameters
    ----------
    in_fns : str or list of str
        If str, path to the file listing BAM files to be merged;
        If list of str, a list of BAM files to be merged.
    out_fn : str
        Path to the output merged BAM file.
    ncores : int
        Number of cores.
    
    Returns
    -------
    Void.
    """
    is_input_list = is_vector(in_fns)
    
    list_fn = None
    if is_input_list:
        list_fn = out_fn + ".lst"
        list2file(in_fns, list_fn)
    else:
        list_fn = in_fns
        
    sam_merge_from_file(list_fn, out_fn, ncores = ncores)
    
    if is_input_list:
        os.remove(list_fn)


        
def sam_merge_from_file(list_fn, out_fn, ncores = 1):
    """Merge BAM files listed in a file, one per line.
    
    Parameters
    ----------
    list_fn : str
        The file listing BAM files to be merged.
    out_fn : str
        Path to the output merged BAM file.
    ncores : int
        Number of cores.
    
    Returns
    -------
    Void.
    """
    pysam.merge(
        "-b", list_fn, 
        "-f", "-c", "-p", "--no-PG", 
        "-@", str(ncores-1), 
        out_fn
    )


    
def __sam_merge_batch(
    list_fn,
    n_sam,
    out_fn, 
    tmp_dir, 
    max_per_batch = 100, ncores = 1,
    depth = 0
):
    if n_sam <= max_per_batch:
        sam_merge(list_fn, out_fn, ncores = ncores)
        return
    
    # Note, here
    # - max_n_batch: to account for the max allowed files and subfolders in
    #   one folder.
    # - max_per_batch: to account for the max number of BAM files that
    #   pysam.merge() allows.
    bd_m, bd_n, bd_indices = split_n2batch(
            n_sam, ncores, max_n_batch = 10000, max_per_batch = max_per_batch)
    
    in_fn_list = file2list(list_fn)
    
    res_dir = os.path.join(tmp_dir, str(depth))
    os.makedirs(res_dir, exist_ok = True)
    
    bd_out_fn_list = []
    if ncores <= 1:
        for idx, (b, e) in enumerate(bd_indices):
            bd_list_fn = os.path.join(res_dir, "%d.%d.bam.lst" % (depth, idx))
            list2file(in_fn_list[b:e], bd_list_fn)
            bd_out_fn = os.path.join(res_dir, "%d.%d.bam" % (depth, idx))
            bd_out_fn_list.append(bd_out_fn)

            sam_merge(bd_list_fn, bd_out_fn, ncores = 1)
    else:
        mp_res = []
        pool = multiprocessing.Pool(processes = min(ncores, bd_m))
        for idx, (b, e) in enumerate(bd_indices):
            bd_list_fn = os.path.join(res_dir, "%d.%d.bam.lst" % (depth, idx))
            list2file(in_fn_list[b:e], bd_list_fn)
            bd_out_fn = os.path.join(res_dir, "%d.%d.bam" % (depth, idx))
            bd_out_fn_list.append(bd_out_fn)
            
            mp_res.append(pool.apply_async(
                func = sam_merge,
                kwds = dict(
                    list_fn = bd_list_fn,
                    out_fn = bd_out_fn,
                    ncores = 1
                ),
                callback = None,
                error_callback = mp_error_handler
            ))
        pool.close()
        pool.join()
        
    for fn in bd_out_fn_list:
        pysam.index(fn)
        
    new_list_fn = os.path.join(tmp_dir, "%d.bam.lst" % (depth + 1, ))
    list2file(bd_out_fn_list, new_list_fn)
    
    __sam_merge_batch(
        list_fn = new_list_fn,
        n_sam = len(bd_out_fn_list),
        out_fn = out_fn,
        tmp_dir = tmp_dir, 
        max_per_batch = max_per_batch,
        ncores = ncores,
        depth = depth + 1
    )
    

def sam_merge_batch(
    in_fn_list, out_fn, 
    tmp_dir, 
    max_per_batch = 100, ncores = 1,
    remove_tmp = False
):
    """Merge large number of BAM files by splitting them into batches.
    
    Parameters
    ----------
    in_fn_list : list of str
        A list of BAM files to be merged.
    out_fn : str
        Path to the output merged BAM file.
    tmp_dir : str
        Path to folder for storing temporary files.
    max_per_batch : int, default 100
        Maximum number of BAM files to be merged in one batch.
    ncores : int, default 1
        Number of cores.
    remove_tmp : bool, default False
        Whether to remove tmp files and `tmp_dir` after merging.
        
    Returns
    -------
    Void.    
    """
    depth = -1
    list_fn = os.path.join(tmp_dir, "%d.bam.lst" % (depth + 1, ))
    list2file(in_fn_list, list_fn)
            
    __sam_merge_batch(
        list_fn = list_fn,
        n_sam = len(in_fn_list),
        out_fn = out_fn,
        tmp_dir = tmp_dir, 
        max_per_batch = max_per_batch,
        ncores = ncores,
        depth = depth + 1
    )
    if remove_tmp:
        shutil.rmtree(tmp_dir)
    


def sam_sort(in_bam, out_bam = None, tag = None, max_mem = "4G", ncores = 1):
    """Sort BAM file by genomic position or specific tag.

    Parameters
    ----------
    in_bam : str
        The BAM file to be sorted.
    out_bam : str or None, default None
        The output BAM file.
        If `None`, sort in place.
    tag : str or None, default None
        The tag by which the BAM alignments will be sorted, e.g., "CB" to sort
        by cell barcodes.
        Set to None if do not use `tag` for sorting.
    max_mem : str
        The maximum memory to be used per thread.
        This value will be passed to the "-m" option of "samtools/pysam sort".
    ncores : int, default 1
        Number of cores.
    
    Returns
    -------
    int
        Return code. 0 if success, negative otherwise.
    """
    inplace = False
    if out_bam is None or out_bam == in_bam:
        inplace = True
        out_bam = in_bam + ".tmp.bam"

    args = ["-o", out_bam]
    if tag is not None:
        args.extend(["-t", tag])
    args.extend([
        "-m", str(max_mem),
        "-@", str(ncores - 1),
        "--no-PG",
        in_bam
    ])
    
    pysam.sort(*args)

    if inplace:
        os.replace(out_bam, in_bam)
    return(0)
