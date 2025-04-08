# sam.py


import pysam
from ..utils.sam import sam_cat, sam_sort


def sam_cat_and_sort(in_fns, out_fn, max_mem = "4G", ncores = 1, index = True):
    """Concatenate a list of BAM files and then sort by genomic position.

    Parameters
    ----------
    in_fns : str or list of str
        If str, path to the file listing BAM files to be concatenated;
        If list of str, a list of BAM files to be concatenated.
    out_fn : str
        Path to the output concatenated BAM file.
    max_mem : str
        The maximum memory to be used per thread.
        This value will be passed to the "-m" option of "samtools/pysam sort".
    ncores : int
        Number of cores.
    index : bool, default True
        Whether to index the `out_fn` after sorting.
    
    Returns
    -------
    Void.
    """
    sam_cat(
        in_fns = in_fns,
        out_fn = out_fn,
        ncores = ncores
    )
    
    sam_sort(
        in_bam = out_fn,
        out_bam = out_fn,
        tag = None,
        max_mem = max_mem,
        ncores = ncores
    )
    
    if index:
        pysam.index(out_fn)
        
        

def calc_max_mem(ncores, total_mem = 40, min_mem = 1, max_mem = 4):
    if max_mem * ncores <= total_mem:
        return(("%sG" % max_mem, ncores))
    
    if min_mem * ncores > total_mem:
        ncores = max(1, total_mem // min_mem)
        return(("%sG" % min_mem, ncores))
    
    mem = int(total_mem / ncores)
    return(("%sG" % mem, ncores)) 
