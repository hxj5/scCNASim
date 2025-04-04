# io.py


from ..afc.io import merge_tsv as afc_merge_tsv
from ..utils.zfile import ZF_F_PLAIN



# internal use only!
def merge_tsv(in_fn_list, out_fn, remove = False):
    """Merge TSV files.

    Parameters
    ----------
    in_fn_list : list of str
        Pathes to the TSV files outputed by each (sub-)process.
        The TSV files list the output features, each per line.
    out_fn : str
        Path to the output merged feature TSV file.
    remove : bool, default False
        Whether to remove the files in `in_fn_list` after merging.

    Returns
    -------
    int
        Return code. 0 if success, negative otherwise.
    """
    return afc_merge_tsv(
        in_fn_list = in_fn_list,
        in_format = ZF_F_PLAIN, 
        out_fn = out_fn,
        out_fmode = "w",
        out_format = ZF_F_PLAIN, 
        remove = remove
    )
