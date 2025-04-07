# gcna.py - genomic CNA pre-processing.


import functools
from logging import error
from .utils import cmp_two_intervals
from ..io.base import load_cnas



def merge_cna_profile(in_fn, out_fn, max_gap = 1):
    """Merge adjacent regions with the same CNA profiles.

    Merge adjacent regions with the same allele-specific copy number
    profile in each CNA clone.

    Parameters
    ----------
    in_fn : str
        Path to input file.
    out_fn : str
        Path to output file.
    max_gap : int, default 1
        The maximum gap length that is allowed between two adjacent regions.
        `1` for strict adjacence.

    Returns
    -------
    int
        The return code. 0 if success, negative if error.
    int
        Number of records before merging.
    int
        Number of records after merging.
    """
    sep = "\t"
    n_old, n_new = -1, -1

    # load data
    try:
        df = load_cnas(in_fn, sep = sep)
    except Exception as e:
        error("load CNA profile file failed '%s'." % str(e))
        return((-3, n_old, n_new))
    n_old = df.shape[0]

    dat = {}
    for i in range(df.shape[0]):
        rec = df.iloc[i, ]
        chrom = rec["chrom"]
        clone = rec["clone"]
        if clone not in dat:
            dat[clone] = {}
        if chrom not in dat[clone]:
            dat[clone][chrom] = {}
        cns = "%d_%d" % (rec["cn_ale0"], rec["cn_ale1"])
        if cns not in dat[clone][chrom]:
            dat[clone][chrom][cns] = []
        dat[clone][chrom][cns].append((rec["start"], rec["end"]))


    # merge (clone-specific) adjacent CNAs.
    for clone, cl_dat in dat.items():
        for chrom, ch_dat in cl_dat.items():
            for cns in ch_dat.keys():
                iv_list = sorted(
                    ch_dat[cns], 
                    key = functools.cmp_to_key(cmp_two_intervals)
                )
                s1, e1 = iv_list[0]
                new_list = []
                for s2, e2 in iv_list[1:]:
                    if s2 <= e1 + max_gap:    # overlap adjacent region
                        e1 = max(e1, e2)
                    else:                     # otherwise
                        new_list.append((s1, e1))
                        s1, e1 = s2, e2
                new_list.append((s1, e1))
                ch_dat[cns] = new_list


    # check whether there are (strictly) overlapping regions with 
    # distinct profiles.
    for clone, cl_dat in dat.items():
        for chrom, ch_dat in cl_dat.items():
            iv_list = []
            for cns in ch_dat.keys():
                cn_ale0, cn_ale1 = [int(x) for x in cns.split("_")]
                iv_list.extend(
                    [(s, e, cn_ale0, cn_ale1) for s, e in ch_dat[cns]])
            iv_list = sorted(
                iv_list, 
                key = functools.cmp_to_key(cmp_two_intervals)
            )
            s1, e1 = iv_list[0][:2]
            for iv in iv_list[1:]:
                s2, e2 = iv[:2]
                if s2 <= e1:    # overlap adjacent region
                    error("distinct CNA profiles '%s', (%d, %d) and (%d, %d)." % 
                        (chrom, s1, e1, s2, e2))
                    return((-5, n_old, n_new))
            cl_dat[chrom] = iv_list


    # save profile
    n_new = 0
    fp = open(out_fn, "w")
    for clone in sorted(dat.keys()):
        cl_dat = dat[clone]
        for chrom in sorted(cl_dat.keys()):
            ch_dat = cl_dat[chrom]
            for s, e, cn_ale0, cn_ale1 in ch_dat:
                fp.write("\t".join([chrom, str(s), str(e), \
                    clone, str(cn_ale0), str(cn_ale1)]) + "\n")
                n_new += 1
    fp.close()
    return((0, n_old, n_new))
