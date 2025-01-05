# io.py - input and output.

import functools
from logging import error
from ..io.base import load_cnas, load_features, save_features
from ..utils.grange import format_chrom, reg2str


def __cmp_two_intervals(x1, x2):
    """Compare two intervals.
    
    Parameters
    ----------
    x1 : list of (int, int)
        The begin and end coordinates of the first interval.
    x2 : list of (int, int)
        The begin and end coordinates of the second interval.

    Returns
    -------
    int
        The return code.
        0 if equal; positive if `x1` is larger; negative if `x2` is larger.
    """
    s1, e1 = x1[:2]
    s2, e2 = x2[:2]
    if s1 == s2:
        if e1 == e2:
            return(0)
        else:
            return e1 - e2
    else:
        return s1 - s2
    

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
        rec = df.loc[i, ]
        chrom = rec["chrom"]
        clone_id = rec["clone"]
        if clone_id not in dat:
            dat[clone_id] = {}
        if chrom not in dat[clone_id]:
            dat[clone_id][chrom] = {}
        ale_key = "%d_%d" % (rec["cn_ale0"], rec["cn_ale1"])
        if ale_key not in dat[clone_id][chrom]:
            dat[clone_id][chrom][ale_key] = []
        dat[clone_id][chrom][ale_key].append((rec["start"], rec["end"]))


    # merge (clone-specific) adjacent CNAs.
    for clone_id, cl_dat in dat.items():
        for chrom, ch_dat in cl_dat.items():
            for ale_key in ch_dat.keys():
                iv_list = sorted(
                    ch_dat[ale_key], 
                    key = functools.cmp_to_key(__cmp_two_intervals)
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
                ch_dat[ale_key] = new_list


    # check whether there are (strictly) overlapping regions with 
    # distinct profiles.
    for clone_id, cl_dat in dat.items():
        for chrom, ch_dat in cl_dat.items():
            iv_list = []
            for ale_key in ch_dat.keys():
                cn_ale0, cn_ale1 = [int(x) for x in ale_key.split("_")]
                iv_list.extend(
                    [(s, e, cn_ale0, cn_ale1) for s, e in ch_dat[ale_key]])
            iv_list = sorted(
                iv_list, 
                key = functools.cmp_to_key(__cmp_two_intervals)
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
    for clone_id in sorted(dat.keys()):
        cl_dat = dat[clone_id]
        for chrom in sorted(cl_dat.keys()):
            ch_dat = cl_dat[chrom]
            for s, e, cn_ale0, cn_ale1 in ch_dat:
                region_id = reg2str(chrom, s, e)
                fp.write("\t".join([chrom, str(s), str(e), region_id, \
                    clone_id, str(cn_ale0), str(cn_ale1)]) + "\n")
                n_new += 1
    fp.close()
    return((0, n_old, n_new))


def filter_features_by_chroms(in_fn, out_fn, chrom_list):
    """Filter features by chromosomes.
    
    This function filters features that are not in the input chromosomes.

    Parameters
    ----------
    in_fn : str
        Path to the input file.
    out_fn : str
        Path to the output file.
    chrom_list: list of str
        A list of chromosome names.
    
    Returns
    -------
    int
        The return code. 0 if success, negative if error.
    int
        Number of records before filtering.
    int
        Number of records after filtering.
    """
    sep = "\t"
    
    df = load_features(in_fn, sep = sep)
    n_old = df.shape[0]
    
    chrom_list = [format_chrom(c) for c in chrom_list]
    df_new = df[df["chrom"].isin(chrom_list)]
    n_new = df_new.shape[0]
    
    save_features(df_new, out_fn, sep = sep)
    
    return((0, n_old, n_new))


def merge_features_bidel(in_fn, out_fn, max_gap = 1, max_frac = 0.1):
    """Merge adjacent features and remove overlapping bi-features.
    
    Here, bi-features means two features overlapping with each other.

    Parameters
    ----------
    in_fn : str
        Path to the input file.
    out_fn : str
        Path to the output file.
    max_gap : int, default 1
        The maximum gap length that is allowed between two adjacent regions.
        `1` for strict adjacence.
    max_frac : float, default 0.1
        The maximum fraction of overlapping part of either one of bi-features.
        If exceeds the value, both features will be removed.
    
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
    max_fn_len = 127
    n_old, n_new = -1, -1


    # load data
    dat = {}
    df = load_features(in_fn, sep = sep)
    n_old = df.shape[0]
    for i in range(df.shape[0]):
        rec = df.loc[i, ]
        chrom = rec["chrom"]
        if chrom not in dat:
            dat[chrom] = []
        dat[chrom].append((rec["start"], rec["end"], rec["feature"]))


    # merge adjacent features
    bi_dat = set()      # features to be removed.
    for chrom, ch_dat in dat.items():
        iv_list = sorted(ch_dat, 
                    key = functools.cmp_to_key(__cmp_two_intervals))
        n = len(iv_list)
        for i in range(n - 1):
            s1, e1, f1 = iv_list[i]
            for j in range(i + 1, n):
                s2, e2, f2 = iv_list[j]
                if s2 <= e1 + max_gap:    # overlap adjacent region
                    len1 = e1 - s1 + 1
                    len2 = e2 - s2 + 1
                    len_overlap = min(e1, e2) - s2 + 1
                    if len_overlap / float(min(len1, len2)) > max_frac:
                        bi_dat.add(f1)
                        bi_dat.add(f2)


    # save features
    n_new = 0
    fp = open(out_fn, "w")
    for chrom in sorted(dat.keys()):
        ch_dat = dat[chrom]
        for s, e, f in ch_dat:
            if f not in bi_dat:
                fp.write("\t".join([chrom, str(s), str(e), f]) + "\n")
                n_new += 1
    fp.close()
    return((0, n_old, n_new))


def merge_features_first(in_fn, out_fn, max_gap = 1, new_name_how = "join"):
    """Merge adjacent features and only keep the first feature.
    
    Keep the first feature and remove features overlapping with it.

    Parameters
    ----------
    in_fn : str
        Path to the input file.
    out_fn : str
        Path to the output file.
    max_gap : int, default 1
        The maximum gap length that is allowed between two adjacent regions.
        `1` for strict adjacence.
    new_name_how : str, default "join"
        How to name the merged features.
        "join": join the names of the two features with string "__".
    
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
    max_fn_len = 127
    n_old, n_new = -1, -1


    # load data
    dat = {}
    df = load_features(in_fn, sep = sep)
    n_old = df.shape[0]
    for i in range(df.shape[0]):
        rec = df.loc[i, ]
        chrom = rec["chrom"]
        if chrom not in dat:
            dat[chrom] = []
        dat[chrom].append((rec["start"], rec["end"], rec["feature"]))


    # merge adjacent features
    for chrom, ch_dat in dat.items():
        iv_list = sorted(ch_dat, 
                    key = functools.cmp_to_key(__cmp_two_intervals))
        s1, e1, f1 = iv_list[0]
        new_list = []
        for s2, e2, f2 in iv_list[1:]:
            if s2 <= e1 + max_gap:    # overlap adjacent region
                if new_name_how == "join":
                    if len(f1) < max_fn_len:
                        f1 = f1 + "__" + f2
                        if len(f1) >= max_fn_len:
                            f1 += "__"      # as a marker of truncated string.
            else:                     # otherwise
                new_list.append((s1, e1, f1))
                s1, e1, f1 = s2, e2, f2
        new_list.append((s1, e1, f1))
        dat[chrom] = new_list


    # save features
    n_new = 0
    fp = open(out_fn, "w")
    for chrom in sorted(dat.keys()):
        ch_dat = dat[chrom]
        for s, e, f in ch_dat:
            fp.write("\t".join([chrom, str(s), str(e), f]) + "\n")
            n_new += 1
    fp.close()
    return((0, n_old, n_new))


def merge_features_union(in_fn, out_fn, max_gap = 1, new_name_how = "join"):
    """Merge adjacent features and keep their union genomic range.

    Parameters
    ----------
    in_fn : str
        Path to the input file.
    out_fn : str
        Path to the output file.
    max_gap : int, default 1
        The maximum gap length that is allowed between two adjacent regions.
        `1` for strict adjacence.
    new_name_how : str, default "join"
        How to name the merged features.
        "join": join the names of the two features with string "__".
    
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
    max_fn_len = 127
    n_old, n_new = -1, -1


    # load data
    dat = {}
    df = load_features(in_fn, sep = sep)
    n_old = df.shape[0]
    for i in range(df.shape[0]):
        rec = df.loc[i, ]
        chrom = rec["chrom"]
        if chrom not in dat:
            dat[chrom] = []
        dat[chrom].append((rec["start"], rec["end"], rec["feature"]))


    # merge adjacent features
    for chrom, ch_dat in dat.items():
        iv_list = sorted(ch_dat, 
                    key = functools.cmp_to_key(__cmp_two_intervals))
        s1, e1, f1 = iv_list[0]
        new_list = []
        for s2, e2, f2 in iv_list[1:]:
            if s2 <= e1 + max_gap:    # overlap adjacent region
                e1 = max(e1, e2)
                if new_name_how == "join":
                    if len(f1) < max_fn_len:
                        f1 = f1 + "__" + f2
                        if len(f1) >= max_fn_len:
                            f1 += "__"      # as a marker of truncated string.
            else:                     # otherwise
                new_list.append((s1, e1, f1))
                s1, e1, f1 = s2, e2, f2
        new_list.append((s1, e1, f1))
        dat[chrom] = new_list


    # save features
    n_new = 0
    fp = open(out_fn, "w")
    for chrom in sorted(dat.keys()):
        ch_dat = dat[chrom]
        for s, e, f in ch_dat:
            fp.write("\t".join([chrom, str(s), str(e), f]) + "\n")
            n_new += 1
    fp.close()
    return((0, n_old, n_new))
