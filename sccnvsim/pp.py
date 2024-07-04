# pp.py - preprocessing

import functools

from logging import error

from .utils.grange import format_chrom, format_start, format_end
from .utils.zfile import zopen


def __cmp_two_intervals(x1, x2):
    s1, e1 = x1[:2]
    s2, e2 = x2[:2]
    if s1 == s2:
        if e1 == e2:
            return(0)
        else:
            return e1 - e2
    else:
        return s1 - s2


def merge_features(in_fn, out_fn, max_gap = 1, new_name_how = "join"):
    """Merge adjacent features.

    Parameters
    ----------
    in_fn : str
        Path to the input file.
    out_fn : str
        Path to the output file.
    max_gap : int
        The maximum gap length that is allowed between two adjacent regions. 
        `1` for strict adjacence.
    new_name_how : str
        How to name the merged features. `join`: join the names of the two 
        features with char '>'.
    
    Returns
    -------
    int
        0 if success, negative if error.
    """
    sep = "\t"

    # load data
    fp = zopen(in_fn, "rt")
    dat = {}
    nl = 0
    for line in fp:
        nl += 1
        items = line.strip().split(sep)
        if len(items) < 4:
            error("too few columns of line %d.\n" % nl)
            return(-3)
        chrom, start, end, feature = items[:4]
        start, end = format_start(start), format_end(end)
        feature = feature.strip('"')
        if chrom not in dat:
            dat[chrom] = []
        dat[chrom].append((start, end, feature))
    fp.close()

    # merge adjacent features
    for chrom, ch_dat in dat.items():
        iv_list = sorted(ch_dat, key = functools.cmp_to_key(__cmp_two_intervals))
        s1, e1, f1 = iv_list[0]
        new_list = []
        for s2, e2, f2 in iv_list[1:]:
            if s2 <= e1 + max_gap:    # overlap adjacent region
                e1 = max(e1, e2)
                if new_name_how == "join":
                    f1 = f1 + ">" + f2
            else:                     # otherwise
                new_list.append((s1, e1, f1))
                s1, e1, f1 = s2, e2, f2
        new_list.append((s1, e1, f1))
        dat[chrom] = new_list

    # save features
    fp = open(out_fn, "w")
    for chrom in sorted(dat.keys()):
        ch_dat = dat[chrom]
        for s, e, f in ch_dat:
            fp.write("\t".join([chrom, str(s), str(e), f]) + "\n")
    fp.close()
    return(0)
