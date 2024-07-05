# pp.py - preprocessing

import functools
import os
import time

from logging import info, error

from ..utils.grange import format_chrom, format_start, format_end, reg2str
from ..utils.zfile import zopen


COMMAND 



# conf: config::GlobalConfig
def pp_core(conf):
    os.makedirs(conf.out_dir, exist_ok = True)

    merged_feature_fn = os.path.join(conf.out_dir, "features.tsv")
    r = merge_features(
        in_fn = conf.g.feature_fn,
        out_fn = merged_feature_fn,
        max_gap = 1,
        new_name_how = "join"
    )
    if r < 0:
        error("merge features failed (%d)." % r)
        raise ValueError
    
    merged_feature_fn = os.path.join(conf.out_dir, "features.merged.tsv")
    r = merge_features(
        in_fn = conf.g.feature_fn,
        out_fn = merged_feature_fn,
        max_gap = 1,
        new_name_how = "join"
    )
    if r < 0:
        error("merge features failed (%d)." % r)
        raise ValueError    
    
    res = {
        "merged_feature_fn": merged_feature_fn
    }
    
    return(res)


def pp_run(conf):
    ret = -1
    res = None

    start_time = time.time()
    time_str = time.strftime(
        "%Y-%m-%d %H:%M:%S", time.localtime(start_time))
    info("start time: %s." % time_str)

    try:
        res = pp_core(conf)
    except ValueError as e:
        error(str(e))
        error("Running program failed.")
        error("Quiting ...")
        ret = -1
    else:
        info("All Done!")
        ret = 0
    finally:
        end_time = time.time()
        time_str = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime(end_time))
        info("end time: %s" % time_str)
        info("time spent: %.2fs" % (end_time - start_time, ))

    return((ret, res))


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
    

def merge_cnv_profile(in_fn, out_fn, max_gap = 1):
    """Merge adjacent regions with the same CNV profiles

    Merge adjacent regions with the same allele-specific copy number
    profile in each CNV clone.

    Parameters
    ----------
    in_fn : str
        Path to input file.
    out_fn : str
        Path to output file.
    max_gap : int
        The maximum gap length that is allowed between two adjacent regions.
        `1` for strict adjacence.

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
        if len(items) < 7:
            error("too few columns of line %d." % nl)
            return(-3)
        chrom, start, end, region_id, clone_id, cn_ale0, cn_ale1 = items[:7]
        chrom = format_chrom(chrom)
        start, end = format_start(start), format_end(end)
        cn_ale0, cn_ale1 = int(cn_ale0), int(cn_ale1)
        region_id = region_id.strip('"')
        clone_id = clone_id.strip('"')
        if clone_id not in dat:
            dat[clone_id] = {}
        if chrom not in dat[clone_id]:
            dat[clone_id][chrom] = {}
        ale_key = "%d_%d" % (cn_ale0, cn_ale1)
        if ale_key not in dat[clone_id][chrom]:
            dat[clone_id][chrom][ale_key] = []
        dat[clone_id][chrom][ale_key].append((start, end))
    fp.close()

    # merge (clone-specific) adjacent CNVs.
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
                    error("distinct CNV profiles '%s', (%d, %d) and (%d, %d)." % 
                        (chrom, s1, e1, s2, e2))
                    return(-5)
            cl_dat[chrom] = iv_list

    # save profile
    fp = open(out_fn, "w")
    for clone_id in sorted(dat.keys()):
        cl_dat = dat[clone_id]
        for chrom in sorted(cl_dat.keys()):
            ch_dat = cl_dat[chrom]
            for s, e, cn_ale0, cn_ale1 in ch_dat:
                region_id = reg2str(chrom, s, e)
                fp.write("\t".join([chrom, str(s), str(e), region_id, \
                    clone_id, str(cn_ale0), str(cn_ale1)]) + "\n")
    fp.close()


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
            error("too few columns of line %d." % nl)
            return(-3)
        chrom, start, end, feature = items[:4]
        chrom = format_chrom(chrom)
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
