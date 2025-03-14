# io.py - input and output.

import os

from logging import error
from logging import warning as warn
from .gfeature import SNP, SNPSet, BlockRegion
from ..io.base import load_features, load_snps
from ..utils.vcf import vcf_load
from ..utils.zfile import zopen


def load_feature_from_txt(fn, sep = "\t"):
    """Load features from plain file.

    Parameters
    ----------
    fn : str
        Path to header-free plain file listing features, each per line.
        The first 5 columns should be
        chrom : str
            The chromosome name of the feature.
        start : int
            The start genomic position of the feature, 1-based and inclusive.
        end : int
            The end genomic position of the feature, 1-based and inclusive.
        name : str
            The name of the feature.
        strand : str
            DNA strand orientation of the feature, "+" (positive) or 
            "-" (negative).
    sep : str, default '\t'
        The delimiter of the `fn`.

    Returns
    -------
    list of afc.gfeature.BlockRegion or None
        A list of :class:`~afc.gfeature.BlockRegion` objects if success,
        `None` otherwise.
    """
    reg_list = []
    df = load_features(fn, sep = sep)
    for i in range(df.shape[0]):
        rec = df.loc[i, ]
        reg = BlockRegion(
            rec["chrom"],
            rec["start"],
            rec["end"] + 1,
            rec["feature"],
            rec["strand"]
        )
        reg_list.append(reg)
    return reg_list


def load_snp_from_tsv(fn, sep = "\t"):
    """Load phased SNPs from TSV file.

    Parameters
    ----------
    fn : str
        Path to TSV file containing 6 columns without header:
        chrom : str
            Chromosome name.
        pos : int
            1-based genomic position.
        ref : str
            The reference (REF) base of the SNP.
        alt : str
            The alternative (ALT) base of the SNP.
        ref_hap : {0, 1}
            The haplotype index for the `ref` base.
        alt_hap : {1, 0}
            The haplotype index for the `alt` base.
    sep : str, default '\t'
        The delimiter of the `fn`.

    Returns
    -------
    afc.gfeature.SNPSet or None
        A :class:`~afc.gfeature.SNPSet` object if success, `None` otherwise.
    """
    snp_set = SNPSet()
    df = load_snps(fn, sep = sep)
    for i in range(df.shape[0]):
        nl = i + 1
        rec = df.loc[i, ]
        ref, alt = rec["ref"].upper(), rec["alt"].upper()
        if len(ref) != 1 or ref not in "ACGTN":
            warn("invalid REF base of line %d." % nl)
            continue
        if len(alt) != 1 or alt not in "ACGTN":
            warn("invalid ALT base of line %d." % nl)
            continue
        a1, a2 = rec["ref_hap"], rec["alt_hap"]
        if (a1 == 0 and a2 == 1) or (a1 == 1 and a2 == 0):
            snp = SNP(
                chrom = rec["chrom"], 
                pos = rec["pos"],
                ref = ref, 
                alt = alt, 
                ref_hap = a1, 
                alt_hap = a2
            )
            if snp_set.add(snp) < 0:
                error("failed to add SNP of line %d." % nl)
                return None
        else:
            warn("invalid GT of line %d." % nl)
            continue
    return snp_set


def load_snp_from_vcf(fn):
    """Load phased SNPs from VCF file.

    Parameters
    ----------
    fn : str
        Path to VCF file.

    Returns
    -------
    afc.gfeature.SNPSet or None
        A :class:`~afc.gfeature.SNPSet` object if success, `None` otherwise.
    """
    snp_set = SNPSet()
    df, header = vcf_load(fn)
    for i in range(df.shape[0]):
        nl = i + 1
        rec = df.loc[i, ]
        ref, alt = rec["REF"].upper(), rec["ALT"].upper()
        if len(ref) != 1 or ref not in "ACGTN":
            warn("invalid REF base of line %d." % nl)
            continue
        if len(alt) != 1 or alt not in "ACGTN":
            warn("invalid ALT base of line %d." % nl)
            continue          
        fields = rec["FORMAT"].split(":")
        if "GT" not in fields:
            warn("GT not in line %d." % nl)
            continue
        idx = fields.index("GT")
        values = rec.iloc[9].split(":")
        if len(values) != len(fields):
            warn("len(fields) != len(values) in line %d." % nl)
            continue
        gt = values[idx]
        sep = ""
        if "|" in gt: 
            sep = "|"
        elif "/" in gt: 
            sep = "/"
        else:
            warn("invalid delimiter of line %d." % nl)
            continue
        a1, a2 = [int(x) for x in gt.split(sep)[:2]]
        if (a1 == 0 and a2 == 1) or (a1 == 1 and a2 == 0):
            snp = SNP(
                chrom = rec["CHROM"],
                pos = rec["POS"], 
                ref = ref, 
                alt = alt, 
                ref_hap = a1, 
                alt_hap = a2
            )
            if snp_set.add(snp) < 0:
                error("failed to add SNP of line %d." % nl)
                return None
        else:
            warn("invalid GT of line %d." % nl)
            continue
    return snp_set


def _fmt_line(ln, k):
    items = ln.split("\t")
    items[0] = str(int(items[0]) + k)
    return("\t".join(items))


# internal use only!
def merge_mtx(in_fn_list, in_format,
              out_fn, out_fmode, out_format,
              nrow_list, ncol, nrecord, remove = False):
    """Merge count matrix files generated by multi-processing.

    In each (sub-)process, the feature index in the output count matrix file
    always starts from 0 (due to feature filtering, the transcriptomics-scale
    index should not be assigned beforehand).
    This function merges the count matrices outputed by each (sub-)processes,
    and update the feature index to 0-based, transcriptomics-scale.

    Parameters
    ----------
    in_fn_list : list of str
        Pathes to the count matrix files outputed by each (sub-)process.
    in_format : int
        The format of each file in `in_fn_list`.
        Its value should be compatible with the "file_type" option in
        :func:`~utils.zfile.zopen`.
    out_fn : str
        Path to the output merged count matrix file.
    out_fmode : str
        The file mode of the `out_fn`.
        Its value should be compatible with the "mode" option in
        :func:`~utils.zfile.zopen`.
    out_format : int
        The file format of `out_fn`.
        Its value should be compatible with the "file_type" option in
        :func:`~utils.zfile.zopen`.
    nrow_list : list of int
        Number of unique features (row indexes) in each file within
        `in_fn_list`.
    ncol : int
        Number of unique cells (column indexes).
        Every file in `in_fn_list` should have the same `ncol`.
    nrecord : int
        Total number of records in all count matrices.
    remove : bool, default False
        Whether to remove the files in `in_fn_list` after merging.

    Returns
    -------
    int
        Return code. 0 if success, negative otherwise.
    """
    if len(in_fn_list) != len(nrow_list):
        return(-1)
    bufsize = 1048576    # 1M

    is_bytes = "b" in out_fmode
    out_fp = zopen(out_fn, out_fmode, out_format, is_bytes = is_bytes)

    nrow_total = sum(nrow_list)
    header  = "%%MatrixMarket matrix coordinate integer general\n"
    header += "%%\n"
    header += "%d\t%d\t%d\n" % (nrow_total, ncol, nrecord)
    if is_bytes:
        header = bytes(header, "utf8")
    out_fp.write(header)

    nline = 0
    k = 0
    for in_fn, nrow in zip(in_fn_list, nrow_list):
        with zopen(in_fn, "rt", in_format) as in_fp:
            while True:
                lines = in_fp.readlines(bufsize)
                if not lines:
                    break
                nline += len(lines)
                lines = [_fmt_line(ln, k) for ln in lines]
                s = "".join(lines)
                if is_bytes:
                    s = bytes(s, "utf8")
                out_fp.write(s)
        k += nrow
    out_fp.close()
    if nline != nrecord:
        return(-1)
    if remove:
        for in_fn in in_fn_list:
            os.remove(in_fn)
    return(0) 


# internal use only!
def merge_tsv(in_fn_list, in_format, 
              out_fn, out_fmode, out_format, 
              remove = False):
    """Merge feature TSV files generated by multi-processing.

    Parameters
    ----------
    in_fn_list : list of str
        Pathes to the TSV files outputed by each (sub-)process.
        The TSV files list the output features, each per line.
    in_format : int
        The format of each file in `in_fn_list`.
        Its value should be compatible with the "file_type" option in
        :func:`~utils.zfile.zopen`.
    out_fn : str
        Path to the output merged feature TSV file.
    out_fmode : str
        The file mode of the `out_fn`.
        Its value should be compatible with the "mode" option in
        :func:`~utils.zfile.zopen`.
    out_format : int
        The file format of `out_fn`.
        Its value should be compatible with the "file_type" option in
        :func:`~utils.zfile.zopen`.
    remove : bool, default False
        Whether to remove the files in `in_fn_list` after merging.

    Returns
    -------
    int
        Return code. 0 if success, negative otherwise.
    """
    bufsize = 1048576   # 1M
    is_bytes = "b" in out_fmode
    out_fp = zopen(out_fn, out_fmode, out_format, is_bytes)
    in_fmode = "rb" if is_bytes else "rt"
    for in_fn in in_fn_list:
        with zopen(in_fn, in_fmode, in_format) as in_fp:
            while True:
                dat = in_fp.read(bufsize)
                if not dat:
                    break
                out_fp.write(dat)
    out_fp.close()
    if remove:
        for in_fn in in_fn_list:
            os.remove(in_fn)
    return(0)


# internal use only!
def rewrite_mtx(in_fn, in_format, 
                out_fn, out_fmode, out_format, 
                nrow, ncol, nrecord, remove = False):
    """Add file header to a temporary count matrix file.

    The temporary count matrix file does not have the file header, i.e., the
    file header of the matrix market exchange format.

    Parameters
    ----------
    in_fn : str
        Path to the temporary count matrix file.
    in_format : int
        The file format of `in_fn`.
        Its value should be compatible with the "file_type" option in
        :func:`~utils.zfile.zopen`.
    out_fn : str
        Path to the updated output count matrix file.
    out_fmode : str
        The file mode of the `out_fn`.
        Its value should be compatible with the "mode" option in
        :func:`~utils.zfile.zopen`.
    out_format : int
        The file format of `out_fn`.
        Its value should be compatible with the "file_type" option in
        :func:`~utils.zfile.zopen`.
    nrow : int
        Number of unique features (row indexes) in `in_fn`.
    ncol : int
        Number of unique cells (column indexes) in `in_fn`.
    nrecord : int
        Number of records in `in_fn`.
    remove : bool, default False
        Whether to remove the original `in_fn` after adding file header.

    Returns
    -------
    int
        Return code. 0 if success, negative otherwise.
    """
    if in_fn == out_fn:
        return(-1)
    bufsize = 1048576   # 1M
  
    is_bytes = "b" in out_fmode
    out_fp = zopen(out_fn, out_fmode, out_format, is_bytes = is_bytes)
    header  = "%%MatrixMarket matrix coordinate integer general\n"
    header += "%%\n"
    header += "%d\t%d\t%d\n" % (nrow, ncol, nrecord)
    if is_bytes:
        header = bytes(header, "utf8")
    out_fp.write(header)

    nline = 0
    in_fmode = "rb" if is_bytes else "rt"
    sep = b"" if is_bytes else ""
    with zopen(in_fn, in_fmode, in_format) as in_fp:
        while True:
            lines = in_fp.readlines(bufsize)
            if not lines:
                break
            nline += len(lines)
            s = sep.join(lines)
            out_fp.write(s)
    out_fp.close()
    if nline != nrecord:
        return(-1)
    if remove:
        os.remove(in_fn)
    return(0)
