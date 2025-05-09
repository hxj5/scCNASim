# xfile.py - wrapper of file object, supporting GZIP/BGZIP. 

# Note that there's a bug when reading BGZIP file with pysam (issue #438):
#   it will stop reading when encounting a blank line, this could be caused 
#   by the issue that `bgzf_getline()` function always drops the tailing '\n'.
#   This bug exists in pysam v0.18.0 and before(?).
# So we need to use gzip to read GZIP and BGZIP files
# it's probably ok to use pysam to write BGZIP files.


import gzip
import os
import pysam



class ZFile:
    """Simple wrapper of file object that supports plain/GZIP/BGZF formats."""

    def __init__(self, file_name, mode, file_type, 
                 is_bytes = False, encoding = None):
        """
        Parameters
        ----------
        file_name : str
            Path to the file.
        mode : str
            File mode.
        file_type : int
            File type / format, one of `ZF_F_XXX`.
        is_bytes : bool, default False
            Whether the data is of type `bytes`.
        encoding : str or None, default None
            Encoding for `bytes` data.
            If None, set to "utf8".
        """
        self.file_name = file_name
        self.file_type = file_type
        self.mode = mode
        self.is_bytes = is_bytes
        self.encoding = encoding if encoding else "utf8"

        # fp
        #   The file object.
        self.fp = None

        # buf : str or bytes
        #   The buffer.
        self.buf = self.__reset_buf()

        if file_type == ZF_F_AUTO:
            fn = file_name.lower()
            if fn.endswith(".gz") or fn.endswith(".gzip"):
                file_type = ZF_F_GZIP
            else:
                file_type = ZF_F_PLAIN

        if file_type == ZF_F_PLAIN:
            self.fp = open(self.file_name, self.mode)
        elif file_type == ZF_F_GZIP:
            self.fp = gzip.open(self.file_name, self.mode)
        elif file_type == ZF_F_BGZIP:
            if "r" in self.mode:
                self.fp = gzip.open(self.file_name, self.mode)
            else:
                self.fp = pysam.BGZFile(self.file_name, self.mode)
        else:
            raise ValueError("invalid file type")

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

    def __iter__(self):
        return self.fp

    def __next__(self):
        line = self.readline()
        if not line:
            raise StopIteration()
        return(line)

    def __reset_buf(self):
        return(bytes("", self.encoding) if self.is_bytes else "")

    def close(self):
        if self.fp:
            if self.buf:
                self.fp.write(self.buf)
                self.buf = None
            self.fp.close()
            self.fp = None

    def read(self, size = None):
        if not self.fp:
            raise OSError()
        if size is None:
            return(self.fp.read())
        else:
            return(self.fp.read(size))

    def readline(self, size = None):
        if not self.fp:
            raise OSError()
        if size is None:
            return(self.fp.readline())
        else:
            return(self.fp.readline(size))

    def readlines(self, size = None):
        if not self.fp:
            raise OSError()
        if size is None:
            return(self.fp.readlines())
        else:
            return(self.fp.readlines(size))

    def write(self, data):
        if not self.fp:
            raise OSError()
        self.buf += data
        if len(self.buf) >= ZF_BUFSIZE:
            ret = self.fp.write(self.buf)
            self.buf = self.__reset_buf()
            return(ret)
        return(len(data))

    

def zopen(file_name, mode, file_type = None, is_bytes = False, encoding = None):
    """Open a file.
    
    Parameters
    ----------
    file_name : str
        Path to the file.
    mode : str
        File mode.
    file_type : int or None, default None
        File type / format, one of ZF_F_XXX.
        If None, set to ZF_F_AUTO.
    is_bytes : bool, default False
        Whether the data is of type `bytes`.
    encoding : str or None, default None
        Encoding for `bytes` data.
        If None, set to "utf8".
    
    Returns
    -------
    utils.zfile.ZFile
        The file object.
    """
    if not file_name:
        raise OSError()
    if file_type is None:
        file_type = ZF_F_AUTO
    return ZFile(file_name, mode, file_type, is_bytes, encoding)



def zconcat(in_fn_list, in_format, 
              out_fn, out_fmode, out_format, 
              remove = False):
    """Concatenate a list of files.

    Parameters
    ----------
    in_fn_list : list of str
        Pathes to the files to be concatenated.
    in_format : int
        The format of each file in `in_fn_list`.
        It should be compatible with the `file_type` in :func:`zopen`.
    out_fn : str
        Path to the output file.
    out_fmode : str
        The file mode of the `out_fn`.
        It should be compatible with the `mode` in :func:`zopen`.
    out_format : int
        The file format of `out_fn`.
        It should be compatible with the `file_type` in :func:`zopen`.
    remove : bool, default False
        Whether to remove the files in `in_fn_list` after concatenation.

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



# file type / format.
ZF_F_PLAIN = 0
ZF_F_GZIP = 1
ZF_F_BGZIP = 2
ZF_F_AUTO = 3

ZF_BUFSIZE = 1048576   # 1M



# TODO: update the debugging codes below
if __name__ == "__main__":
    pass
