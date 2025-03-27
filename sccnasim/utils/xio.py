# xio.py


def file2list(fn):
    with open(fn, "r") as fp:
        lst = [line.strip().strip('"') for line in fp]
    return(lst)


def list2file(lst, fn):
    with open(fn, "w") as fp:
        for item in lst:
            fp.write("%s\n" % item)
