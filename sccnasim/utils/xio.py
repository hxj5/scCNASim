# xio.py


import pickle



def file2list(fn):
    with open(fn, "r") as fp:
        lst = [line.strip().strip('"') for line in fp]
    return(lst)


def list2file(lst, fn):
    with open(fn, "w") as fp:
        for item in lst:
            fp.write("%s\n" % item)

            
            
def load_pickle(fn):
    with open(fn, "rb") as fp:
        obj = pickle.load(fp)
    return(obj)


def save_pickle(obj, fn):
    with open(fn, "wb") as fp:
        pickle.dump(obj, fp)
