# grange.py - genomic range



def cmp_two_intervals(x1, x2):
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
