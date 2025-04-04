# gcna.py



import numpy as np



def calc_cn_fold(cna_profile, cna_features, p, normal = 2, loss_allele_freq = 0.001):
    """Calculate the feature-specific copy number (CN) fold in each clone.
    
    Parameters
    ----------
    pandas.DataFrame
        The clonal CNA profile.
        It should contain a column "cn".
    cna_features : dict of {str : numpy.ndarray}
        The indices of overlapping features of each CNA region.
    p : int
        Number of all features.
    normal : int
        The copy number in normal cells.
    loss_allele_freq : float, default 0.001
        The frequency of the lost allele, to mimic real error rate, i.e., 
        sometimes we observe reads from the lost allele.
    
    Returns
    -------
    dict of {str : numpy.ndarray}
        Keys are clones, values are feature-specific CN folds of corresponding
        clone.
        The copy number (CN) fold, e.g., 1.0 for copy neutral; >1.0 for copy
        gain; and <1.0 for copy loss.
    """
    assert normal > 0
    normal = float(normal)
    
    cn_fold = {}        # clone x feature copy number fold.
    for i in range(cna_profile.shape[0]):
        rec = cna_profile.iloc[i, ]
        clone, region = rec["clone"], rec["region"]
        if clone not in cn_fold:
            cn_fold[clone] = np.repeat(1.0, p)
        assert region in cna_features
        feature_idx = cna_features[region]
        r = float(max(rec["cn"] / normal, loss_allele_freq))
        cn_fold[clone][feature_idx] = r
    return(cn_fold)
