# mf_umi.py - multi-feature UMIs.



import os
from ..io.base import load_feature_objects, load_samples, load_features



def mfu_main(
    fet_obj_fn,
    sample_fn,
    feature_fn,
    count_dir,
    tmp_dir
):
    """Main function for processing multi-feature UMIs.
    
    Parameters
    ----------
    fet_obj_fn : str
        The python pickle file storing `~..utils.gfeature.Feature` objects.
    sample_fn : str
        File storing cell IDs or barcodes.
    feature_fn : str
        File storing feature annotations.
    count_dir : str
        The output folder to store the count matrices.
    tmp_dir : str
        Folder to store temporary data.

    Returns
    -------
    Dict
        The results.
    """
    # check args.
    reg_list = load_feature_objects(fet_obj_fn)
    samples = load_samples(sample_fn)         # ndarray
    features = load_features(feature_fn)      # DataFrame
    features = features["feature"].to_numpy()
    
    assert len(reg_list) == features.shape[0]
    for i, reg in enumerate(reg_list):
        assert reg.name == features[i]
    
    os.makedirs(count_dir, exist_ok = True)
    os.makedirs(tmp_dir, exist_ok = True)
    
    
    # 