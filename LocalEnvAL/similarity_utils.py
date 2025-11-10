from .LEF import LEF
import numpy as np
from .const.PeriodicTable import PerTab

def pseudo_JS_div(x:LEF, y:LEF):
    """
    based on JS gradient, but x, y  are not normalized
    JS gradient' range is [0, ln2], smaller is more similar
    
    Args:
        x (LEF): _description_
        y (LEF): _description_

    Returns:
        float: JS gradient
    """
    rdf_x = x.rdf
    rdf_y = y.rdf
    cut = min(rdf_x.shape[0], rdf_y.shape[0])
    rdf_x = rdf_x[:cut]
    rdf_y = rdf_y[:cut]
    tar_x = rdf_x.flatten()
    tar_y = rdf_y.flatten()
    def JS_div(x,y):
        coeff = max(np.trapezoid(x), np.trapezoid(y))
        x /= coeff 
        y /= coeff 
        x += 1e-6
        y += 1e-6
        m = (x+y)/2
        KL_xm = np.sum(x * np.log(x/m))
        KL_ym = np.sum(y * np.log(y/m))
        return (KL_xm + KL_ym) * 0.5
    
    return JS_div(tar_x, tar_y)

def similar_loc_env_search(target, slef_queue, threshold:float= 0.1):
    similar_lef_list: list = []
    similarities:list = []
    ln2 = np.log(2)
    for slef in slef_queue:
        for lef in slef:
            tmp_similarity = 1- pseudo_JS_div(target, lef)/ln2
            if tmp_similarity > 0.9:
                similar_lef_list.append(lef)
                similarities.append(tmp_similarity)
    return similar_lef_list, similarities