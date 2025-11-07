from .LEF import LEF
import numpy as np
from .const.PeriodicTable import PerTab

def pseudo_JS_grad(x:LEF, y:LEF):
    """
    based on JS gradient, but x, y  are not normalized
    JS gradient' range is [0, ln2], smaller is more similar
    
    Args:
        x (LEF): _description_
        y (LEF): _description_

    Returns:
        float: JS gradient
    """
    tar_x = x.rdf.flatten()
    tar_y = y.rdf.flatten()
    def JS_grad(x,y):
        coeff = max(np.trapezoid(x), np.trapezoid(y))
        x /= coeff 
        y /= coeff 
        x += 1e-6
        y += 1e-6
        m = (x+y)/2
        KL_xm = np.sum(x * np.log(x/m))
        KL_ym = np.sum(y * np.log(y/m))
        return (KL_xm + KL_ym) * 0.5
    
    return JS_grad(tar_x, tar_y)