try:
    from scipy.ndimage import gaussian_filter1d
    from scipy.signal import find_peaks
    _HAS_SCIPY = True
except Exception:
    _HAS_SCIPY = False

# third part functions
from pymatgen.core import Structure

# local class and functions    
from .RDF import RDF
from .visualize import *
from .soap_utils import *
from .LEF import LEF, SLEF
from .file_save import SLEF_read, SLEF_save
