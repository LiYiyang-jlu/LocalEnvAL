from random import gauss
import numpy as np

from pymatgen.core import Structure
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

from typing import Optional, Tuple, List, Dict
from types import FunctionType
from functools import partial
from const.TableOfElement import toe

from math import floor

class RDF : 
    def __init__(self, 
                 structure: Structure,
                 symmetrize: bool = True,
                 r_max: float = 6.0,
                 dr: float = 0.1,
                 sigma: float = 0.0) -> None:
        if symmetrize:
            finder = SpacegroupAnalyzer(structure, symprec=0.006, angle_tolerance=5)
            crystal = finder.get_conventional_standard_structure()
        
        
        self.S: Structure = structure
        self.r_max = float(r_max)
        self.dr = float(dr)
        self.sigma = float(sigma)
        
        self.edges = np.arange(0.0, self.r_max + self.dr, self.dr)
        self.r = 0.5 * (self.edges[:-1] * self.edges[1:])
        self.shell_vol = _shell_volumes(self.edges)
        
        self.N = len(self.S)
        self.V = self.S.lattice.volume
        self.rho = self.N /self.V 
        
        self.species : List[str] = sorted(list({str(s.specie) for s in self.S}) )
        self._species_indices : Dict [str, np.ndarray] = {
            sp: np.array([i for i, site in enumerate(self.S) if str(site.specie) == sp ], dtype= int) for sp in self.species
        }
        
        self.g_total: Optional[np.ndarray] = None
        self.g_partial: Dict[str, np.ndarray] = {}
        self.g_single_atom: Dict[str, Dict[np.int64, np.ndarray]] = {}
    
    def _accumulate_counts(self,
                           center_indice: int) -> np.ndarray: 
        """
        counts center_atom's neighbor
        
        1.0:remove g_tot and g_partial support, 
        only keep single atom neighbor counts

        Args:
            center_indice (int): _description_

        Returns:
            np.ndarray: 2-dim matrix
        """
        #counts = np.zeros(len(self.r), dtype=float)
        neighbor_all = self.S.get_all_neighbors(r= self.r_max, include_image= True)
        neighbor = neighbor_all[center_indice]
        counts:np.ndarray = np.zeros(( len(self.edges), len(toe.keys()) ))
        
        for nb in neighbor:
            index_r = float(nb[1])/self.dr
            index_r = floor(index_r)+1 if (index_r - int(index_r) > 0.5) else floor(index_r)
            index_s = nb.label
            counts[index_r][index_s] += 1
            
        return counts        
                
        
    def compute_single_atom(self, atom_species: Optional[List[str]]= None) -> Dict:
        atom_dict:dict = {}
        element_tmp = {}
        if atom_species is None:
            atom_dict = self._species_indices
        else: 
            for key in atom_species:
                atom_dict[key] = self._species_indices[key]
        
        for key in atom_dict:
            flag = 0
            tmp_g_dict = {}
            tmp_meta_dict = {}
            
            for indice in atom_dict[key]:
                tmp = int(indice)
                counts = self._accumulate_counts(center_indice= tmp)

                expected = len(atom_dict[key])
                with np.errstate(divide='ignore', invalid= 'ignore'):
                    g = counts / expected
                    g[~np.isfinite(g)] = 0.0
                
                g = _smooth(g, self.sigma)
                subkey = indice 
                tmp_g_dict[subkey] = g 
                # tmp_meta_dict[subkey]
            self.g_single_atom[key] = tmp_g_dict
            element_tmp[key] = tmp_meta_dict 
        #self.meta["single"]
        return self.g_single_atom
# TODO:finish it

# === Local function === # 

from typing import Dict, List, Tuple, Iterable, Optional
from . import _HAS_SCIPY

def _shell_volumes(edges: np.ndarray) -> np.ndarray:
    """$\Delta V = \frac{4\pi}{3} (r_{i+1}^3-r_{i}^3)$ for each radial shell given bin edges.
 
    Args:
        edges (np.ndarray): all edges

    Returns:
        np.ndarray: volumes
    """
    return (4.0 * np.pi / 3.0) * (edges[1: ] ** 3 - edges[:-1]**3)

def _smooth(y:np.ndarray, sigma_bins: float) -> np.ndarray:
    g_new = y.copy()
    gaussian: FunctionType
    if sigma_bins and sigma_bins > 0:
        if _HAS_SCIPY: 
            from . import gaussian_filter1d            
            gaussian = gaussian_filter1d
        
        else:
            
            win = int (6 * sigma_bins + 1)
            x = np.arange(win) - win//2
            kernel = np.exp(-0.5 * (x / sigma_bins) ** 2)
            kernel /= kernel.sum()    
            gaussian = lambda xx: np.convolve(xx, kernel, mode="same")

        for i in range(y.shape()[1]):        
            g_new[:,i] = gaussian(y[:,i])
    
    return g_new

""" 
def _find_first_two_peaks_and_valleys(
    r: np.ndarray,
    g: np.ndarray,
    prominece: float =0.02,
    distance_bins: int = 5,
    drop_frac: float = 0.05,
    window_to_second: float = 0.6,
    window_abs_A: float = 0.8,
    min_run: int = 3,
    abs_floor: float = 1e-8
) -> Dict[str, Optional[float]]:
    out = {"first_peak_r": None, "second_peak_r": None, "rcut1": None, "rcut2":None}
    y = g 
    
    # find peaks
    if _HAS_SCIPY:
        from . import find_peaks
        pidx, _ = find_peaks(y, prominence=prominece, distance= distance_bins)
    else: 
        # rewrite the peaks_find algorithm with Sliding Window Method -- LiYiyang
        pidx = []
        for i in range(distance_bins, len(y)-distance_bins): 
            if y[i] == np.max(y[i-distance_bins:i+distance_bins]):
                pidx.append(i)
        pidx = np.array(pidx, dtype=int)
    
    if pidx.size == 0:
        return out 
    
    out["first_peak_r"] = float(r[pidx[0]])
    if pidx.size >= 2:
        out["second_peak_r"] = float(r[pidx[1]])
    
    def _bounded_valley_after(p0: int, p1: Optional[int]) -> Optional[float]:
        # 1. set search boundary, keep r_cut1 close to first peak
        j_start = p0 + 1
        j_end = len(y) - 2
        # 1.1 limiting from second peak
        if p1 is not None:
            j_end = min(j_end, int(p0 + window_to_second * (p1 - p0)))
        # 1.2 limiting from absolute distance
        j_abs = np.searchsorted(r, r[p0] + window_abs_A)
        j_end = min(j_end, max(j_start + 1, j_abs))
        
        if j_start >= j_end:
            return None
        seg = y[j_start : j_end + 1]
        if seg.size == 0:
            return None 
        
        peak_h = max(y[p0], abs_floor)
        thr = max(abs_floor, drop_frac * peak_h)
        
        # 2. find first local minimum
        for j in range(j_start + 1, j_end):
            if (y[j] <= thr) and (y[j] < y[j-1]) and y[j] < y[j+1]:
                return float(r[j])
        
        # 3. find first starting indice of the series which all bins less than min_run
        # it means that return when finding a flat site 
        below = seg <= thr 
        run = 0 
        j_rel = None
        for t, flag in enumerate(below):
            run = run + 1 if flag else 0
            if run >= min_run:
                j_rel = t - min_run + 1
                break
        if j_rel is not None: 
            return float(r[j_start + j_rel])
        
        # 4. find local minimum in window
        j_loc = int(np.argmin(seg))
        return float(r[j_start + j_loc])
    
    out["rcut1"] = _bounded_valley_after(int(pidx[0]), int(pidx[1]) if pidx.size >= 2 else None)
    if pidx.size >= 2:
        out["rcut2"] = _bounded_valley_after(int(pidx[1]), int(pidx[2]) if pidx.size >= 3 else None)
        
    return out

"""