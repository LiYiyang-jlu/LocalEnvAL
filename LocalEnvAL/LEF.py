from pymatgen.core import Structure
from .RDF import RDF
from .soap_utils import *

import numpy as np
class SLEF:
    """Structure Local Environment Frame
    """
    def __init__(self, structure:Structure,
                 symmetrize: bool = True,
                 r_max = 6,
                 dr: float = 0.05,
                 sigma: float = 0.03,
                 n_max: int =6,
                 l_max: int =6) -> None:
        self.rdf = RDF(
            structure,
            symmetrize= symmetrize,
            r_max= r_max,
            dr= dr,
            sigma= sigma
        )
        self.S = self.rdf.S 
        self.species = sorted(list(set([("".join(filter(str.isalpha, str(s)))) for s in self.S.species])))
        self.results = self.rdf.rdf_cut()
        self.rcut = self.rdf.rcut
        self.r = self.rdf.r
        
        soap = compute_soap(
            rcut= self.rcut,
            n_max= n_max,
            l_max= l_max,
            sigma= 1,
            struct= self.S,
            species_list= self.species
        )
        
        self.atom_local_env_list = creat_local_env_list(self.S, 
                                                        self.rdf, 
                                                        soap, 
                                                        len(self.r),
                                                        self)

    def __getitem__(self, indice) -> "LEF":
        return self.atom_local_env_list[indice]
    def __str__(self):
        return str(self.S)
        
    def __iter__(self):
        for le in self.atom_local_env_list:
            yield le

class LEF:
    def __init__(self, structure, 
                 indice, 
                 rdf: np.ndarray, 
                 soap, 
                 parent: SLEF) -> None:
        self.S = structure
        self.indice = indice
        self.rdf = rdf
        self.soap = soap
        self.specie = str(self.S.species[indice])        
        self.parent: SLEF = parent
    
    def __str__(self) -> str:
       # log = '--- Local Enviroment ---\n'
       # log += f'rdf= {self.rdf}\n'
       # log += f'soap= {}'   
        log = f"{self.parent.S.formula}: {self.specie}{self.indice}" 
        return log    


        

def creat_local_env_list(structure: Structure, rdf : RDF, soap, x_lim, obj) -> list[LEF]:
    local_env_list: list[LEF] = []

    rdf.compute_single_atom()
    for key in rdf.g_single_atom.keys():
        for sub_key in rdf.g_single_atom[key].keys():
            radial_dist_func = rdf.g_single_atom[key][sub_key]
            radial_dist_func : np.ndarray = np.array(radial_dist_func[:x_lim])
            local_env = LEF(structure= structure,
                            indice= int(sub_key), 
                            rdf= radial_dist_func,#rdf.g_single_atom[key][sub_key],#[:r_len], 
                            soap= soap[sub_key],
                            parent= obj)

            for le in local_env_list:
                if le == local_env:
                    continue
            local_env_list.append(local_env)        
            
    return local_env_list