from pymatgen.core import Structure
from dscribe.descriptors import SOAP
from pymatgen.io.ase import AseAtomsAdaptor

def load_struct(filepath: str):
    struct = Structure.from_file(filepath)
    structure_prim = struct.get_primitive_structure(tolerance= 0.1)
    return structure_prim

def get_soap(rcut, species_list, n_max=8, l_max=6, sigma=1.0):
    soap = SOAP(
        species= species_list,
        r_cut= rcut,
        n_max= n_max,
        l_max= l_max,
        sigma= sigma,
        periodic= True,
        sparse= False
    )
    return soap 

def to_ase_atoms(struct):
    ase_atoms = AseAtomsAdaptor.get_atoms(struct)
    return ase_atoms

def compute_soap(rcut, n_max, l_max, sigma, struct, species_list):
    soap = get_soap(rcut, species_list, n_max, l_max, sigma)
    atoms = to_ase_atoms(struct)
    soap_matrix = soap.create(atoms, n_jobs=1)
    return soap_matrix