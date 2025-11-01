import pickle as pk 
import os
from .LEF import SLEF
def SLEF_save(SLEF: SLEF, path:str = "./output/") -> None:
    name = f"{SLEF.S.formula}({SLEF.S.get_space_group_info()[1]}).SLEF"
    origin_path = os.getcwd()
    if not os.path.exists(path):
        os.makedirs(path)
    os.chdir(path)     
    
    with open(name, "wb") as file:
        pk.dump(SLEF, file)
        os.chdir(origin_path)

def SLEF_read(path) -> SLEF|None:
    with open(path, "rb") as file:
        data = pk.load(file)
        try:
            if type(data) != SLEF:
                raise TypeError("This file is not a SLEF, please check file path")
            return data 
        except TypeError as error:
            print("Type Error:",error)
            return None
            