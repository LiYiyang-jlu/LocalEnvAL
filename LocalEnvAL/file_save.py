from math import ceil
import pickle as pk 
import os
import re
from .LEF import SLEF
def SLEF_save(slef: SLEF, path:str = "./output/") -> None:
    volume = ceil(slef.S.volume)
    space_group_info = slef.S.get_space_group_info()[0]
    name = f"{slef.S.formula}({space_group_info}){volume}.SLEF"
    illegal_chars = r'[<>:"/\|?*]'
    name = re.sub(illegal_chars, '_', name)
    
    origin_path = os.getcwd()
    if not os.path.exists(path):
        os.makedirs(path)
         
    file_path = os.path.join(path, name)
    print(f"file saved at {file_path}")
    with open(file_path, "wb") as file:
        pk.dump(slef, file)
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
            