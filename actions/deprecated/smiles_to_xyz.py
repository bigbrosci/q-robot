#!/usr/bin/env python3
'''Convert the SMILES TO xyz and gjf files '''
from smiles2xyz import * 
from gpm import *
from xyz2mol import *
import sys 

smiles = sys.argv[1]
coords = get_coords_from_rdkit(smiles)

name = smiles_file_name(smiles)

save_gjf(coords, name)
save_xyz(coords, name)
save_poscar(coords, name)

