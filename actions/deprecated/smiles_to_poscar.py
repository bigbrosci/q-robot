#!/usr/bin/env python3
from gpm import * 
from  smiles2xyz import * 
import sys 
smiles = sys.argv[1]

#get_xyz_from_chemml(smiles)
coords = get_coords_from_chemml(smiles)
save_xyz(coords, smiles)
save_poscar(coords, 'POSCAR')
