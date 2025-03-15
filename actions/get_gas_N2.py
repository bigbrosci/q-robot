# -*- coding: utf-8 -*-
"""
Created on Fri May 31 11:03:17 2024

@author: lqlhz
"""

#!/usr/bin/env python3
from ase.io import read, write
#from ase.geometry import get_distances
import numpy as np



def find_ru_bonded_to_n(poscar_path, cutoff=2.3):
    # Read the POSCAR file
    structure = read(poscar_path)
    symbols = structure.symbols
    n_indexs = [idx for idx, symbol in enumerate(symbols) if symbol == 'N']
    ru_indices = [idx for idx, symbol in enumerate(symbols) if symbol == 'Ru']
    # Calculate distances from N to all Ru atoms
    n_coordinates = [structure.get_positions()[i] for i in n_indexs]
    N0 = np.mean([n_coordinates[0]],axis=0)
    d_NN = structure.get_distances(n_indexs[0], n_indexs[1], mic=True)
    
    
    distances = structure.get_distances(n_indexs[0], ru_indices, mic=True)
    short_distances = [dist for dist in distances if dist <= cutoff]
    # Find Ru atoms that are within the cutoff distance
    bonded_ru_indices = [idx  for idx, dist in enumerate(distances) if dist <= cutoff]
    C1 = np.mean([structure.positions[i] for i in bonded_ru_indices], axis=0)
    
    distances = structure.get_distances(n_indexs[1], ru_indices, mic=True)
    short_distances = [dist for dist in distances if dist <= cutoff]
    # Find Ru atoms that are within the cutoff distance
    bonded_ru_indices = [idx  for idx, dist in enumerate(distances) if dist <= cutoff]
    C2 = np.mean([structure.positions[i] for i in bonded_ru_indices], axis=0)
    C0 = np.mean([C1, C2],axis = 0)
    
    v1 = (N0 - C0 ) / np.linalg.norm(N0 - C0)


    a0 = N0 + v1 * 1.2
    print(a0)
    
    v2 = (n_coordinates[1] - n_coordinates[0]) / np.linalg.norm(n_coordinates[1] - n_coordinates[0])

    N1 = a0 - 0.55 * v2 
    N2 = a0 + 0.55 * v2

    print(N1, N2)
    
    structure.positions[n_indexs[0]] = N1 
    structure.positions[n_indexs[1]] = N2
    write('POSCAR_N2',structure)
# Example usage of the function

find_ru_bonded_to_n('POSCAR', cutoff=2.3)
