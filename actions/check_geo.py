#!/usr/bin/env python3
from ase.io import read
#from ase.geometry import get_distances
import numpy as np
import sys 


metal, anchor, cutoff = sys.argv[1:]
def find_ru_bonded_to_n(poscar_path, cutoff=2.7):
    # Read the POSCAR file
    structure = read(poscar_path)
    symbols = structure.symbols
    ru_indices = [idx for idx, symbol in enumerate(symbols) if symbol == metal ]
    try:
        n_index = [idx for idx, symbol in enumerate(symbols) if symbol == anchor][0]
    except:
        print('check anchor atom')
#        n_index = [idx for idx, symbol in enumerate(symbols) if symbol == 'H'][0] # single H
    # Calculate distances from N to all Ru atoms
    distances = structure.get_distances(n_index, ru_indices, mic=True)
    short_distances = [dist for dist in distances if dist <= cutoff]
    # Find Ru atoms that are within the cutoff distance
    bonded_ru_indices = [idx  for idx, dist in enumerate(distances) if dist <= cutoff]

    # Sort the indices and format the output
    bonded_ru_indices.sort()
    formatted_indices = '_'.join(f'{idx + 1}' for idx in bonded_ru_indices)  # +1 to convert from 0-based to 1-based indexing

    return formatted_indices, short_distances  

# Ex*ample usage of the function
poscar_path = 'POSCAR'  # Specify the path to your POSCAR file
cutoff = float(cutoff)  # M-H: 2.1, M-N: 2.7 


try: 
    bonded_ru, short_distances  = find_ru_bonded_to_n(poscar_path, cutoff)
except:
    poscar_path = 'CONTCAR'  # Specify the path to your POSCAR file
    bonded_ru, short_distances  = find_ru_bonded_to_n(poscar_path, cutoff)

#print("Bonded Ru atoms to N:", bonded_ru)
print(bonded_ru)#, short_distances)
