#!/usr/bin/env python3 
from ase.io import read, write
from ase.geometry import get_distances

def modify_structure(poscar_path):
    # Read the POSCAR file
    structure = read(poscar_path)
    
    # Find indices of all H atoms and the single N atom
    h_indices = [atom.index for atom in structure if atom.symbol == 'H']
    n_index = [atom.index for atom in structure if atom.symbol == 'N'][0]  # Assuming there is only one N atom
    
    # Get the position of the N atom
    n_position = structure.positions[n_index:n_index+1]  # Ensure 2D array for distance calculation
    
    # Initialize a list to collect H atoms to be deleted
    h_to_delete = []
    distances = structure.get_distances(n_index, h_indices, mic=True)
    print(sorted(distances))
poscar_path = 'POSCAR'  # Specify the path to your POSCAR file
modify_structure(poscar_path)

