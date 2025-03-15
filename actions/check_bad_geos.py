#!/usr/bin/env python3
import numpy as np
from ase.io import read

# Read the POSCAR using ASE
cluster = read('POSCAR')

# Define the metal element
metal = 'Ru'
bond_threshold = 1.2  # in Angstroms

# Function to calculate and print short bonds
def check_short_bonds(cluster, metal, bond_threshold):
    metal_indices = [atom.index for atom in cluster if atom.symbol == metal]
    non_metal_indices = [atom.index for atom in cluster if atom.symbol != metal]

#    print(f"Checking bonds shorter than {bond_threshold} Å between non-{metal} atoms and {metal} atoms:")

    count = 0 
    for non_metal_index in non_metal_indices:
        for metal_index in metal_indices:
            distance = np.linalg.norm(cluster[non_metal_index].position - cluster[metal_index].position)
            if distance < bond_threshold:
                count +=1 
#                print(f"Non-{metal} atom {non_metal_index} ({cluster[non_metal_index].symbol}) and {metal} atom {metal_index} ({cluster[metal_index].symbol}) have a bond length of {distance:.2f} Å")
    if count >=1: 
        print('Bad')
    else:
        print('Good')

# Run the function
check_short_bonds(cluster, metal, bond_threshold)

