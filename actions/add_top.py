#!/usr/bin/env python3 

from ase.io import read, write
from ase.neighborlist import neighbor_list
from ase import Atoms
import numpy as np
import copy, sys

'''
Command to use it: 

python3 add_top.py  Num 

Note: Num is the atom index in the POSCAR, i.e. the atom you want to place CO on.

'''

def calculate_normal_vector(atom_index, cluster, cutoff):
    vector_sum = np.zeros(3)
    atom = cluster[atom_index]

    for neighbor in cluster:
        if neighbor != atom:
            vector = atom.position - neighbor.position
            distance = np.linalg.norm(vector)
            if distance < cutoff and distance > 0:  # Avoid zero division
                vector_sum += vector / distance

    if np.linalg.norm(vector_sum) > 0:
        return vector_sum / np.linalg.norm(vector_sum)
    else:
        return np.array([0, 0, 1])  # Default to z-direction

def add_N_to_single_site(poscar_path,exposed_index,output_prefix='POSCAR'):
    cluster = read(poscar_path)
    cutoff = 3.0  # Adjust based on your system
    modified_cluster = copy.deepcopy(cluster)
    ni_atom = modified_cluster[exposed_index - 1]  # Adjust index for 0-based indexing

    # Calculate the normal vector for positioning the N atom
    normal_vector = calculate_normal_vector(exposed_index - 1, modified_cluster, cutoff)

    # Position for N atom
    c_position = ni_atom.position + normal_vector * 2.0  # 2.0 Å away from Ni atom
    o_position = ni_atom.position + normal_vector * 3.128  # 3.128 Å away from Ni atom

    # Add N atom to the structure
    modified_cluster += Atoms('C', positions=[c_position])
    modified_cluster += Atoms('O', positions=[o_position])

    # Save the modified structure
    output_file = f'{output_prefix}_{exposed_index}'
    write(output_file, modified_cluster, format='vasp', vasp5=True)

def add_nitrogen_to_exposed_atoms(poscar_path, exposed_indices, output_prefix='POSCAR'):
    # cluster = read(poscar_path)
    # cutoff = 3.0  # Adjust based on your system

    for exposed_index in exposed_indices:
        add_N_to_single_site(poscar_path,exposed_index,output_prefix='POSCAR')
        # modified_cluster = copy.deepcopy(cluster)
        # ni_atom = modified_cluster[exposed_index - 1]  # Adjust index for 0-based indexing

        # # Calculate the normal vector for positioning the N atom
        # normal_vector = calculate_normal_vector(exposed_index - 1, modified_cluster, cutoff)

        # # Position for N atom
        # n_position = ni_atom.position + normal_vector * 2.0  # 2.0 Å away from Ni atom

        # # Add N atom to the structure
        # modified_cluster += Atoms('N', positions=[n_position])

        # # Save the modified structure
        # output_file = f'{output_prefix}_{exposed_index}'
        # write(output_file, modified_cluster, format='vasp', vasp5=True)

# Exposed Ni atoms indices (1-based)
# exposed_indices = [14, 15, 2, 4, 3, 19, 18, 20, 21, 29, 30, 31, 39, 38, 28, 42, 43, 40, 32, 22, 16, 23, 17, 24, 33, 34, 41, 36, 35, 37, 7, 6, 25, 26, 27]

# Path to the POSCAR file
# poscar_path = './POSCAR'

# Add nitrogen to exposed Ni atoms and save the structures
# add_nitrogen_to_exposed_atoms(poscar_path, exposed_indices)



poscar_path = './POSCAR'
exposed_index = int(sys.argv[1])
add_N_to_single_site(poscar_path,exposed_index,output_prefix='POSCAR')
