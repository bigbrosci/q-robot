from ase.io import read, write
from ase.neighborlist import neighbor_list
from ase import Atoms
import numpy as np
import copy, sys
import itertools
from ase.io import read, write
from scipy.spatial import ConvexHull
import numpy as np
import copy, os
from numpy.linalg import norm
from scipy.spatial.transform import Rotation as R


def calculate_triangle_center(atoms):
    return np.mean(atoms, axis=0)


def calculate_normal_vector_of_plane(atoms):
    vector1 = atoms[1] - atoms[0]
    vector2 = atoms[2] - atoms[0]
    normal_vector = np.cross(vector1, vector2)
    normal_vector /= np.linalg.norm(normal_vector)
    if np.dot(normal_vector, np.mean(atoms, axis=0)) > 0:
        normal_vector = -normal_vector
    return normal_vector

def is_valid_triangle(atoms_positions, max_length=3.0):
    for i in range(3):
        for j in range(i + 1, 3):
            if np.linalg.norm(atoms_positions[i] - atoms_positions[j]) > max_length:
                return False
    return True

def count_short_nitrogen_metal_bonds(cluster, nitrogen_index, max_distance=2.0):
    count = 0
    for i, atom in enumerate(cluster):
        if i != nitrogen_index and atom.symbol != 'N':
            distance = np.linalg.norm(cluster[nitrogen_index].position - atom.position)
            if distance < max_distance:
                count += 1
    return count

def add_one_atom(poscar_path, exposed_indices, output_prefix='POSCAR'):
    cluster = read(poscar_path)
    
    for combination in itertools.combinations(exposed_indices, 3):
        ## get the atoms that form the hollow site 
        sorted_combination = sorted(combination)
        triangle_atom_indices = [i - 1 for i in sorted_combination] # indice count from 1, ase from 0
        triangle_atoms_positions = [cluster[i].position for i in triangle_atom_indices]

        if not is_valid_triangle(triangle_atoms_positions):
            continue

        center = calculate_triangle_center(triangle_atoms_positions)
        normal_vector = calculate_normal_vector_of_plane(triangle_atoms_positions) # there should be two directions

        # Check geometry with N placed in one direction
        modified_cluster_1 = copy.deepcopy(cluster)
        n_position_1 = center + normal_vector * 1.5 # 1.5 \AA above the hollow site
        modified_cluster_1 += Atoms('N', positions=[n_position_1])
        nitrogen_index_1 = len(modified_cluster_1) - 1
        num_bonds_1 = count_short_nitrogen_metal_bonds(modified_cluster_1, nitrogen_index_1)

        # Check geometry with N placed in the opposite direction
        modified_cluster_2 = copy.deepcopy(cluster)
        n_position_2 = center - normal_vector * 1.5
        modified_cluster_2 += Atoms('N', positions=[n_position_2])
        nitrogen_index_2 = len(modified_cluster_2) - 1
        num_bonds_2 = count_short_nitrogen_metal_bonds(modified_cluster_2, nitrogen_index_2)

        # Choose the geometry with fewer short M-N bonds
        best_cluster = modified_cluster_1 if num_bonds_1 < num_bonds_2 else modified_cluster_2
        destination_path = os.path.join('hollow', '_'.join(map(str, sorted_combination))) 
        os.makedirs(destination_path, exist_ok=True)
        output_file = destination_path + '/POSCAR'
        write(output_file, best_cluster, format='vasp', vasp5=True)

def add_one():
    ## add one atom
    exposed_indices  = [45,41,26,9,44,15,6,37,43,7,38,29,32,1,12,39,19,42,13,18,25,4,33,21,30,5,36,11,8,17,28,10,20,40,35,2,34,16,3,27]
    #list_file = './list'
    #with open(list_file, 'r') as f_in:
    #    exposed_indices = f_in.readline().rstrip().split()
    # Convert the exposed indices to integers
    exposed_indices = [int(index) for index in exposed_indices]
    # Path to the POSCAR file
    poscar_path = './POSCAR_bare'
    os.makedirs('hollow', exist_ok=True)
    # Run the function to choose the best geometry
    add_one_atom(poscar_path, exposed_indices)
    print('WARNING:'*3 + '\n' + 'Check structures with your EYEs before running jobs!')

def get_new_coordiates(A1, A2, B1, list_H):
    # Calculate B2 as the center (midpoint) of C and D
    B2 = calculate_triangle_center(list_H)
    
    translation_vector = A2 - B1
    B2_translated = B2 + translation_vector 
    
    # Step 2: Calculate vectors a and b
    a = A2 - A1
    b = B2_translated - A2  # New b after moving B1 to A2
    
    # Normalize vectors a and b
    a_normalized = a / norm(a)
    b_normalized = b / norm(b)
    
    # Calculate rotation axis and angle
    axis = np.cross(b_normalized, a_normalized)
    axis_normalized = axis / norm(axis) if norm(axis) != 0 else b_normalized
    angle = np.arccos(np.clip(np.dot(b_normalized, a_normalized), -1.0, 1.0))
    
    # Apply rotation around A2
    rotation = R.from_rotvec(axis_normalized * angle)
    
    new_list_H = [rotation.apply(i - B1) + A2 for i  in list_H]
    
    return new_list_H# C_final, D_final, E_final

def add_more_atoms(site):
    sites = [int(i) -1 for i in site.split('_')]
    atoms = read(site + '/POSCAR')
    coordinates = [atoms[i].position for i in sites]
    A1 = calculate_triangle_center(coordinates)  # fcc site 
    A2 = atoms[-1].position     # N site, add N first and then put the NH,NH2,NH3 et al using the N site as reference.

    atoms_temp = read('../POSCAR_temp')
    B1 = atoms_temp[0].position

    list_H = atoms_temp.get_positions()[1:]
    
    new_list_H  = get_new_coordiates(A1, A2, B1, list_H)
    list_ele = atoms_temp.get_chemical_symbols()[1:]
    
    out_geo = copy.deepcopy(atoms)
    
    for num, coord in enumerate(new_list_H):
        out_geo += Atoms(list_ele[num], positions=[coord])
    write(site+'/POSCAR_temp', out_geo, format='vasp', vasp5=True)

def add_more():
    os.chdir('hollow')
    all_items = os.listdir('.')
    folders = [item for item in all_items if os.path.isdir(item)]
    folders = [folder for folder in folders if folder.count('_') == 2]
        
    for site in folders:     
        add_more_atoms(site)


add_one()
add_more() 


