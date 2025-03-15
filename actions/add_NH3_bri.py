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

def is_valid_bri(atoms_positions, max_length=3.0):
    if np.linalg.norm(atoms_positions[0] - atoms_positions[1]) > max_length:
        return False
    return True

def calculate_normal_vector_bri(atom_index, cluster, cutoff):
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

def add_N_to_single_site_bri(poscar_path,exposed_indexs,output_prefix='POSCAR'):
    cluster = read(poscar_path)
    cutoff = 3.0  # Adjust based on your system

    
    def get_geo(exposed_indexs):
        modified_cluster = copy.deepcopy(cluster)
        A = modified_cluster[exposed_indexs[0] - 1]  # Adjust index for 0-based indexing
        B = modified_cluster[exposed_indexs[1] - 1]  
        
        translation_vector = B.position - A.position
        # Calculate the normal vector for positioning the N atom
        normal_vector = calculate_normal_vector_bri(exposed_indexs[0] - 1, modified_cluster, cutoff)
    
        # Position for N atom
        N_position = A.position + normal_vector * 1.6 +  0.5 * translation_vector  # 2.0 Å away from Ni atom
        
        # Add N atom to the structure
        modified_cluster += Atoms('N', positions=[N_position])
        return modified_cluster

    cluster_1 = get_geo(exposed_indexs)
    cluster_2 = get_geo(exposed_indexs[::-1])
    num_bonds_1 = count_short_nitrogen_metal_bonds(cluster_1, max_distance=1.2) 
    num_bonds_2 = count_short_nitrogen_metal_bonds(cluster_2, max_distance=1.2) 

    final_cluster = cluster_1
    if num_bonds_1 > num_bonds_2: 
        final_cluster = cluster_2 

    destination_path = os.path.join('bri', '_'.join([str(i) for i in exposed_indexs])) 
    os.makedirs(destination_path, exist_ok=True)
    output_file = destination_path + '/POSCAR'

    write(output_file, final_cluster, format='vasp', vasp5=True)
    # return cluster 

def count_short_nitrogen_metal_bonds(cluster, max_distance=1.2):
    nitrogen_index = 1
    for i, atom in enumerate(cluster):
        if atom.symbol == 'N':
            nitrogen_index = i 
            
    count = 0
    for i, atom in enumerate(cluster):
        if i != nitrogen_index and atom.symbol != 'N':
            distance = np.linalg.norm(cluster[nitrogen_index].position - atom.position)
            if distance < max_distance:
                count += 1
    return count


def add_one():
    ## add single one atom 'N' to the bri sites.
    # exposed_indices  = [45,41,26,9,44,15,6,37,43,7,38,29,32,1,12,39,19,42,13,18,25,4,33,21,30,5,36,11,8,17,28,10,20,40,35,2,34,16,3,27]
    list_file = './list'  ### list all the bri sites in one line
    with open(list_file, 'r') as f_in:
        exposed_indices = f_in.readline().rstrip().split()
    print(exposed_indices)
    # Convert the exposed indices to integers
    exposed_indices = [int(index) for index in exposed_indices]
    
        # Path to the POSCAR file
    poscar_path = './POSCAR_bare'
    cluster = read(poscar_path)
    os.makedirs('bri', exist_ok=True)
    for combination in itertools.combinations(exposed_indices, 2):
        ## get the atoms that form the hollow site 
        sorted_combination = sorted(combination)
        bri_atom_indices = [i - 1 for i in sorted_combination] # indice count from 1, ase from 0
        bri_atoms_positions = [cluster[i].position for i in bri_atom_indices]

        if not is_valid_bri(bri_atoms_positions):
            continue
        
        add_N_to_single_site_bri(poscar_path,sorted_combination,output_prefix='POSCAR')
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
    ''' Add species with more than three atoms '''
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
        
    write(site+'/POSCAR_bri', out_geo, format='vasp', vasp5=True)

def add_more():
    os.chdir('bri')
    all_items = os.listdir('.')
    folders = [item for item in all_items if os.path.isdir(item)]
    folders = [folder for folder in folders if folder.count('_') == 1]
        
    for site in folders:     
        add_more_atoms(site)

### Run script, add N to the bri sites, then move the NH3 to the N site.
add_one()
add_more() 


