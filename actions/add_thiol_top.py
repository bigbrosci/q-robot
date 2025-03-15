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
from scipy.spatial.distance import euclidean
from sklearn.decomposition import PCA

anchor_atom = 'N'
metal = 'Ru'
def calculate_triangle_center(list_coordinates):
    return np.mean(list_coordinates, axis=0)

#def calculate_normal_vector_TOP(anchor_atom_index, cluster, cutoff):
#    vector_sum = np.zeros(3)
#    atom = cluster[anchor_atom_index]
#
#    for neighbor in cluster:
#        if neighbor == metal:
#            vector = atom.position - neighbor.position
#            distance = np.linalg.norm(vector)
#            if distance < cutoff and distance > 0:  # Avoid zero division
#                vector_sum += vector / distance
#
#    if np.linalg.norm(vector_sum) > 0:
#        return vector_sum / np.linalg.norm(vector_sum)
#    else:
#        return np.array([0, 0, 1])  # Default to z-direction

def calculate_normal_vector_TOP(anchor_atom_index, cluster, cutoff):
    # Extract the position of the anchor atom
    anchor_atom = cluster[anchor_atom_index]
    anchor_position = anchor_atom.position
    
    # Find neighbors within the cutoff distance
    neighbors = []
    for neighbor in cluster:
        if not np.array_equal(neighbor.position, anchor_position):  # Ensure the neighbor is not the anchor atom
            vector = neighbor.position - anchor_position
            distance = np.linalg.norm(vector)
            if distance < cutoff:
                neighbors.append(neighbor.position)
    
    # Ensure we have enough neighbors to define a plane
    if len(neighbors) < 3:
        raise ValueError("Not enough neighbors to define a plane.")
    
    # Perform PCA on the neighbors to find the plane
    neighbors = np.array(neighbors)
    pca = PCA(n_components=3)
    pca.fit(neighbors)
    
    # The normal vector to the plane is the third principal component
    normal_vector = pca.components_[2]
    
    # Normalize the normal vector
    normal_vector /= np.linalg.norm(normal_vector)
    
    return normal_vector
    

def count_short_metal_nonmetal_bonds(cluster, metal_symbol=metal, max_distance=1.0):
    """
    Count the number of bonds shorter than a specific distance between metal atoms and non-metal atoms in the cluster.

    Parameters:
    - cluster (ASE Atoms object): The cluster of atoms.
    - metal_symbol (str): Symbol of the metal atoms in the cluster.
    - max_distance (float): The maximum distance threshold for bond counting.

    Returns:
    - count (int): Number of bonds shorter than max_distance between metal atoms and non-metal atoms.
    """
    # Find indices of metal and non-metal atoms
    metal_indices = [i for i, atom in enumerate(cluster) if atom.symbol == metal_symbol]
    nonmetal_indices = [i for i, atom in enumerate(cluster) if atom.symbol != metal_symbol]

    count = 0

    # Calculate distances between metal and non-metal atoms
    for i in metal_indices:
        for j in nonmetal_indices:
            distance = cluster.get_distance(i, j, mic=True)
            if distance < max_distance:
                count += 1

    return count


def add_anchor_to_single_site(poscar_path,exposed_index,output_prefix='POSCAR'):
    cluster = read(poscar_path)
    cutoff = 3.0  # Adjust based on your system
    modified_cluster = copy.deepcopy(cluster)
    metal_atom = modified_cluster[exposed_index - 1]  # Adjust index for 0-based indexing
    print(exposed_index)
    # Calculate the normal vector for positioning the N atom
    normal_vector = calculate_normal_vector_TOP(exposed_index - 1, modified_cluster, cutoff)

    # Position for N atom
    A_position1 = metal_atom.position + normal_vector * 2.0  # 2.0 Å away from Ni atom
    A_position2 = metal_atom.position - normal_vector * 2.0  # 2.0 Å away from Ni atom


    cluster_1 = modified_cluster + Atoms(anchor_atom, positions=[A_position1])
    n_bad_bonds1 =  count_short_metal_nonmetal_bonds(cluster_1, metal, 1.7)
    
    modified_cluster = copy.deepcopy(cluster)
    cluster_2 = modified_cluster + Atoms(anchor_atom, positions=[A_position2])
    n_bad_bonds2 =  count_short_metal_nonmetal_bonds(cluster_2, metal, 1.7)
   
    final_cluster = cluster_1
    if n_bad_bonds1 > n_bad_bonds2:
        final_cluster = cluster_2
        
    destination_path = os.path.join('top', str(exposed_index)) 
    os.makedirs(destination_path, exist_ok=True)
    output_file = destination_path + '/POSCAR'

    write(output_file, final_cluster, format='vasp', vasp5=True)


def count_short_anchor_metal_bonds(cluster, anchor_index, max_distance=2.0):
    count = 0
    for i, atom in enumerate(cluster):
        #if i != anchor_index and atom.symbol != anchor_atom:
        if atom.symbol == metal:
            distance = np.linalg.norm(cluster[anchor_index].position - atom.position)
            if distance < max_distance:
                count += 1
    return count


def add_one():
    ## add single one anchor atom to the top sites. create the folders that contain the anchor atoms on the top site.
    # exposed_indices  = [45,41,26,9,44,15,6,37,43,7,38,29,32,1,12,39,19,42,13,18,25,4,33,21,30,5,36,11,8,17,28,10,20,40,35,2,34,16,3,27]
    list_file = './list'  ### list all the top sites in one line
    with open(list_file, 'r') as f_in:
        exposed_indices = f_in.readline().rstrip().split()
    # Convert the exposed indices to integers
    exposed_indices = [int(index) for index in exposed_indices]
    # Path to the POSCAR file
    poscar_path = './POSCAR_bare'
    os.makedirs('top', exist_ok=True)
    # Run the function to choose the best geometry
    for indice in exposed_indices:
        add_anchor_to_single_site(poscar_path, indice, output_prefix='POSCAR')

    #print('WARNING:'*3 + '\n' + 'Check structures with your EYEs before running jobs!')


 
def get_new_coordiates(cluster, A1, A2, B1, list_H):
    '''A1, A2, B1 are the coordinates, list_H is the list of coordinates '''
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
    
    new_list_H = [rotation.apply(i - B1) + A2 for i  in cluster.get_positions()[1:]]
    
    return new_list_H# C_final, D_final, E_final

def add_more_atoms(site):
    '''add a molecules that contains more than 2 atoms to one top site '''
    sites = [int(i) -1 for i in site.split('_')] ### the site is counted from 1, not 0 from ase
    atoms = read(site + '/POSCAR')
    coordinates = [atoms[i].position for i in sites]
    A1 = calculate_triangle_center(coordinates)  # fcc site 
    A2 = atoms[-1].position     # N site, add N first and then put the NH,NH2,NH3 et al using the N site as reference.

    atoms_temp = read('../POSCAR_temp')  ### POSCAR_temp contains the molecules coordinates
    B1 = atoms_temp[0].position  ## The first atom in the POSCAR_temp should be the anchor atom


    # Initialize list_H to store atoms forming bonds with atoms_temp[0]
    list_H = []
    # Define bond distance threshold
    bond_distance = 1.9  # in Ångstroms
    positions = atoms_temp.get_positions()

    # Iterate over all atoms except the first one (atoms_temp[0])
    for i in range(1, len(positions)):
        atom_position = positions[i]
        distance = euclidean(B1, atom_position)
        
        if distance <= bond_distance:
            list_H.append(atom_position)
    
    # Convert list_H to numpy array if needed
    list_H = np.array(list_H)
        
    new_list_H  = get_new_coordiates(atoms_temp, A1, A2, B1, list_H)
    list_ele = atoms_temp.get_chemical_symbols()[1:]
    
    out_geo = copy.deepcopy(atoms)
    
    for num, coord in enumerate(new_list_H):
        out_geo += Atoms(list_ele[num], positions=[coord])
    write(site+'/POSCAR_top', out_geo, format='vasp', vasp5=True)

def add_more():
    os.chdir('top')
    all_items = os.listdir('.')
    folders = [item for item in all_items if os.path.isdir(item)]
    folders = [folder for folder in folders if folder.count('_') == 0]
        
    for site in folders:     
        add_more_atoms(site)

### Run script, add N to the top sites, then move the NH3 to the N site.
add_one()
add_more() 


