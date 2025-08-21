import os
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from ase.io import read
from ase.neighborlist import NeighborList, natural_cutoffs

# Set default font properties
# plt.rc('font', size=18)  # Default text font size
# plt.rc('axes', titlesize=18)  # Title font size
# plt.rc('axes', labelsize=18)  # X and Y label font size
# plt.rc('xtick', labelsize=16)  # X tick label font size
# plt.rc('ytick', labelsize=16)  # Y tick label font size
# plt.rc('legend', fontsize=16)  # Legend font size
# plt.rc('figure', titlesize=20)  # Figure title font size
#plt.rcParams['font.family'] = 'Times New Roman'
# from xgboost import XGBRegressor

from ase.data import covalent_radii as ase_covalent_radii, atomic_numbers

def round_scientific(number, decimals):
    # Format the number in scientific notation with the desired number of decimals
    formatted_number = f"{number:.{decimals}e}"
    return formatted_number

covalent_radii = np.array([1,
                           0.32, 0.46, 1.20, 0.94, 0.77, 0.75, 0.71, 0.63, 0.64, 0.67,
                           1.40, 1.25, 1.13, 1.04, 1.10, 1.02, 0.99, 0.96, 1.76, 1.54,
                           1.33, 1.22, 1.21, 1.10, 1.07, 1.04, 1.00, 0.99, 1.01, 1.09,
                           1.12, 1.09, 1.15, 1.10, 1.14, 1.17, 1.89, 1.67, 1.47, 1.39,
                           1.32, 1.24, 1.15, 1.13, 1.13, 1.08, 1.15, 1.23, 1.28, 1.26,
                           1.26, 1.23, 1.32, 1.31, 2.09, 1.76, 1.62, 1.47, 1.58, 1.57,
                           1.56, 1.55, 1.51, 1.52, 1.51, 1.50, 1.49, 1.49, 1.48, 1.53,
                           1.46, 1.37, 1.31, 1.23, 1.18, 1.16, 1.11, 1.12, 1.13, 1.32,
                           1.30, 1.30, 1.36, 1.31, 1.38, 1.42, 2.01, 1.81, 1.67, 1.58,
                           1.52, 1.53, 1.54, 1.55])

cnmax_values = {
    'top': 12,
    'bridge': 18,
    'three_fold': 22,
    'four_fold': 26
}

def calculate_triangle_center(list_coordinates):
    return np.mean(list_coordinates, axis=0)

def is_valid_bri(atoms_positions, max_length=3.0):
    if np.linalg.norm(atoms_positions[0] - atoms_positions[1]) > max_length:
        return False
    return True

def calculate_normal(v1, v2, v3):
    """Calculate the unit normal vector of a triangle defined by three points"""
    normal = np.cross(v2 - v1, v3 - v1)
    return normal / np.linalg.norm(normal)

def calculate_cn(structure, metal, cutoff=3.0):  # Note: Cutoff is metal-specific
    """Calculate the coordination number for each atom in the structure, considering only metal atoms."""
    distances = structure.get_all_distances(mic=True)
    cn = np.zeros(len(structure))
    for i, atom in enumerate(structure):
        if atom.symbol == metal:
            neighbors = np.where((distances[i] < cutoff) & (distances[i] > 0))[0]
            cn[i] = sum(1 for n in neighbors if structure[n].symbol == metal)
    return cn

def get_CN(atoms, mult=1.0):
    radii = natural_cutoffs(atoms, mult=mult)
    nl = NeighborList(radii, bothways=True, self_interaction=False)
    nl.update(atoms)
    CN_matrix = nl.get_connectivity_matrix(sparse=False)
    CN = CN_matrix.sum(axis=0)
    return CN, CN_matrix

def calculate_gcn(structure, cn, cnmax, metal='Ru'):
    """Calculate the GCN for each atom in the structure, considering only metal atoms."""
    distances = structure.get_all_distances(mic=True)
    gcn = np.zeros(len(structure))
    for i, atom in enumerate(structure):
        if atom.symbol == metal:
            neighbors = np.where((distances[i] < 3) & (distances[i] > 0))[0]
            cn_neighbors = [cn[n] for n in neighbors if structure[n].symbol == metal]
            gcn[i] = np.sum(cn_neighbors) / cnmax
    return np.round(gcn, 2)  # 保留两位小数

def get_GCN(CN, CN_matrix, cn_max: float=None):
    if cn_max is None:
        cn_max = CN_matrix.sum(axis=0).max()
    num_atoms = CN_matrix.shape[0]
    GCN = np.zeros(num_atoms)
    for A in range(num_atoms):
        for B in range(num_atoms):
            GCN[A] += CN_matrix[A, B] * CN[B]
        GCN[A] /=  cn_max
    return GCN

def calculate_gcn_for_sites(structure, site_indices, gcn):
    """Calculate GCN values for specified sites in a structure, considering only metal atoms."""
    gcn_values = []
    for sites in site_indices:
        total_gcn = 0
        for site in sites:
            total_gcn += gcn[site]
        gcn_values.append(round(total_gcn / len(sites), 2))  # 保留两位小数
    return gcn_values

def calculate_cluster_size(structure, metal='Ru'):
    """Calculate the cluster size based on the maximum distance between metal atoms."""
    metal_positions = [atom.position for atom in structure if atom.symbol == metal]
    max_distance = 0
    for i in range(len(metal_positions)):
        for j in range(i + 1, len(metal_positions)):
            distance = np.linalg.norm(metal_positions[i] - metal_positions[j])
            if distance > max_distance:
                max_distance = distance
    return max_distance



# Set the global multiplier for bond cutoff
GLOBAL_MULT = 1.0  # <-- Set global mult here

def get_connection(atoms_in, metal='Cu', mult=1.0):
    from itertools import combinations  # Only need combinations now
    """
    Build connection info for metal atoms using direct distance and covalent radii sum as bond criterion.
    Handles periodic boundary conditions.
    """
    if mult is None:
        mult = GLOBAL_MULT  # Use global mult if not specified

    # Filter only metal atoms and keep their indices in the original atoms object
    metal_indices = [i for i, atom in enumerate(atoms_in) if atom.symbol == metal]
    atoms = atoms_in[metal_indices]  # ASE Atoms object of only metal atoms

    # Get covalent radii for all metal atoms
    from ase.data import covalent_radii as ase_covalent_radii, atomic_numbers
    radii = [ase_covalent_radii[atomic_numbers[metal]]] * len(metal_indices)
    radii = np.array(radii) * mult

    n = len(metal_indices)
    CN_matrix = np.zeros((n, n), dtype=int)
    connections = {i: [] for i in range(n)}

    for i in range(n):
        for j in range(i + 1, n):
            d = atoms_in.get_distance(metal_indices[i], metal_indices[j], mic=True)
            cutoff = radii[i] + radii[j]
            if d < cutoff:
                CN_matrix[i, j] = 1
                CN_matrix[j, i] = 1
                connections[i].append(j)
                connections[j].append(i)

    CN = CN_matrix.sum(axis=0)
    exposed_top_sites = [i for i, cn in enumerate(CN) if cn < 12]

    cn_of_connected_atoms = {}
    for i in range(n):
        connected_indices = connections[i]
        cn_of_connected_atoms[i] = [CN[j] for j in connected_indices]

    exposed_connections = {i: [j for j in connections[i] if j in exposed_top_sites] for i in exposed_top_sites}
    # Bridge sites (pairs of connected exposed top sites)
    bridge_sites = []
    long_bridge_sites = []
    checked_pairs = set()
    checked_long_pairs = set()
    for i in exposed_top_sites:
        for j in exposed_top_sites:
            if i < j and j in exposed_connections.get(i, []) and i in exposed_connections.get(j, []):
                pair = tuple(sorted([i, j]))
                if pair not in checked_pairs:
                    idx_a = metal_indices[i]
                    idx_b = metal_indices[j]
                    d_ab = atoms_in.get_distance(idx_a, idx_b, mic=True)
                    # Get covalent radii for both atoms
                    r_a = ase_covalent_radii[atomic_numbers[metal]]
                    r_b = ase_covalent_radii[atomic_numbers[metal]]
                    radii_sum = r_a + r_b
                    # Normal bridge: distance <= radii_sum
                    if d_ab <= radii_sum:
                        bridge_sites.append(sorted([i, j]))
                        checked_pairs.add(pair)
                    # Long bridge: radii_sum < distance <= radii_sum * 1.5
                    elif radii_sum < d_ab <= radii_sum * 1.5:
                        if pair not in checked_long_pairs:
                            long_bridge_sites.append(sorted([i, j]))
                            checked_long_pairs.add(pair)
                
    # Hollow sites (triplets of closely connected exposed top sites forming a closed loop or nearly closed loop)
    hollow_sites = []
    checked_hollow = set()
    for a, b, c in combinations(exposed_top_sites, 3):
        idx_a = metal_indices[a]
        idx_b = metal_indices[b]
        idx_c = metal_indices[c]
        d_ab = atoms_in.get_distance(idx_a, idx_b, mic=True)
        d_bc = atoms_in.get_distance(idx_b, idx_c, mic=True)
        d_ca = atoms_in.get_distance(idx_c, idx_a, mic=True)
        dists = [d_ab, d_bc, d_ca]
        num_short = sum(d < 3.0 for d in dists)
        num_mid = sum(3.0 <= d < 3.6 for d in dists)
        # Consider triangles with either all three bonds < 3.0, or two < 3.0 and one < 3.6
        if (num_short == 3) or (num_short == 2 and num_mid == 1):
            hollow = tuple(sorted([a, b, c]))
            if hollow not in checked_hollow:
                hollow_sites.append(sorted([a, b, c]))
                checked_hollow.add(hollow)

    # Square sites (quartets of connected exposed top sites forming a closed loop)
    square_sites = []
    checked = set()
    for i, j, k, l in combinations(exposed_top_sites, 4):
        # Only check one cyclic order (i->j->k->l->i)
        if (j in exposed_connections[i] and
            k in exposed_connections[j] and
            l in exposed_connections[k] and
            i in exposed_connections[l]):
            square = tuple(sorted([i, j, k, l]))
            if square not in checked:
                square_sites.append(sorted([i, j, k, l]))
                checked.add(square)
    print(connections)
    return connections, cn_of_connected_atoms, exposed_top_sites, bridge_sites, hollow_sites, square_sites

def number_to_letter(num):
    """Conver the CN to letters"""
    if num == 0:
        return '0'
    return chr(ord('a') + num - 1)


def get_CN_GA(path, mult=1.0, metal='Cu'):
    
    '''
    Important: it is the path for the clean cluster. 
    1) the exposed surface of the cluster will be scanned and all the top, bridge, and hollow sites will be stored.
    2) the sites will be stored as a dictionray 
    3) all the groups will be stored too to  fix the order of the  groups in the matrix
    '''    
    file_in = os.path.join(path,'POSCAR')
    atoms = read(file_in)
    print(metal)
    connections, cn_of_connected_atoms, exposed_top_sites, bridge_sites, hollow_sites, square_sites  = get_connection(atoms, metal=metal)

    # Save exposed_top_sites to a file as a list (1-based indices)
    exposed_top_sites_list = list(exposed_top_sites)  # already 0-based
    exposed_top_sites_file = os.path.join(path, 'exposed_top_sites.txt')
    with open(exposed_top_sites_file, 'w') as f:
        json.dump(exposed_top_sites_list, f)
    
    def get_one_site_string(site):
        """Obtain the text string representation for single site """
        cn_values = sorted(cn_of_connected_atoms[site])
        print(cn_values)
        # print(cn_values)
        if len(cn_values) <= 13:
            cn_values += [0] * (13 - len(cn_values))  # Pad with zeros if less than 12 elements
        else: 
            cn_values += [0] * (13 - 12)  # Pad with zeros if less than 12 elements
        
        cn_string = ''.join(number_to_letter(num) for num in cn_values)
            
        return cn_string
    
    def get_groups(site_indices, site_type):
        """
        Generic function to get group strings for a site (top, bridge, hollow, square).
        site_indices: int for top, list of ints for others.
        site_type: 'top', 'bridge', 'hollow', or 'square'
        """
        if site_type == 'top':
            site = site_indices
            groups = [get_one_site_string(site)]
            groups += [get_one_site_string(s) for s in connections[site]]
            return groups
        else:
            sites = site_indices
            strings = [get_one_site_string(i) for i in sites]
            groups = ['-'.join(sorted(strings))]
            # Collect all unique surrounding atoms except the site atoms themselves
            all_neighbors = set()
            for i in sites:
                all_neighbors.update(connections[i])
            surrounding_atoms = list(all_neighbors - set(sites))
            groups += [get_one_site_string(s) for s in surrounding_atoms]
            return groups

    GA_dict = {}
    for site in exposed_top_sites:
        key = site + 1
        GA_dict[key] = get_groups(site, 'top')

    for sites in bridge_sites:
        key = '_'.join(str(i + 1) for i in sites)
        GA_dict[key] = get_groups(sites, 'bridge')

    for sites in hollow_sites:
        key = '_'.join(str(i + 1) for i in sites)
        GA_dict[key] = get_groups(sites, 'hollow')

    for sites in square_sites:
        key = '_'.join(str(i + 1) for i in sites)
        GA_dict[key] = get_groups(sites, 'square')

    GA_file = os.path.join(path, 'GA_dict.txt')

    # Writing JSON data
    with open(GA_file, 'w') as f:
        json.dump(GA_dict, f)
        
    # Summarize and extract all unique groups from GA_dict
    groups = {group for sublist in GA_dict.values() for group in sublist}
    
    # Save groups to a text file, one per line
    groups_file = os.path.join(path, 'groups.txt')
    with open(groups_file, 'w') as f:
        for group in groups:
            f.write(f"{group}\n")


def get_full_GA_matrix(cluster_path):
    """
    Generate the GA (Group Additivity) matrix for all the surface sites
    """
    GA_dict, groups = load_ga_data(cluster_path)
    # Prepare header for the GA matrix DataFrame
    header = ['site'] + groups
    GA_matrix = []

    for key,values in GA_dict.items():
        groups_geo = values
        row = [key]
        # Create the GA matrix row: first column is the site label, then counts for each group
        for group in groups:
            count = groups_geo.count(group)
            row.append(count)
        GA_matrix.append(row)

    # Create a DataFrame from the GA matrix
    df_GA = pd.DataFrame(GA_matrix, columns=header)

    
    # Save the GA matrix to CSV
    output_csv = cluster_path + "/GA_matrix_full.csv"
    
    df_GA.to_csv(output_csv, index=False)
    print("GA matrix saved as", output_csv)
    
    return df_GA

def load_ga_data(cluster_path):
    """
    Load GA_dict and groups from the specified cluster path.
    """
    ga_file = os.path.join(cluster_path, 'GA_dict.txt')
    groups_file = os.path.join(cluster_path, 'groups.txt')

    # Load GA_dict
    with open(ga_file, 'r') as f:
        GA_dict = json.load(f)

    # Load groups
    with open(groups_file, 'r') as f:
        groups = [line.strip() for line in f if line.strip()]

    return GA_dict, groups
