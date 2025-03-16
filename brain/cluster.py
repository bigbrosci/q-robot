import os, sys
import csv
import json
import numpy as np
import pandas as pd
import itertools
import copy
import math
import matplotlib.pyplot as plt
import joblib  # import the joblib library
from ase import Atoms, io
import platform
import re



# Set default font properties
# plt.rc('font', size=18)  # Default text font size
# plt.rc('axes', titlesize=18)  # Title font size
# plt.rc('axes', labelsize=18)  # X and Y label font size
# plt.rc('xtick', labelsize=16)  # X tick label font size
# plt.rc('ytick', labelsize=16)  # Y tick label font size
# plt.rc('legend', fontsize=16)  # Legend font size
# plt.rc('figure', titlesize=20)  # Figure title font size
plt.rcParams['font.family'] = 'Times New Roman'


from matplotlib.patches import Patch
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap
from matplotlib import cm
from matplotlib.colors import ListedColormap


from scipy.spatial import KDTree, ConvexHull
from scipy.spatial.transform import Rotation as R
from scipy.spatial.distance import euclidean
from scipy.integrate import odeint
from scipy.integrate import solve_ivp
from scipy.constants import h,  eV

from ase import Atoms
from ase.io import read, write
from ase.neighborlist import NeighborList, natural_cutoffs, neighbor_list
from ase.io.vasp import read_vasp, write_vasp
from ase.visualize import view

from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
from sklearn.model_selection import cross_val_score, KFold
from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet, BayesianRidge, LogisticRegression
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.svm import SVR
from sklearn.metrics import make_scorer, mean_absolute_error, mean_squared_error, r2_score
from sklearn.neural_network import MLPRegressor
from sklearn.decomposition import PCA
from sklearn.cross_decomposition import PLSRegression
from sklearn.pipeline import make_pipeline
from xgboost import XGBRegressor
from scipy.optimize import curve_fit
from scipy.stats import pearsonr
from sklearn.metrics import mean_absolute_error, mean_squared_error



# def convert_unix_to_windows_path(unix_path):
#     # Replace forward slashes with backslashes
#     windows_path = unix_path.replace('/', '\\')
    
#     # Check if the path starts with '/mnt/', which indicates a mounted drive in WSL
#     if windows_path.startswith('\\mnt\\'):
#         # Extract the drive letter and the rest of the path
#         drive_letter = windows_path[5]
#         rest_of_path = windows_path[6:]
#         # Construct the Windows path
#         windows_path = f'{drive_letter.upper()}:\\{rest_of_path}'
    
#     return windows_path


# def convert_linux_path_to_relative(path):
#     # Only adjust the path on Windows
#     if platform.system() == 'Windows':
#         # This regex looks for any folder ending with _ads, followed by a path
#         match = re.search(r'[^/\\]+_ads([/\\].*)', path)
#         if match:
#             # Return the relative path starting with '.' 
#             # (e.g., './1/POSCAR' from '/.../H_ads/1/POSCAR')
#             return '.' + match.group(1)
#     return path




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

# 设置不同吸附位点的 cnmax 值
cnmax_values = {
    'top': 12,
    'bridge': 18,
    'three_fold': 22,
    'four_fold': 26
}

gas_dict = {
    "CH": -6.05202894,
    "CH2": -11.78901352,
    "CH3": -18.19317746,
    "CO": -14.78130054,
    "H": -1.11671722,
    "H2": -6.76668776,
    "N": -3.12298738,
    "N2": -16.62922486,
    "NH": -8.10061060,
    "NH2": -13.53307777,
    "NH3": -19.53586573,
    "O": -1.90384114,
    "O2": -9.85600535,
    "OH": -7.72794651,
    "OH2": -14.21849996,
    "H2O": -14.21849996,
    "S": -1.08009627,
    "SH": -6.01226074,
    "SH2": -11.19981098,
    "H2S": -11.19981098,
    "thiol": -55.96808919,
    'thiolate_lateral':-50.67527044,
    "thiolate_upright": -50.67527044,
    "slab":-354.42635463
}

def calculate_triangle_center(list_coordinates):
    return np.mean(list_coordinates, axis=0)

def is_valid_bri(atoms_positions, max_length=3.0):
    if np.linalg.norm(atoms_positions[0] - atoms_positions[1]) > max_length:
        return False
    return True

def calculate_normal(p1, p2, p3):
    """Calculate the normal vector of the triangle formed by points p1, p2, p3"""
    v1 = p2 - p1
    v2 = p3 - p1
    normal = np.cross(v1, v2)
    return normal / np.linalg.norm(normal)

def calculate_cn(structure, metal, cutoff=2.7):
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


def get_top_sites(path, metal = 'Ru', mult=0.9):
    """Obtain the exposed  """
    try:
        atoms_in = read(path + '/CONTCAR')
    except:
        atoms_in = read(path + '/POSCAR')
    filtered_atoms = [atom for atom in atoms_in if atom.symbol  in [metal]]
    atoms = Atoms(filtered_atoms)
    radii = natural_cutoffs(atoms, mult=mult)
    nl = NeighborList(radii, bothways=True, self_interaction=False)
    nl.update(atoms)
    CN_matrix = nl.get_connectivity_matrix(sparse=False)
    CN = CN_matrix.sum(axis=0)
    exposed_top_sites = [i for i, cn in enumerate(CN) if cn <= 9]
    return exposed_top_sites

def get_connection(path, metal='Ru', mult=0.9):
    
    # if platform.system() == 'Windows':
    #     print('Windows')
    #     path = convert_linux_path_to_relative(path)
    #     print(path)
    atoms_in = read(path + '/POSCAR')
    
    filtered_atoms = [atom for atom in atoms_in if atom.symbol in [metal]]
    atoms = Atoms(filtered_atoms)
    radii = natural_cutoffs(atoms, mult=mult)
    nl = NeighborList(radii, bothways=True, self_interaction=False)
    nl.update(atoms)
    CN_matrix = nl.get_connectivity_matrix(sparse=False)
    CN = CN_matrix.sum(axis=0)
    
    connections = {}     # Get the connections of all atoms
    for i in range(len(atoms)):
        connections[i] = list(np.where(CN_matrix[i])[0])
        
    cn_of_connected_atoms = {}      # Request 2: Calculate the CN of the connected atoms
    for i in range(len(atoms)):
        connected_indices = connections[i]
        cn_of_connected_atoms[i] = [CN[j] for j in connected_indices]

    exposed_top_sites = [i for i, cn in enumerate(CN) if cn <= 9]       # Exposed top sites (CN <= 9)

    return atoms, connections, cn_of_connected_atoms, exposed_top_sites

def number_to_letter(num):
    """Conver the CN to letters"""
    if num == 0:
        return '0'
    return chr(ord('a') + num - 1)


def get_CN_GA(path, mult=0.9):
    
    '''
    Important: it is the path for the clean cluster. 
    1) the exposed surface of the cluster will be scanned and all the top, bridge, and hollow sites will be stored.
    2) the sites will be stored as a dictionray 
    3) all the groups will be stored too to  fix the order of the  groups in the matrix
    '''    
    
    atoms, connections, cn_of_connected_atoms, exposed_top_sites = get_connection(path, metal='Ru', mult=0.9)
    
    def get_one_site_string(site):
        """Obtain the text string representation for single site """
        cn_values = sorted(cn_of_connected_atoms[site])
        if len(cn_values) < 13:
            cn_values += [0] * (13 - len(cn_values))  # Pad with zeros if less than 12 elements
            cn_string = ''.join(number_to_letter(num) for num in cn_values)
        return cn_string
    
    def get_groups_one_top_site(site):
        groups = []
        groups.append(get_one_site_string(site))
        surrounding_atoms = connections[site]
        for s_atom in surrounding_atoms:
            groups.append(get_one_site_string(s_atom))
        return groups
        
    def get_groups_one_bridge_site(sites) :
        groups = []
        strings = [get_one_site_string(i) for i in sites]
        sorted_strings = sorted(strings)
        bri_site = '-'.join(sorted_strings)
        groups.append(bri_site)
        
        # Step 1: Extend all_atoms with atoms from connections at indices specified by sites
        all_atoms = []
        all_atoms.extend(connections[i] for i in sites)
        all_atoms = [item for sublist in all_atoms for item in sublist]

        surrounding_atoms = list(set(all_atoms))
        for i in sites: 
            surrounding_atoms.remove(i)
        
        for s_atom in surrounding_atoms:
            groups.append(get_one_site_string(s_atom))
        return groups
  
    def get_groups_one_hollow_site(sites) :
        groups = []
        strings = [get_one_site_string(i) for i in sites]
        sorted_strings = sorted(strings)
        hollow_site = '-'.join(sorted_strings)
        groups.append(hollow_site)
        
        # Step 1: Extend all_atoms with atoms from connections at indices specified by sites
        all_atoms = []
        all_atoms.extend(connections[i] for i in sites)
        
        # Step 2: Flatten the list of lists into a single list
        all_atoms = [item for sublist in all_atoms for item in sublist]
        
        # Step 3: Remove duplicates to get surrounding_atoms
        surrounding_atoms = list(set(all_atoms))

        for i in sites: 
            surrounding_atoms.remove(i)
        for s_atom in surrounding_atoms:
            groups.append(get_one_site_string(s_atom))
        return groups    

    
    # Filter connections to only include exposed top sites
    exposed_connections = {i: [j for j in connections[i] if j in exposed_top_sites] for i in exposed_top_sites}
    
    # Bridge sites (pairs of connected exposed top sites)
    bridge_sites = []
    for i in exposed_top_sites:
        for j in exposed_connections[i]:
            if i < j:  # Ensure each pair is considered only once
                bridge_sites.append([i, j])
    
    # Hollow sites (triplets of closely connected exposed top sites)
    hollow_sites = []
    for i in exposed_top_sites:
        for j in exposed_connections[i]:
            for k in exposed_connections[j]:
                if i < j and j < k and i in exposed_connections[k]:  # Ensure each triplet is considered only once and forms a triangle
                    hollow_sites.append([i, j, k])

    GA_dict = {}
    for site in exposed_top_sites:
        key = site + 1 
        GA_dict[key] = get_groups_one_top_site(site)

    for sites in bridge_sites:
        key = '_'.join([str(i+1) for i in sites])
        GA_dict[key] = get_groups_one_bridge_site(sites)

    for sites in hollow_sites:
        key = '_'.join([str(i+1) for i in sites])
        GA_dict[key] = get_groups_one_bridge_site(sites)

    
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

# Define bonding criteria (in Angstrom)
R_N_BOND_MAX = 2.4  # Max Ru-N bond distance
R_H_BOND_MAX = 2.2  # Max Ru-H bond distance
N_N_BOND_MAX = 1.3  # Max N-N bond distance for N2 detection

def get_bonded_atoms(atoms, atom_index, element, cutoff):
    """Find all atoms of a given element bonded to a specified atom."""
    bonded = []
    for i, atom in enumerate(atoms):
        if atom.symbol == element and i != atom_index:
            dist = atoms.get_distance(atom_index, i)
            if dist < cutoff:
                bonded.append(i)
    return bonded

def format_ru_sites(ru_indices):
    """Convert a list of Ru indices into a sorted string format."""
    ru_indices = [i + 1 for i in ru_indices]  # Convert to 1-based index
    ru_indices.sort()  # Sort in ascending order
    return "_".join(map(str, ru_indices)) if ru_indices else "N/A"

def classify_n2_adsorption(atoms):
    """Classify N2 adsorption sites based on Ru coordination."""
    nitrogen_indices = [i for i, atom in enumerate(atoms) if atom.symbol == "N"]
    n2_molecules = []
    visited = set()

    # Identify N2 pairs based on N-N bond distance
    for i in nitrogen_indices:
        if i in visited:
            continue
        for j in nitrogen_indices:
            if i != j and atoms.get_distance(i, j) < N_N_BOND_MAX:
                n2_molecules.append((i, j))
                visited.add(i)
                visited.add(j)
                break

    results = []
    for n1, n2 in n2_molecules:
        ru_bonded_n1 = get_bonded_atoms(atoms, n1, "Ru", R_N_BOND_MAX)
        ru_bonded_n2 = get_bonded_atoms(atoms, n2, "Ru", R_N_BOND_MAX)
        bond_counts = tuple(sorted([len(ru_bonded_n1), len(ru_bonded_n2)]))

        # Classify adsorption type
        if bond_counts == (0, 1):
            adsorption_type = "top1"
            # Identify which N is bonded to Ru
            if len(ru_bonded_n1) == 1:  # If N1 is bonded to Ru, use its Ru index
                ru_site_str = format_ru_sites(ru_bonded_n1)
            else:  # Otherwise, N2 must be bonded, so use its Ru index
                ru_site_str = format_ru_sites(ru_bonded_n2)
        elif bond_counts == (1, 1) and ru_bonded_n1 == ru_bonded_n2:
            adsorption_type = "top2"
            ru_site_str = format_ru_sites(ru_bonded_n1)
        elif bond_counts == (1, 1):
            adsorption_type = "bridge-1"
            ru_site_str = format_ru_sites(set(ru_bonded_n1 + ru_bonded_n2))
        elif bond_counts == (1, 2):
            adsorption_type = "bridge-2"
            ru_site_str = format_ru_sites(set(ru_bonded_n1 + ru_bonded_n2))
        elif bond_counts == (2, 2) and len(set(ru_bonded_n1 + ru_bonded_n2)) == 3:
            adsorption_type = "fcc-1"
            ru_site_str = format_ru_sites(set(ru_bonded_n1 + ru_bonded_n2))
        elif bond_counts == (2, 3):
            adsorption_type = "fcc-2"
            ru_site_str = format_ru_sites(set(ru_bonded_n1 + ru_bonded_n2))
        elif bond_counts == (3, 3):
            adsorption_type = "fcc-3"
            ru_site_str = format_ru_sites(set(ru_bonded_n1 + ru_bonded_n2))
        elif bond_counts == (1, 3):
            adsorption_type = "fcc-4"
            ru_site_str = format_ru_sites(set(ru_bonded_n1 + ru_bonded_n2))
        else:
            adsorption_type = "unknown"
            ru_site_str = "N/A"

        results.append((n1, n2, adsorption_type, ru_site_str))
    # return results
    return ru_site_str

def classify_single_atom_adsorption(atoms, element, bond_cutoff):
    """Classify adsorption of single atoms (N or H) based on Ru coordination."""
    atom_indices = [i for i, atom in enumerate(atoms) if atom.symbol == element]
    results = []

    for idx in atom_indices:
        ru_bonded = get_bonded_atoms(atoms, idx, "Ru", bond_cutoff)

        # Determine adsorption type
        if len(ru_bonded) == 1:
            adsorption_type = "top"
        elif len(ru_bonded) == 2:
            adsorption_type = "bridge"
        elif len(ru_bonded) == 3:
            adsorption_type = "hollow"
        else:
            adsorption_type = "unknown"

        ru_site_str = format_ru_sites(ru_bonded)
        results.append((idx, adsorption_type, ru_site_str))

    # return results
    return ru_site_str

def load_ga_data(cluster_path):
    """
    Loads the GA dictionary and the groups list from files.

    Parameters:
        load_path (str): The directory where the files are stored.
                         Files are expected to be named 'GA_dict.txt' and 'groups.txt'.
    
    Returns:
        tuple: A tuple containing the GA dictionary and the groups list.
    """
    # Load GA_dict from JSON file
    ga_file = os.path.join(cluster_path, 'GA_dict.txt')
    with open(ga_file, 'r') as f:
        GA_dict = json.load(f)
    
    # Load groups from the text file
    groups_file = os.path.join(cluster_path, 'groups.txt')
    with open(groups_file, 'r') as f:
        groups = [line.strip() for line in f if line.strip()]
    
    return GA_dict, groups



def get_dipole_polarization(GA_matrix):
    """
    Perform quadratic curve fitting on the Eads columns of the GA_matrix DataFrame for each row.
    The x-values are [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6] and the y-values are taken from the 
    corresponding columns: 'Eads_-0.6', 'Eads_-0.4', 'Eads_-0.2', 'Eads_0.0', 'Eads_0.2', 'Eads_0.4', 'Eads_0.6'.
    
    For each row, a quadratic function is fitted:
         func(x, a, b, c) = -0.5 * a * x**2 - b * x + c
    where the parameter 'a' is interpreted as the polarizability and 'b' as the dipole.
    In addition, the goodness-of-fit metrics R², MAE, and RMSE are computed using scikit-learn's functions.
    
    Returns:
      A DataFrame with the following columns:
         'site', 'polarizability', 'dipole', 'c', 'R2', 'MAE', 'RMSE'
    """
    
    # Define the quadratic function.
    def func(x, a, b, c):
        return -0.5 * a * x**2 - b * x + c

    # Define the x values corresponding to the Eads columns.
    x_values = np.array([-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6])
    
    # Determine the correct column for the site label.
    site_col = 'Site' if 'Site' in GA_matrix.columns else 'site'
    
    # Define the Eads columns exactly as in GA_matrix.csv.
    eads_cols = ["Eads_-0.6", "Eads_-0.4", "Eads_-0.2", "Eads_0.0", "Eads_0.2", "Eads_0.4", "Eads_0.6"]
    
    results_list = []
    
    # Loop through each row of the GA_matrix.
    for idx, row in GA_matrix.iterrows():
        site = row[site_col]
        # Extract y values for the current row.
        y_values = row[eads_cols].values.astype(float)
        
        # Remove any NaN values.
        valid_mask = ~np.isnan(y_values)
        x_fit = x_values[valid_mask]
        y_fit = y_values[valid_mask]
        
        # Ensure there are at least 3 data points to perform a quadratic fit.
        if len(x_fit) < 3:
            print(f"Skipping site '{site}' (index {idx}) due to insufficient data points.")
            continue
        
        try:
            # Perform curve fitting.
            params, covariance = curve_fit(func, x_fit, y_fit)
            a, b, c = params

            # Predicted values for the valid x data.
            y_pred = func(x_fit, *params)
            
            # Calculate R² using scikit-learn's r2_score.
            R2 = r2_score(y_fit, y_pred)
            
            # Calculate MAE and RMSE.
            mae = mean_absolute_error(y_fit, y_pred)
            rmse = np.sqrt(mean_squared_error(y_fit, y_pred))
            
            # Append the results using descriptive names for a and b.
            results_list.append({
                'site': site,
                'polarizability': a,
                'dipole': b,
                'c': c,
                'R2': R2,
                'MAE': mae,
                'RMSE': rmse
            })
        except Exception as e:
            print(f"Error fitting curve for site '{site}' (index {idx}): {e}")
            continue
    
    # Convert the list to a DataFrame and return.
    results_df = pd.DataFrame(results_list)
    return results_df

def update_GA_matrix_with_dipole_polarization(csv_filepath='./GA_matrix.csv'):
    
    """
    Update the GA_matrix by adding 'polarizability' and 'dipole' columns based on curve fitting.
    
    This function reads the GA_matrix CSV file, computes the dipole polarization fits for each row 
    (using get_dipole_polarization, which returns columns 'polarizability' and 'dipole'), merges these 
    values into the original GA_matrix based on the 'site' column, and then overwrites the CSV file.
    
    Args:
        csv_filepath (str): Path to the GA_matrix CSV file (default: './GA_matrix.csv').
    
    Returns:
        updated_df (pandas.DataFrame): The updated GA_matrix DataFrame with the new columns.
        
            
    """
    # Read the current GA_matrix.
    df_GA = pd.read_csv(csv_filepath)
    
    # Get the dipole polarization results.
    # This function must be defined and should return a DataFrame with columns: 
    # 'site', 'polarizability', 'dipole', 'c', 'R2', 'MAE', 'RMSE'
    dipole_results = get_dipole_polarization(df_GA)
    
    # Merge the dipole polarization results into the GA_matrix based on the 'site' column.
    updated_df = df_GA.merge(dipole_results[['site', 'polarizability', 'dipole']], on='site', how='left')
    
    # Overwrite the original GA_matrix CSV file with the updated DataFrame.
    updated_df.to_csv(csv_filepath, index=False)
    print(f"Updated GA_matrix saved to {csv_filepath}.")
    
    return updated_df


def generate_GA_matrix(cluster_path, list_EF=[str(x) for x in [-0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6]]):
    """
    Generate the GA (Group Additivity) matrix and compute adsorption energies.
    
    This function reads a sorted data CSV file, computes a GA matrix based on adsorption sites 
    and group additivity data, and then computes adsorption energies for each specified EF value.
    Rows corresponding to 'slab' or undefined ('NA') sites are dropped from the final matrix.
    
    Parameters:
        cluster_path (str): Path to the cluster data.
        list_EF (list of str): List of EF values (as strings) for which adsorption energies 
                               will be computed.
    
    Returns:
        df_GA (pandas.DataFrame): DataFrame containing the GA matrix and the computed adsorption energies.
    """
    # Ensure results directory exists
    if not os.path.exists('results'):
        os.makedirs('results')
    
    data_file = './sorted_data.csv'
    df = pd.read_csv(data_file)
    list_path = df['path'].tolist()
    # Set the index to 'Species' for later energy retrieval
    df.set_index('Species', inplace=True)
    
    # Load or generate GA data
    try: 
        GA_dict, groups = load_ga_data(cluster_path)
    except Exception as e:
        print("Error loading GA data:", e)
        get_CN_GA(cluster_path, mult=0.9)
        GA_dict, groups = load_ga_data(cluster_path)
    
    # Prepare header for the GA matrix DataFrame
    header = ['site'] + groups
    GA_matrix = []

    # Loop over each path to build the GA matrix
    for path in list_path: 
        poscar_file = os.path.join(path, 'POSCAR')
        try:
            atoms = read(poscar_file)
        except Exception as e:
            print(f"Error reading POSCAR in path {path}: {e}")
            continue
        
        full_path = os.path.abspath(path)        
        try:
            # Determine the adsorption site type based on the file path
            if 'N2_ads' in full_path:
                ru_site_str = classify_n2_adsorption(atoms)
                print(path, ru_site_str)
            elif any(keyword in full_path for keyword in ['N_ads', 'NH_ads', 'NH2_ads', 'NH3_ads']):    
                ru_site_str = classify_single_atom_adsorption(atoms, 'N', 2.4)
            elif 'H_ads' in full_path:
                ru_site_str = classify_single_atom_adsorption(atoms, 'H', 2.2)
            else:
                ru_site_str = 'NA'
            # Get GA groups associated with the adsorption site; default to empty list if not found.
            groups_geo = GA_dict.get(ru_site_str, [])
        except Exception as e:
            print(f"{path}: Error classifying adsorption site for {full_path}: {e}")
            ru_site_str = 'NA'
            groups_geo = ['NA']
    
        # Create the GA matrix row: first column is the site label, then counts for each group
        row = [ru_site_str]
        for group in groups:
            count = groups_geo.count(group)
            row.append(count)
        GA_matrix.append(row)
    
    # Create a DataFrame from the GA matrix
    df_GA = pd.DataFrame(GA_matrix, columns=header)
    
    # Reset the original dataframe index to retrieve energy data (Species becomes a column)
    df = df.reset_index()
    # Retrieve the species type (e.g., 'H', 'N', 'O', etc.) from the current working directory or elsewhere
    species = get_species_from_cwd()  

    # For each EF value, compute the adsorption energy and add it as a column
    for ef_val in list_EF: 
        col_name = 'Eads_' + ef_val
        # Retrieve the slab energy from df where Species equals 'slab'
        slab_rows = df[df['Species'] == 'slab']
        if slab_rows.empty:
            raise KeyError("'slab' entry is missing in the Species column.")
        try:
            slab_energy = slab_rows.iloc[0][ef_val]
        except KeyError:
            raise KeyError(f"EF value '{ef_val}' not found in the data.")
    
        # Compute the adsorption energy: subtract slab energy and a gas-phase correction
        if species in ['H', 'N', 'O']:
            gas_species = species + '2'
            correction = float(gas_dict[gas_species]) / 2
        else:
            gas_species = species
            correction = float(gas_dict[gas_species])
            
        try:
            df[col_name] = df[ef_val] - slab_energy - correction
        except KeyError:
            raise KeyError(f"Column '{ef_val}' not found in the data.")
        
        # Add the computed energy column to the GA matrix DataFrame.
        # (Assumes the row order of df and df_GA is the same.)
        df_GA[col_name] = df[col_name]
        
    # Remove rows where the site is 'slab' or 'NA'
    df_GA = df_GA[(df_GA['site'] != 'slab') & (df_GA['site'] != 'NA')]
    
    # Save the GA matrix to CSV
    output_csv = "GA_matrix.csv"
    
    df_GA.to_csv(output_csv, index=False)
    print("GA matrix saved as", output_csv)
    
    return df_GA


def get_GA_matrix(cluster_path, EF):
    """
    Retrieve the GA matrix from the CSV file if it exists; otherwise, generate it.
    Then extract the feature matrix X and target vector Y.

    Args:
        EF (str): The EF value (as a string) used to label the adsorption energy column.
                  For example, if EF='0.6', the energy column is assumed to be 'E_ads_0.6'.

    Returns:
        tuple: (df_matrix, X, Y, n_group)
            - df_matrix: The full GA matrix DataFrame (including 'site' and adsorption energy columns).
            - X: The features as a NumPy array (with 'site' and all adsorption energy columns removed).
            - Y: The target adsorption energies (from the corresponding column).
            - n_group: The number of feature columns in X.
    """
    csv_filepath = './GA_matrix.csv'
    
    # Check if the GA_matrix.csv exists
    if os.path.exists(csv_filepath):
        print("GA_matrix.csv found. Reading file...")
        df_matrix = pd.read_csv(csv_filepath)
    else:
        print("GA_matrix.csv not found. Generating GA matrix...")
        # Assumes that generate_GA_matrix is defined and cluster_path is valid.        
        df_matrix = generate_GA_matrix(cluster_path)
    
    # Determine the adsorption energy prefix.
    if any(col.startswith('Eads_') for col in df_matrix.columns):
        prefix = 'Eads_'
    else:
        raise KeyError("No adsorption energy columns found in GA_matrix.")
    
    # Determine the target adsorption energy column name based on EF.
    # If EF is empty, we remove the trailing underscore from the prefix.
    col_name = prefix[:-1] if EF == "" else prefix + EF
    
    # Extract Y values (adsorption energies)
    try:
        Y = df_matrix[col_name].values
    except KeyError:
        raise KeyError(f"The expected adsorption energy column '{col_name}' is not found in the GA matrix.")
    
    # Drop the 'site' column and all adsorption energy columns (those starting with the prefix) to get X.
    cols_to_drop = ['site'] + [col for col in df_matrix.columns if col.startswith(prefix)]
    df_features = df_matrix.drop(columns=cols_to_drop)
    X = df_features.values
    n_group = df_features.shape[1]
    
    return df_matrix, X, Y, n_group
    
def get_species_from_cwd():
    # cwd = os.getcwd()
    cwd = os.path.realpath(os.getcwd())
    # print(cwd, 'cwd')
    species = None
    try:
        # path_parts = cwd.replace(' ', '').split(os.path.sep)
        part = os.path.basename(cwd).strip()
        # print('part', part)
        if '_ads' in part:
            species = part.split('_')[0]
        
    except Exception as e:
        print(f"An error occurred: {e}")
    
    if species:
        return species
    else:
        print("No species found in the current working directory.")

def GA_prediction(X_training, Y_training):
    '''Run conventional GA approach '''
 
    model = PseudoInverseLinearRegression() 
    model.fit(X_training,Y_training)

    Y_pred = model.predict(X_training)
    
    # Calculate the MAE, RMSE, and R²
    mae = mean_absolute_error(Y_training, Y_pred)
    rmse = np.sqrt(mean_squared_error(Y_training, Y_pred))
    r2 = r2_score(Y_training, Y_pred)
    
    # Print the errors and R²
    print(f'Mean Absolute Error (MAE): {mae}')
    print(f'Root Mean Squared Error (RMSE): {rmse}')
    print(f'R-squared (R²): {r2}')
         
    # # Plot the data and regression line
    plt.figure(figsize=(10, 6))
    plt.scatter(Y_training, Y_pred, color='blue', label='Model')
    
    # Add text box for metrics
    textstr = f'MAE: {mae:.2f}\nRMSE: {rmse:.2f}\nR2: {r2:.2f}'
    plt.gca().text(0.65, 0.25, textstr, transform=plt.gca().transAxes,
                    fontsize=16, verticalalignment='top', bbox=dict(facecolor='white', alpha=0))
    
    # Customize the plot
    plt.xlabel('X / Unit', fontsize=16)
    plt.ylabel('Y / Unit', fontsize = 16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(loc='upper left', frameon=False, fontsize=18)
    #plt.title('Linear Regression of De vs Ea')
    plt.tight_layout()
    # Save the figure
    plt.savefig('XY.png')
    return Y_training, Y_pred

    
def save_one_model_split(species, df, EF,  X, Y, models, model_name, ratio=0.95):
    model = models[model_name]
    n_samples = int(len(X) * ratio)  # 80% of the data for training
    indices = np.random.choice(len(X), size=len(X), replace=True)
    X_train, Y_train = X[indices[:n_samples]], Y[indices[:n_samples]]   
    model.fit(X_train,Y_train)
    
    joblib.dump(model, f'{species}_{model_name}_{EF}.pkl')


def test_one_model_split(df, EF,  X, Y, models, model_name, ratio):
    '''Run conventional GA approach '''
    model = models[model_name]
    
    n_samples = int(len(X) * ratio)  # 80% of the data for training
    
    indices = np.random.choice(len(X), size=len(X), replace=True)
    X_train, Y_train = X[indices[:n_samples]], Y[indices[:n_samples]]
    X_test, Y_test = X[indices[n_samples:]], Y[indices[n_samples:]]
    
    model.fit(X_train,Y_train)
    Y_pred = model.predict(X_test)
    
    # Calculate the MAE, RMSE, and R²
    mae = mean_absolute_error(Y_test, Y_pred)
    rmse = np.sqrt(mean_squared_error(Y_test, Y_pred))
    r2 = r2_score(Y_test, Y_pred)

    fig, ax  = plt.subplots(figsize=(10,8))  
    
    species = get_species_from_cwd()
    text = '%s adsorption' %(species)
    # ax.set_title(text)
    ax.scatter(Y_test, Y_pred, edgecolor='k', alpha=0.7, s=200, label = text, color='#1f77b4')
    ax.plot([min(Y), max(Y)], [min(Y), max(Y)], 'r--', linewidth=3)
    ax.set_xlabel(f'DFT Eads (eV) at EF = {EF} eV/A', fontsize=36)
    ax.set_ylabel('Pred. Eads(eV)', fontsize=36)

    textstr = f'MAE: {mae:.2f}  \nRMSE: {rmse:.2f}  \nR2: {r2:.2f}'
    plt.gca().text(0.65, 0.25, textstr, transform=plt.gca().transAxes,
                    fontsize=20, verticalalignment='top', bbox=dict(facecolor='white', alpha=0))


    for i in range(len(Y_test)):
        if abs(Y_test[i] - Y_pred[i]) > 0.15: # 0.15 eV difference
            ax.annotate(df['site'][i], (Y_test[i], Y_pred[i]), textcoords="offset points", xytext=(0,10), ha='center')
            ax.scatter(Y_test[i], Y_pred[i], color='red', s=100, edgecolor='k')
    
    
    plt.legend(loc='upper left', frameon=False, fontsize=40)

    plt.grid(False)
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(3)
        
    #plt.title('Linear Regression of De vs Ea')
    ax.tick_params(axis='both', which='major', width=2, length=8, labelsize=28)
    plt.tight_layout()
    plt.savefig(f'results/XY_%s_testing_%s_{EF}.png' %(model_name, str(1-ratio)), dpi=600)
    
    return mae, rmse, r2 
        
 
class PseudoInverseLinearRegression:
    def __init__(self, fit_intercept=True):
        self.beta_ = None
        self.fit_intercept = fit_intercept

    def fit(self, X, Y):
        if self.fit_intercept:
            # Add a column of ones to X to account for the intercept
            X = np.hstack([np.ones((X.shape[0], 1)), X])
        
        # Compute the beta coefficients using the Moore-Penrose pseudoinverse
        self.beta_ = np.linalg.pinv(X).dot(Y)

    def predict(self, X):
        if self.fit_intercept:
            # Add a column of ones to X for predictions
            X = np.hstack([np.ones((X.shape[0], 1)), X])
        
        # Make predictions using the computed beta coefficients
        return X.dot(self.beta_)
    

# # Now include this custom model in the models dictionary
def refine_models(n_group):
    models = {
        "Ridge": Ridge(alpha=1.0),
        "Lasso": Lasso(alpha=1/n_group),
        "EN": ElasticNet(alpha=0.1, l1_ratio=0.5),
        "DT": DecisionTreeRegressor(max_depth=n_group),
        "RF": RandomForestRegressor(n_estimators=n_group*5, max_depth=n_group*5),
        "GB": GradientBoostingRegressor(n_estimators=n_group*5, max_depth=n_group),
        "SVR": SVR(kernel='linear', C=1.0),
        "GA": PseudoInverseLinearRegression(),
        "XGB": XGBRegressor(n_estimators=11*5, max_depth=11, learning_rate=0.1),
        }
    return models


def bootstrap_cross_validation_one_model(models, model_name, EF, X, Y, ratio, n_bootstrap=100):
    """
    Perform bootstrapped cross-validation for a single model and save results.
    
    Args:
        models (dict): Dictionary of models with their names as keys.
        model_name (str): The name of the model to use for cross-validation.
        X (np.array): Features dataset.
        Y (np.array): Target values.
        ratio (float): Proportion of data to use for training (e.g., 0.8 for 80%).
        n_bootstrap (int): Number of bootstrap iterations.
    
    Returns:
        dict: Dictionary containing metric results for the model.
    """
    model = models[model_name]
    results = {metric: [] for metric in ['MAE', 'MSE', 'RMSE', 'R2', 'MPD', 'MND']}
    n_samples = int(len(X) * ratio)  # Number of training samples

    for i in range(n_bootstrap):
        print(f'Iteration: {i + 1}/{n_bootstrap}')
        
        # Generate bootstrap indices
        indices = np.random.choice(len(X), size=len(X), replace=True)
        X_train, Y_train = X[indices[:n_samples]], Y[indices[:n_samples]]
        X_test, Y_test = X[indices[n_samples:]], Y[indices[n_samples:]]
        
        # Train the model and predict
        model.fit(X_train, Y_train)
        predictions = model.predict(X_test)
        
        # Calculate metrics
        results['MAE'].append(mean_absolute_error(Y_test, predictions))
        results['MSE'].append(mean_squared_error(Y_test, predictions))
        results['RMSE'].append(np.sqrt(mean_squared_error(Y_test, predictions)))
        results['R2'].append(r2_score(Y_test, predictions))
        results['MPD'].append(np.max(predictions - Y_test))
        results['MND'].append(np.min(predictions - Y_test))
    
    # Save results to CSV files
    results_dir = 'results'
    os.makedirs(results_dir, exist_ok=True)  # Ensure the results directory exists
    
    for metric, values in results.items():
        metric_file = f'{results_dir}/{metric.lower()}_{ratio}_{EF}.csv'
        averages_file = f'{results_dir}/{metric.lower()}_averages_{ratio}_{EF}.csv'
        
        # Save all values
        pd.DataFrame({model_name: values}).to_csv(metric_file, index=False)
        print(f"Saved all {metric} values to {metric_file}_{EF}")
        
        # Save averages
        pd.DataFrame({model_name: [np.mean(values)]}).to_csv(averages_file, index=False)
        print(f"Saved average {metric} value to {averages_file}_{EF}")

    return results

def active_learning_one_model(models, model_name, EF, X, Y, sites, 
                              initial_ratio=0.4, increment=0.05, max_ratio=0.95,
                              error_threshold=0.5, results_dir='results'):
    """
    Perform an active learning process for a single model.
    
    Initially, a fraction of the data (default 40%) is used for training and the rest for testing.
    In each iteration, the model is trained on the current training set and evaluated on the testing set.
    Then, the top outliers (largest absolute errors) are selected from the testing set and added 
    to the training set until the training set reaches the maximum ratio (default 95%).
    
    For each iteration, the following metrics are computed:
      - Mean Absolute Error (MAE)
      - Mean Squared Error (MSE)
      - Root Mean Squared Error (RMSE)
      - R² score (R2)
      - Maximum Positive Deviation (MPD): max(predictions - Y_test)
      - Maximum Negative Deviation (MND): min(predictions - Y_test)
    
    Additionally, the function prints out the 'site' labels for test samples with an absolute 
    prediction error larger than error_threshold (default 0.5 eV).
    
    Args:
        models (dict): Dictionary of models.
        model_name (str): Name of the model to use.
        EF: An identifier for saving results.
        X (np.array): Features dataset.
        Y (np.array): Target values.
        sites (np.array or list): Labels corresponding to each sample (e.g., the 'site' column).
        initial_ratio (float): Initial fraction of data for training (default 0.4).
        increment (float): Fraction of total data to add each iteration (default 0.05).
        max_ratio (float): Maximum fraction of data for training (default 0.95).
        error_threshold (float): Threshold in eV to flag outliers (default 0.5 eV).
        results_dir (str): Directory where results will be saved.
    
    Returns:
        list of dict: A list of dictionaries containing performance metrics for each iteration.
    """
    total_samples = len(X)
    indices = np.arange(total_samples)
    np.random.shuffle(indices)  # Randomize the data order

    # Initial split: use first initial_ratio of samples for training.
    n_train = int(total_samples * initial_ratio)
    train_indices = indices[:n_train].tolist()
    test_indices = indices[n_train:].tolist()

    results = []
    current_ratio = len(train_indices) / total_samples

    # Ensure the results directory exists.
    os.makedirs(results_dir, exist_ok=True)
    iteration = 0

    while current_ratio < max_ratio and len(test_indices) > 0:
        iteration += 1

        # Build training and testing sets.
        X_train = X[train_indices]
        Y_train = Y[train_indices]
        X_test = X[test_indices]
        Y_test = Y[test_indices]

        # Train the model on the current training set.
        model = models[model_name]
        model.fit(X_train, Y_train)
        predictions = model.predict(X_test)
        
        
        # Save the model for future prediction 
        species  = get_species_from_cwd() 
        joblib.dump(model, f'{species}_{model_name}_{EF}.pkl')

        # Calculate performance metrics.
        mae = mean_absolute_error(Y_test, predictions)
        mse = mean_squared_error(Y_test, predictions)
        rmse = np.sqrt(mse)
        r2 = r2_score(Y_test, predictions)

        # Compute deviations.
        deviations = predictions - Y_test
        mpd = np.max(deviations)  # Maximum positive deviation.
        mnd = np.min(deviations)  # Maximum negative deviation.

        # Record the results for this iteration.
        iter_result = {
            'iteration': iteration,
            'training_ratio': current_ratio,
            'n_train': len(train_indices),
            'n_test': len(test_indices),
            'MAE': mae,
            'MSE': mse,
            'RMSE': rmse,
            'R2': r2,
            'MPD': mpd,
            'MND': mnd
        }
        results.append(iter_result)

        # Print out 'site' labels for test samples with error > error_threshold.
        errors = np.abs(Y_test - predictions)
        test_sites = np.array(sites)[test_indices]
        outlier_mask = errors > error_threshold
        if np.any(outlier_mask):
            outlier_indices = np.where(outlier_mask)[0]
            sorted_indices = outlier_indices[np.argsort(-errors[outlier_indices])]
            print(f"Iteration {iteration}: Outlier site labels with error > {error_threshold} eV:")
            for idx in sorted_indices:
                print(f"  Site: {test_sites[idx]}, Error: {errors[idx]:.3f} eV")
        else:
            print(f"Iteration {iteration}: No outliers with error > {error_threshold} eV.")

        # Determine the number of samples to add.
        n_to_add = int(total_samples * increment)
        # Safeguard: ensure at least one sample is added per iteration.
        if n_to_add < 1:
            n_to_add = 1
        if n_to_add > len(test_indices):
            n_to_add = len(test_indices)

        # Select the top n_to_add samples (by error ranking).
        sorted_error_indices = np.argsort(-errors)  # Largest errors first.
        selected_indices = [test_indices[i] for i in sorted_error_indices[:n_to_add]]

        # Move selected indices from testing to training.
        train_indices.extend(selected_indices)
        test_indices = [idx for idx in test_indices if idx not in selected_indices]

        # Update current training ratio.
        current_ratio = len(train_indices) / total_samples

    # Save overall results to a single CSV file.
    overall_csv_file = os.path.join(results_dir, f'active_learning_overall_results_{EF}.csv')
    results_df = pd.DataFrame(results)
    results_df.to_csv(overall_csv_file, index=False)
    print(f"Saved overall active learning results to {overall_csv_file}")
    
    return results


def plot_active_learning_from_csv(csv_filepath,EF):
    """
    Read the active learning results from a CSV file and plot performance metrics vs. training ratio.
    
    The CSV file is expected to have columns such as:
    'training_ratio', 'MAE', 'MSE', 'RMSE', 'R2', 'MPD', 'MND'
    
    Args:
        csv_filepath (str): The path to the CSV file containing the active learning results.
    """
    df = pd.read_csv(csv_filepath)
    df.sort_values('training_ratio', inplace=True)
    plt.figure(figsize=(10, 6))
    plt.plot(df['training_ratio'], df['MAE'], marker='o', label='MAE')
    # plt.plot(df['training_ratio'], df['MSE'], marker='s', label='MSE')
    plt.plot(df['training_ratio'], df['RMSE'], marker='^', label='RMSE')
    # plt.plot(df['training_ratio'], df['R2'], marker='v', label='R2')
    plt.plot(df['training_ratio'], df['MPD'], marker='D', label='MPD')
    plt.plot(df['training_ratio'], df['MND'], marker='x', label='MND')
    plt.tick_params(axis='both', labelsize=18)
    plt.xlabel('Training Ratio', fontsize = 20)
    plt.ylabel('Metric Value', fontsize = 20)
    plt.title(f'Active Learning Performance Metrics vs. Training Ratio at {EF} eV/A')
    plt.legend(framealpha=0, fontsize=18)
    plt.tight_layout()
    plt.grid(True)
    plt.savefig(f'results/mae_{EF}.png', dpi=300)
    plt.close()
    
def run_or_plot_active_learning(models, model_name, EF, X, Y, sites,
                                initial_ratio=0.4, increment=0.05, max_ratio=0.95,
                                results_dir='results'):
    """
    Check for an existing overall CSV file with active learning results.
    If found, read and plot the results. Otherwise, run the active learning process,
    save the merged results to one CSV file, and then plot the results.
    
    Args:
        models (dict): Dictionary of models.
        model_name (str): Name of the model to use.
        EF: An identifier for saving results.
        X (np.array): Features dataset.
        Y (np.array): Target values.
        sites (np.array or list): Labels corresponding to each sample (e.g., the 'site' column).
        initial_ratio (float): Initial fraction of data for training.
        increment (float): Fraction of total data to add each iteration.
        max_ratio (float): Maximum fraction of data for training.
        results_dir (str): Directory where results will be saved.
    """
    overall_csv_file = os.path.join(results_dir, f'active_learning_overall_results_{EF}.csv')
    
    if os.path.exists(overall_csv_file):
        print(f"CSV file found at {overall_csv_file}. Reading and plotting results...")
        plot_active_learning_from_csv(overall_csv_file, EF)
    else:
        print("CSV file not found, running active learning process...")
        active_learning_one_model(models, model_name, EF, X, Y, sites,
                                  initial_ratio=initial_ratio, increment=increment, 
                                  max_ratio=max_ratio, results_dir=results_dir)
        print("Active learning process completed. Plotting results...")
        plot_active_learning_from_csv(overall_csv_file, EF)



### GA Applilcations
    
def find_all_rhombuses(atoms, connections, surface_indices, bond_length_threshold):
    # atoms, connections, cn_of_connected_atoms, exposed_top_sites = get_connection(path, metal='Ru', mult=0.9)
    """Find all possible rhombuses based on atom connections"""
    rhombuses = []
    for A_index in surface_indices:
        A_neighbors = connections[A_index]
        for B_index in A_neighbors:
            if B_index not in surface_indices or A_index == B_index:
                continue
            if np.linalg.norm(atoms[A_index].position - atoms[B_index].position) > bond_length_threshold:
                continue

            B_neighbors = connections[B_index]
            for C_index in B_neighbors:
                if C_index not in surface_indices or C_index in [A_index, B_index]:
                    continue
                if np.linalg.norm(atoms[B_index].position - atoms[C_index].position) > bond_length_threshold:
                    continue

                C_neighbors = connections[C_index]
                for D_index in C_neighbors:
                    if D_index not in surface_indices or D_index in [A_index, B_index, C_index]:
                        continue
                    if np.linalg.norm(atoms[C_index].position - atoms[D_index].position) > bond_length_threshold:
                        continue
                    if np.linalg.norm(atoms[D_index].position - atoms[A_index].position) > bond_length_threshold:
                        continue

                    # Save rhombus indices while maintaining the original order
                    rhombus_indices = [A_index, B_index, C_index, D_index]
                    rhombuses.append(rhombus_indices)

    return rhombuses


def get_matrix_to_be_predicted(cluster_path, site):
    # Load or generate GA data
    try: 
        GA_dict, groups = load_ga_data(cluster_path)
    except Exception as e:
        print("Error loading GA data:", e)
        get_CN_GA(cluster_path, mult=0.9)
        GA_dict, groups = load_ga_data(cluster_path)
        
    """ Get the group matrix of one site"""
    try:
        groups_site = GA_dict[site]
    except:
        site = str(site)
        groups_site = GA_dict[site]
    matrix_site = []
    for group in groups:
        num_group = 0
        for group_site in groups_site:
            if group_site == group:
                num_group += 1
        matrix_site.append(num_group)
    matrix_site = np.array([matrix_site])
    return matrix_site

def predict_Eads_site(cluster_path, species, site, Prop):
    '''Print the adsorption energy at the specific site, 
    the data in the data_path is applied to train the GA model to predict the 
    same species at the provided site.
    Prop values: -0.6, -0.4, -0.2, 0.0, 0.2, 0.4, 0.6, polarizability, dipole, all are strings.
    '''
    GA_model = joblib.load(f'{species}_GA_{Prop}.pkl')  
    preidct_matrix =  get_matrix_to_be_predicted(cluster_path, site)
    ### predict the energies of species at the site    
    E_species = GA_model.predict(preidct_matrix)

    return E_species


# ------------------ Site Conversion ------------------
def convert_sites(sites_list):
    """
    Given a list of four site numbers corresponding to atoms A, B, C, and D (in clockwise order),
    return a dictionary mapping standard site labels to their numeric string representations.
    
    For example, if sites_list = [20, 23, 24, 42]:
      - "top_A": "20"
      - "top_B": "23"
      - "top_C": "24"
      - "top_D": "42"
      - "bridge_A-B": "20-23"
      - "bridge_B-C": "23-24"
      - "bridge_C-D": "24-42"
      - "bridge_D-A": "20-42"  (always smaller-first)
      - "hollow_ABD": "20-23-42"  (sorted ascending among A, B, D)
      - "hollow_BCD": "23-24-42"  (sorted ascending among B, C, D)
    """
    if len(sites_list) != 4:
        raise ValueError("sites_list must contain exactly four elements corresponding to A, B, C, and D.")
    A, B, C, D = sites_list
    site_dict = {}
    site_dict["top_A"] = str(A)
    site_dict["top_B"] = str(B)
    site_dict["top_C"] = str(C)
    site_dict["top_D"] = str(D)
    site_dict["bridge_A-B"] = f"{min(A, B)}-{max(A, B)}"
    site_dict["bridge_B-C"] = f"{min(B, C)}-{max(B, C)}"
    site_dict["bridge_C-D"] = f"{min(C, D)}-{max(C, D)}"
    site_dict["bridge_D-A"] = f"{min(D, A)}-{max(D, A)}"
    site_dict["hollow_ABD"] = "-".join(str(x) for x in sorted([A, B, D]))
    site_dict["hollow_BCD"] = "-".join(str(x) for x in sorted([B, C, D]))
    return site_dict

# ------------------ Modified Helper Function ------------------
def select_best_site(cluster_path, species, sites, Prop, site_mapping=None):
    """
    Evaluate the predicted adsorption energy for a given species at each candidate site.
    If site_mapping is provided, convert the candidate label to its numeric string.
    
    Args:
        cluster_path (str): Path to the cluster data.
        species (str): The adsorbate species (e.g., "NH3", "NH2", "NH", "N", "H", "N2").
        sites (list): List of candidate site labels (e.g., "top_A", "bridge_A-B", etc.).
        Prop: Additional properties required by predict_Eads_site.
        site_mapping (dict, optional): Mapping of candidate labels to numeric strings.
    
    Returns:
        tuple: (best_site, energy_dict) where energy_dict is keyed by the candidate (symbolic label).
    """
    energy_dict = {}
    for site in sites:
        candidate = site_mapping[site] if (site_mapping is not None and site in site_mapping) else site
        Eads = predict_Eads_site(cluster_path, species, candidate, Prop)
        energy_dict[site] = Eads
        print(f"{species} at site {site} ({candidate}): Eads = {Eads:.3f} eV")
    best_site = min(energy_dict, key=energy_dict.get)
    print(f"Best {species} site: {best_site} with Eads = {energy_dict[best_site]:.3f} eV\n")
    return best_site, energy_dict

# ------------------ Candidate Mappings ------------------

# (1) For NH3 adsorption:
nh3_candidates = ["A", "B", "C", "D"]

# (2) For NH3 → NH2 + H:
def get_nh2_candidates(nh3_site):
    mapping = {
        "A": ["top_A", "bridge_A-B", "bridge_D-A", "hollow_ABD"],
        "B": ["top_B", "bridge_A-B", "bridge_B-C", "hollow_ABD", "hollow_BCD"],
        "C": ["top_C", "bridge_B-C", "bridge_C-D", "hollow_BCD"],
        "D": ["top_D", "bridge_C-D", "bridge_D-A", "hollow_ABD", "hollow_BCD"]
    }
    return mapping.get(nh3_site, [])

# For H in Step 2 (H determined by NH2):
nh2_to_h_candidates_step2 = {
    "top_A":     ["top_B", "top_D", "hollow_ABD", "hollow_BCD"],
    "top_B":     ["top_A", "top_C", "hollow_ABD", "hollow_BCD"],
    "top_C":     ["top_B", "top_D", "hollow_ABD", "hollow_BCD"],
    "top_D":     ["top_A", "top_C", "hollow_ABD", "hollow_BCD"],
    "bridge_A-B": ["top_A", "top_B", "bridge_A-B", "hollow_ABD"],
    "bridge_B-C": ["top_B", "top_C", "bridge_B-C", "hollow_BCD"],
    "bridge_C-D": ["top_C", "top_D", "bridge_C-D", "hollow_BCD"],
    "bridge_D-A": ["top_D", "top_A", "bridge_D-A", "hollow_ABD"],
    "hollow_ABD": ["hollow_ABD", "top_A", "top_B", "top_D", "bridge_A-B", "bridge_D-A"],
    "hollow_BCD": ["hollow_BCD", "top_B", "top_C", "top_D", "bridge_B-C", "bridge_C-D"]
}

# (3) For NH2 → NH + H:
nh2_to_nh_candidates = {
    "top_A":     ["top_A", "bridge_A-B", "bridge_D-A", "hollow_ABD"],
    "top_B":     ["top_B", "bridge_A-B", "bridge_B-C", "hollow_ABD", "hollow_BCD"],
    "top_C":     ["top_C", "bridge_B-C", "bridge_C-D", "hollow_BCD"],
    "top_D":     ["top_D", "bridge_C-D", "bridge_D-A", "hollow_ABD", "hollow_BCD"],
    "bridge_A-B": ["top_A", "top_B", "bridge_A-B", "hollow_ABD"],
    "bridge_B-C": ["top_B", "top_C", "bridge_B-C", "bridge_A-B", "hollow_BCD"],
    "bridge_C-D": ["top_C", "top_D", "bridge_C-D", "hollow_BCD"],
    "bridge_D-A": ["top_D", "top_A", "bridge_D-A", "hollow_ABD"],
    "hollow_ABD": ["hollow_ABD", "top_A", "top_B", "top_D", "bridge_A-B", "bridge_D-A"],
    "hollow_BCD": ["hollow_BCD", "top_B", "top_C", "top_D", "bridge_B-C", "bridge_C-D"]
}

# For H in Step 3 (H determined by NH):
nh_to_h_candidates = {
    "top_A":     ["top_B", "top_D", "top_C", "hollow_ABD", "hollow_BCD"],
    "top_B":     ["top_A", "top_C", "top_D", "bridge_A-B", "bridge_B-C", "bridge_A-D", "bridge_C-D", "hollow_ABD", "hollow_BCD"],
    "top_C":     ["top_B", "top_D", "top_A", "bridge_B-C", "bridge_C-D", "bridge_A-B", "bridge_D-A", "hollow_BCD", "hollow_ABD"],
    "top_D":     ["top_A", "top_B", "top_C", "bridge_C-D", "bridge_D-A", "bridge_A-B", "bridge_B-C", "hollow_ABD", "hollow_BCD"],
    "bridge_A-B": ["top_C", "top_D", "bridge_B-C", "hollow_ABD", "hollow_BCD"],
    "bridge_B-C": ["top_A", "top_D", "bridge_D-A", "bridge_A-B", "bridge_C-D", "hollow_ABD", "hollow_BCD"],
    "bridge_C-D": ["top_A", "top_B", "bridge_A-B", "bridge_B-C", "bridge_D-A", "hollow_ABD", "hollow_BCD"],
    "bridge_D-A": ["top_B", "top_C", "bridge_B-C", "bridge_A-B", "bridge_C-D", "hollow_ABD", "hollow_BCD"],
    "hollow_ABD": ["hollow_BCD", "top_A", "top_B", "top_C", "top_D"],
    "hollow_BCD": ["hollow_ABD", "top_A", "top_B", "top_C", "top_D"]
}

# (4) For NH → N + H:
nh_to_n_candidates = {
    "top_A":     ["top_A", "bridge_A-B", "bridge_D-A", "hollow_ABD"],
    "top_B":     ["top_B", "bridge_A-B", "bridge_B-C", "hollow_ABD", "hollow_BCD"],
    "top_C":     ["top_C", "bridge_B-C", "bridge_C-D", "hollow_BCD"],
    "top_D":     ["top_D", "bridge_C-D", "bridge_D-A", "hollow_ABD", "hollow_BCD"],
    "bridge_A-B": ["top_A", "top_B", "bridge_A-B", "hollow_ABD"],
    "bridge_B-C": ["top_B", "top_C", "bridge_B-C", "bridge_A-B", "hollow_BCD"],
    "bridge_C-D": ["top_C", "top_D", "bridge_C-D", "hollow_BCD"],
    "bridge_D-A": ["top_D", "top_A", "bridge_D-A", "hollow_ABD"],
    "hollow_ABD": ["hollow_ABD", "top_A", "top_B", "top_D", "bridge_A-B", "bridge_D-A"],
    "hollow_BCD": ["hollow_BCD", "top_B", "top_C", "top_D", "bridge_B-C", "bridge_C-D"]
}
# For H in Step 4, reuse nh_to_h_candidates.

# ------------------ Reaction Chain Functions ------------------

# Step 1: NH₃ Adsorption
def determine_NH3_configuration(cluster_path, Prop, site_mapping=None):
    best_NH3_site, energies = select_best_site(cluster_path, "NH3", nh3_candidates, Prop, site_mapping)
    return {"site": best_NH3_site, "energy": energies[best_NH3_site]}

# Step 2: NH₃ → NH₂ + H (H determined by NH2)
def determine_NH2_configuration(cluster_path, Prop, best_NH3_site, site_mapping=None):
    nh2_candidates = get_nh2_candidates(best_NH3_site)
    best_NH2_site, energies_NH2 = select_best_site(cluster_path, "NH2", nh2_candidates, Prop, site_mapping)
    candidate_H_sites = nh2_to_h_candidates_step2.get(best_NH2_site, [])
    if candidate_H_sites:
        best_H_site, energies_H = select_best_site(cluster_path, "H", candidate_H_sites, Prop, site_mapping)
    else:
        best_H_site, energies_H = None, {}
    return {
        "NH2": {"site": best_NH2_site, "energy": energies_NH2[best_NH2_site]},
        "H1": {"site": best_H_site, "energy": energies_H.get(best_H_site, None)}
    }

# Step 3: NH₂ → NH + H (H determined by NH)
def determine_NH_configuration(cluster_path, Prop, best_NH2_site, site_mapping=None):
    candidate_NH_sites = nh2_to_nh_candidates.get(best_NH2_site, [])
    best_NH_site, energies_NH = select_best_site(cluster_path, "NH", candidate_NH_sites, Prop, site_mapping)
    candidate_H_sites = nh_to_h_candidates.get(best_NH_site, [])
    if candidate_H_sites:
        best_H_site, energies_H = select_best_site(cluster_path, "H", candidate_H_sites, Prop, site_mapping)
    else:
        best_H_site, energies_H = None, {}
    return {
        "NH": {"site": best_NH_site, "energy": energies_NH[best_NH_site]},
        "H2": {"site": best_H_site, "energy": energies_H.get(best_H_site, None)}
    }

# Step 4: NH → N + H (H determined by N)
def determine_N_configuration(cluster_path, Prop, best_NH_site, site_mapping=None):
    candidate_N_sites = nh_to_n_candidates.get(best_NH_site, [])
    if candidate_N_sites:
        best_N_site, energies_N = select_best_site(cluster_path, "N", candidate_N_sites, Prop, site_mapping)
    else:
        best_N_site, energies_N = None, {}
    candidate_H_sites = nh_to_h_candidates.get(best_N_site, [])
    if candidate_H_sites:
        best_H_site, energies_H = select_best_site(cluster_path, "H", candidate_H_sites, Prop, site_mapping)
    else:
        best_H_site, energies_H = None, {}
    return {
        "N": {"site": best_N_site, "energy": energies_N.get(best_N_site, None)},
        "H3": {"site": best_H_site, "energy": energies_H.get(best_H_site, None)}
    }

# N₂ Adsorption
def determine_N2_configuration(cluster_path, Prop, site_mapping=None):
    candidate_sites = [
        "top_A", "top_B", "top_C", "top_D",
        "bridge_A-B", "bridge_B-C", "bridge_C-D", "bridge_D-A",
        "hollow_ABD", "hollow_BCD"
    ]
    best_N2_site, energies_N2 = select_best_site(cluster_path, "N2", candidate_sites, Prop, site_mapping)
    return {"site": best_N2_site, "energy": energies_N2[best_N2_site]}


# ------------------ Top-Level Full Configuration ------------------
def determine_full_configuration(cluster_path, Prop, sites_list=None):
    """
    Determine the full stable configuration for the reaction chain and return a dictionary
    where each species is mapped to a tuple (site, energy):
      1. NH3 → NH2 + H      (H determined by NH2)
      2. NH2 → NH + H       (H determined by NH)
      3. NH  → N + H        (H determined by N)
      4. Also determine the N2 adsorption site.
    
    If sites_list is provided (a list of four numbers for atoms A, B, C, D in clockwise order),
    it is converted to a mapping for energy predictions.
    
    Returns a dictionary with keys:
      "NH3", "NH2", "NH", "N", "N2", "H1", "H2", and "H3".
    """
    site_mapping = convert_sites(sites_list) if sites_list is not None else None

    config_NH3 = determine_NH3_configuration(cluster_path, Prop, site_mapping)
    config_NH2 = determine_NH2_configuration(cluster_path, Prop, config_NH3["site"], site_mapping)
    config_NH = determine_NH_configuration(cluster_path, Prop, config_NH2["NH2"]["site"], site_mapping)
    config_N = determine_N_configuration(cluster_path, Prop, config_NH["NH"]["site"], site_mapping)
    config_N2 = determine_N2_configuration(cluster_path, Prop, site_mapping)

    final_config = {
        "NH3": (config_NH3["site"], config_NH3["energy"]),
        "NH2": (config_NH2["NH2"]["site"], config_NH2["NH2"]["energy"]),
        "NH":  (config_NH["NH"]["site"], config_NH["NH"]["energy"]),
        "N":   (config_N["N"]["site"], config_N["N"]["energy"]),
        "N2":  (config_N2["site"], config_N2["energy"]),
        "H1":  (config_NH2["H1"]["site"], config_NH2["H1"]["energy"]),
        "H2":  (config_NH["H2"]["site"], config_NH["H2"]["energy"]),
        "H3":  (config_N["H3"]["site"], config_N["H3"]["energy"])
    }
    return final_config


def check_coplanar(points, threshold=0.1):
    """Check if four points are coplanar, unit is Å"""
    p1, p2, p3, p4 = points
    # Vectors
    v1 = p2 - p1
    v2 = p3 - p1
    v3 = p4 - p1
    # Normal vector
    normal = np.cross(v1, v2)
    normal /= np.linalg.norm(normal)
    # Calculate distance of p4 to the plane
    distance = np.abs(np.dot(v3, normal))
    return distance < threshold

def filter_coplanar_rhombuses(rhombuses, atoms, coplanar_threshold=0.5):
    """Filter coplanar rhombuses"""
    coplanar_rhombuses = []
    for indices in rhombuses:
        points = np.array([atoms[idx].position for idx in indices])
        if check_coplanar(points, threshold=coplanar_threshold):
            coplanar_rhombuses.append(indices)
    return coplanar_rhombuses

def filter_rhombuses_by_dihedral(rhombuses, atoms, dihedral_min, dihedral_max):
    def calculate_dihedral(p1, p2, p3, p4):
        """Calculate the dihedral angle of four points"""
        b1 = p2 - p1
        b2 = p3 - p2
        b3 = p4 - p3

        # Normal vectors
        n1 = np.cross(b1, b2)
        n1 /= np.linalg.norm(n1)
        n2 = np.cross(b2, b3)
        n2 /= np.linalg.norm(n2)

        # Dihedral angle
        m1 = np.cross(n1, b2 / np.linalg.norm(b2))
        x = np.dot(n1, n2)
        y = np.dot(m1, n2)
        angle = np.arctan2(y, x)
        return np.degrees(angle)
    
    """Filter rhombuses by dihedral angle"""
    filtered_rhombuses = []
    for indices in rhombuses:
        points = np.array([atoms[idx].position for idx in  indices])
        dihedral_angle = calculate_dihedral(points[0], points[1], points[2], points[3])
        if dihedral_angle < 0:
            dihedral_angle += 180
        # print(dihedral_angle)            
        if dihedral_min <= dihedral_angle <= dihedral_max:
            filtered_rhombuses.append(indices)
            # print(f"Rhombus {indices} with dihedral angle: {dihedral_angle:.2f} degrees")
    return filtered_rhombuses

def filter_sites(all_rhombuses, atoms):
    neat_sites = []
    for indices in all_rhombuses:
        if len(indices) != 4:
            raise ValueError("List must contain exactly 4 elements")
        
        # Get positions using indices
        positions = atoms.get_positions()
        
        max_distance = 0
        max_pair = (None, None)
        for i in range(len(indices)):
            for j in range(i + 1, len(indices)):
                distance = np.linalg.norm(positions[indices[i]] - positions[indices[j]])
                if distance > max_distance:
                    max_distance = distance
                    max_pair = (indices[i], indices[j])
                    
        # Find the pair with the longest distance
        first, fourth = max_pair
        
        # Get the remaining two indices
        remaining = [index for index in indices if index not in [first, fourth]]
        second, third = sorted(remaining)
        
        # Combine the sorted parts
        sorted_indices = [min(first, fourth), second, third, max(first, fourth)]        
        
        if sorted_indices not in neat_sites:
            neat_sites.append(sorted_indices)
    return neat_sites    

def point_in_prism(point, base_vertices, top_vertices):
    """Check if a point is within the prism defined by base and top triangles"""
    def is_point_in_triangle(p, v1, v2, v3):
        d1 = np.dot(np.cross(v2 - v1, p - v1), np.cross(v2 - v1, v3 - v1))
        d2 = np.dot(np.cross(v3 - v2, p - v2), np.cross(v3 - v2, v1 - v2))
        d3 = np.dot(np.cross(v1 - v3, p - v3), np.cross(v1 - v3, v2 - v3))
        return d1 >= 0 and d2 >= 0 and d3 >= 0
    
    # Check if the point is between the two parallel triangles
    prism_height = np.linalg.norm(top_vertices[0] - base_vertices[0])
    normal = calculate_normal(base_vertices[0], base_vertices[1], base_vertices[2])
    distance_to_base = np.dot(point - base_vertices[0], normal)
    
    if abs(distance_to_base) > prism_height:
        return False
    
    # Project the point onto the plane of the base triangle
    projected_point = point - distance_to_base * normal
    
    # Check if the projected point is inside the base triangle
    return is_point_in_triangle(projected_point, *base_vertices)

def count_atoms_in_prism(atoms, indices, height):
    """Count the number of atoms within the prism formed by translating the base triangle along the normal vector"""
    pos = atoms.get_positions()
    
    # Get positions of the vertices
    p1, p2, p3 = pos[indices[0]], pos[indices[1]], pos[indices[2]]
    
    # Calculate the normal vector
    normal = calculate_normal(p1, p2, p3)
    
    # Translate base vertices along the normal vector
    top_vertices = [v + height * normal for v in [p1, p2, p3]]
    bottom_vertices = [v - height * normal for v in [p1, p2, p3]]
    
    # Check other atoms
    atom_count_top = 0
    atom_count_bottom = 0
    for i in range(len(pos)):
        if i in indices:
            continue
        p = pos[i]
        
        if point_in_prism(p, [p1, p2, p3], top_vertices):
            atom_count_top += 1
        if point_in_prism(p, [p1, p2, p3], bottom_vertices):
            atom_count_bottom += 1
    
    return atom_count_top, atom_count_bottom

def is_exposed_triangle(atoms, indices, height):
    """Check if a triangle formed by three indices is exposed by analyzing prisms"""
    atom_count_top, atom_count_bottom = count_atoms_in_prism(atoms, indices, height)
    return atom_count_top, atom_count_bottom

def is_exposed_rhombus(atoms, rhombus_indices, height=2.5):
    """Check if a rhombus formed by four indices is exposed"""
    triangles = [
        [rhombus_indices[0], rhombus_indices[1], rhombus_indices[2]],
        # [rhombus_indices[0], rhombus_indices[2], rhombus_indices[3]],
        # [rhombus_indices[0], rhombus_indices[1], rhombus_indices[3]],
        [rhombus_indices[1], rhombus_indices[2], rhombus_indices[3]]
    ]
    
    false_or_true = []
    for triangle in triangles:
        atom_count_top, atom_count_bottom = is_exposed_triangle(atoms, triangle, height)
        # print(f"Triangle {triangle} has {atom_count_top} atoms in the top prism and {atom_count_bottom} atoms in the bottom prism.")
        if atom_count_top > 0 and atom_count_bottom > 0:
            false_or_true.append(False)
        else:
            false_or_true.append(True)
    if all(false_or_true) : 
        return True 
    else:
        return False


def find_rhombus(atoms, top_list, bond_threshold=2.6):
    """
    Given a list of exactly four atom indices (top_list), determine their ordering as a rhombus.
    
    The procedure is:
      1) Find the two atoms that are farthest apart; assign the one with the smaller index as A and the other as C.
      2) Of the two remaining atoms, assign the one with the smaller index as B and the other as D.
      3) Check that the bonds A-B, B-C, C-D, and D-A are all shorter than or equal to the bond_threshold.
         If any bond exceeds the threshold, the configuration is considered invalid and an empty list is returned.
    
    Parameters:
        atoms: An object with a method get_distance(i, j) to compute the distance between atoms.
        top_list: A list of exactly four atom indices.
        bond_threshold: The maximum allowed bond length.
    
    Returns:
        A list [A, B, C, D] representing the ordered atoms if valid, or an empty list if invalid.
    """
    if len(top_list) != 4:
        raise ValueError("top_list must contain exactly 4 atoms.")
    
    max_distance = 0.0
    pair = (None, None)
    # Find the pair of atoms with the maximum distance.
    for i in range(4):
        for j in range(i + 1, 4):
            d = atoms.get_distance(top_list[i], top_list[j])
            if d > max_distance:
                max_distance = d
                pair = (top_list[i], top_list[j])
    
    # Assign A and C: the atom with the smaller index becomes A.
    A, C = pair if pair[0] < pair[1] else (pair[1], pair[0])
    
    # The remaining two atoms become B and D; the one with the smaller index is B.
    remaining = [atom for atom in top_list if atom not in (A, C)]
    B, D = remaining if remaining[0] < remaining[1] else (remaining[1], remaining[0])
    
    # Check the bonds: A-B, B-C, C-D, and D-A.
    if (atoms.get_distance(A, B) > bond_threshold or
        atoms.get_distance(B, C) > bond_threshold or
        atoms.get_distance(C, D) > bond_threshold or
        atoms.get_distance(D, A) > bond_threshold):
        return []  # Invalid rhombus configuration.
    
    return [A, B, C, D]

def write_indices_to_file(file_path, indices_list):
    with open(file_path, 'w') as file:
        for indices in indices_list:
            text = ' '.join(map(str, indices))
            file.write(text + '\n')

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

def count_short_metal_nonmetal_bonds(cluster, metal='Ru', max_distance=1.6):
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
    metal_indices = [i for i, atom in enumerate(cluster) if atom.symbol == metal]
    nonmetal_indices = [i for i, atom in enumerate(cluster) if atom.symbol != metal]

    count = 0

    # Calculate distances between metal and non-metal atoms
    for i in metal_indices:
        for j in nonmetal_indices:
            distance = cluster.get_distance(i, j, mic=True)
            if distance < max_distance:
                count += 1
    return count

def find_shortest_bond(atoms, cutoff=2.5):
    # read POSCAR
    # atoms = read(path)
    
    # get the index for N
    n_index = None
    ru_indices = []
    
    for i, atom in enumerate(atoms):
        if atom.symbol == 'N':
            n_index = i
        elif atom.symbol == 'Ru':
            ru_indices.append(i)

    if n_index is None:
        raise ValueError("No nitrogen (N) atom found in the POSCAR file.")  

    #Calculate the distance between N and each Ru atom, and filter out the distances that are less than the cutoff.
    distances = []

    for ru_index in ru_indices:
        distance = atoms.get_distance(n_index, ru_index, mic=True)
        if distance < cutoff:
            distances.append((ru_index, distance))

    
    ##"Sort by distance and find the shortest distance."
    if distances:
        distances.sort(key=lambda x: x[1])
        shortest_distance = distances[0]
        num_short_bonds = len(distances)
    else:
        shortest_distance = None
        num_short_bonds = 0

    return shortest_distance, num_short_bonds


##################=========================

def add_anchor_to_single_site_bri(poscar_path,exposed_indexs, metal='Ru', anchor_atom='N', output_prefix='POSCAR'):
    cluster = read(poscar_path)
    cutoff = 3.0  # Adjust based on your system
    modified_cluster = copy.deepcopy(cluster)
    def get_geo(modified_cluster, exposed_indexs):

        A = modified_cluster[exposed_indexs[0] ]  # Adjust index for 0-based indexing
        B = modified_cluster[exposed_indexs[1] ]  
        
        translation_vector = B.position - A.position
        # Calculate the normal vector for positioning the N atom
        normal_vector = calculate_normal_vector_bri(exposed_indexs[0] , modified_cluster, cutoff)
    
        # Position for N atom
        anchor_position = A.position + normal_vector * 1.6 +  0.5 * translation_vector  
        
        # Add N atom to the structure
        modified_cluster += Atoms(anchor_atom, positions=[anchor_position])
        return modified_cluster

    modified_cluster = copy.deepcopy(cluster)
    cluster_1 = get_geo(modified_cluster, exposed_indexs)
    num_bonds_1 = count_short_metal_nonmetal_bonds(cluster_1, metal, max_distance=1.6) 

    modified_cluster = copy.deepcopy(cluster)
    cluster_2 = get_geo(modified_cluster, exposed_indexs[::-1])
    
    num_bonds_2 = count_short_metal_nonmetal_bonds(cluster_2, metal, max_distance=1.6) 

    final_cluster = cluster_1
    if num_bonds_1 > num_bonds_2: 
        final_cluster = cluster_2 

    shortest_distance, num_short_bonds = find_shortest_bond(final_cluster, 2.5)
    # print(shortest_distance, num_short_bonds)
    if shortest_distance[1] > 1.5 and num_short_bonds < 3:     
        destination_path = os.path.join('bri', '_'.join([str(i) for i in exposed_indexs])) 
        os.makedirs(destination_path, exist_ok=True)
        output_file = destination_path + '/POSCAR'
    
        write(output_file, final_cluster, format='vasp', vasp5=True)

def add_one_bri(path, exposed_indices):
    poscar_path = path +'/POSCAR_bare'
    cluster = read(poscar_path)
    os.makedirs(path + '/bri', exist_ok=True)
    for combination in itertools.combinations(exposed_indices, 2):
        ## get the atoms that form the hollow site 
        sorted_combination = sorted(combination)
        bri_atom_indices = [i  for i in sorted_combination] # indice count from 0
        bri_atoms_positions = [cluster[i].position for i in bri_atom_indices]

        if not is_valid_bri(bri_atoms_positions):
            continue
        
        add_anchor_to_single_site_bri(poscar_path,sorted_combination,output_prefix='POSCAR')
 
def get_new_coordiates(cluster, A1, A2, B1, list_H):
    '''A1, A2, B1 are the coordinates, list_H is the list of coordinates '''
    # Calculate B2 as the center (midpoint) of C and D
    B2 = calculate_triangle_center(list_H) 
    translation_vector = A2 - B1
    B2_translated = B2 + translation_vector 
    
    # Step 2: Calculate vectors a and b
    a = A2 - A1
    a_normalized = a / np.linalg.norm(a)
    b = B2_translated - A2  # New b after moving B1 to A2
    b_normalized = b / np.linalg.norm(a)
    
    # Calculate rotation axis and angle
    axis = np.cross(b_normalized, a_normalized)
    
    axis_normalized = axis / np.linalg.norm(axis) if np.linalg.norm(axis) != 0 else b_normalized
    
    angle = np.arccos(np.clip(np.dot(b_normalized, a_normalized), -1.0, 1.0))     # Angle of rotation
    
    # Apply rotation around A2
    rotation = R.from_rotvec(axis_normalized * angle)
    
    new_list_H = [rotation.apply(i - B1) + A2 for i  in cluster.get_positions()[1:]]
    
    return new_list_H# C_final, D_final, E_final

def add_more_atoms(site):
    ''' Add species with more than three atoms to the anchor atoms. This function works for top, bri, and hollows sites '''
    sites = [int(i)  for i in site.split('_')]
    atoms = read(site + '/POSCAR')
    coordinates = [atoms[i].position for i in sites]
    A1 = calculate_triangle_center(coordinates)  # fcc site 
    A2 = atoms[-1].position     # N site, add N first and then put the NH,NH2,NH3 et al using the N site as reference.

    atoms_temp = read('../POSCAR_temp')
    B1 = atoms_temp[0].position

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
        
    write(site+'/POSCAR_bri', out_geo, format='vasp', vasp5=True)

def add_more_bri(path):
    """Add larger moelcules to all the bridge sites. The sites are the folders generated by the add_one_bri function"""
    os.chdir(path + '/bri')
    all_items = os.listdir('.')
    folders = [item for item in all_items if os.path.isdir(item)]
    folders = [folder for folder in folders if folder.count('_') == 1]
        
    for site in folders:     
        add_more_atoms(site)
        
def get_active_sites(path):   # Path: the cluster model directory
    atoms = read(path + '/POSCAR')
    mult = 0.9 
    bond_length_threshold = 2.7
    metal = 'Ru'

    ############## Step1: obtain the bridge sites
    list_file = list_file = os.path.join(path, 'list')   ### list all the top sites in one line
    if not os.path.exists(list_file):
        # print("Warning: No list file found. Examine the exposed sites using Coordination Numbers.")
        surface_indices = get_top_sites(path, metal, mult=0.9)  ## Seems that this method is not ideal
    else:
        with open(list_file, 'r') as f_in:
            top_indices = f_in.readline().rstrip().split()
            top_indices = [int(i) for i in top_indices]
            surface_indices = [i-1 for i in top_indices]
 
    add_one_bri(path,surface_indices)
    ############## Step2: obtain the rhombus sites
    bri_path = os.path.abspath('./bri')
    all_bri_folders = os.listdir(bri_path)
    bri_sites = [item for item in all_bri_folders if os.path.isdir(os.path.join(bri_path, item))]

    l_rhombus  = []
    for sites in bri_sites: 
        l_sites = sites.split('_')  
        # Bridge site configuration folder example: 43_44, 43 and 44 are the index of atoms in the POSCAR, counting from 0
        B, C  = [int(i) for i in l_sites]
        indices = find_rhombus(atoms, surface_indices, B, C, bond_threshold=3.0)
        l_rhombus.append(indices)
       
    #Check the rhombus sites cannot be formed from some bridge sites 
    for num, site in enumerate(all_bri_folders):
        if len(l_rhombus[num]) == 0:
            sites = [int(i) for i in site.split('_') ]
            # print('No active site is found for: ', site)        
            
    l_rhombus = [item for item in l_rhombus if item] ## remove the empty 

    ############## Step3: categorize  rhombus sites: planar or edge
    coplanar_threshold = 0.5  # Adjust threshold as needed, unit is Å
    coplanar_rhombus = filter_coplanar_rhombuses(l_rhombus, atoms, coplanar_threshold)
    edge_rhombus = filter_rhombuses_by_dihedral(l_rhombus, atoms, 40, 120)

    write_indices_to_file(os.path.join(path, 'all_sites.txt'), l_rhombus)
    write_indices_to_file(os.path.join(path, 'planar_sites.txt'), coplanar_rhombus)
    write_indices_to_file(os.path.join(path, 'edge_sites.txt'), edge_rhombus)

    return l_rhombus, coplanar_rhombus, edge_rhombus


################## Section MKM

def get_k(Ea, T=673):
    kb = 8.617333262145E-5  # Boltzmann's constant in eV/K
    kbT = kb * T  # Boltzmann's constant times Temperature (in eV)
    A = kbT / (h / eV)  # Pre-exponential factor in s^-1, converting h to eV·s
    kf = A * math.exp(-Ea / kbT)  # Activation energy in eV, so no conversion needed
    return kf

def get_rate(kf0, kr0):
    # Given constants and conversion factors
    J2eV = 6.24150974E18            # eV/J
    Na = 6.0221415E23               # mol-1
    T = 823                         # Temperature in Kelvin
    h = 6.626068E-34 * J2eV         # Planck's constant in eV*s
    kb = 1.3806503E-23 * J2eV       # Boltzmann's constant in eV/K
    kbT = kb * T                    # Boltzmann's constant times Temperature
    kJ_mol2eV = 0.01036427
    wave2eV = 0.00012398
    FNH30 = 1                       # Inlet molar flow rate (assumed in mol/s)
    p0 = 1                          # Total pressure in atm
    # Function to calculate the partial pressures
    def get_pressures(FNH30, X, p0):
        # Calculating flow rates of NH3, N2, and H2
        FNH3 = FNH30 * (1 - X)
        FN2 = FNH30 * 0.5 * X
        FH2 = FNH30 * 1.5 * X
        
        # Total flow rate
        FT = FNH3 + FN2 + FH2
        
        # Calculating partial pressures
        pNH3 = (FNH3 / FT) * p0
        pN2 = (FN2 / FT) * p0
        pH2 = (FH2 / FT) * p0
        
        return pNH3, pN2, pH2

    # kf0 = [5631186144.388563, 50.01374140854966, 167.3703043694318, 35375.898901179375, 1.4277669985218864e-08, 10257.81417760215]
    # kr0 = [133466475.11618523, 78037.92020523513, 51662039.04911003, 526671805.90369725, 1.6185464638531345e-08, 66391314202.84105]
    
    # Function to get the rates
    def get_rates(theta, kf, kr, X, FNH30, p0):
        pNH3, pN2, pH2 = get_pressures(FNH30, X, p0)
        
        tNH3, tNH2, tNH, tN, tH = theta
        tstar = 1.0 - tNH3 - tNH2 - tNH - tN - tH  # Site balance for tstar
        
        rate = np.zeros(6)
        rate[0] = kf[0] * pNH3 * tstar - kr[0] * tNH3
        rate[1] = kf[1] * tNH3 * tstar - kr[1] * tNH2 * tH
        rate[2] = kf[2] * tNH2 * tstar - kr[2] * tNH * tH
        rate[3] = kf[3] * tNH * tstar - kr[3] * tN * tH
        rate[4] = kf[4] * tN ** 2 - kr[4] * pN2 * tstar
        rate[5] = kf[5] * tH ** 2 - kr[5] * pH2 * tstar
        return rate
    
    # Function to get the ODEs
    def get_odes(rate):
        dt = [0] * 5
        dt[0] = rate[0] - rate[1]                          # d(tNH3)/dt
        dt[1] = rate[1] - rate[2]                          # d(tNH2)/dt
        dt[2] = rate[2] - rate[3]                          # d(tNH)/dt
        dt[3] = rate[3] - 2 * rate[4]                      # d(tN)/dt
        dt[4] = rate[1] + rate[2] + rate[3] - 2 * rate[5]  # d(tH)/dt
        return dt
    
    def combined_odes(initial_state, t, kf, kr, FNH30, p0):
        X = initial_state[-1]
        theta = initial_state[:-1]
        
        rates = get_rates(theta, kf, kr, X, FNH30, p0)
        
        ds_dt = get_odes(rates)
        
        # Example improvement: Calculate dx_dt based on a specific rate
        # Assuming rate[5] corresponds to a key product formation step
        dx_dt = rates[4] / FNH30  # Normalize by inlet flow to represent conversion rate
        
        comb_odes = np.append(ds_dt, dx_dt)
        
        return comb_odes
    
    def solve_combined_odes(kf, kr, X0, FNH30, p0):
        thetas0 = np.zeros(5)
        initial_state = np.append(thetas0, X0)
        start_time = 0
        end_time = 1e6
        
        def combined_odes_wrapper(t, y):
            return combined_odes(y, t, kf, kr, FNH30, p0)
            
        # Solving the ODEs
        solution = solve_ivp(combined_odes_wrapper, [start_time, end_time], initial_state, method='LSODA', dense_output=True)
    
        t = np.linspace(start_time, end_time, 10000) 
        sol = solution.sol(t)
        thetas_solution = sol[:-1, :]
        X_solution = sol[-1, :]
            
        return thetas_solution, X_solution, t 
    
    # Solving for theta0 and rates
    X0 = 0.0
    # print('Running')
    thetas_t, X_t, t = solve_combined_odes(kf0, kr0, X0, FNH30, p0)
    # print('Done')
    rates = get_rates(thetas_t[:,-1], kf0, kr0, X_t[-1],  FNH30 , p0)
    tof = p0 *  X_t[-1]  / 10000 
    print('tof:', tof, 'lg(tof):', math.log10(tof), 'Conversion:', X_t[-1], 'rate_rds:',  min(rates)  )
    
    # print(thetas_t[:,-1],  X_t[-1])
    # print(math.log(p0 /100000 * X_t[-1]  / 1000000 ))
    # print("XNH3:", X_t[-1], )

    return rates
