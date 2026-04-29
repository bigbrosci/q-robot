import os
import sys
import csv
import json
import numpy as np
import pandas as pd
import itertools
import copy
import math
import matplotlib.pyplot as plt

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


def convert_unix_to_windows_path(unix_path):
    # Replace forward slashes with backslashes
    windows_path = unix_path.replace('/', '\\')
    
    # Check if the path starts with '/mnt/', which indicates a mounted drive in WSL
    if windows_path.startswith('\\mnt\\'):
        # Extract the drive letter and the rest of the path
        drive_letter = windows_path[5]
        rest_of_path = windows_path[6:]
        # Construct the Windows path
        windows_path = f'{drive_letter.upper()}:\\{rest_of_path}'
    
    return windows_path


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

def get_CN_GA(path, mult=0.9):
    """Get the matrix for group additivity"""
    atoms, connections, cn_of_connected_atoms, exposed_top_sites = get_connection(path, metal='Ru', mult=0.9)
    
    def number_to_letter(num):
        """Conver the CN to letters"""
        if num == 0:
            return '0'
        return chr(ord('a') + num - 1)
    
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

    groups = []
    GA_file = path + '/GA_dict.txt'
    # Writing JSON data
    with open(GA_file, 'w') as f:
        json.dump(GA_dict, f)

    ## summarize and get all the groups
    for key, value in GA_dict.items():
        for group in value:
            groups.append(group)
    
    groups = list(set(groups))
    # Write the groups to a txt file
    groups_file = path + '/groups.txt'
    with open(groups_file, 'w') as f:
        for group in groups:
            f.write(f"{group}\n")

    ### SAVE using csv
    header = ['site'] + groups
    
    # Prepare data for CSV
    csv_data = [header]
    for key, value in GA_dict.items():
        GA_num = [str(value.count(group)) for group in groups]
        row = [str(key)] + GA_num
        csv_data.append(row)
    
    # Save to CSV
    csv_filename = path + '/GA_matrix.csv'
    with open(csv_filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(csv_data)
    return groups


def get_CN_GA_cluster(path,metal='Ru', mult=0.9):
    
#    path = convert_unix_to_windows_path(path)
    atoms_in = read(path + '/POSCAR')
    filtered_atoms = [atom for atom in atoms_in if atom.symbol  in [metal]]
    # filtered_atoms = [atom for atom in atoms_in if atom.symbol not in ['N', 'H']]
    # Create a new Atoms object with only 'Ru'
    atoms = Atoms(filtered_atoms)
    radii = natural_cutoffs(atoms, mult=mult)
    nl = NeighborList(radii, bothways=True, self_interaction=False)
    nl.update(atoms)
    CN_matrix = nl.get_connectivity_matrix(sparse=False)
    CN = CN_matrix.sum(axis=0)
    
    # Request 1: Get the connections of all atoms
    connections = {}
    for i in range(len(atoms)):
        connections[i] = list(np.where(CN_matrix[i])[0])
    
    # Request 2: Calculate the CN of the connected atoms
    cn_of_connected_atoms = {}
    for i in range(len(atoms)):
        connected_indices = connections[i]
        cn_of_connected_atoms[i] = [CN[j] for j in connected_indices]
    
    def number_to_letter(num):
        if num == 0:
            return '0'
        return chr(ord('a') + num - 1)
    
    def get_one_site_string(site):
        cn_values = sorted(cn_of_connected_atoms[site])
        if len(cn_values) < 13:
            cn_values += [0] * (13 - len(cn_values))  # Pad with zeros if less than 12 elements
            cn_string = ''.join(number_to_letter(num) for num in cn_values)
        return cn_string
    
    def get_one_top_site(site):
        groups = []
        groups.append(get_one_site_string(site))
        surrounding_atoms = connections[site]
        for s_atom in surrounding_atoms:
            groups.append(get_one_site_string(s_atom))
        return groups
        
    def get_one_bridge_site(sites) :
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
  
    def get_one_hollow_site(sites) :
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

    # Exposed top sites (CN <= 9)
    exposed_top_sites = [i for i, cn in enumerate(CN) if cn <= 9]
    
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
        # key = site + 1 
        key = site
        GA_dict[key] = get_one_top_site(site)

    for sites in bridge_sites:
        # key = '_'.join([str(i+1) for i in sites])
        key = '_'.join([str(i) for i in sites])
        GA_dict[key] = get_one_bridge_site(sites)

    for sites in hollow_sites:
        key = '_'.join([str(i) for i in sites])
        GA_dict[key] = get_one_bridge_site(sites)

    groups = []
    GA_file = path + '/GA_dict.txt'
    # Writing JSON data
    with open(GA_file, 'w') as f:
        json.dump(GA_dict, f)

    ## summarize and get all the groups
    for key, value in GA_dict.items():
        for group in value:
            groups.append(group)
    
    groups = list(set(groups))
    # Write the groups to a txt file
    groups_file = path + '/groups.txt'
    with open(groups_file, 'w') as f:
        for group in groups:
            f.write(f"{group}\n")
            
    ### SAVE using csv
    header = ['site'] + groups
    
    # Prepare data for CSV
    csv_data = [header]
    for key, value in GA_dict.items():
        GA_num = [str(value.count(group)) for group in groups]
        row = [str(key)] + GA_num
        csv_data.append(row)
    
    # Save to CSV
    csv_filename = path + '/GA_matrix.csv'
    with open(csv_filename, mode='w', newline='') as file:
        writer = csv.writer(file)
        writer.writerows(csv_data)
    return groups

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

def GA_prediction(X,Y):
    '''Run conventional GA approach '''
    # # df = pd.read_csv('group_matrix.csv')
    # Y = df_matrix['Eads'].values
    # df_matrix = df_matrix.drop(columns=['site', 'Eads'])
    # X = df_matrix.values # Convert to numpy arrays
    
    model = PseudoInverseLinearRegression() 
    model.fit(X,Y)
    Y_pred = model.predict(X)
    
    # Calculate the MAE, RMSE, and R²
    mae = mean_absolute_error(Y, Y_pred)
    rmse = np.sqrt(mean_squared_error(Y, Y_pred))
    r2 = r2_score(Y, Y_pred)
    
    # Print the errors and R²
    print(f'Mean Absolute Error (MAE): {mae}')
    print(f'Root Mean Squared Error (RMSE): {rmse}')
    print(f'R-squared (R²): {r2}')
         
    # # Plot the data and regression line
    plt.figure(figsize=(10, 6))
    plt.scatter(Y, Y_pred, color='blue', label='Model')
    
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
    return Y, Y_pred
 

def XY_linear(X, Y):
    # Reshape X and Y
    X = X.reshape(-1, 1)
    Y = Y.reshape(-1, 1)

    model = LinearRegression()
    model.fit(X, Y)
    Y_pred = model.predict(X)
    
    # Calculate metrics
    mae = mean_absolute_error(Y, Y_pred)
    rmse = np.sqrt(mean_squared_error(Y, Y_pred))
    r2 = r2_score(Y, Y_pred)
    
    equation = r'Ea = {:.2f} $\Delta$E + {:.2f}'.format(model.coef_[0, 0], model.intercept_[0])
    
    # Plot the data and regression line
    plt.figure(figsize=(10, 6))
    plt.scatter(X, Y, color='blue', label='DFT')
    plt.plot(X, Y_pred, color='red', label='Reg.')
    
    # Add text box for metrics
    textstr = f'{equation}\nMAE: {mae:.2f}\nRMSE: {rmse:.2f}\nR2: {r2:.2f}'
    plt.gca().text(0.65, 0.25, textstr, transform=plt.gca().transAxes,
                    fontsize=16, verticalalignment='top', bbox=dict(facecolor='white', alpha=0.5))
    
    # Customize the plot
    plt.xlabel(r'X /eV', fontsize=16)
    plt.ylabel('Y /eV', fontsize=16)
    plt.xticks(fontsize=16)
    plt.yticks(fontsize=16)
    plt.legend(loc='upper left', frameon=False, fontsize=18)
    plt.tight_layout()
    
    # Save the figure
    plt.savefig('XY_linear.png')
    # plt.show()
    

def test_one_model(X, Y, models, model_name):
    '''Run conventional GA approach '''
    model = models[model_name]
    model.fit(X,Y)
    Y_pred = model.predict(X)
    
    # Calculate the MAE, RMSE, and R²
    mae = mean_absolute_error(Y, Y_pred)
    rmse = np.sqrt(mean_squared_error(Y, Y_pred))
    r2 = r2_score(Y, Y_pred)
    
    # Print the errors and R²
    # print(f'Mean Absolute Error (MAE): {mae}')
    # print(f'Root Mean Squared Error (RMSE): {rmse}')
    # print(f'R-squared (R²): {r2}')
    
    fig, ax  = plt.subplots(figsize=(10,8))  
    
    species = get_species_from_cwd()
    text = '%s adsorption' %(species)
    # ax.set_title(text)
    ax.scatter(Y, Y_pred, edgecolor='k', alpha=0.7, s=200, label = text, color='#1f77b4')
    ax.plot([min(Y), max(Y)], [min(Y), max(Y)], 'r--', linewidth=3)
    ax.set_xlabel('DFT Eads (eV)', fontsize=36)
    ax.set_ylabel('Predicted Eads (eV)', fontsize=36)

    textstr = f'MAE: {mae:.2f}  \nRMSE: {rmse:.2f}  \nR2: {r2:.2f}'
    # plt.gca().text(0.65, 0.25, textstr, transform=plt.gca().transAxes,
    #                 fontsize=20, verticalalignment='top', bbox=dict(facecolor='white', alpha=0))
    # ax.tick_params(axis='both', which='major', labelsize=30)
# 
    # ax.xticks(fontsize=30)
    # ax.yticks(fontsize=30)
    # plt.legend(loc='upper left', frameon=False, fontsize=26)
 
    # plt.legend(framealpha=0, fontsize = 18)
    plt.legend(loc='upper left', frameon=False, fontsize=40)

    plt.grid(False)
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(3)
        
    #plt.title('Linear Regression of De vs Ea')
    ax.tick_params(axis='both', which='major', width=2, length=8, labelsize=28)
    plt.tight_layout()
    # Save the figure
    plt.savefig('results/XY_%s.png' %(model_name), dpi=600)
    return textstr, Y, Y_pred
 
    
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
    
 # plt.plot(EF, E_Hollow, label='E_Hollow', color='#ff7f0e',linewidth=3,linestyle=":", marker='.', markersize=14)
 # plt.plot(EF, E_Edge, label='E_Edge', color='#1f77b4', linewidth=3, marker='*',linestyle=":", markersize=14)# marker='s')

 # plt.xlabel(r'EF (V/$\mathrm{\AA}$)', fontsize=18)
 # plt.ylabel('Eads$^N$ (eV)', fontsize=18)
 # # plt.title('Energy vs EF Plot')
 # plt.legend()
 # plt.grid(False)
 # ax = plt.gca()
 # for spine in ax.spines.values():
 #     spine.set_linewidth(3)
 # ax.tick_params(axis='both', which='major', width=2, length=8, labelsize=18)
 # plt.legend(loc='upper right', frameon=False, fontsize=18)
 
    

def ads_training(test_path):
    # test_path = '/mnt/c/Users/lqlhz/OneDrive - UMass Lowell/Projects/P1_cluster/Ru_clusters/run_3/Ru30/sum/NH2_ads'
    df_GA, X, Y, n_group = get_GA_X_Y(test_path)
    # print(n_group)
    models  = refine_models(n_group)
    model_name = 'GA'
    species = get_species_from_cwd(test_path)
    
    GA_model = models[model_name]
    GA_model.fit(X,Y)

    #### Perform the bootstrap cross-validation
    
    results = bootstrap_cross_validation_one_model(models, model_name, X, Y, n_bootstrap=100)
    
    try:
        df_mae = pd.read_csv(test_path + '/results/mae.csv')
    except: 
#        test_path =  convert_unix_to_windows_path(test_path)
        df_mae = pd.read_csv(test_path + '/results/mae.csv')
    
    average_values = df_mae.mean()
    plt.figure(figsize=(10, 6))
    average_values.plot(kind='bar', label = species)
         
    
    plt.grid(False)
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(3)
    #plt.title('Linear Regression of De vs Ea')
    ax.tick_params(axis='both', which='major', width=2, length=8, labelsize=18)
    plt.legend(loc='upper left', frameon=False, fontsize=26)
    
    plt.xlabel('Model')
    plt.ylabel('MAE / eV')
    # plt.title('Bar Chart of MAE Values')
    # plt.xticks(rotation=45)
    plt.tight_layout()
    plt.legend(framealpha=0, fontsize = 18)
    plt.savefig('results/mae.png', dpi=300)
 
    
 
def get_GA_X_Y_self(EF):
    ''' Gas phase energy of itself is the reference energy'''
    if not os.path.exists('results'):
        os.makedirs('results')
    data_file = './sorted_data.csv'
    df = pd.read_csv(data_file)
    list_path = df['path'].tolist()
    df.set_index('Species', inplace=True)

#    list_path = [convert_unix_to_windows_path(i) for i in list_path]
    
    ### Collect all the groups in the dataset. 
    groups_data = []
    for path in list_path:
        try: 
            groups = get_CN_GA(path, mult=0.9)
        except: 
            #path = convert_unix_to_windows_path(path)
            groups = get_CN_GA(path, mult=0.9)
        groups_data.extend(groups)
        
    groups_data = list(set(groups_data))   # The sequence of the groups are fixed here.
    n_group = len(groups_data)    # number of groups

    ### analyze the individual configurations and get the Group and its numbers 
    list_site = df['site'].tolist()
    data = []    
    
    for num, path in enumerate(list_path):
        GA_file = os.path.join(path, 'GA_dict.txt')
        try:
            with open(GA_file, 'r') as f:
                GA_dict = json.load(f)
            sites = str(list_site[num]) ## Top, Bri, and Hollow only
            local_groups = []
            local_groups.extend(GA_dict[sites]) 
        except:
            print('Error', path, sites)     
            local_groups = []
                
        GA_num = [local_groups.count(group) for group in groups_data]        
        
        row = [list_site[num]] + GA_num
        data.append(row)

    # Create the DataFrame
    columns = ['site'] + groups_data  
    df_GA = pd.DataFrame(data, columns=columns)  # Group Matrix has been stored in the dataframe
    
    # Get the Y value for GA model, adsorption at EF
    species  = get_species_from_cwd()  
    
        # Get adsorption energy
    try:
        slab_energy = df.loc['slab', EF]  # Ensure 'slab' exists
    except KeyError:
        raise KeyError(f"'slab' entry is missing in the Species column.")


    if species in ['H', 'N', 'O']:
        gas_species = species + '2'
        df['Eads'] = df[EF] - slab_energy -  float(gas_dict[gas_species]/2)
    else:
        gas_species = species
        df['Eads'] = df[EF] - slab_energy - float(gas_dict[gas_species] )

    # print(df.index)
    # print(df_GA.index)
    
    # Save to CSV
    df = df.reset_index()  # Convert Species index back to a column
    df_GA['Eads'] = df['Eads']
    # print(df['Eads'])
    
    df_GA = df_GA[df_GA['site'] != 'slab']
    print(df_GA['Eads'])
    csv_filename = './results/group_matrix.csv'
    df_GA.to_csv(csv_filename, index=False)
    df_matrix = df_GA.copy()
    
    # df_GA = df.read_csv(csv_filename)
    Y = df_matrix['Eads'].values
    df_matrix = df_matrix.drop(columns=['site', 'Eads'])
    X = df_matrix.values # Convert to numpy arrays
    
    return df_GA, X, Y, n_group


def get_GA_X_Y(data_path, cluster_path):
    """NH3 gas phase energy is the reference""" 
    
#    data_path = convert_unix_to_windows_path(data_path)
    if not os.path.exists('results'):
        os.makedirs('results')

    data_file = data_path + '/sorted_data.csv'
    df = pd.read_csv(data_file)

        
    list_path = df['path'].tolist()
    list_path = [element.replace('/OUTCAR', '') for element in  list_path]
    
    # groups_data = []
    # for path in list_path:
    #     path = convert_unix_to_windows_path(path)
    #     groups = get_CN_GA(path, mult=0.9)
    #     groups_data.extend(groups)
    # groups_data = list(set(groups_data))   
    # n_group = len(groups_data)   
  
    
    GA_file =  cluster_path + '/GA_dict.txt'
#    GA_file = convert_unix_to_windows_path(GA_file)
    with open(GA_file, 'r') as f:
        GA_dict = json.load(f)
    
    groups_file = cluster_path + '/groups.txt'
    
    # Read the groups from the text file
    with open(groups_file, 'r') as f:
        groups_data = f.read().splitlines()
        # n_group = len(groups_data)  

    list_site = df['site'].tolist()
    data = []
    
    for num, path in enumerate(list_path):

        # print(path,'test')
        try:
            sites = str(list_site[num]-1)
        except:
            sites = list_site[num].split('_')
            sites = [int(i)-1 for i in sites]
            sites = '_'.join([str(i) for i in sites])
        local_groups = []
        try:
            local_groups.extend(GA_dict[sites]) 
        except:
            print('Error', path, sites) 

        GA_num = [local_groups.count(group) for group in groups_data]        
        
        row = [list_site[num]] + GA_num
        data.append(row)
 
    # Create the DataFrame
    columns = ['site'] + groups_data  
    df_GA = pd.DataFrame(data, columns=columns)
    species  = get_species_from_cwd(data_path)

    # print(species, 'test')
    if species == 'H' :# 'N', 'O']:
        # gas_species = species + '2'
        df_GA['Eads'] = df['E'] - gas_dict['H2'] / 2  -  gas_dict['slab']
        # print(df_GA['Eads'])
    else:
        gas_species = 'NH3'
        if species == 'NH3':
            df_GA['Eads'] = df['E'] + gas_dict['H2'] * 0.0  - gas_dict[gas_species]  -  gas_dict['slab']
        elif species == 'NH2':
            df_GA['Eads'] = df['E'] + gas_dict['H2'] * 0.5  - gas_dict[gas_species]  -  gas_dict['slab']
        elif species == 'NH':
            df_GA['Eads'] = df['E'] + gas_dict['H2'] * 1.0  - gas_dict[gas_species]  -  gas_dict['slab']
        elif species == 'N':
            df_GA['Eads'] = df['E'] + gas_dict['H2'] * 1.5  - gas_dict[gas_species]  -  gas_dict['slab']
                           
        # gas_species = species
        # df_GA['Eads'] = df['E'] - gas_dict[gas_species]  -  gas_dict['slab']

    # Save to CSV
    csv_filename = 'results/group_matrix.csv'
    df_GA.to_csv(csv_filename, index=False)
    df_matrix = df_GA.copy()
    # print(len(df_matrix.columns.tolist()))
    # df_GA = df.read_csv(csv_filename)
    Y = df_matrix['Eads'].values
    df_matrix = df_matrix.drop(columns=['site', 'Eads'])
    # print(len(df_matrix.columns.tolist()))
    # X = df_matrix.values # Convert to numpy arrays
    X = np.array(df_matrix)

    return df_GA, X, Y, groups_data


def get_beta(df, E) :
    # # Read the energy data and get the GAVs  
    # Ensure that the DataFrame df and the energy vector E have compatible dimensions
    if len(E) != len(df):
        raise ValueError("The length of the energy vector E does not match the number of rows in the matrix df.")
    
    # Compute beta using the pseudoinverse of the matrix df
    beta_E = np.linalg.pinv(df).dot(E)
    
    # Save the beta values to a CSV file
    beta_out = pd.DataFrame({'GAV': beta_E})
    beta_out.to_csv('beta.csv', index=False) 
    
    return beta_E 

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

    
    
# def bootstrap_cross_validation_one_model(models, model_name, X, Y, ratio, n_bootstrap=100):
    
#     model = models[model_name]
    
#     results = {metric: {model_name: [] } for metric in ['MAE', 'MSE', 'RMSE', 'R2', 'MPD', 'MND']}

#     n_samples = int(len(X) * ratio)  # 80% of the data for training
    
#     for i in range(n_bootstrap):
#         print('Iteration: ', i)
#         indices = np.random.choice(len(X), size=len(X), replace=True)
#         X_train, Y_train = X[indices[:n_samples]], Y[indices[:n_samples]]
#         X_test, Y_test = X[indices[n_samples:]], Y[indices[n_samples:]]

#         # print(name, 'running')
#         model.fit(X_train, Y_train)
#         predictions = model.predict(X_test)
        
#         # Store each metric for this iteration
#         results['MAE'][model_name].append(mean_absolute_error(Y_test, predictions))
#         results['MSE'][model_name].append(mean_squared_error(Y_test, predictions))
#         results['RMSE'][model_name].append(np.sqrt(mean_squared_error(Y_test, predictions)))
#         results['R2'][model_name].append(r2_score(Y_test, predictions))
#         results['MPD'][model_name].append(np.max(predictions - Y_test))
#         results['MND'][model_name].append(np.min(predictions - Y_test))
   

#     ### Step 2: Save results to CSV files
#     for metric, model_results in results.items():
#         df = pd.DataFrame(model_results)
#         averages = df.mean()
#         df.to_csv(f'results/{metric.lower()}_{ratio}.csv', index=False)
#         averages.to_csv(f'results/{metric.lower()}_averages_{ratio}.csv')

#     return results

def bootstrap_cross_validation(models, X, Y, n_bootstrap=100):
    print('Run Cross Validation >>> \n')
    # This will store the full list of scores for each metric for each model
    # results = {metric: {name: [] for name in models} for metric in ['MAE', 'MSE', 'RMSE', 'R2', 'Max Positive Deviation', 'Min Negative Deviation']}
    results = {metric: {name: [] for name in models} for metric in ['MAE', 'MSE', 'RMSE', 'R2', 'MPD', 'MND']}

    n_samples = int(len(X) * 0.70)  # 80% of the data for training

    for i in range(n_bootstrap):
        print('Iteration: ', i)
        indices = np.random.choice(len(X), size=len(X), replace=True)
        X_train, Y_train = X[indices[:n_samples]], Y[indices[:n_samples]]
        X_test, Y_test = X[indices[n_samples:]], Y[indices[n_samples:]]

        for name, model in models.items():
            # print(name, 'running')
            model.fit(X_train, Y_train)
            predictions = model.predict(X_test)
            
            # Store each metric for this iteration
            results['MAE'][name].append(mean_absolute_error(Y_test, predictions))
            results['MSE'][name].append(mean_squared_error(Y_test, predictions))
            results['RMSE'][name].append(np.sqrt(mean_squared_error(Y_test, predictions)))
            results['R2'][name].append(r2_score(Y_test, predictions))
            results['MPD'][name].append(np.max(predictions - Y_test))
            results['MND'][name].append(np.min(predictions - Y_test))
   

    ### Step 2: Save results to CSV files
    for metric, model_results in results.items():
        df = pd.DataFrame(model_results)
        averages = df.mean()
        df.to_csv(f'results/{metric.lower()}.csv', index=False)
        averages.to_csv(f'results/{metric.lower()}_averages.csv')

    return results

# # Now include this custom model in the models dictionary
# def refine_models(n_group):
#     models = {
#         "Ridge": Ridge(alpha=1.0),
#         "Lasso": Lasso(alpha=1/n_group),
#         "ElasticNet": ElasticNet(alpha=0.1, l1_ratio=0.5),
#         "DecisionTree": DecisionTreeRegressor(max_depth=n_group),
#         "RandomForest": RandomForestRegressor(n_estimators=n_group*5, max_depth=n_group*5),
#         "GradientBoosting": GradientBoostingRegressor(n_estimators=n_group*5, max_depth=n_group),
#         "SVR": SVR(kernel='linear', C=1.0),
#         "GroupAdditivity": PseudoInverseLinearRegression(),
#         "XGBoost": XGBRegressor(n_estimators=11*5, max_depth=11, learning_rate=0.1),
#         }
#     return models

def refine_models(n_group):
    models = {
        "Ridge": Ridge(alpha=1.0),
        "DT": DecisionTreeRegressor(max_depth=n_group),
        "RF": RandomForestRegressor(n_estimators=n_group*5, max_depth=n_group*5),
        "GB": GradientBoostingRegressor(n_estimators=n_group*5, max_depth=n_group),
        "GA": PseudoInverseLinearRegression(),
        "XGB": XGBRegressor(n_estimators=11*5, max_depth=11, learning_rate=0.1),
        }
    return models


def models_test(models, X, Y, df):
    # Fit each model and calculate metrics
    results = {}
    # print(df)
    
    for name, model in models.items():
        model.fit(X, Y)
        Y_pred = model.predict(X)
        
        mae = mean_absolute_error(Y, Y_pred)
        rmse = np.sqrt(mean_squared_error(Y, Y_pred))
        r2 = r2_score(Y, Y_pred)
        
        results[name] = {"MAE": mae, "RMSE": rmse, "R²": r2}
        df[f'E_ads_predicted_{name}'] = Y_pred
    
    # Print the results
    # for name, metrics in results.items():
    #     print(f'{name}: MAE = {metrics["MAE"]}, RMSE = {metrics["RMSE"]}, R² = {metrics["R²"]}')
    
    # Save the results to a new CSV file
    output_path = './results/thermo_with_predicted_E_ads_comparison.csv'
    df.to_csv(output_path, index=False)
    
    # Plotting
    fig, axes = plt.subplots(nrows=2, ncols=3, figsize=(30,18))
    axes = axes.flatten()
    
    # cwd = os.getcwd()
    species  = get_species_from_cwd()
    
    for idx, (name, metrics) in enumerate(results.items()):
        text = f'{name}' + ' ' + species
        ax = axes[idx]
        Y_pred = df[f'E_ads_predicted_{name}']
        ax.scatter(Y, Y_pred, edgecolor='k', alpha=0.7, s=100, label = text)
        ax.plot([min(Y), max(Y)], [min(Y), max(Y)], 'r--')
        ax.set_xlabel('DFT Eads / eV', fontsize=22)
        ax.set_ylabel('Predicted Eads / eV', fontsize=22)
        # ax.set_title(text)
        ax.text(0.05, 0.95, f'MAE: {metrics["MAE"]:.2f}\nRMSE: {metrics["RMSE"]:.2f}\nR²: {metrics["R²"]:.2f}', 
                transform=ax.transAxes, fontsize=22, verticalalignment='top')
        legend_entries = [text]
        handles = [plt.Line2D([], [], color='none', label=entry) for entry in legend_entries]
        
        ax.legend(handles=handles, loc='lower right', frameon=False, fontsize=26)
        # Identify outliers
        for i in range(len(Y)):
            if abs(Y[i] - Y_pred[i]) > 0.2:
                ax.annotate(df['site'][i], (Y[i], Y_pred[i]), textcoords="offset points", xytext=(0,10), ha='center')
                ax.scatter(Y[i], Y_pred[i], color='red', s=100, edgecolor='k')
    
    plt.tight_layout()
    plt.savefig('results/ML_results.png', dpi=300)    
 
    
 
### GA Applilcation 

    
    
def find_all_rhombuses(atoms, connections, surface_indices, bond_length_threshold):
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

def find_rhombus(atoms, top_list, B, C, bond_threshold=2.6):
    # 计算 BC 的距离
    BC_distance = atoms.get_distance(B, C)

    # 对 top_list 进行排序
    top_list_sorted = sorted(top_list)

    # 用于存储结果的列表
    rhombus_list = []

    for A in top_list_sorted:
        if A == B or A == C:
            continue

        # 判断 A 与 B 和 C 是否成键
        if atoms.get_distance(A, B) <= bond_threshold and atoms.get_distance(A, C) <= bond_threshold:
            # print('A',A)
            for D in top_list_sorted:
                if D == A or D == B or D == C:
                    continue

                # 判断 D 与 B 和 C 是否成键
                if atoms.get_distance(D, B) <= bond_threshold and atoms.get_distance(D, C) <= bond_threshold:
                    # print('D', D)
                    AD_distance = atoms.get_distance(A, D)
                    if AD_distance > BC_distance:
                        rhombus_list = [A, B, C, D]
                        return rhombus_list  # 

    return rhombus_list

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
    # 读取 POSCAR 文件
    # atoms = read(path)
    
    # 提取 N 和 Ru 原子的索引
    n_index = None
    ru_indices = []
    
    for i, atom in enumerate(atoms):
        if atom.symbol == 'N':
            n_index = i
        elif atom.symbol == 'Ru':
            ru_indices.append(i)

    if n_index is None:
        raise ValueError("No nitrogen (N) atom found in the POSCAR file.")
    
    # 计算 N 和每个 Ru 原子之间的距离，并筛选小于 cutoff 的距离
    distances = []

    for ru_index in ru_indices:
        distance = atoms.get_distance(n_index, ru_index, mic=True)
        if distance < cutoff:
            distances.append((ru_index, distance))

    # 按距离排序并找出最短的距离
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
        bri_atom_indices = [i   for i in sorted_combination] # indice count from 0
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
    # path = os.getcwd()
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
