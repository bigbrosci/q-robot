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
import shutil
import subprocess
from openpyxl import load_workbook
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
from sklearn.metrics import r2_score, mean_absolute_error, mean_squared_error


import numpy as np
from ase import Atoms
from ase.io import read, write

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
# from xgboost import XGBRegressor
from scipy.optimize import curve_fit
from scipy.stats import pearsonr
from sklearn.metrics import mean_absolute_error, mean_squared_error

from mkm import * 

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



def get_connection(atoms_in, metal='Cu', mult=1.0):
    # atoms_in = read(path + '/POSCAR')

    filtered_atoms = [atom for atom in atoms_in if atom.symbol in [metal]]
    atoms = Atoms(filtered_atoms)
    radii = natural_cutoffs(atoms, mult=mult)
    nl = NeighborList(radii, bothways=True, self_interaction=False)
    nl.update(atoms)
    CN_matrix = nl.get_connectivity_matrix(sparse=False)
    CN = CN_matrix.sum(axis=0)
    exposed_top_sites = [i for i, cn in enumerate(CN) if cn <= 12]       # Exposed top sites (CN <= 9)

    connections = {}     # Get the connections of all atoms
    for i in range(len(atoms)):
        connections[i] = list(np.where(CN_matrix[i])[0])

    cn_of_connected_atoms = {}      # Request 2: Calculate the CN of the connected atoms
    for i in range(len(atoms)):
        connected_indices = connections[i]
        cn_of_connected_atoms[i] = [CN[j] for j in connected_indices]

    exposed_connections = {i: [j for j in connections[i] if j in exposed_top_sites] for i in exposed_top_sites}

    # Bridge sites (pairs of connected exposed top sites)
    bridge_sites = []
    for i in exposed_top_sites:
        for j in exposed_connections[i]:
            if i < j:  # Ensure each pair is considered only once
                # Add constraint: at least one atom has CN <= 10
                #if CN[i] <= 10 or CN[j] <= 10:
                bridge_sites.append([i, j])

    # Hollow sites (triplets of closely connected exposed top sites)
    hollow_sites = []
    for i in exposed_top_sites:
        for j in exposed_connections[i]:
            for k in exposed_connections[j]:
                if i < j and j < k and i in exposed_connections[k]:  # Ensure each triplet is considered only once and forms a triangle
                    # Add constraint: at least one atom has CN <= 9
                    #if CN[i] <= 9 or CN[j] <= 9 or CN[k] <= 9:
                    hollow_sites.append([i, j, k])
    # Square sites (quartets of connected exposed top sites forming a closed loop)
    square_sites = []
    for i in exposed_top_sites:
        for j in exposed_connections[i]:
            if j <= i:
                continue
            for k in exposed_connections[j]:
                if k in (i, j) or k <= j:
                    continue
                for l in exposed_connections[k]:
                    if l in (i, j, k) or l <= k:
                        continue
                    # Check if l connects back to i to form a closed loop
                    if i in exposed_connections[l]:
                        square = [i, j, k, l]
                        # Add constraint: at least one atom has CN <= 9
                        #if CN[i] <= 9 or CN[j] <= 9 or CN[k] <= 9 or CN[l] <= 9:
                        square_sites.append(square)

    return connections, cn_of_connected_atoms, exposed_top_sites, bridge_sites, hollow_sites, square_sites

def number_to_letter(num):
    """Conver the CN to letters"""
    if num == 0:
        return '0'
    return chr(ord('a') + num - 1)


def get_CN_GA(path, mult=1.0, metal='Ru'):
    
    '''
    Important: it is the path for the clean cluster. 
    1) the exposed surface of the cluster will be scanned and all the top, bridge, and hollow sites will be stored.
    2) the sites will be stored as a dictionray 
    3) all the groups will be stored too to  fix the order of the  groups in the matrix
    '''    
    file_in = os.path.join(path,'POSCAR')
    atoms = read(file_in)
    print(metal)
    connections, cn_of_connected_atoms, exposed_top_sites, bridge_sites, hollow_sites, square_sites  = get_connection(atoms, metal=metal, mult=0.9)
    
    def get_one_site_string(site):
        """Obtain the text string representation for single site """
        cn_values = sorted(cn_of_connected_atoms[site])
        # print(cn_values)
        if len(cn_values) < 13:
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
