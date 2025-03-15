import csv
import numpy as np
import pandas as pd
from sys import argv
import os
import json

from ase.io import read, write
from ase.neighborlist import NeighborList, natural_cutoffs
from ase.visualize import view
from ase import Atoms
from ase.io.vasp import read_vasp, write_vasp

import matplotlib.pyplot as plt
plt.rc('font', size=18)
import matplotlib.patches as patches
from matplotlib.collections import PatchCollection
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm
from matplotlib.colors import ListedColormap

from sklearn.model_selection import cross_val_score, KFold
from sklearn.linear_model import LinearRegression, Ridge, Lasso, ElasticNet, BayesianRidge, LogisticRegression
from sklearn.tree import DecisionTreeRegressor
from sklearn.ensemble import RandomForestRegressor, GradientBoostingRegressor
from sklearn.svm import SVR
from sklearn.metrics import make_scorer, mean_absolute_error, mean_squared_error, r2_score
from sklearn.neural_network import MLPRegressor
from xgboost import XGBRegressor
from sklearn.decomposition import PCA
from sklearn.cross_decomposition import PLSRegression
from sklearn.pipeline import make_pipeline


metal = 'Ru'

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

def calculate_cn(structure, metal, cutoff=2.7):
    """Calculate the coordination number for each atom in the structure, considering only metal atoms."""
    distances = structure.get_all_distances(mic=True)
    cn = np.zeros(len(structure))
    for i, atom in enumerate(structure):
        if atom.symbol == metal:
            neighbors = np.where((distances[i] < cutoff) & (distances[i] > 0))[0]
            cn[i] = sum(1 for n in neighbors if structure[n].symbol == metal)
    return cn

def get_CN(atoms, mult=0.9):
    radii = natural_cutoffs(atoms, mult=mult)
    nl = NeighborList(radii, bothways=True, self_interaction=False)
    nl.update(atoms)
    CN_matrix = nl.get_connectivity_matrix(sparse=False)
    CN = CN_matrix.sum(axis=0)
    return CN, CN_matrix


def calculate_gcn(structure, cn, cnmax):
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

def calculate_cluster_size(structure):
    """Calculate the cluster size based on the maximum distance between metal atoms."""
    metal_positions = [atom.position for atom in structure if atom.symbol == metal]
    max_distance = 0
    for i in range(len(metal_positions)):
        for j in range(i + 1, len(metal_positions)):
            distance = np.linalg.norm(metal_positions[i] - metal_positions[j])
            if distance > max_distance:
                max_distance = distance
    return max_distance


def get_CN_GA(path, mult=0.9):
    atoms_in = read(path + '/POSCAR')
    filtered_atoms = [atom for atom in atoms_in if atom.symbol not in ['N', 'H']]
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
        key = site + 1 
        GA_dict[key] = get_one_top_site(site)

    for sites in bridge_sites:
        key = '_'.join([str(i+1) for i in sites])
        GA_dict[key] = get_one_bridge_site(sites)

    for sites in hollow_sites:
        key = '_'.join([str(i+1) for i in sites])
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
        
    exposed_top_sites = [i + 1 for i in exposed_top_sites]
    return exposed_top_sites

path = '.'
top_sites = get_CN_GA(path, mult=0.9)
print(top_sites)