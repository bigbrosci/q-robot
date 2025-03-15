import os, sys
import csv
import json
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib.patches import Patch
import matplotlib.patches as patches
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.collections import PatchCollection
from matplotlib.colors import ListedColormap
from matplotlib import cm
from matplotlib.colors import ListedColormap


from ase import Atoms
from ase.io import read, write
from ase.neighborlist import NeighborList, natural_cutoffs, neighbor_list

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

######### Section 0: Set the basic parameters.
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

"""Energy of Gas phase species """
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


######### Section-1: Analyze the geometry and generate the Matrix for GA and other ML models.

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


def get_CN_GA_cluster(path,metal='Ru', mult=0.9):
    
    path = convert_unix_to_windows_path(path)
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

    
def get_species_from_cwd(cwd):
    # cwd = os.getcwd()
    # print(cwd)
    species = None
    try:
        path_parts = cwd.replace(' ', '').split(os.path.sep)
        for part in path_parts:
            if 'ads' in part:
                species = part.replace('_ads', '')
                break
    except Exception as e:
        print(f"An error occurred: {e}")
    
    if species:
        return species
    else:
        print("No species found in the current working directory.")
    
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

   

def get_GA_X_Y(data_path, cluster_path):
    """NH3 gas phase energy is the reference""" 
    
    data_path = convert_unix_to_windows_path(data_path)
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
    GA_file = convert_unix_to_windows_path(GA_file)
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



######## Section 2 ML prediction
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
    '''Update the hyperameters '''
    models = {
        "Ridge": Ridge(alpha=1.0),
        "DT": DecisionTreeRegressor(max_depth=n_group),
        "RF": RandomForestRegressor(n_estimators=n_group*5, max_depth=n_group*5),
        "GB": GradientBoostingRegressor(n_estimators=n_group*5, max_depth=n_group),
        "GA": PseudoInverseLinearRegression(),
        "XGB": XGBRegressor(n_estimators=11*5, max_depth=11, learning_rate=0.1),
        }
    return models


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
    
    cwd = os.getcwd() 
    species = get_species_from_cwd(cwd)
    text = '%s adsorption' %(species)
    # ax.set_title(text)
    ax.scatter(Y, Y_pred, edgecolor='k', alpha=0.7, s=200, label = text, color='#1f77b4')
    ax.plot([min(Y), max(Y)], [min(Y), max(Y)], 'r--', linewidth=3)
    ax.set_xlabel('DFT Eads (eV)', fontsize=36)
    ax.set_ylabel('Predicted Eads (eV)', fontsize=36)

    textstr = f'MAE: {mae:.2f}  \nRMSE: {rmse:.2f}  \nR2: {r2:.2f}'
    # plt.gca().text(0.65, 0.25, textstr, transform=plt.gca().transAxes,
    #                 fontsize=20, verticalalignment='top', bbox=dict(facecolor='white', alpha=0))
    plt.legend(loc='upper left', frameon=False, fontsize=28)
    plt.grid(False)
    ax = plt.gca()
    for spine in ax.spines.values():
        spine.set_linewidth(3)
        
    #plt.title('Linear Regression of De vs Ea')
    ax.tick_params(axis='both', which='major', width=2, length=8, labelsize=28)
    # ax.xticks(fontsize=30)
    # ax.yticks(fontsize=30)
    plt.tight_layout()
    # Save the figure
    plt.savefig('results/XY_%s.png' %(model_name), dpi=600)
    return textstr, Y, Y_pred
    
def bootstrap_cross_validation_one_model(models, model_name, X, Y, n_bootstrap=100):
    
    model = models[model_name]
    
    results = {metric: {model_name: [] } for metric in ['MAE', 'MSE', 'RMSE', 'R2', 'MPD', 'MND']}

    n_samples = int(len(X) * 0.80)  # 80% of the data for training
    
    for i in range(n_bootstrap):
        print('Iteration: ', i)
        indices = np.random.choice(len(X), size=len(X), replace=True)
        X_train, Y_train = X[indices[:n_samples]], Y[indices[:n_samples]]
        X_test, Y_test = X[indices[n_samples:]], Y[indices[n_samples:]]

        # print(name, 'running')
        model.fit(X_train, Y_train)
        predictions = model.predict(X_test)
        
        # Store each metric for this iteration
        results['MAE'][model_name].append(mean_absolute_error(Y_test, predictions))
        results['MSE'][model_name].append(mean_squared_error(Y_test, predictions))
        results['RMSE'][model_name].append(np.sqrt(mean_squared_error(Y_test, predictions)))
        results['R2'][model_name].append(r2_score(Y_test, predictions))
        results['MPD'][model_name].append(np.max(predictions - Y_test))
        results['MND'][model_name].append(np.min(predictions - Y_test))
   

    ### Step 2: Save results to CSV files
    for metric, model_results in results.items():
        df = pd.DataFrame(model_results)
        averages = df.mean()
        df.to_csv(f'results/{metric.lower()}.csv', index=False)
        averages.to_csv(f'results/{metric.lower()}_averages.csv')

    return results



df_GA, X, Y,  n_group  = get_GA_X_Y()
models  = refine_models(n_group)
test_one_model(X, Y, models, 'GA')


#### Perform the bootstrap cross-validation for One ML approaches

results = bootstrap_cross_validation_one_model(models, 'GA', X, Y, n_bootstrap=101)
# print (results)
df_mae = pd.read_csv('results/mae.csv')
# plt.figure(figsize=(10, 8))
fig, ax  = plt.subplots(figsize=(10,8))  
average_values = df_mae.mean()
average_values = df_mae
cwd = os.getcwd()
species  = get_species_from_cwd(cwd)
average_values.plot(kind='bar', label = species)

plt.xlabel('Number of CV test')
plt.ylabel('MAE / eV')
# plt.title('Bar Chart of MAE Values')
# plt.xticks(rotation=45)
plt.xticks(ticks=range(0, len(df_mae) , 10), rotation=0)  # Set xticks every 10th label
plt.xticks(fontsize=22)
plt.yticks(fontsize=22)
plt.tight_layout()
plt.legend(framealpha=0, fontsize = 18)
plt.grid(False)
ax = plt.gca()
for spine in ax.spines.values():
    spine.set_linewidth(1.5)
#plt.title('Linear Regression of De vs Ea')
ax.tick_params(axis='both', which='major', width=2, length=8, labelsize=18)
plt.legend(loc='upper right', frameon=False, fontsize=22)
plt.savefig('results/mae.png', dpi=600)