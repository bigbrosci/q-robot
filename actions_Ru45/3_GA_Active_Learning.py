from sys import argv
import os,sys
import platform
import re

if platform.system() == 'Windows':
    q_robot_path = r'C:\Users\lqlhz\OneDrive - UMass Lowell\bin\q-robot\brain'
elif platform.system() == 'Linux':
    q_robot_path = '/home/qli/bin/q-robot/brain'
else:
    raise EnvironmentError("Unsupported operating system.")

sys.path.append(q_robot_path)

from cluster import *

metal = 'Ru'
cluster_path = '../slab'
species = get_species_from_cwd()  # Customize this function as 
#list_EF  = [str(x) for x in [-0.7, -0.5,  -0.3,  -0.1, 0.0, 0.1, 0.3, 0.5, 0.7]]
#list_EF  = [str(x) for x in [-0.7,-0.6, -0.5, -0.4, -0.3, -0.2, -0.1, 0.0, 0.1, 0.2, 0.3, 0.4, 0.5,0.6, 0.7]]
list_EF  = [str(x) for x in [-0.6,-0.4,-0.2, 0.0,  0.2, 0.4, 0.6]]
generate_GA_matrix_species(cluster_path, list_EF)
csv_file='GA_matrix.csv'
update_GA_matrix_with_dipole_polarization(csv_file, list_EF)


for EF in list_EF:
    df_GA, X, Y, n_group = get_GA_matrix(cluster_path, EF, list_EF)
    models  = refine_models(n_group)
    sites = df_GA['site'].values
    active_learning_one_model(models, 'GA', EF, X, Y, sites,
                                    initial_ratio=0.45, increment=0.1, max_ratio=0.95,
                                    results_dir='results')

# update_GA_matrix_with_dipole_polarization()
# list_EF.append('polarizability')
# list_EF.append('dipole')
# list_EF.append('c')

