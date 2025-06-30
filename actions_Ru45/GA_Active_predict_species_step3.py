from sys import argv
import os,sys
import platform
import re

if platform.system() == 'Windows':
    q_robot_path = r'C:\Users\lqlhz\OneDrive - UMass Lowell\bin\q-robot\brain'
elif platform.system() == 'Linux':
    q_robot_path = '/home/qiangli/bin/q-robot/brain'
else:
    raise EnvironmentError("Unsupported operating system.")
sys.path.append(q_robot_path)

from cluster import *

metal = 'Ru'
path = '.'
species = get_species_from_cwd()  # Customize this function as 
cluster_path = './slab'

list_EF  = [str(x) for x in [-0.6,-0.4,-0.2,0.0,0.2,0.4,0.6]]
### Scan all the EF
for EF in list_EF:
    df_GA, X, Y, n_group = get_GA_matrix(cluster_path, EF)
    models  = refine_models(n_group)
    sites = df_GA['site'].values
    run_or_plot_active_learning(models, 'GA', EF, X, Y, sites,
                                    initial_ratio=0.4, increment=0.05, max_ratio=0.95,
                                    results_dir='results')


df_GA, X, Y, n_group = get_GA_matrix(cluster_path, '0.0') 
models  = refine_models(n_group)
sites = df_GA['site'].values

EF = '0.0'
run_or_plot_active_learning(models, 'GA', EF, X, Y, sites,
                                initial_ratio=0.4, increment=0.05, max_ratio=0.95,
                                results_dir='results')


df_EF = get_dipole_polarization(df_GA)

Y_polar = df_EF['polarizability'].values
Y_dipole = df_EF['dipole'].values

run_or_plot_active_learning(models, 'GA', 'polarizability', X, Y_polar, sites,
                                initial_ratio=0.4, increment=0.05, max_ratio=0.95,
                                results_dir='results')

run_or_plot_active_learning(models, 'GA', 'dipole', X, Y_dipole, sites,
                                initial_ratio=0.4, increment=0.05, max_ratio=0.95,
                                results_dir='results')
