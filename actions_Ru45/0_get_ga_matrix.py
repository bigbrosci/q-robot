
from sys import argv
import os,sys

#q_robot_path = r'C:\Users\lqlhz\OneDrive - UMass Lowell\bin\q-robot\brain'
q_robot_path='/home/qli/bin/q-robot/brain'

if q_robot_path not in sys.path:
    sys.path.append(q_robot_path)
from cluster import *

metal = 'Ru'

cluster_path = './slab'
atoms_in = read(cluster_path+'/POSCAR')
connections, cn_of_connected_atoms, exposed_top_sites, bridge_sites, hollow_sites, square_sites = get_connection(atoms_in, metal='Ru', mult=0.9)
print(connections)
get_CN_GA(cluster_path, mult=0.9, metal='Ru')
get_full_GA_matrix(cluster_path)
