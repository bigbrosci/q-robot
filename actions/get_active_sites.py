#!/usr/bin/env python3
""" 
Core ides in determining the active sites:
    1) identify the exposed sites, can be achieved either via CN or manual check
    2) add the N to the bridge sites formed by the exposed sites (B,C)
    3) Identify the rhombus site (A,B,C,D) starting from the available bridge sites (B,C) 
        i) dis(BC)< dis(AD)
        ii) the rhombus is formed by ABC and BCD triangles.
    4) there are two types of rhombus sites: planar on the facet and edge on the intersection of two facets.
"""
import sys
import os 
q_robot_path = r'C:\Users\lqlhz\OneDrive - UMass Lowell\bin\q-robot\brain'
if q_robot_path not in sys.path:
    sys.path.append(q_robot_path)

from cluster import *


#########################

path = os.getcwd()
atoms = read(path + '/POSCAR')
mult = 0.9 
bond_length_threshold = 2.7
metal = 'Ru'

############## Step1: obtain the bridge sites
list_file = list_file = os.path.join(path, 'list')   ### list all the top sites in one line
if not os.path.exists(list_file):
    print("Warning: No list file found. Examine the exposed sites using Coordination Numbers.")
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
        print('No active site is found for: ', site)        
        
l_rhombus = [item for item in l_rhombus if item] ## remove the empty 


############## Step3: categorize  rhombus sites: planar or edge
coplanar_threshold = 0.5  # Adjust threshold as needed, unit is Å
coplanar_rhombus = filter_coplanar_rhombuses(l_rhombus, atoms, coplanar_threshold)
edge_rhombus = filter_rhombuses_by_dihedral(l_rhombus, atoms, 40, 120)

write_indices_to_file(os.path.join(path, 'all_sites.txt'), l_rhombus)
write_indices_to_file(os.path.join(path, 'planar_sites.txt'), coplanar_rhombus)
write_indices_to_file(os.path.join(path, 'edge_sites.txt'), edge_rhombus)            
