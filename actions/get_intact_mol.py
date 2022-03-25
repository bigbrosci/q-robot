#!/usr/bin/env python3 

from lattice import * 
import sys 
from lignin_ga import *

file_in = sys.argv[1]

lines_f, dict_1, dict_2 = read_car(file_in)

atom_s = ['C', 'H', 'O']
atom_list = get_atom_list(lines_f, atom_s)

for key, val in dict_2.items():
    if 'O' in key:
        anchor_atom = int(val[0])

if len(sys.argv[:]) == 3: 
    anchor_atom = int(sys.argv[2])

'''Correct the periodict problem and return the coordination of the atom list ''' 
coords_intact = get_intact_molecule(lines_f, atom_list, anchor_atom)

coord_list = []
mass_list = []

for num, i in enumerate(atom_list):
    ele = get_ele_name(lines_f, i)
    mass_list.append(atomic_mass.get(ele))
    coord_list.append(list(coords_intact[num]))
    
get_mass_center = lambda c,m:list(map(sum,zip(*[[i[0]*j/sum(m),i[1]*j/sum(m), i[2]*j/sum(m)]for i,j in zip(*([c,m]))])))
mass_center = np.array(get_mass_center(coord_list, mass_list))


la1 = np.array([ float(i) for i in lines_f[2].strip().split() ])
la2 = np.array([ float(i) for i in lines_f[3].strip().split() ])
la3 = np.array([ float(i) for i in lines_f[4].strip().split() ])
box_center = (la1+la2+la3)/2

trans_vector = np.array([10,10,10]) - mass_center
#trans_vector = box_center + np.array([20,20,20]) - mass_center

coords = [] 
for num, i in enumerate(atom_list):
    ele = get_ele_name(lines_f, i)
    coord_ele = coords_intact[num] + trans_vector
    infor_ele = ele + ' ' + ' '.join([str(k) for k in list(coord_ele)])
    coords.append(infor_ele)
    

dict_ele_coord = xyz_analyzer(coords)[0]
out_name = 'POSCAR_periodic' 
file_out = open(out_name, 'w')
file_out.write('Periodic_corrected\n1.0 \n')

dire_a = np.array([[20,0,0],[0, 20,0], [0,0,20]])

lattice_a = ' '.join([str(i) for i in dire_a[0]]) + '\n'
lattice_b = ' '.join([str(i) for i in dire_a[1]]) + '\n'
lattice_c = ' '.join([str(i) for i in dire_a[2]]) + '\n'

file_out.write('%s %s %s' %(lattice_a, lattice_b, lattice_c))

file_out.write(' '.join(dict_ele_coord.keys()) + '\n')
file_out.write(' '.join([str(len(i)) for i in dict_ele_coord.values()]) + '\n')
file_out.write('Selective dynamics\nCartesian\n')
for i in dict_ele_coord.values():
    for j in i:
        file_out.write(' '.join([str(k) for k in j]) + ' T  T  T\n')
file_out.close()
