#!/usr/bin/env python3 

from lattice import * 
import sys 
from gpm import *

#file_in = sys.argv[1]
file_in = 'POSCAR_deleted'

lines_f, dict_1, dict_2 = read_car(file_in)

atom_s = ['C', 'H', 'O']
atom_list = get_atom_list(lines_f, atom_s)

for key, val in dict_2.items():
    if 'O' in key:
        anchor_atom = int(val[0])

'''Correct the periodict problem and return the coordination of the atom list ''' 
coords_intact = get_intact_molecule(lines_f, atom_list, anchor_atom)

coords = [] 

for num, i in enumerate(atom_list):
    ele = get_ele_name(lines_f, i)
    infor_ele = ele + ' ' + ' '.join([str(k) for k in list(coords_intact[num])])
    coords.append(infor_ele)
#
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
