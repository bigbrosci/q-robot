#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
This script is used to copy the atoms (molecule) from file_from (template) to the file_to.
You can use this script for the following cases:
1) optimize the adsorbates on the smaller slab with 2 layers, and then move the atoms to the larger slab with 4 layers.    
2) optimize the adsorbates on Pt(111) surface, and copy the atoms to Pd(111) surface. 
3) the atoms are copied to the end of the file_to, No merging for duplicated elements.
4) if you want to merge the duplicate elements, use add.py action instead.
    
Remember:
The coordinates of the copied atoms need to be modified when the surface in two slabs are different.    
'''

import sys 
from lattice import *

if len(sys.argv[:]) <= 3:
    print('Command: move_atoms.py file_from, file_to, atoms')
    print('Example: to move C H and 12th atom from POSCAR1 to POSCAR2: \n')
    print('Command: copy.py POSCAR1 POSCAR2 C H 12')
    exit()
else:     
    script, file_from, file_to = sys.argv[0:3]
    atom_s = sys.argv[3:]

lines_f = read_car(file_from)[0]
atom_list = get_atom_list(lines_f, atom_s)
#print(atom_list)
ele_list, ele_num, line_atoms, coord_s = get_selected_lines(lines_f, atom_list)
#
lines_t = read_car(file_to)[0]
total_lines_t = sum([int(i) for i in lines_t[6].strip().split()]) + 9
#
out_name = 'POSCAR_move'
f_out = open(out_name, 'w')
f_out.writelines(lines_t[0:5])  ## Write head part
### Write element 
f_out.write('%s  %s \n'  %(' '.join(lines_t[5].split()), ' '.join(ele_list)))  
### Write Numbers 
f_out.write('%s  %s \n'  %(' '.join(lines_t[6].split()), ' '.join([str(i) for i in ele_num])))  
### Write Coordinates part of file_to
f_out.writelines(lines_t[7: total_lines_t])
f_out.writelines(line_atoms[:])
#
f_out.close()

print('\nThe output file is:\t%s' %(out_name))
