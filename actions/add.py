#!/usr/bin/env python3 
# -*- coding: utf-8 -*-
'''
This script is used to copy atoms from one POSCAR to the  other POSCAR file.
1) the action is similar to move_atom.py
2) move_atoms.py only copy the atoms to the end of the POSCAR, this script can merge the duplicated atoms.
C H O N 
1 2 3 4 
If you add another O to the POSCAR, the output will be 
C H O N 
1 2 4 4 
The output of move_atoms.py will be :
C H O N O 
1 2 3 4 1    
'''
import sys
from lattice import *


if len(sys.argv[:]) < 3:
    print('\nCommand Usage: copy.py file_from file_to atoms_to_be_copied\n')
    print('Example: to Copy C H and 12th atom from POSCAR1 to POSCAR2: \n')
    print('Command: copy.py POSCAR1 POSCAR2 C H 12')
    exit()
else:    
    script, file_from, file_to = sys.argv[0:3]
    atom_s = sys.argv[3:]

lines_f = read_car(file_from)[0]
atom_list = get_atom_list(lines_f, atom_s)
ele_list, ele_num, line_atoms, coord_s = get_selected_lines(lines_f, atom_list)
#print(ele_list)
lines_t = read_car(file_to)[0]
#
out_name = 'POSCAR_add'
for i in coord_s:
    lines_t = add_one_atom(lines_t, i)
file_out = open(out_name, 'w')
file_out.writelines(lines_t)
file_out.close()

print('\nThe output file is:\t%s' %(out_name) )
