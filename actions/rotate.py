#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# Author : Qiang 
# 14-01-2019
import sys 
from lattice import * 

'''
step1: use two atom to definle the axis
step2: get the list for rotated atoms.
step3: find anchor atom to make the molecule intact, in other words, solving the 
periodic boundary condition problem.
step4: rotate the atoms in the list

Along the axis(atom_A, atom_B), rotate atoms by theta angle (in degree).
Command: rotate.py atom_A atom_B atom_1 atom_2 ... Angle    
For example: to rotate the C H O atoms 15 degree along the axis formed by No.2 and No.10 atoms.
rotate.py 2 10 C H O 15 
'''
raw_infor = sys.argv[1:]
if len(raw_infor) < 4:
    print_rotate()
    exit()
    
lines = read_car('POSCAR')[0]
infor = get_rotate_infor(lines, raw_infor)
lines_rotated = get_atoms_pbc_rot(lines, infor)
f_out = open('POSCAR_rotated', 'w')
for i in lines_rotated: 
    f_out.writelines(i) 
f_out.close()
print('\nThe output file is named as: POSCAR_rotated\n')
