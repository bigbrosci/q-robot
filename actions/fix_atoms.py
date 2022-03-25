#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Written By Qiang on 17th, August, 2018
This script is used to fix and relax the atoms in POSCAR

for example:

1) fix the atoms 1-6 with T T F 
python fix_atoms.py 1-6 TTF
2) fix all C atoms in POSCAR with T F F 
python fix_atoms.py C TFF 
3) fix all C and O atoms with F F F:
python fix_atoms.py C O  FFF
4) relax the atoms: 1 3 4 9-12 with T T T
python fix_atoms.py 1 3 4 9-12 TTT

'''

from lattice import *
import sys
import os


if len(sys.argv) <= 2:
    print('\nCommand to use it: fix_atoms file atoms FFF \n')
    exit()
elif len(sys.argv) == 3:
    print("\nWe are going to fix %s......" %sys.argv[1]) 
    script, file_in, atom_s = sys.argv[:]
    tf_in = 'FFF'
else:
    print("\nWe are going to fix %s......" %sys.argv[1]) 
    script, file_in = sys.argv[:2]
    tf_in = sys.argv[-1]
    if len(tf_in) != 3 and 'T' not in tf_in and 'F' not in tf_in:
        atom_s = sys.argv[2:]
        tf_in  = 'FFF'
    else:
        atom_s = sys.argv[2:-1]
        
tf = list(tf_in)
      
lines = read_car(file_in)[0]
atom_list = get_atom_list(lines, atom_s)
line_list = [i+8 for i in atom_list]
print(atom_list)

def fix_one_line(line, tf_in):
    infor = line.strip().split()
    new_line = ' %s  %s\n' %('  '.join(infor[0:3]), ' '.join(tf_in))
    return new_line

out_name = file_in + '_fixed'
file_out = open(out_name, 'w')
for num, line in enumerate(lines):    
    if num in line_list:
        new_line = fix_one_line(line, tf_in)
        file_out.write(new_line)
    else:
        file_out.write(line)         
file_out.close()
print('\nThe output file is:\t %s' %(out_name))
