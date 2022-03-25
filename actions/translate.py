#!/usr/bin/env python3

import numpy as np  
import sys 
import math
from lattice import * 
import argparse

PARSER = argparse.ArgumentParser(description='Apply a vector to a molecule.')
PARSER.add_argument("-x", type=float, default=None,
                    help="Distance vector (in x axis) to apply to molecule.")
PARSER.add_argument("-y", type=float, default=None,
                    help="Distance vector (in y axis) to apply to molecule.")
PARSER.add_argument("-z", type=float, default=None,
                    help="Distance vector (in z axis) to apply to molecule.")
PARSER.add_argument("-a", nargs = 2, default=None,
                    help="Point_A in A-->B vector: center of the two selected atoms")
PARSER.add_argument("-b", nargs = 2, default=None,
                    help="Point_B in A-->B vector: center of the two selected atoms")
PARSER.add_argument("-s", nargs = '+',
                    help="Selected Atoms to move")
ARGS = PARSER.parse_args()
#

lines = read_car('POSCAR')[0]

T = [0,0,0]
if ARGS.a and ARGS.b: 
    A1, A2 = ARGS.a
    A3, A4 = ARGS.b
    T = get_vector_T(lines, A1, A2, A3, A4)
    if ARGS.z:
        T[2] = ARGS.z
    
if ARGS.x:
    print(ARGS.x)
    T[0] = ARGS.x
if ARGS.y:
    print(ARGS.y)
    T[1] = ARGS.y
if ARGS.z:
    print(ARGS.z)
    T[2] = ARGS.z

atom_s = ARGS.s

is_direct = is_direct_or_not(lines)[0]
if is_direct:
    print('Error: Find Direct Coordinates. \n\n Please use command: dire2cart.py to convert the coordinates in POSCAR \n\n and then run this again.')
else:           
    atom_list = get_atom_list(lines, atom_s)
    lines_translated = shift_atoms(lines, atom_list, T)
     
    f_out = open('POSCAR_translated', 'w')
    for i in lines_translated:
        f_out.writelines(i)
    f_out.close()

print('\nThe translated file is named as:  POSCAR_translated\n')
