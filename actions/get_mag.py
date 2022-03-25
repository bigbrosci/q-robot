#!/usr/bin/env python3
import sys 
import os
from outcar import *
from lattice import * 

atom_s = sys.argv[1:]

if os.path.isfile('OUTCAR'):
    f = open('OUTCAR')
    lines_o = f.readlines()
    f.close()
else:
    print('No OUTCAR in Current Path. Bye!')
    exit()    

lines = read_car('POSCAR')[0]
atom_list = get_atom_list(lines, atom_s)

dict_mag = get_mag(lines_o)
for atom in atom_list:
    print(atom, '\t', '\t'.join([str(i) for i in dict_mag.get(atom)]))

f_out = open('Magnetization.txt', 'w')    
for k, v in dict_mag.items():
    f_out.write('%s:\t %s\n' %(k, '\t'.join([str(i) for i in v])))
f_out.close()    

