#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''
Modified on 5th, March, 2019
Written By Qiang on 8th, Nov, 2018
This script is used to prepare the POSCAR for XPS calculations.
More Information, please see: 
1) http://cms.mpi.univie.ac.at/vasp/vasp/ICORELEVEL_tag_core_level_shifts.html
2) our wiki: http://aliga.iciq.es/wiki/index.php/Xps.py 
'''

import sys 
from lattice import * 
from mouth import *

if len(sys.argv[:]) == 2:
    script, atom_s = sys.argv[:]
    atom_s = int(atom_s)  # _s stands for selected 
    line_num_s = atom_s + 8
else:
    print('Command To use this script: xps.py  atom_numer')

lines = read_car('POSCAR')[0]
line_total = sum([int(i) for i in lines[6].strip().split()]) + 9  # Total line of POSCAR 
ele_s = get_ele_name(lines, atom_s)    

lines_deleted = delete_one_atom(lines, atom_s)

### Save POSCAR_XPS file 
f_out = open('POSCAR_xps', 'w')
f_out.writelines(lines[0:5]) # Write head part
f_out.write(ele_s + ' ' + ' '.join(lines_deleted[5].strip().split()) + '\n') # Write the Element line 
# Write the Element Numbers line
f_out.write('1 ' + ' ')
f_out.write(' '.join(lines_deleted[6].strip().split()) + '\n')
f_out.writelines(lines[7:9]) # Write Selective, Cartesian part
f_out.writelines(lines[line_num_s])  # Write the selected atom coordinate
f_out.writelines(lines_deleted[9:])  # Write the rest coordinates   
f_out.close()

print_xps()
