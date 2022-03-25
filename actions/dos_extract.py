#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np 
import sys 
from lattice import *
from read_doscar import * 

arg_l = sys.argv[:]
input_list = arg_l[1:-1]
if '.dat' not in arg_l[-1].strip():
    out_name = 'dos_temp.dat'
else:
    out_name = arg_l[-1].strip()

lines = read_car('POSCAR')[0]
atom_list_raw, orbital_list_raw = get_atom_orbital_list(lines, input_list)
atom_list = get_atom_list(lines, atom_list_raw)
orbital_list = get_orbital_list(orbital_list_raw)

lines_dos, information = read_doscar()
natoms, E_max, E_min, nedos, e_fermi, ISPIN  = information[:]
energy = get_energy_list(lines_dos, information)

#### Generate the sum dos data of orbitals for the atoms we selected ###
spin_1 = np.zeros(nedos)
spin_2 = np.zeros(nedos)
for i in atom_list:
    for j in orbital_list:
        if ISPIN == 2 : 
            spin_1 = spin_1  +  get_single_orbital(i,j, lines_dos, information)[0]
            spin_2 = spin_2  +  get_single_orbital(i,j, lines_dos, information)[1]
        elif ISPIN == 1: 
            spin_1 = spin_1  +  get_single_orbital(i,j, lines_dos, information)[0]

#### Write the out-put file ###
f_out = open(out_name, 'w')
for i in range(0, len(spin_1)):
    if ISPIN == 2 :
        f_out.write('%10.5f %10.5f %10.5f \n' %(energy[i], spin_1[i], spin_2[i]))
    elif  ISPIN == 1 :
        f_out.write('%10.5f %10.5f \n' %(energy[i], spin_1[i]))
f_out.close()
#
