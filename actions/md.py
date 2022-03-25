#!/usr/bin/env python3
# -*- coding: utf-8 -*-
''' This script is used to 
1) check the steps that O--O bond breaks in the MD simulations of OOH on the metal surfaces
2) save the xy coordinates of the O atoms in trajectory.
3) vasprun.xml and POSCAR files will be read
# Written By Qiang on 9th-Jan-2019
'''
import sys, os
import numpy as np 
from  lattice  import *

script, atom1, atom2 = sys.argv
lines_pos = read_car('POSCAR')[0]
vector, la1, la2 = get_vectors(lines_pos)

def save_atoms(atom_list):
    atom_list = sorted(atom_list)
    file_out = []
    file_name = []
    break_steps = []
    for i in atom_list:
        out = 'file_out_' + str(i)
        file = 'data_' + str(i)
        file_out.append(out)  
        file_name.append(file)  
    for i in range(0,len(file_out)):
        file_out[i] = open(file_name[i], 'w')

    with open('vasprun.xml') as infile:
        lines = infile.readlines()
        infile.seek(0)
        count = 0
        for i, l in enumerate(infile):
            if '<varray name="positions" >' in l:
                count +=1 
                xyz = []
                def get_xyz(line): 
                    line_atom = lines[line]   ## the atom number you want to focus
                    atom_xyz   = np.array([ float(i) for i in line_atom.strip().split()[1:4] ])
                    x, y, z = [sum(k) for k in atom_xyz * vector ]   
                    return x, y, z
                for num, atom in enumerate(atom_list):
                    x, y, z = get_xyz(i + atom)
                    xyz.append(get_xyz(i + atom))
                    file_out[num].write('%s\t%s\t%s\n' %(x,y,z)) 
                a = np.array(xyz[0])
                b = np.array(xyz[1])
                if get_distance(lines_pos, a, b) > 1.7:
                    break_steps.append(count)
    for i in file_out:
        i.close()
    if len(break_steps) >= 1:
        print(break_steps[1])
    else:
        print('No')
    
atoms_list = [int(atom1),int(atom2)]
save_atoms(atoms_list)
