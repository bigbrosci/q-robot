#!/usr/bin/env python3
# -*- coding: utf-8 -*-
from lattice import *
from data import bcc, fcc, hcp
from data import *
from ase import Atoms
import ase.io 
from ase.build import molecule
from ase.build import bulk 
from ase.build import surface
from ase.build import add_vacuum 
from ase.build import fcc111, bcc110, hcp0001
from ase.constraints import FixAtoms
import subprocess
import sys, os

#data_dict = dict_metals_novdW   ## {element:[E, Natom, lattice_a, lattice_c]

#vdW_type = sys.argv[1]
vdW_list=['vdWNo','vdWD2','vdWD3zero','vdWD3BJ','vdWoptB86b','vdWoptB88','vdWoptPBE','vdWDF','vdWDF2','vdWrevDF2','vdWSCAN','vdWTS2','vdWTS21','vdWTS4']
print('\n' * 2, 'vdW Types: ', '  '.join(vdW_list), '\n'*2)
vdW_type = input("Choos one vdW type from above >>>  ")

dict_type = 'dict_metals_' + vdW_type
data_dict = eval(dict_type)   ## {element:[E, Natom, lattice_a, lattice_c]

def cssm(metal, data_dict):  # cleave_stable_surfaces_from_metals 
    name = 'POSCAR_' + metal
    if metal in bcc:   # For bcc metals, cleave 110 surface 
        lattice_a = float(data_dict.get(metal)[2])
        for i in range(1,6):
            name_out = name + '_' + str(i)
            slab = bcc110(metal, a = lattice_a, size = (i, i, 4), vacuum = 7.5)
            constraint_l = FixAtoms(indices=[atom.index for atom in slab if atom.index < i*i*2])
            slab.set_constraint(constraint_l)
            ase.io.write(name_out, slab, format='vasp')    
            ### Add the element line to the POSCAR file ###
            subprocess.call(['sed -i ' + '\'5a' + metal + '\'  ' + name_out], shell = True  )
            bottom(name_out)            
    elif metal in hcp:   # For hcp metals, cleave 0001 surface 
        lattice_a, lattice_c = [float(i) for i in data_dict.get(metal)[2:]]
        for i in range(1,6):
            name_out = name + '_' + str(i)
            slab = hcp0001(metal, a = lattice_a, c = lattice_c, size = (i, i, 4), vacuum = 7.5)
            constraint_l = FixAtoms(indices=[atom.index for atom in slab if atom.index < i*i*2])
            slab.set_constraint(constraint_l)
            ase.io.write(name_out, slab, format='vasp')
            subprocess.call(['sed -i ' + '\'5a' + metal + '\'  ' + name_out], shell = True  )
            bottom(name_out)            
            
    elif metal in fcc:   # For fcc metals, cleave 111 surface
        lattice_a = float(data_dict.get(metal)[2])
        for i in range(1,6):
            name_out = name + '_' + str(i)
            slab = fcc111(metal, a = lattice_a, size = (i, i, 4), vacuum = 7.5)
#            slab.center(vacuum=7.5, axis = 2)
            constraint_l = FixAtoms(indices=[atom.index for atom in slab if atom.index < i*i*2])
            slab.set_constraint(constraint_l)
            ase.io.write(name_out, slab, format='vasp')
            subprocess.call(['sed -i ' + '\'5a' + metal + '\'  ' + name_out], shell = True  )
            bottom(name_out)            
    else: 
        print(metal)
        print('Please add your element in the crystal structure lists: bcc, hcp, and fcc')  
        
for metal in data_dict.keys():
    cssm(metal, data_dict)
