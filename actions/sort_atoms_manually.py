#!/usr/bin/env python3
''' 
Usage Command: python3 sort_atoms_manually.py file_read  ele1 ele2  ele3  2.0 
file_read: POSCAR or CONTCAR 
ele1, ele2, ele3: the elements in the POSCAR
2.0: All atoms in Z directions below 2.0 A will be fixed. 
'''
import sys
import ase
from ase.build.tools import sort
from  ase.io import read
from  ase.io import write
from ase import Atoms
from ase.constraints import FixAtoms
file_in = sys.argv[1]
sort_list = sys.argv[2:-1]
fix_z = float(sys.argv[-1])


#file_in = 'POSCAR'
#sort_list = ['Ni', 'C', 'H', 'O']
#fix_z = 1.0

model = ase.io.read(file_in, format='vasp')

l_total = [] # the sorted list of the atoms
for ele in sort_list:
    l_total.extend([i.index for i in model if i.symbol == ele])

model_sorted = Atoms(cell= model.cell, pbc=True)

for num in l_total:
    model_sorted.append(model[num])

c = FixAtoms(mask=model_sorted.positions[:, 2] < fix_z)
model_sorted.set_constraint(c)

## Sort atoms by the z coordinates.
# model_sorted = sort(model, tags=model.positions[:, 2])
## Sort atoms by the elements by default
# model_sorted = sort(model, tags=model.get_chemical_symbols())

op = 'POSCAR_sorted'
ase.io.write(op, model_sorted, format='vasp', vasp5=True)

print('Output file is:\t', op)
