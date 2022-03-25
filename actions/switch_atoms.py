#!/usr/bin/env python3 
import sys
import ase
from  ase.io import read
from  ase.io import write

file_in = sys.argv[1]
atom_1 = sys.argv[2]
atom_2 = sys.argv[3]

model = ase.io.read(file_in, format='vasp')
positions = model.get_positions()
new_positions = positions.copy()

new_positions[atom_1] =  positions[atom_2] 
new_positions[atom_2] =  positions[atom_1] 

model.positions = new_positions
# model[114].positions = positions_2
# model[115].positions = positions_1 

ase.io.write('POSCAR_switched', model, format='vasp', vasp5=True)
