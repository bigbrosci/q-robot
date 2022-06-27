#!/usr/bin/env python3 
import sys
import ase
from  ase.io import read
from  ase.io import write
import numpy as np


A, B = [int(i)-1  for i in sys.argv[1:3]]

model = ase.io.read('POSCAR', format='vasp')

positions = model.get_positions().tolist()

positions[A], positions[B] = positions[B], positions[A]

model.positions = np.array(positions)

ase.io.write('POSCAR_swaped', model, format='vasp', vasp5=True)
