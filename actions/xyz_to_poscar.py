#!/usr/bin/env python3 
# -*- coding: utf-8 -*-
"""
Use ase to convert xyz files to POSCAR
"""

import sys 
from ase.io import read, write
from ase import Atoms


xyz_in = sys.argv[1]


mol = read(xyz_in)
atoms = Atoms(mol,
              cell = [(40, 0.0, 0.),(0, 40.0, 0.),(0.0, 0.0, 40.0)],
              pbc = True)

## Another way
# atoms  = Atoms(mol)
# atoms.set_cell([(18, 0.0, 0.),(0, 19.0, 0.),(0.0, 0.0, 20.0)])
# atoms.set_pbc((True, True, True))

write('POSCAR', atoms, format='vasp')





