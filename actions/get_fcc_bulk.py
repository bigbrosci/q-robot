#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 18:49:32 2020

@author: qli
"""

""" 
https://wiki.fysik.dtu.dk/ase/ase/build/build.html?highlight=bulk
"""

import sys 
from ase.io import read, write
from ase import Atoms
from ase.build import bulk

metal, lattice_constant = sys.argv[1:3]

atoms = bulk(str(metal), 'fcc', a=float(lattice_constant), cubic = True)

write('POSCAR', atoms, format='vasp')