#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''https://wiki.fysik.dtu.dk/ase/ase/build/surface.html '''

import ase.io 
from ase.build import surface
from ase.build import fcc111, bcc110, hcp0001
from ase.constraints import FixAtoms
import sys



bcc = ['V',  'Cr', 'Mn', 'Fe', 'Nb', 'Pb']
hcp = ['Mg', 'Sc', 'Ti', 'Co', 'Zn', 'Y', 'Zr', 'Tc', 'Ru', 'Cd', 'Hf', 'Re', 'Os']
fcc = ['Al', 'Ca', 'Ni', 'Cu', 'Rh', 'Pd', 'Ag', 'Ir', 'Pt', 'Au']


bulk = ase.io.read("CONTCAR")
lattice_a = bulk.cell[0][0]
lattice_c = bulk.cell[2][2]


def cleave_surfaces(metal):
    '''Cleave the most stable surfaces for different close-packed metals '''
    if metal in bcc:
        slab = bcc110(metal, a = lattice_a, size = (1, 1, 4), vacuum = 7.5)
    elif metal in fcc:
        slab = fcc111(metal, a = lattice_a, size = (1, 1, 4), vacuum = 7.5)
    elif metal in hcp:
        slab = hcp0001(metal, a = lattice_a, c = lattice_c, size = (1, 1, 4), vacuum = 7.5)
    else:
        '''try to cleave the surfaces for alloy metals '''
        slab = surface(bulk, (1,1,1), 4)
        slab.center(vacuum=7.5, axis=2)

    '''Shif the atoms from center to bottom '''
    positions = slab.get_positions()
    positions[:,-1] = positions[:,-1] -  min(positions[:,-1])
    slab.positions = positions

    '''fix the bottom two layer metals to mimic the bulk '''
    highest_z = max(positions[:,-1])
    constraint_z = FixAtoms(indices=[atom.index for atom in slab if atom.z <= highest_z / 2])
    slab.set_constraint(constraint_z)

    ase.io.write('POSCAR', slab, format='vasp')  

metal = sys.argv[1]
cleave_surfaces(metal)
