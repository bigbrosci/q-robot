#!/usr/bin/env python3 
import sys
import ase
from  ase.io import read
from  ase.io import write

from optparse import OptionParser
parser = OptionParser()


import argparse

PARSER = argparse.ArgumentParser(description='Apply a vector to a molecule.')
PARSER.add_argument("-A", type=str, default=None,
                    help="atoms in File_1 will be replaced")
PARSER.add_argument("-B", type=str, default=None,
                    help="Atoms in File_2 will be used for replacement")
PARSER.add_argument("-s", nargs = '+',
                    help="Selected Atoms to be replaced")
PARSER.add_argument("-f", nargs = '+',
                    help="Selected Atoms to be used for replacement")
ARGS = PARSER.parse_args()




file_1 = ARGS.A
file_2 = ARGS.B

atoms_1 = [int(i)-1 for i in ARGS.s]
atoms_2 = [int(i)-1 for i in ARGS.f]


model_1 = ase.io.read(file_1, format='vasp')
model_2 = ase.io.read(file_2, format='vasp')

positions_1 = model_1.get_positions()
positions_2 = model_2.get_positions()


for num, atom in enumerate(atoms_1):
    A = atom 
    B = atoms_2[num]
    positions_1[A] = positions_2[B]

model_1.positions = positions_1

ase.io.write('POSCAR', model_1, format='vasp', vasp5=True)
