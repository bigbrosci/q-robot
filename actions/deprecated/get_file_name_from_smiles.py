#!/usr/bin/env python3
import os, sys
from smiles2xyz import * 
from os import path

xyz_file = sys.argv[1]
name = xyz_file.split('.')[0]
smiles, formula  = get_smiles_m1(xyz_file)
if '=' in smiles:
    smiles, formula = get_smiles_m2(xyz_file)

file_name = smiles.replace('(', 'L').replace(')', 'R').replace('[', 'l').replace(']','r').replace('=', 'e')
#print(xyz_file, '\t', formula, '\t', smiles)
print(file_name)


os.rename(xyz_file, file_name+'.xyz')

if path.exists(name+'.log'):
    os.rename(name+'.log', file_name+'.log')
if path.exists(name+'.gjf'):
    os.rename(name+'.gjf', file_name+'.gjf')
if path.exists(name+'.chk'):
    os.rename(name+'.chk', file_name+'.chk')
if path.exists(name+'_hg.txt'):
    os.rename(name+'_hg.txt', file_name+'_hg.txt')
