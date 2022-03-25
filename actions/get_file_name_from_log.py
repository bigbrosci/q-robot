#!/usr/bin/env python3
import os, sys
from smiles2xyz import *
from  xyz2mol import *
from os import path
from gpm import *
import openbabel, pybel
from pybel import *
from lattice import *

log_file = sys.argv[1]
name = log_file.split('.')[0]


smiles, formula  = get_smiles_m1_log(log_file)
if '=' in smiles:
    smiles, formula  = get_smiles_m2_log(log_file)   
print(log_file, smiles, formula)
#
#file_name = smiles.replace('(', 'L').replace(')', 'R').replace('[', 'l').replace(']','r').replace('=', 'e')
#
#if path.exists(name+'.log'):
#    os.rename(name+'.log', file_name+'.log')
#if path.exists(name+'.gjf'):
#    os.rename(name+'.gjf', file_name+'.gjf')
#if path.exists(name+'.chk'):
#    os.rename(name+'.chk', file_name+'.chk')
#if path.exists(name+'_hg.txt'):
#    os.rename(name+'_hg.txt', file_name+'_hg.txt')
#
#print("%s has been converted to %s" %(log_file, file_name+'.log'))
