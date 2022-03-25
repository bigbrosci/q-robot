#!/usr/bin/env python3 
# -*- coding: utf-8 -*-
"""
"""

import os, sys
import shutil 
import pandas as pd
import subprocess
import csv
from data import *
from smiles2xyz import get_smiles_file


data=pd.read_csv('/home/qli/ccei/database/DFT_results/log/thermo_m062x.csv')
mols = list(data['smiles'])
logs = list(data['log_file'])


xyz_file = sys.argv[1]
smiles = get_smiles_file(xyz_file)

in_or_not = False
if smiles in mols:
    num = mols.index(smiles)
    in_or_not = True
    print('%s is alreay in the database as %s' %(xyz_file, logs[num])) 

if in_or_not == False:
    print('%s is not in the database' %(xyz_file))
