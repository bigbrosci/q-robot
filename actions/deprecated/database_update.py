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

import rdkit  
from rdkit.Chem import AllChem as Chem
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import rdMolDescriptors
from rdkit.Chem import Lipinski

from smiles2xyz import get_smiles_file


def append_data(dict_one):
    list_T = list(range(50,1001,20))
    item_list = ['log_file', 'smiles', 'formula', 'E_zpe', 'Hf_0', 'Hf_T', 'Cp', 'S']
    item_list_values = [str(dict_one.get(i)) for i in item_list]
   
    f = open('/home/qli/ccei/database/DFT_results/log/thermo_m062x.csv', 'a+') 
    f.write('%s' %(','.join(item_list_values)))    
    for T in list_T:
        cp = dict_one.get(T).get('Cp')
        s = dict_one.get(T).get('S')
        f.write(',%s,%s' %(str(cp), str(s)))
    f.write('\n')  
    f.close()

data=pd.read_csv('/home/qli/ccei/database/DFT_results/log/thermo_m062x.csv')
mols = list(data['smiles'])
logs = list(data['log_file'])
last_log = logs[-1]
last_num = int(last_log.replace('.log', '').replace('S', ''))


log_file = sys.argv[1]
command = 'get_log_infor.py  ' + log_file
subprocess.call(command, shell=True)

thermal_file = log_file.replace('.log', '_thermo.txt')
xyz_file = log_file.replace('log', 'xyz')
smiles = get_smiles_file(xyz_file)


in_or_not = False
if smiles in mols:
    num = mols.index(smiles)
    in_or_not = True
    print('%s is alreay in the database as %s' %(log_file, logs[num])) 

if in_or_not == False:
   
   dict_one = eval(open(thermal_file).read())
   if dict_one.get('Stability') == 'Stable' :  
       print('No Imaginary Frequencies were found.\n Adding it to the database \n')
       new_name = 'S' + str(last_num+1)  
       new_log = new_name + '.log'
       new_xyz = new_name + '.xyz'
       new_thermal = new_name + '_thermo.txt' 
       dict_one['log_file'] = new_log
       dict_one['name'] = new_name
       save_dict_txt(dict_one, new_name+'_thermo')
       shutil.move(log_file, new_log )
       shutil.move(xyz_file, new_xyz )
       shutil.move(thermal_file, new_thermal)
       append_data(dict_one)
       print('Done')
