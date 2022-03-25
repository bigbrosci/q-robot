#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  9 15:22:21 2020

@author: qli
"""

#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May  9 15:06:03 2020

@author: qli
"""
import pandas as pd
import sys 

cas_dict = {}
smiles_dict = {}
#f_in = open('/home/qli/ccei/database/data_02_04/thermo_m062x.csv', 'r')


data_base = pd.read_csv('/home/qli/ccei/database/data_02_04/thermo_m062x.csv') #,index_col=(0))

log_base = data_base['log_file']
smiles_base = data_base['smiles']
E_zpe = data_base['E_zpe']
cas_dict = dict(zip(log_base,smiles_base))
smiles_dict = dict(zip(smiles_base, log_base))
E_zpe_dict=dict(zip(smiles_base,E_zpe))

#for line in f_in.readlines():
#    key, val = line.rstrip().split(',')
#    cas_dict[key] = val
#    smiles_dict[val] = key

look_for = sys.argv[1].strip()

if look_for in cas_dict.keys():
    print(cas_dict[look_for])
elif look_for in smiles_dict.keys():
    print(smiles_dict[look_for], E_zpe_dict[look_for])
else:
    print('%s is not in the database' %(look_for))   
