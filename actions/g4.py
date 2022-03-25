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

import sys 

cas_dict = {}
smiles_dict = {}
f_in = open('/home/qli/ccei/database/G4/g4.dict', 'r')
for line in f_in.readlines():
    key, val = line.rstrip().split(',')
    cas_dict[key] = val
    smiles_dict[val] = key

look_for = sys.argv[1].strip()

if look_for in cas_dict.keys():
    print(cas_dict[look_for])
elif look_for in smiles_dict.keys():
    print(smiles_dict[look_for])
else:
    print('%s is not in the database' %(look_for))   
