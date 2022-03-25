#!/usr/bin/env python3

import os, sys
from data import save_dict_txt, eval_dict_txt


data_file = '/home/qli/ccei/database/data_11_04/04Nov2019.out'

#### Generate the simple dictionary
#dict_simple = '/home/qli/ccei/database/data_11_04/dict_simple'
#dict_smiles = {}
#with open(data_file) as f:
#    lines = f.readlines()
#    for num, line in enumerate(lines):
#        infor = line.rstrip().split(',')
#        smiles = infor[1]
#        dict_smiles.update({smiles:num})
#        
#save_dict_txt(dict_smiles, dict_simple) 
       

#### Eval the simple dictionary
dict_simple = '/home/qli/ccei/database/data_11_04/dict_simple.txt'
dict_simple = eval_dict_txt(dict_simple)
#print(dict_simple)

#### get the detail information of one SMILES
f = open(data_file, 'r')
lines = f.readlines()
f.close()

def is_in_or_not(smiles):
    '''To check is a smiles in the database or not '''
    smiles = smiles.strip()
    if smiles in dict_simple.keys():
        print('Yes')
    else:
        print('No')
        
def get_data_infor_from_smiles(smiles):
    dict_smiles = {}
    try:
        smiles_line_num = dict_simple.get(smiles)
        infor = lines[smiles_line_num].rstrip().split(',')
        dict_smiles.update({'log_file': infor[0]})
        dict_smiles.update({'smiles': infor[1]})
        dict_smiles.update({'formula': infor[2]})
        dict_smiles.update({'hbond': infor[3]})
        dict_smiles.update({'Stablity': infor[7]})
        dict_smiles.update({'E_ele': infor[8]})
        dict_smiles.update({'E_zpe': infor[9]})
        dict_smiles.update({'H_298': infor[10]})
        dict_smiles.update({'G_298': infor[11]})
        freq_list = infor[12:]
        freq_list[0]= freq_list[0].replace('[', '')
        freq_list[-1]= freq_list[-1].replace(']', '')
        freq_list = [float(i) for i in freq_list]
        dict_smiles.update({'frequencies': freq_list})
        print('%s is found in the database.' %(smiles))
        return dict_smiles
    except:
        print('%s is not in the database, please check it again.' %(smiles))

#get_data_infor_from_smiles('CCC([O])C(O)c1ccc(O)cc1')   
        
smiles = sys.argv[1]
is_in_or_not(smiles)
#get_data_infor_from_smiles(smiles)
 
