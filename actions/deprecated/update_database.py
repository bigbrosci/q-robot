#!/usr/bin/env python3
import os, sys, shutil
from smiles2xyz import *
from  xyz2mol import *
from os import path
from gpm import *
from openbabel import openbabel, pybel
from pybel import *
from data import save_dict_txt, eval_dict_txt

#log_file = sys.argv[1]
#name = log_file.split('.')[0]
#
#smiles, formula  = get_smiles_m1_log(log_file)
#if '=' in smiles:
#    smiles, formula  = get_smiles_m2_log(log_file)   
#
#file_name = smiles.replace('(', 'L').replace(')', 'R').replace('[', 'l').replace(']','r').replace('=', 'e')
#if path.exists(name+'.log'):
#    os.rename(name+'.log', file_name+'.log')
#if path.exists(name+'.gjf'):
#    os.rename(name+'.gjf', file_name+'.gjf')
#if path.exists(name+'.chk'):
#    os.rename(name+'.chk', file_name+'.chk')
#if path.exists(name+'_hg.txt'):
#    os.rename(name+'_hg.txt', file_name+'_hg.txt')

#print("%s has been converted to %s" %(log_file, file_name+'.log'))

def get_simple_dict():
    ### Generate the simple dictionary
    data_file = '/home/qli/ccei/database/data_11_04/04Nov2019.out'
    dict_simple_file = '/home/qli/ccei/database/data_11_04/dict_simple'
    dict_simple = {}
    with open(data_file) as f:
        lines = f.readlines()
        for num, line in enumerate(lines):
            infor = line.rstrip().split(',')
            smiles = infor[1]
            dict_simple.update({smiles:num})
    save_dict_txt(dict_simple, dict_simple_file) 
    return dict_simple
       

def get_simple_dict_from_file():
    #### Eval the simple dictionary
    dict_simple = '/home/qli/ccei/database/data_11_04/dict_simple.txt'
    dict_simple = eval_dict_txt(dict_simple)
    return dict_simple

def get_lines_from_db():
    #### get the detail information of one SMILES
    data_file = '/home/qli/ccei/database/data_11_04/04Nov2019.out'
    f = open(data_file, 'r')
    lines = f.readlines()
    f.close()
    return lines

def is_in_or_not(smiles):
    '''To check is a smiles in the database or not '''
    smiles = smiles.strip()
    dict_simple = get_simple_dict()
    if smiles in dict_simple.keys():
        return  True
    else:
        return  False
        
def get_data_infor_from_smiles(smiles):
    dict_smiles = {}
    dict_simple = get_simple_dict()
    lines = get_lines_from_db()
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
#        print('%s is found in the database.' %(smiles))
        return dict_smiles
    except:
        print('%s is not in the database, please check it again.' %(smiles))

db_path = '/home/qli/ccei/database/data_11_26/'
data_file = db_path+'data.txt'


def update_data_file(lines):
    f = open(data_file, 'w')
    f.writelines(lines)
    f.close()

lines_db = get_lines_from_db()
dict_simple = get_simple_dict()

dict_check_file = sys.argv[1]    
f = open(dict_check_file, 'r')
lines_check = f.readlines()
f.close()

for line in lines_check: 
    infor = line.rstrip().split(',')
    smiles = infor[1]
    E_zpe = float(infor[9])
    if is_in_or_not(smiles):
        dict_db = get_data_infor_from_smiles(smiles)
        E_zpe_db = float(dict_db.get('E_zpe'))
        if E_zpe < E_zpe_db:
            num_db = dict_simple.get(smiles)
            lines_db[num_db] = line  # update the database summary file
            update_data_file(lines_db)            
            file_db = db_path + dict_db.get('log_file')
            os.rename(file_db, file_db+'deleted') 
            shutil.move(infor[0], file_db) 
        else:
            print('%s is not updated to the database' %(infor[0]))
    else: 
        print(smiles)
        lines_db.append(line)
        update_data_file(lines_db)            
        name = smiles.replace('(', 'L').replace(')', 'R').replace('[', 'l').replace(']','r').replace('=', 'e')
        print(db_path+name+'.log')
        shutil.move(infor[0], db_path+name+'.log') 
