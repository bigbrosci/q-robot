#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import random
import string
import os, sys
from gpm import *
from smiles2xyz import *

data_file = '/home/qli/ccei/database/data/summary.out'
path_data_run = '/home/qli/ccei/database/data_add/mol/'


ring_out = sys.argv[1]
get_smiles_from_ring(ring_out)
analyze_ring_radical_smiles('rad_ring.out')
analyze_ring_molecule_smiles('mol_ring.out')

def check_lines(data_file):
    f_in = open(data_file, 'r')
    lines = f_in.readlines()
    f_in.close()
    done_list = []
    new_lines = lines.copy()
    for line in lines:
        infor = line.strip().split()
        if len(infor) != 5:
            print(line)
            new_lines.remove(line)
        else:
            smiles = infor[4].strip()
            done_list.append(smiles)
    return new_lines, done_list

lines, done_list = check_lines(data_file)

file_in = open('mol_final_rdkit.out', 'r')        
lines_mol = file_in.readlines()
file_in.close()
#
lack_out = open('lack.out', 'w')
done_out = open('done.out', 'w')
for line in lines_mol:
    smiles = line.rstrip()
    if smiles not in done_list:
        try: 
            smiles_to_gjf(smiles, path_data_run)
            lack_out.write('%s\n' %(smiles))
        except  ValueError:
            print(smiles)
    else:
        done_out.write('%s\n' %(smiles))
    
lack_out.close()
done_out.close()
