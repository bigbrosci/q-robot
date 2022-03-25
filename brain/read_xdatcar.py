#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Read XDATCAR files
@author: qli
"""
import numpy as np
import os 
def split_l(line):
    list_line = line.strip().split()
    return list_line

if os.path.isfile('./XDATCAR'):
    file_in = open('XDATCAR', 'r')
    lines_xd = file_in.readlines()
    file_in.close()
    system = lines_xd[0]    
    scale  = float(split_l(lines_xd[1])[0])
    la1 = np.array([float(i) for i in split_l(lines_xd[2])])
    la2 = np.array([float(i) for i in split_l(lines_xd[3])])
    la3 = np.array([float(i) for i in split_l(lines_xd[4])])
    vector   = np.transpose(np.array([la1, la2, la3]))
    elements = split_l(lines_xd[5])
    ele_num  = [int(i) for i in split_l(lines_xd[6])]
    ele_sum  = sum(ele_num)
    num_block = (len(lines_xd) - 7 ) / (ele_sum + 1)  # how many steps optimized in the calculation
else: 
    print('There is no XDATCAR file in path: %s' %(os.path.abspath('.')))
    exit

def get_block_num():
    num = 0
    look_up = 'Direct'
    for line in lines_xd:
        if look_up in line:
            num += 1 
    return num
        
def get_block_from_xdatcar(num):
    n_start = (num - 1) * (ele_sum + 1) + 7
    n_end   = (num) * (ele_sum + 1) + 7 
    list_block = lines_xd[n_start : n_end]
    return list_block

def get_atom_from_xdatcar(num_block, num_atom):
    list_block = get_block_from_xdatcar(num_block)
    return list_block[num_atom]

def get_atoms_from_xdatcar(num_block, num_list):
    lines_xd_atoms = []
    list_block = get_block_from_xdatcar(num_block)
    for i in num_list:
        lines_xd_atoms.append(list_block[i])        
    return lines_xd_atoms

