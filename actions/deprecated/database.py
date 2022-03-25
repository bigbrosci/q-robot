#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
Created on Sun Sep  1 12:04:15 2019

@author: qli
"""
import sys, os
import random, string


def save_out(name, lines):
    out = open(name, 'w')
    out.writelines(lines)
    out.close()

def analyze_summary_out(data_file='/home/qli/ccei/database/data/summary.out'): 
    ''' Step 1 Analyze the summary.out file and separate the good and bad data: remove the duplicate lines, uncomplete data, and initially separate the data'''
    f_in = open(data_file, 'r')
    lines = f_in.readlines()
    f_in.close()
    lines = list(dict.fromkeys(lines))  
    lines_normal = []
    lines_nonstable = []
    lines_problems  = []
    for line in lines:
        infor = line.strip().split()
        if len(infor) != 5:
            lines_problems.append(line)
        elif 'NonStable' in line:
            lines_nonstable.append(line)
        else:
            lines_normal.append(line)

#    name_normal = data_path + 'step1_normal.txt'
#    name_nonstable = data_path + 'step1_nonstable.txt'
#    name_problems  = data_path + 'step1_problems.txt'
#    
#    save_out(name_normal, lines_normal)
#    save_out(name_nonstable, lines_nonstable)
#    save_out(name_problems, lines_problems)
    
    '''Step2: Analyze the normal data, find the multiple calculations, e.g one structure was calculated several times due to the various of isomers '''
    l_smiles = []
    dict_ns  = {}  # {Smiles: [(energy_1, line number_1), (energy_2, line number_2), ...]}
   
    for num, line in enumerate(lines_normal):
        '''Get two dictionaries '''
        infor = line.strip().split()
        energy = float(infor[2])
        smiles = infor[4]
        if smiles not in dict_ns.keys():
            dict_ns[smiles] = [[energy, num]]
#            dict_ns.update({smiles:[[energy, num]]})
        else:
            dict_ns[smiles].append([energy, num])

    lines_mul = []
    lines_uni = []
    lines_mul_s = []
    for key, val in dict_ns.items():
        if len(val) == 1:
            lines_uni.append(lines_normal[val[0][1]])        
        else:
            for i in val:
                lines_mul.append(lines_normal[i[1]])
            val.sort()
            lines_mul_s.append(lines_normal[val[0][1]]) 
    
    '''Save the output files'''
    lines_2 = lines_uni + lines_mul_s

#    name_uni = data_path + 'step2_uni.out'
#    name_mul = data_path + 'step2_mul.out'
#    name_mul_s = data_path + 'step2_mul_s.out'
#    name_final = data_path + 'final.out'
#    save_out(name_mul, lines_mul)
#    save_out(name_uni, lines_uni)
#    save_out(name_mul_s, lines_mul_s)
    save_out('step2_final.out', lines_2)
    

    '''Step3: Analyze the Alpha-C radicals '''
    lines = lines_2[:]
    lines_d = []
    lines_f = []
    smiles_f = []
    for line in lines:
        if '=' in line:
            lines_d.append(line)
        else:
            lines_f.append(line)
            smiles_f.append(line.strip().split()[4] + '\n')

    save_out('alpha_radical.txt', lines_d)
    save_out('step3_final.out', lines_f)
    save_out('step3_smiles.out', smiles_f)

