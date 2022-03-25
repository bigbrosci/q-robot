#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This script is used to stwitch the two atoms in POSCAR for seggregation calculations.
@author: Qiang
"""
#from itertools import izip as zip # works in python2
from lattice import * 
import sys 

threshold = 0.5
try:
    lines, dict_car1, dict_car2 = read_car('POSCAR')
    layers = determinelayers(lines, threshold)
    ele_name  = lines[5].strip().split()
    ele_num = [int(i) for i in lines[6].strip().split()] 
    total_num = sum(ele_num)
    line_total = total_num + 9
except IOError:
    print('Can not find POSCAR, Read Check!')
    exit()    

####Let the user choose the layers that he/she intersted
def get_selected_layers():
    layer_max =  max(layers.keys())
    layer_min =  min(layers.keys())
    s_layer = [] 
    if len(sys.argv[:]) == 1:
        print('\nFind %s layers in your POSCAR by using %s \AA as creteria.\n' %(layer_max, threshold) )
        print( '\nPlese follow the rules below and input the TWO layers you intersted:'    )
        print( '1) The Bottom and top layers are: %s %s' %(layer_min, layer_max))
        print( '2) Select the 1st and 3rd layers: %s %s ' %(layer_max, layer_max-3) )
        print( '3) Select 1st layer twice, input_1 : %s %s' %(layer_max, layer_max) )    
        print( '4) Select 1st layer twice, input_2 : %s' %(layer_max) )
        print( '5) 1st layer will be selected by default if no inputs.\n'  )
        
        #layer_s = raw_input('Iuput the two layers you want to focus on >>>\t')     ## Python2
        layer_s = input('Iuput the two layers you want to focus on >>>\t')    
        
        if len(layer_s.split()) == 0 :
            print( '\nYou did not select any layers.' )
            s_layer = [layer_max, layer_max]
        elif len(layer_s.split()) == 1 :
            s_layer = [int(layer_s), int(layer_s)]
        elif len(layer_s.split()) == 2 :
            s_layer = [int(i) for i in layer_s.split()]
        else:
            print('\nERROR: Read the instructions above carefully and rerun this script!\n')
            exit()
    elif len(sys.argv[:]) == 2:
        layer_s = int(sys.argv[1])        
        s_layer = [int(layer_s), int(layer_s)]
    elif len(sys.argv[:]) == 3: 
        s_layer = [int(i) for i in sys.argv[1:]]
    else:
        print('Command: switch.py layer1, layer2')
        exit()
    return s_layer      # return a list containing the layers selected   

layer1, layer2 = get_selected_layers()
atom_list_layer1 = [int(i)  for i in layers.get(layer1) ]
atom_list_layer2 = [int(i)  for i in layers.get(layer2) ]
list_pairs_raw = get_pairs(atom_list_layer1, atom_list_layer2)

def reduce_pairs(pairs):
    pairs_reduced = pairs[:]
    for pair in pairs:
        num_1, num_2 = pair[:]
        atom_1, atom_2 = [get_ele_name(lines, i) for i in pair]
        if num_1 == num_2 or atom_1 == atom_2 :
            pairs_reduced.remove(pair)
    return pairs_reduced

list_pairs = reduce_pairs(list_pairs_raw)       
#print(list_pairs)
 
for pair in list_pairs:
    '''pair = [atom_1, atom_2], you are going to swithc their coordinates ''' 
    num_1, num_2 = pair[:]
    atom_1, atom_2 = [get_ele_name(lines, i) for i in pair]
    print('Switch \t %s-%s \t  %s-%s' %(atom_1, num_1, atom_2, num_2))
    name_pair = '-'.join([str(i) for i in pair])
    lines_pair = switch_atoms(lines, pair[0], pair[1])
    file_out = open('switch_'+name_pair, 'w')
    file_out.writelines(lines_pair)
    file_out.close()

