#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Convert direc coordiation to cartesian Writen By Qiang 
#version 1.3_Qiang add parser functions

import sys, os
import numpy as np 
#from mouth import print_dire_to_cart
from lattice import * 
import sys 

from optparse import OptionParser
parser = OptionParser()  

parser.add_option("-i",   
                  type = 'string', dest="file", default="POSCAR",  
                  help="File to be converted or fixed")  

parser.add_option("-s", 
                  type = "int", dest="selected", default=0,
                  help="select how many bottom layers to fix") 
parser.add_option("-t", 
                  type = 'float',dest="threshold", default=0.5,  
                  help="Threshold of distance to separate layers in z direction!") 
(options, args) = parser.parse_args()

threshold = options.threshold ## The distance in z direction between two layers   
fixedlayer = options.selected 

if os.path.isfile(options.file):
    file_to_be_converted = options.file
else:
    if os.path.isfile("POSCAR"):    
        file_to_be_converted = "POSCAR"
    elif os.path.isfile("CONTCAR"):
        file_to_be_converted = "CONTCAR"
    else: 
        print(u"\n Error, %s, or POSCAR or CONTCAR does not exist!" %(options.file))
        exit()

def get_infor(file_to_be_converted):
    lines = read_car(file_to_be_converted)[0]
    num_atoms = sum([int(x) for x in lines[6].split()])
    is_direct, is_select = is_direct_or_not(lines)
    ### To get the vector for conversion from direct to Cartesian
    if is_direct:  
        vector = get_vectors(lines)[0]
    else: 
        vector = np.array([[1, 0 , 0], [0, 1, 0], [0, 0, 1]])
    
    start_num = 9
    if not is_select:
        start_num = 8
        
    return  lines, vector, start_num, num_atoms, is_direct


def determinelayers(z_cartesian):
    threshold = 1.5
    seq = sorted(z_cartesian)
    n_layers = int((seq[-1]-seq[0])/threshold) + 1 
    
    layerscount = {}
    for layer in range(1, n_layers+1):
        list_layer = []
        for i, z in enumerate(z_cartesian):
            if  (layer - 1) * threshold <=  z - seq[0] < layer * threshold:
                list_layer.append(i)
        layerscount.update({layer:list_layer})
    return layerscount
#
    
def convert():
    x_cartesian = []
    y_cartesian = []
    z_cartesian = []
    tf = []
    for i in range(start_num, num_atoms + start_num):
        line_data =  [float(ele) for ele in lines[i].split()[0:3]]
        line_data = np.array([line_data])
        x, y, z = [sum(k) for k in line_data * vector ]     
        x_cartesian.append(x)
        y_cartesian.append(y)
        z_cartesian.append(z)
        tf.append(' T T T')
        
    layerscount =determinelayers(z_cartesian)
    print( '\n Find %s layers!'  %(len(layerscount)-1))
    
    file_out = open(file_to_be_converted, 'w')
    file_out.writelines(lines[0:7])  # first 7 lines are kept the same 
    
    file_out.write('Selective\n')  
    file_out.write('Cartesian' + '\n') 
    for i in range(1,len(layerscount)+1):
        for j in layerscount.get(i):
            if i <= fixedlayer: 
                tf[j] = ' F F F'
            else:
                tf[j] = ' T T T'
                
    for i in range(0,len(x_cartesian)):
        file_out.write("\t%+-3.10f   %+-3.10f   %+-3.10f  %s\n" %(x_cartesian[i], y_cartesian[i], z_cartesian[i], tf[i]))
        
    file_out.close()
    
lines, vector, start_num, num_atoms, is_direct = get_infor(file_to_be_converted)
file_back = open(file_to_be_converted+'_back', 'w')
file_back.writelines(lines)
file_back.close

if is_direct : 
    print("\n%s has Direct Coordinates, Contersion starts.... "  %(file_to_be_converted))
    convert()
else:
    print("\n%s has Cartesian Coordinates Already! We are going to fix layers only."  %(file_to_be_converted))
    convert()
    
print( '-----------------------------------------------------\n')
print( '\n %s Now has Cartesian Coordiates\n' %(file_to_be_converted))
print( '-----------------------------------------------------\n')

