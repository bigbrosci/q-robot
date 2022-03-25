#!/usr/bin/env python3
import sys, os
from potcar import *

is_poscar_here = False
if os.path.isfile('POSCAR'):
    is_poscar_here = True
else:    
    print('%s\n' %('%' * 40))
    print('\nWARNING!!! No POSCAR in current folder!!\n' * 2 )
    print('%s\n' %('%' * 40))  
    
if len(sys.argv[:]) == 1: 
    if is_poscar_here:
        poscar =  open('POSCAR', 'r')
        lines = poscar.readlines()
        poscar.close()
        ele_list =  lines[5].rstrip().split()  ## Read the element line
        if ele_list[0].isdigit():  
            '''For VASP4 and ASE type POSCAR without elements in line 6'''
            print('\nASE type POSCAR\n')
            ele_list = lines[0].rstrip().split()
        concatenate(ele_list)
        if os.path.isfile('POTCAR'):
            read_potcar('POTCAR')
    else: 
        if os.path.isfile('POTCAR'):
            read_potcar('POTCAR')          
        else:     
            print('\nBut you can generate POTCAR by using command: pp.py element1 element2 ....\n')
if len(sys.argv)  > 1 :
    ele_list = sys.argv[1:] 
    concatenate(ele_list)
    if os.path.isfile('POTCAR'):
        read_potcar('POTCAR')          