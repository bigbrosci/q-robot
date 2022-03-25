#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Sep 30 07:42:44 2019

@author: qli
"""
from gpm import *
import sys 
coord_file = sys.argv[1]

coords = file_analyzer(coord_file)
#cc_single, cc_double, co_single, ch_single, oh_single, oo_hbond, oh_long_bond = get_dis_all(coords) 
calc_hbond_num1(coords)  
calc_hbond_num2(coords)        
#print(oo_hbond)    
        
    
    


    
        
    

