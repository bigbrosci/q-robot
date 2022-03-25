#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 16:12:54 2019

@author: qiang
"""

import numpy as np
from gpm import *
import sys

'''
A will be fixed
shift, double rotate B groups to match its structure to A
'''
file_in = sys.argv[1]
coords_A = file_analyzer('template.xyz')
coords_B = file_analyzer(file_in)

PA1 = get_coord_atom(coords_A, 22)
PA2 = get_coord_atom(coords_A, 27)
PA3 = get_coord_atom(coords_A, 24)

PB1 = get_coord_atom(coords_B, 1)

###Step1 shift rings
vector_shift = PA1 - PB1    
coords_B_S = shift_coords(coords_B, vector_shift)

# summarize the rotation information
'''Use p1, p2, p3 to define the plane '''
p1 = PA1
p2 = PA2
p3 = get_coord_atom(coords_B_S, 6)
theta_1 = get_angle((p3-p1), (p2-p1)) ## Rotate angle
nv = get_normal_vector(p1, p2, p3)    ## Rotate axis 

atom_list = get_atom_list_g09(coords_B, ['1-']) # Rotate atom list

coords_B_R1 = rotate_atoms(coords_B_S, PA1, PA1+nv, atom_list, theta_1)
#save_xyz(coords_B_R1, 'B1')
#
PB2  = get_coord_atom(coords_B_R1, 6) 
PB3 = get_coord_atom(coords_B_R1, 2)
theta_2 = get_angle(PA3-PA1, PB3-PA1)

if theta_2 > 1.5707963268:
    PB3 = get_coord_atom(coords_B_R1, 3)
    theta_2 = get_angle(PA3-PA1, PB3-PA1)

print(theta_2)    
coords_B_R2 = rotate_atoms(coords_B_R1, PA1, PB2, atom_list, theta_2)

final_coords=coords_A[:21] + coords_B_R2

save_xyz(final_coords, 'B2')
