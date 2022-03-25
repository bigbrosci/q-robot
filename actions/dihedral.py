#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 17 15:06:37 2019

@author: qiang
"""

from gpm import *
import sys 
file_in = sys.argv[1]
points = sys.argv[2:6]

coords = file_analyzer(file_in)

p = [get_coord_atom(coords, int(i)) for i in points]

dihedral = get_dihedral(p)

dihedral = abs(dihedral)
if dihedral > 90:
    dihedral = 180 - dihedral
    
print(dihedral)
