#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul  5 15:47:16 2019
@author: qiang
This script is used to merge the coordinates from different gjf or xyz files
"""

import sys 
from gpm import * 
files_in_list = sys.argv[1:]

coords = []
for file_in in files_in_list:
    coords_file_in = file_analyzer(file_in)
    coords += coords_file_in

save_xyz(coords, 'merged')    
