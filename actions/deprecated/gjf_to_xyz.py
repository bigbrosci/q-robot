#!/usr/bin/env python3
# -*- coding: utf-8 -*-
'''
Written by Qiang in UD, 07-02-2019
This script will Update the gjf file, the old one will be covered. (be careful)
''' 
from gpm import * 
import sys 

script, file_gjf = sys.argv[:]
name, lines = read_file(file_gjf)
coords = get_coord_gjf(lines)
save_xyz(coords, name)
