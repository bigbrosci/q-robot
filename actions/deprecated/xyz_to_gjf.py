#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul  3 21:49:16 2019

@author: qiang
"""
import sys
from gpm import * 
from lattice import *

script , file_in = sys.argv[:]
name, lines = read_file(file_in)
coords = get_coord_xyz(lines)
save_gjf(coords, name)
