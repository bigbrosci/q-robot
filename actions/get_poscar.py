#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 14:09:12 2019
@author: qiang
This script can convert the gjf, com, and xyz files to POSCAR
"""

import sys 
from gpm import * 

file_in, name = sys.argv[1:]

coords = file_analyzer(file_in)
save_poscar(coords, name)
