#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  9 14:09:12 2019
@author: qiang
This script can convert the gjf, com, and xyz files to POSCAR
"""

import sys 
from gpm import * 

file_in = sys.argv[1]

save_poscar_from_extended_xyz(file_in)

