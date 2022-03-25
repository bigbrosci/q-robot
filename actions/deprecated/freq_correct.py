#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
To check if there is imaginary frequencies in the log file
"""
import sys 
import numpy as np
from gpm import *

file_in, scale  = sys.argv[1:3]
scale = float(scale)

new_coords = freq_correct(file_in, scale)

new_name = file_in.split('.')[0]+'_new'
save_gjf(new_coords, new_name)
