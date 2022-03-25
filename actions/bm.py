#!/usr/bin/env python3
#Get Lattice paramters from BM state equation
# Written by Qiang 
# To use it : python bm.py data

import sys 
from lattice import bm_fitting
data_file = sys.argv[1]

bm_fitting(data_file)
