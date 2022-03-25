#!/usr/bin/env python3
'''Generate the k_add file to show how the KPOINTS are splitted in the line model.'''
from kpoints import *

lines_k, num_pairs = read_kpoints_band()
lines_k_add = get_k_add_lines(lines_k, num_pairs)

f_out = open('k_add', 'w')
f_out.writelines(lines_k_add)
f_out.close()
       
