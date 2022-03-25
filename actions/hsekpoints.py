#!/usr/bin/env python3
''' 
This script will Generate the KPOINTS file for HSE band calculations.
For HSE band calculation, you need to 
Step1: Perform the GGA-PBE band structure calculation first to generate the IBZKPT file.
Step2: copy the KPOINTS in step1 to KPOINTS_band, cp KPOINTS KPOINTS_band
step3: run this script to generate the KPOINTS for HSE calculation. IBZKPT and KPOINTS_band are both needed.
'''
from kpoints import *
import sys, os

if not os.path.isfile('KPOINTS_band') or  not os.path.isfile('IBZKPT') :
    print('''Error: Can not find KPOINTS_band/IBZKPT file. Do the following and rerun this script:
          0) read the header of this script;
          1) Prepare the KPOINTS file for standard PBE calculation ;
          2) Rename  this file as KPOINTS_band. ''')
    exit() 

lines_k, num_pairs = read_kpoints_band('KPOINTS_band')
lines_k_add = get_k_add_lines(lines_k, num_pairs)
num_k_add   = len(lines_k_add)

lines_ib, num_lines_ib, num_points_ib = read_ibzkpt() 

### Write KPOINTS filr for HSE 
f_k = open('KPOINTS', 'w')
f_k.writelines(lines_ib[0])
f_k.write(' %s \n' %(num_k_add + num_points_ib))
f_k.writelines(lines_ib[2:])
for line in lines_k_add:
    line_weight = line.strip() + '  0.000\n'
    f_k.write(line_weight)
f_k.close()

### Write k_add file 
f_out = open('k_add', 'w')
f_out.writelines(lines_k_add)
f_out.close()