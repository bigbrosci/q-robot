#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
@author: qli
Created on Thu Jun  2 22:09:18 2022
Creat POSCAR for improved dimer calculations
1) run frequency calculation: 
IBRION = 5
POTIM = 0.015
NFREE = 2
NWRITE = 3 ### Must be 3
2) run this script:
get_dimer.py 
3) use the POSCAR for IDM calc.
NSW = 100           
Prec=Normal
IBRION=44           !  use the dimer method as optimization engine
EDIFFG=-0.05
POTIM = 0.05
    
"""

import numpy as np
from ase.io import read, write
import os
from sys import exit


# os.chdir('/home/win.udel.edu/qli/Desktop/freq/')

model = read('POSCAR')
model_positions = model.get_positions()
model.write('POSCAR_dimer', vasp5=True)

# print(model_positions)
# print(len(model))

l_start = 0 # the number of line which contains  'Eigenvectors after division by SQRT(mass)' 
with open('OUTCAR') as f_in:
    lines = f_in.readlines()
    for num, line in enumerate(lines):
        if 'Eigenvectors after division by SQRT(mass)' in line:
            l_start = num

if l_start == 0:
    print('''Check Frequency results and then rerun this script.\n**Remember**: NWRITE must be 3. BYEBYE! ''' )
    exit()
    
freq_infor_block = lines[l_start:]        
l_position = 0
wave_num = 0.0
for num, line in enumerate(freq_infor_block): 
        if 'f/i' in line:           
            wave_tem = float(line.rstrip().split()[6])
            if wave_tem > wave_num:
                wave_num = wave_tem  
                l_position = num+2

pos_dimer = open('POSCAR_dimer', 'a')
pos_dimer.write('  ! Dimer Axis Block\n')

vib_lines = freq_infor_block[l_position:l_position+len(model)]
for line in vib_lines:
    infor = line.rstrip().split()[3:]
    pos_dimer.write(' '.join(infor)+'\n')

pos_dimer.close()
print('''
      DONE!
      Output file is named as: POSCAR_dimer and can be used for dimer calculations.
      Don't forget to rename POSCAR_dimer to POSCAR before you run the dimer jobs.      
      ''')