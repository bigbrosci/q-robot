#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spyder Editor

This is a temporary script file.
"""
import numpy as np
from ase.io import read, write
import os

#os.chdir('/home/win.udel.edu/qli/Desktop/freq/')

## Read POSCAR
print('Read POSCAR >>>\t')
model = read('POSCAR')
model_positions = model.get_positions()
print('Read POSCAR DONE >>>\t')
#Read OUTCAR and get the line number of the largest imaginary frequency
print('Read OUTCAR >>>\t')
l_position = 0
with open('OUTCAR') as f_in:
    lines = f_in.readlines()
    wave_num = 0.0
    for num, line in enumerate(lines):
        if 'f/i' in line:
            wave_tem = float(line.rstrip().split()[6])
            if wave_tem > wave_num:
                wave_num = wave_tem 
                l_position = num+2         
      
vib_lines = lines[l_position:l_position+len(model)]
vib_dis = []
for line in vib_lines:
    infor = [float(i) for i in line.rstrip().split()[3:]]
    vib_dis.append(infor)
vib_dis = np.array(vib_dis)

print('Read OUTCAR DONE >>>\t')         

print("Start Correcting the xyz coordinates with a factor of 0.1")
# 0.3 is the displacement factor to add to the poscar.
new_positions = model_positions + vib_dis * 0.1 
model.positions = new_positions
write('POSCAR_corrected', model, vasp5=True)
print("Done! Output file is named as: POSCAR_corrected")
