#!/usr/bin/env python3
'''Can be used to convert direct to cartesian format  '''
import sys 
import ase 
from  ase.io import read
from  ase.io import write 
import numpy as np


file_in = 'CONTCAR'

model = ase.io.read(file_in)

z_Ru = 4.1829

# list_C = [112, 110, 109, 111, 113, 114]

list_C = [int(i) for i in sys.argv[1:]]

z_C = np.array([model.get_positions()[i-1][2]  for i in list_C])

diss_Ru_C = np.average(z_C - z_Ru)

print(diss_Ru_C)
