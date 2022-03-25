#!/usr/bin/env python3
'''Can be used to convert direct to cartesian format  '''
import sys 
import ase 
from  ase.io import read
from  ase.io import write 
file_in = sys.argv[1]

model = ase.io.read(file_in)
#model.center()
ase.io.write('POSCAR', model, format='vasp', vasp5=True)
