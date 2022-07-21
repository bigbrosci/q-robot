#!/usr/bin/env python3
import sys
import ase
from ase.build.tools import sort
from  ase.io import read
from  ase.io import write
file_in = sys.argv[1]
 
model = ase.io.read('POSCAR', format='vasp')
sort(model, tags=model.get_chemical_symbols())
ase.io.write('POSCAR_sorted_by_ele', model, format='vasp', vasp5=True)
