#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Oct  2 10:54:36 2019
expand the cell to a larger one
@author: qli
"""
import sys, os
import ase
from ase.io import read, write
import ase.io.vasp

file_read = sys.argv[1]
x,y,z     = [int(i) for i in sys.argv[2:5]]
cell = read(file_read, format="vasp")
write("POSCAR_ex",cell*(x, y, z), direct=True,sort=True)

